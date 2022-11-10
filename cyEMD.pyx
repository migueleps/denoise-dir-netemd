
from cpython cimport array
import array
from scipy.optimize import minimize_scalar
from numpy import zeros
import numpy as np
import os
import math

def dis2(s): #Same as dis but now allows count values in R instead of N, which is useful to use the reconstruction from projecting the data using PCA
    vals=sorted(list(set(map(math.floor,s))))
    m=sum(s)/len(s)
    smod=list(map(math.floor,s))
    return [[i-m,i+1-m,smod.count(i)/float(len(smod))] for i in vals]

def dis(s): #maps the empirical values (s) to a normalized histogram with bins of with 1 and centers this to have mean 0. This is then used as the input distribution for NetEmd.
    vals=sorted(list(set(s)))
    m=sum(s)/len(s)
    return [[i-m,i+1-m,s.count(i)/float(len(s))] for i in vals]

def mean(d):
    return sum([((i[0]+i[1])/2.0)*i[2] for i in d])

def variance(d): #computes variance of distribution d. Note that this differs from the variance in the discrete case since d is assumed to be a piecewise uniform/continious pdf. For instance compare X Bernoulli (1/2), so X puts equal mass on 0 and 1. The mean is 1/2, the variance is 1/4; To Y continuous uniform[0,1]: the mean is again 1/2, the variance is 1/12.
    return sum([(i[0]**2+i[1]**2+i[0]*i[1])*i[2]/3.0 for i in d])-mean(d)**2

def rescale(d): #rescale d
    m=variance(d)**(-0.5)
    return [[i[0]*m,i[1]*m,i[2]] for i in d]

def cum(d): #computes the cummulative of d which is a piece-wise linear function which we encode as the coordinates and values at the points where the slope changes.
    cm=[[d[0][0],0.0]]
    dlast=d[0][0]
    c=0
    for i in d:
        if dlast!=i[0]:
            cm.append([i[0],c])
        c+=i[2]
        cm.append([i[1],c])
        dlast=i[1]
    return np.array(cm,dtype='f')

cdef float cP(float [:,:] c, int ind, float p): #computes value of the cummulative c at point p
    if p<=c[0][0]:
        return 0.0
    if p>=c[-1][0]:
        return 1.0
    else:
        return c[ind][1]+((p-c[ind][0])/(c[ind+1][0]-c[ind][0]))*(c[ind+1][1]-c[ind][1])

cdef float KSA(float [:,:] c1, float [:,:] c2): #computes EMD between cdfs c1 and c2
    cdef int ind1=-1
    cdef int ind2=-1
    cdef int i
    cdef float a=0.0
    cdef float x2=0.0
    cdef float x1
    cdef float pp
    cdef int l1=len(c1)
    cdef int l2=len(c2)
    cdef int lt=l1+l2-1
    if c1[0][0]<c2[0][0]:
        ind1+=1
        pp=c1[0][0]
    else:
        ind2+=1
        pp=c2[0][0]
    x2=0.0
    for i in range(lt):
        x1=x2
        if ind1+1==l1 or ind2+1==l2:
            if ind1+1==l1:
                ind2+=1
                x2=1-c2[ind2][1]
                a+=abs(x1+x2)*(c2[ind2][0]-pp)*0.5
                pp=c2[ind2][0]
            else:
                ind1+=1
                x2=c1[ind1][1]-1
                a+=abs(x1+x2)*(c1[ind1][0]-pp)*0.5
                pp=c1[ind1][0]
        else:
            if c1[ind1+1][0]<c2[ind2+1][0]:
                ind1+=1
                x2=c1[ind1][1]-cP(c2,ind2,c1[ind1][0])
                if x1*x2>=0:
                    a+=abs(x1+x2)*(c1[ind1][0]-pp)*0.5
                else:
                    a+=0.5*(c1[ind1][0]-pp)*(x1**2+x2**2)/(abs(x1)+abs(x2))
                pp=c1[ind1][0]
            else:
                ind2+=1
                x2=cP(c1,ind1,c2[ind2][0])-c2[ind2][1]
                if x1*x2>=0:
                    a+=abs(x1+x2)*(c2[ind2][0]-pp)*0.5
                else:
                    a+=0.5*(c2[ind2][0]-pp)*(x1**2+x2**2)/(abs(x1)+abs(x2))
                pp=c2[ind2][0]
    return a

cdef float [:,:] shift(float [:,:] c,float t): #shifts c by t
    cdef int i
    cdef int lc=len(c)
    cdef float[:,:] dd=c.copy()
    for i in range(lc):
        dd[i][0]+=t
    return dd

cdef float KSAt(float [:,:] c1,float [:,:] c2,float t): #computes EMD between c1 and c2+t
    return KSA(shift(c1,t),c2)

def KSAPS(float [:,:] d1, float [:,:] d2): #find t that minimizes EMD and return EMD*(c1,c2)=EMD(c1+t,c2)
    cdef float tt
    cdef float r
    r=minimize_scalar(lambda tt: KSAt(d1,d2,tt),method='brent',tol=0.00001, options={'maxiter': 150}).fun
    return r

def toM(K,queries,n_file): #write K to file
    f = open(n_file,'w')
    f.write('\t')
    for i in range(len(queries)):
            f.write(str(queries[i])+'\t')
    f.write('\n')
    for i in range(len(queries)):
            f.write(str(queries[i])+'\t')
            for j in range(len(queries)):
                    f.write(str(round(K[i][j],6))+'\t')
            f.write('\n')
    f.close()
    return 1

def readcount_undir(fn,orb): #get graphlet degree distribution for orbit orb
    f=open(fn,'r')
    fread= f.readlines()
    cnts=[]
    for line in fread[1:]:
        line=line.split('\t')
        cnts.append(float(line[orb]))
    return cnts

def readcount_dir(fn,orb): #get graphlet degree distribution for orbit orb
    f=open(fn,'r')
    fread= f.readlines()
    cnts=[]
    for line in fread:
        line=line.split(' ')
        cnts.append(float(line[orb]))
    return cnts

def get_queries(indir):
    queries=[]
    for root, dirs, filenames in os.walk(indir):
        for f in filenames:
            if not ".count" in f and not 'temp' in f:
                queries.append(os.path.join(root, f))
    queries.sort()
    return queries

def get_cum_undir_dimred(queries,orb): #get cdfs for orbit orb
    return [cum(rescale(dis2(readcount_undir(queries[i]+'.countsO',orb)))) for i in range(len(queries))]

def get_cum_dir_dimred(queries,orb): #get cdfs for orbit orb
    return [cum(rescale(dis2(readcount_dir(queries[i]+'.countsO',orb)))) for i in range(len(queries))]

def get_cum_undir(queries,orb): #get cdfs for orbit orb
    return [cum(rescale(dis(readcount_undir(queries[i]+'.countsO',orb)))) for i in range(len(queries))]

def get_cum_dir(queries,orb): #get cdfs for orbit orb
    return [cum(rescale(dis(readcount_dir(queries[i]+'.countsO',orb)))) for i in range(len(queries))]

def orb_exists_undir(queries,orb):
    n=len(queries)
    S=zeros((n,n))
    ss = [sum(readcount_undir(queries[i]+'.countsO',orb)) for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if ss[i] > 0 or ss[j] > 0:
                S[i][j] = 1
                S[j][i] = 1
    return S

def orb_exists_dir(queries,orb):
    n=len(queries)
    S=zeros((n,n))
    ss = [sum(readcount_dir(queries[i]+'.countsO',orb)) for i in range(n)]
    for i in range(n):
          for j in range(i+1,n):
              if ss[i] > 0 or ss[j] > 0:
                  S[i][j] = 1
                  S[j][i] = 1
    return S


def MKSA(queries,orb): #computes NetEmd between queries for orbit orb
    n=len(queries)
    C=get_cum_undir(queries,orb)
    K=zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            t=KSAPS(C[i],C[j])
            K[i][j]=t
            K[j][i]=t
    return K

def MKSA_D(queries,orb): #computes NetEmd between queries for orbit orb
    n=len(queries)
    C=get_cum_dir(queries,orb)
    K=zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            t=KSAPS(C[i],C[j])
            K[i][j]=t
            K[j][i]=t
    return K

def MKSA_W(queries,orb): #computes NetEmd between queries for orbit orb
    n=len(queries)
    C=get_cum_undir(queries,orb)
    S=orb_exists_undir(queries,orb)
    K=zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            if S[i][j] == 0:
              continue
            t=KSAPS(C[i],C[j])
            K[i][j]=t
            K[j][i]=t
    return K, S

def MKSA_WD(queries,orb): #computes NetEmd between queries for orbit orb
    n=len(queries)
    C=get_cum_dir(queries,orb)
    S=orb_exists_dir(queries,orb)
    K=zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            if S[i][j] == 0:
              continue
            t=KSAPS(C[i],C[j])
            K[i][j]=t
            K[j][i]=t
    return K, S

def MKSA_dimred(queries,orb): #computes NetEmd between queries for orbit orb
    n=len(queries)
    C=get_cum_undir_dimred(queries,orb)
    K=zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            t=KSAPS(C[i],C[j])
            K[i][j]=t
            K[j][i]=t
    return K

def MKSA_dimredD(queries,orb): #computes NetEmd between queries for orbit orb
    n=len(queries)
    C=get_cum_dir_dimred(queries,orb)
    K=zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            t=KSAPS(C[i],C[j])
            K[i][j]=t
            K[j][i]=t
    return K

def MKSAP_(V):
    return MKSA(V[0],V[1])

def MKSAP_D(V):
    return MKSA_D(V[0],V[1])

def MKSAP_W(V):
    return MKSA_W(V[0],V[1])

def MKSAP_WD(V):
    return MKSA_WD(V[0],V[1])

def MKSAP_dimred(V):
    return MKSA_dimred(V[0],V[1])

def MKSAP_dimredD(V):
    return MKSA_dimredD(V[0],V[1])
