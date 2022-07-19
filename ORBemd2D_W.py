
# coding: utf-8

# In[1]:

import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORBD_W as cy


# In[24]:

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])

def comp_sol(Ms,Ss,n):
    for i in range(n):
        for j in range(n):
            if Ms[i][j] == 0:
                continue
            elif Ss[i][j] == 0:
                print i,j,Ms[i][j]
            else:
                Ms[i][j] = Ms[i][j]/Ss[i][j]
    return Ms

def join_counts(queries):
    for q in queries:
        q2d = q + '.countsO2d'
        qf = q + '.countsO'
        os.system('mv %s %s' % (q2d,qf))

# In[3]:

def KSAQ(indir,outdir,n):
    queries=cy.get_queries(indir)
    join_counts(queries)
    orbs=range(3)
    V=[[queries,orb] for orb in orbs]
    if __name__ == '__main__':
        p = Pool(n)
        c=p.imap(cy.MKSAP,V)
        p.close()
        p.join()
    for i,K in enumerate(c):
        if i==0:
            Ms=K[0]
            Ss=K[1]
        elif i==2:
            Ms=Ms+K[0]
            Ss=Ss+K[1]
            cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_G2D_'+indir)
        else:
            Ms=Ms+K[0]
            Ss=Ss+K[1]
    return 1
        
        
    
    


# In[ ]:

KSAQ(indir,outdir,n)


# In[ ]:



