

# coding: utf-8

# In[1]:

import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORBD_W as cy
import shutil

# In[24]:

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])


def comp_sol(Ms,Ss,len_queries):
    new_Ms = np.copy(Ms)
    for i in range(len_queries):
        for j in range(len_queries):
            if Ss[i][j] == 0 and new_Ms[i][j] != 0:
                print("something wrong!",i,j)
            if new_Ms[i][j] == 0:
                continue
            else:
                new_Ms[i][j] = new_Ms[i][j]/Ss[i][j]
    return new_Ms

def join_counts(queries):
    for q in queries:
        q2d = q + '.countsO2d'
        q3d = q + '.countsO3d'
        q4d = q + '.countsO4d'
        qf = q + '.countsO'
        os.system("paste -d '\\0' %s %s %s > %s" % (q2d,q3d,q4d,qf))


# In[3]:

def KSAQ(indir,outdir,n):
#    if os.path.isdir(outdir):
#        shutil.rmtree(outdir)
#    os.mkdir(outdir)
    queries=cy.get_queries(indir)
    join_counts(queries)
    orbs=range(730)
    V=[[queries,orb] for orb in orbs]
    if __name__ == '__main__':
        p = Pool(n)
        c=p.imap(cy.MKSAP,V)
        p.close()
        p.join()
    for i,K in enumerate(c):
#        cy.toM(K[0],queries,outdir+'/NetEmd_Orb'+str(i)+indir)
        if i==0:
            Ms=K[0]
            Ss=K[1]
        else:
            Ms=Ms+K[0]
            Ss=Ss+K[1]
        if i==2:
            cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_G2D_'+indir)
        elif i==32:
            cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_G3D_'+indir)    
        elif i==729:
            cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_G4D_'+indir)
    return 1
        
        
    
    


# In[ ]:

KSAQ(indir,outdir,n)


# In[ ]:


