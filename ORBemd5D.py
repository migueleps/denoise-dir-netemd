

# coding: utf-8

# In[1]:

import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORBD as cy


# In[24]:

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])


def join_counts(queries):
    for q in queries:
        q2d = q + '.countsO2d'
        q3d = q + '.countsO3d'
        q4d = q + '.countsO4d'
        q5d = q + '.countsO5d'
        qf = q + '.countsO'
        os.system("paste -d '\\0' %s %s %s %s > %s" % (q2d,q3d,q4d,q5d,qf))


# In[3]:

def KSAQ(indir,outdir,n):
    queries=cy.get_queries(indir)
    join_counts(queries)
    orbs=range(45637)
    V=[[queries,orb] for orb in orbs]
    if __name__ == '__main__':
        p = Pool(n)
        c=p.imap(cy.MKSAP,V)
        p.close()
        p.join()
    for i,K in enumerate(c):
        if i==0:
            Ms=K
        else:
            Ms=Ms+K
        if i==45636:
            cy.toM(Ms/(i+1),queries,'NetEmd_G5D_'+indir)
    return 1
        
        
    
    


# In[ ]:

KSAQ(indir,outdir,n)


# In[ ]:



