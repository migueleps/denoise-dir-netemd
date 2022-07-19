
# coding: utf-8

# In[1]:

import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORB4 as cy


# In[24]:

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])


# In[3]:

def KSAQ(indir,outdir,n):
    os.mkdir(outdir)
    queries=cy.get_queries(indir)
    orbs=range(15)
    V=[[queries,orb] for orb in orbs]
    if __name__ == '__main__':
        p = Pool(n)
        c=p.imap(cy.MKSAP,V)
        p.close()
        p.join()
    for i,K in enumerate(c):
        cy.toM(K,queries,outdir+'/NetEmd_Orb'+str(i)+indir)
	if i==0:
            Ms=K
        elif i==3:
            Ms=Ms+K
            cy.toM(Ms/(i+1),queries,'NetEmd_G3_'+indir)
        elif i==14:
            Ms=Ms+K
            cy.toM(Ms/(i+1),queries,'NetEmd_G4_'+indir)
        elif True:
            Ms=Ms+K
    return 1
        
    
    


# In[ ]:

KSAQ(indir,outdir,n)


# In[ ]:



