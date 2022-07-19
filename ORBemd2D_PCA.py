
# coding: utf-8

# In[1]:

import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORBD_PCA as cy
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# In[24]:

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])
norbs=3
perc_var=float(sys.argv[4])

def join_counts(queries):
    for q in queries:
        q2d = q + '.countsO2d'
        qf = q + '.countsO'
        os.system('mv %s %s' % (q2d,qf))


def pca_rewrite(queries):
    for f in queries:
        raw_counts = pd.read_csv(f+".countsO",sep=" ",header=None,index_col=False,usecols=range(norbs))
        data = raw_counts.values
        z_scaler = StandardScaler()
        z_data = z_scaler.fit_transform(data)
        pca_ncomps = PCA().fit(z_data)
        nc = 0
        exp_var = 0
        while exp_var < perc_var:
            exp_var += pca_ncomps.explained_variance_ratio_[nc]
            nc += 1
        pca = PCA(n_components = nc).fit(z_data)
        reduced_data = pca.transform(z_data)
        reconstructed_zdata = pca.inverse_transform(reduced_data)
        reconstructed_data = z_scaler.inverse_transform(reconstructed_zdata)
        reconstructed_raw_counts = pd.DataFrame(data = reconstructed_data,columns = list(raw_counts))
        reconstructed_raw_counts.to_csv(f+".countsO",index=False,sep=" ",header=False)
        
# In[3]:

def KSAQ(indir,outdir,n):
    queries=cy.get_queries(indir)
    join_counts(queries)
    orbs=range(norbs)
    pca_rewrite(queries)
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
        if i==2:
            cy.toM(Ms/(i+1),queries,'NetEmd_G2D_'+indir+"_PCA")
    return 1
        
        
    
    


# In[ ]:

KSAQ(indir,outdir,n)


# In[ ]:



