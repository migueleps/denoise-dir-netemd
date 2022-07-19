
# coding: utf-8

# In[1]:

import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORB_PCA as cy
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import time

# In[24]:

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])
ncomps=int(sys.argv[4])
norbs=15

# In[3]:

def ica_rewrite(queries):
    for f in queries:
        raw_counts = pd.read_csv(f+".countsO",sep="\t",index_col=False,usecols = range(norbs))
        data = raw_counts.values
        colmeans = data.mean(axis=0)
        centered_data = data - colmeans
        ica = FastICA(n_components=ncomps,max_iter=1000).fit(centered_data)
        reduced_data = ica.transform(centered_data)
        reconstructed_data = ica.inverse_transform(reduced_data)
        reconstructed_uncentered_data = reconstructed_data + colmeans
        reconstructed_raw_counts = pd.DataFrame(data = reconstructed_uncentered_data,columns = list(raw_counts))
        reconstructed_raw_counts.to_csv(f+".countsO",index=False,sep="\t")


def KSAQ(indir,outdir,n):
    queries=cy.get_queries(indir)
    t = time.time()
    ica_rewrite(queries)
    print("Type to compute ICA rewrite: {}".format(time.time()-t))
    orbs=range(norbs)
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
        if i==3:
            cy.toM(Ms/(i+1),queries,'NetEmd_G3_4W_{}_ICA{}comps'.format(indir,ncomps))
        elif i==14:
            cy.toM(Ms/(i+1),queries,'NetEmd_G4_{}_ICA{}comps'.format(indir,ncomps))

    return 1
        

KSAQ(indir,outdir,n)




