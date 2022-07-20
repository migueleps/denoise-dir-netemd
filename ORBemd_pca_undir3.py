import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORB_PCA as cy
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import time

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])
norbs=4
expvar=float(sys.argv[4])

def pca_rewrite(queries):
    for f in queries:
        raw_counts = pd.read_csv(f+".countsO",sep="\t",usecols=range(norbs))
        data = raw_counts.values
        z_scaler = StandardScaler()
        z_data = z_scaler.fit_transform(data)
        pca_ncomps = PCA().fit(z_data)
        nc = 0
        exp_var = 0
        while exp_var < expvar:
            exp_var += pca_ncomps.explained_variance_ratio_[nc]
            nc += 1
        pca = PCA(n_components = nc).fit(z_data)
        reduced_data = pca.transform(z_data)
        reconstructed_zdata = pca.inverse_transform(reduced_data)
        reconstructed_data = z_scaler.inverse_transform(reconstructed_zdata)
        reconstructed_raw_counts = pd.DataFrame(data = reconstructed_data,columns = list(raw_counts))
        reconstructed_raw_counts.to_csv(f+".countsO",index=False,sep="\t")


def KSAQ(indir,outdir,n):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    queries=cy.get_queries(indir)
    start = time.time()
    pca_rewrite(queries)
    print("Time taken to rewrite orbits with pca (s):", time.time() - start)
    orbs=range(norbs)
    V=[[queries,orb] for orb in orbs]
    if __name__ == '__main__':
        p = Pool(n)
        c=p.imap(cy.MKSAP,V)
        p.close()
        p.join()
    for i,K in enumerate(c):
        cy.toM(K,queries,"{}/NetEmd_Orb{}{}_PCA{}EXPVAR".format(outdir,i,indir,int(expvar*100)))
        if i==0:
            Ms=K
        else:
            Ms=Ms+K
        if i==3:
            cy.toM(Ms/(i+1),queries,'NetEmd_G3_{}_PCA{}EXPVAR'.format(indir,expvar))

    return 1

KSAQ(indir,outdir,n)
