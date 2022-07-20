import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install(language_level=3)
import cyemdORBD_PCA as cy
import shutil
import time
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import time

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])
ncomps=int(sys.argv[4])
norbs = 730

def join_counts(queries):
    for q in queries:
        q2d = q + '.countsO2d'
        q3d = q + '.countsO3d'
        q4d = q + '.countsO4d'
        qf = q + '.countsO'
        os.system("paste -d '\\0' {} {} {} > {}".format(q2d,q3d,q4d,qf))


def ica_rewrite(queries):
    for f in queries:
        raw_counts = pd.read_csv(f+".countsO",sep=" ",header=None,index_col=False,usecols=range(norbs))
        data = raw_counts.values
        colmeans = data.mean(axis=0)
        data = data - colmeans
        ica = FastICA(n_components=ncomps,max_iter=1000).fit(data)
        reduced_data = ica.transform(data)
        reconstructed_data = ica.inverse_transform(reduced_data)
        reconstructed_uncentered_data = reconstructed_data + colmeans
        reconstructed_raw_counts = pd.DataFrame(data = reconstructed_uncentered_data,columns = list(raw_counts))
        reconstructed_raw_counts.to_csv(f+".countsO",index=False,sep=" ",header=False)


def KSAQ(indir,outdir,n):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    queries=cy.get_queries(indir)
    join_counts(queries)
    t = time.time()
    ica_rewrite(queries)
    print("Type to compute ICA rewrite: {}".format(time.time()-t))
    orbs=range(norbs)
    V=[[queries,orb] for orb in orbs]
    p = Pool(n)
    c=p.imap(cy.MKSAP,V)
    p.close()
    p.join()
    for i,K in enumerate(c):
        cy.toM(K,queries,"{}/NetEmd_Orb{}{}_ICA{}COMPS".format(outdir,i,indir,ncomps))
        if i == 0:
            Ms = K
        else:
            Ms += K
        if i==32:
            cy.toM(Ms/(i+1),queries,'NetEmd_G3D_{}_ICA{}COMPS'.format(indir,ncomps))
        elif i==729:
            cy.toM(Ms/(i+1),queries,'NetEmd_G4D_{}_ICA{}COMPS'.format(indir,ncomps))
    return 1

if __name__ == '__main__':
    KSAQ(indir,outdir,n)
