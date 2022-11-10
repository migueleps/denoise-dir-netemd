import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install(language_level=3)
import cyEMD as cy
import shutil
import time
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import time

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

def join_counts(queries, graphlet_size):
    for q in queries:
        paste_cmd = 'paste -d "\\0" '
        for gsize in range(2,graphlet_size+1):
            paste_cmd = '{}{}.countsO{}d '.format(paste_cmd,q,gsize)
        paste_cmd = '{}> {}.countsO'.format(paste_cmd,q)
        os.system(paste_cmd)

def pca_rewrite(queries, directed_flag, norbs, expvar):
    for f in queries:
        if directed_flag:
            raw_counts = pd.read_csv("{}.countsO".format(f),sep=" ",header=None,index_col=False,usecols=range(norbs))
        else:
            raw_counts = pd.read_csv("{}.countsO".format(f),sep="\t",usecols=range(norbs))
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
        if directed_flag:
            reconstructed_raw_counts.to_csv("{}.countsO".format(f),index=False,header=False,sep=" ")
        else:
            reconstructed_raw_counts.to_csv("{}.countsO".format(f),index=False,sep="\t")

def ica_rewrite(queries, directed_flag, norbs, ncomps):
    for f in queries:
        if directed_flag:
            raw_counts = pd.read_csv("{}.countsO".format(f),sep=" ",header=None,index_col=False,usecols=range(norbs))
        else:
            raw_counts = pd.read_csv("{}.countsO".format(f),sep="\t",usecols=range(norbs))
        data = raw_counts.values
        colmeans = data.mean(axis=0)
        centered_data = data - colmeans
        ica = FastICA(n_components=ncomps,max_iter=1000,whiten="unit-variance").fit(centered_data)
        reduced_data = ica.transform(centered_data)
        reconstructed_data = ica.inverse_transform(reduced_data)
        reconstructed_uncentered_data = reconstructed_data + colmeans
        reconstructed_raw_counts = pd.DataFrame(data = reconstructed_uncentered_data,columns = list(raw_counts))
        if directed_flag:
            reconstructed_raw_counts.to_csv("{}.countsO".format(f),index=False,header=False,sep=" ")
        else:
            reconstructed_raw_counts.to_csv("{}.countsO".format(f),index=False,sep="\t")


graphlet_size_to_norbs = {True: {2: 3, 3: 33, 4: 730, 5: 45637},
                          False: {3: 4, 4: 15, 5: 73}}

print_breakpoints = {True: [2, 32, 729, 45636], False: [3, 14, 72]}

cython_func_choser = {True: {"allorbs": cy.MKSAP_D, "weighted": cy.MKSAP_WD, "ica": cy.MKSAP_dimredD, "pca": cy.MKSAP_dimredD},
                      False: {"allorbs": cy.MKSAP_, "weighted": cy.MKSAP_W, "ica": cy.MKSAP_dimred, "pca": cy.MKSAP_dimred}}

start_labels = {2: "G2D", 32: "G3D", 729: "G4D", 45636: "G5D", 3: "G3", 14: "G4", 72: "G5"}

def compute_dist(indir,outdir,n_threads,directed_flag,graphlet_size,netemd_type,extra_param = None):
    if not os.path.exists(outdir) and outdir != "":
        os.mkdir(outdir)
    queries=cy.get_queries(indir)
    norbs = graphlet_size_to_norbs[directed_flag][graphlet_size]
    end_labels = {"allorbs": "", "weighted": "_weighted", "ica": "_ICA{}comps".format(extra_param), "pca":"_PCA{}expvar".format(extra_param)}
    if directed_flag:
        join_counts(queries,graphlet_size)
    if netemd_type == "pca":
        pca_rewrite(queries,directed_flag,norbs,float(extra_param))
    elif netemd_type == "ica":
        ica_rewrite(queries,directed_flag,norbs,int(extra_param))
    orbs=range(norbs)
    V=[[queries,orb] for orb in orbs]
    p = Pool(n_threads)
    c=p.imap(cython_func_choser[directed_flag][netemd_type],V)
    p.close()
    p.join()
    for i,K in enumerate(c):
        if outdir != "":
            cy.toM(K[0],queries,outdir+'/NetEmd_Orb'+str(i)+indir)
        if i==0:
            if netemd_type == "weighted":
                Ms=K[0]
                Ss=K[1]
            else:
                Ms = K
        else:
            if netemd_type == "weighted":
                Ms=Ms+K[0]
                Ss=Ss+K[1]
            else:
                Ms = Ms+K
        if i in print_breakpoints[directed_flag]:
            if netemd_type == "weighted":
                cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_{}_{}{}'.format(start_labels[i],indir, end_labels[netemd_type]))
            else:
                cy.toM(Ms/(i+1),queries,'NetEmd_{}_{}{}'.format(start_labels[i],indir,end_labels[netemd_type]))
    return 1
