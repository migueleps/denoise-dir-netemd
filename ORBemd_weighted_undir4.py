import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install(language_level=3)
import cyemdORB_W as cy
import shutil

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

def KSAQ(indir,outdir,n):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    queries=cy.get_queries(indir)
    orbs=range(15)
    V=[[queries,orb] for orb in orbs]
    p = Pool(n)
    c=p.imap(cy.MKSAP,V)
    p.close()
    p.join()
    for i,K in enumerate(c):
        cy.toM(K[0],queries,outdir+'/NetEmd_Orb'+str(i)+indir)
        if i==0:
            Ms=K[0]
            Ss=K[1]
        else:
            Ms=Ms+K[0]
            Ss=Ss+K[1]
        if i==3:
            cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_G3_'+indir+'_weighted')
        elif i==14:
            cy.toM(comp_sol(Ms,Ss,len(queries)),queries,'NetEmd_G4_'+indir+'_weighted')
    return 1

if __name__ == '__main__':
    KSAQ(indir,outdir,n)
