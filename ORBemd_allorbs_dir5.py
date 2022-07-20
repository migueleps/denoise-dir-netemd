import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install(language_level=3)
import cyemdORBD as cy

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
        os.system("paste -d '\\0' {} {} {} {} > {}".format(q2d,q3d,q4d,q5d,qf))

def KSAQ(indir,outdir,n):
    queries=cy.get_queries(indir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    join_counts(queries)
    orbs=range(45637)
    V=[[queries,orb] for orb in orbs]
    p = Pool(n)
    c=p.imap(cy.MKSAP,V)
    p.close()
    p.join()
    for i,K in enumerate(c):
        cy.toM(K,queries,outdir+'/NetEmd_Orb'+str(i)+indir)
        if i==0:
            Ms=K
        else:
            Ms=Ms+K
        if i==2:
            cy.toM(Ms/(i+1),queries,'NetEmd_G2D_'+indir)
        elif i==32:
            cy.toM(Ms/(i+1),queries,'NetEmd_G3D_'+indir)
        elif i==729:
            cy.toM(Ms/(i+1),queries,'NetEmd_G4D_'+indir)
        elif i==45636:
            cy.toM(Ms/(i+1),queries,'NetEmd_G5D_'+indir)
    return 1

if __name__ == '__main__':
    KSAQ(indir,outdir,n)
