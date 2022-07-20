import numpy as np
import sys
import os
from multiprocessing import Pool
import pyximport; pyximport.install()
import cyemdORB as cy

indir=sys.argv[1]
outdir=sys.argv[2]
n=int(sys.argv[3])

def KSAQ(indir,outdir,n):
    os.mkdir(outdir)
    queries=cy.get_queries(indir)
    orbs=range(73)
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
        else:
            Ms=Ms+K
        if i==3:
            cy.toM(Ms/(i+1),queries,'NetEmd_G3_'+indir)
        elif i==14:
            cy.toM(Ms/(i+1),queries,'NetEmd_G4_'+indir)
        elif i==72:
            cy.toM(Ms/(i+1),queries,'NetEmd_G5_'+indir)
    return 1

KSAQ(indir,outdir,n)
