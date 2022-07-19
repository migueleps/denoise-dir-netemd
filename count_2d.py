import networkx as nx
import os
import sys

g = nx.read_edgelist(sys.argv[1],data=False,create_using=nx.DiGraph(),nodetype=int)
outf = open(sys.argv[2],"w")

for n in g.nodes():
    pred = set(g.predecessors(n))
    succ = set(g.successors(n))
    recip = pred & succ
    outf.write("%d %d %d \n" % (len(succ) - len(recip), len(pred) - len(recip), len(recip)))

    
