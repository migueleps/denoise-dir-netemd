import os
import sys
import networkx as nx
from multiprocessing import Pool
import shutil

indir=sys.argv[1]
n=int(sys.argv[2])
queries=[]
for root, dirs, filenames in os.walk(indir):
	for f in filenames:
		if not ".count" in f:# and 'mod' not in f:
			queries.append(os.path.join(root, f))


def convertv(fn):
    f=open(fn,'r')
    fout=open(fn+'O1','w')
    for line in f.readlines():
        line=line.split(' ')
        fout.write(str(int(line[0])+1)+'\t'+str(int(line[1])+1)+'\n')
    f.close()
    fout.close()
    return 1
def convert(fn):
    G=nx.read_edgelist(fn,data=False)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    G.remove_nodes_from(list(nx.isolates(G)))
    G = nx.convert_node_labels_to_integers(G,first_label=1)
    nx.write_edgelist(G,fn+'mod',delimiter='\t',data=False)
    return 1

sdir=os.getcwd()+'/'

def count(name):
        n = name.split("/")[1]
        if os.path.isfile("precomputed/"+n+".countsO"):
                shutil.copy2("precomputed/"+n+".countsO",sdir+name+".countsO")
                return
        convert(sdir+name)
        f=open(name+'arg','w')
        f.write(sdir+name+'mod\n')
        f.write(sdir+name+'.countsO\n')
        f.write('orbits\n')
        f.write('five\n')
        f.close()
        os.system('Rscript Subgraphcounts5.R '+name+'arg')
        shutil.copy2(sdir+name+".countsO",sdir+"precomputed/")
        os.system('rm '+name+'arg')
        os.system('rm '+sdir+name+'mod')


if __name__ == '__main__':
        p = Pool(n)
        c=p.map(count, queries)
