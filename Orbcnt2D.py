import os
import sys
import shutil
import networkx as nx

indir=sys.argv[1]
queries=[]
for root, dirs, filenames in os.walk(indir):
	for f in filenames:
		if not ".count" in f and 'mod' not in f:
			queries.append(os.path.join(root, f))

sdir=os.getcwd()+'/'

def prepare_graph(name):
        g = nx.read_edgelist(name,create_using=nx.DiGraph(),data=False,nodetype=int)
        g.remove_edges_from(list(nx.selfloop_edges(g)))
        g.remove_nodes_from(nx.isolates(g))
        g = nx.convert_node_labels_to_integers(g,first_label=1)
        n = g.number_of_edges()
        nx.write_edgelist(g,name+'.mod1',data=False)
        l = ["1" for i in range(n)]
        l.append("")
        open(name+'.mod2',"w").write("\n".join(l))
        os.system('paste -d " " '+name+'.mod1 '+name+'.mod2 > '+name+'.mod')


def count(name):
        n = name.split("/")[1]
        if os.path.isfile("2d_precomputed/"+n+".countsO2d"):
                shutil.copy2("2d_precomputed/"+n+".countsO2d",sdir+name+".countsO2d")
                return

        prepare_graph(name)
        os.system('python count_2d.py '+name+'.mod '+name+'.countsO2d')
        shutil.copy2(sdir+name+".countsO2d",sdir+"2d_precomputed/")
        os.system('rm '+name+'.mod1')
        os.system('rm '+name+'.mod2')
        os.system('rm '+name+'.mod')


for q in queries:
        count(q)
