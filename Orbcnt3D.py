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
        g.remove_edges_from(g.selfloop_edges())
        g.remove_nodes_from(nx.isolates(g))
        g = nx.convert_node_labels_to_integers(g,first_label=1)
        n = g.number_of_edges()
        nx.write_edgelist(g,name+'.mod1',data=False)
        l = ["1" for i in range(n)]
        l.append("")
        open(name+'.mod2',"w").write("\n".join(l))
        os.system('paste -d " " '+name+'.mod1 '+name+'.mod2 > '+name+'.mod3')
        os.system('sort '+name+'.mod3 > '+name+'.mod')
        

def count(name):
        n = name.split("/")[1]

        if os.path.isfile("3d_precomputed/"+n+".countsO3d"):
                shutil.copy2("3d_precomputed/"+n+".countsO3d",sdir+name+".countsO3d")
        else:
                prepare_graph(name)
                os.system('gtscanner/GTScanner -s 3 -d -g '+name+'.mod -m gtrie gtscanner/gtries/or_dir3.gt -or -countsO '+name+'.countsO3d')
                shutil.copy2(sdir+name+".countsO3d",sdir+"3d_precomputed/")    
        if os.path.isfile("2d_precomputed/"+n+".countsO2d"):
                shutil.copy2("2d_precomputed/"+n+".countsO2d",sdir+name+".countsO2d")
        else:
                if not os.path.isfile(name+'.mod'):
                        prepare_graph(name)
                os.system('python count_2d.py '+name+'.mod '+name+'.countsO2d')
                shutil.copy2(sdir+name+".countsO2d",sdir+"2d_precomputed/")    

        if os.path.isfile(name+'.mod'):
                os.system('rm '+name+'.mod1')
                os.system('rm '+name+'.mod2')
                os.system('rm '+name+'.mod')
                os.system('rm '+name+'.mod3')

        
for q in queries:
        count(q)