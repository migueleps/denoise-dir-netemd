import os
import sys
import networkx as nx
from multiprocessing import Pool
import shutil

def get_queries(indir):
    queries=[]
    for root, dirs, filenames in os.walk(indir):
    	for f in filenames:
    		if not ".count" in f and 'mod' not in f:
    			queries.append(os.path.join(root, f))
    return queries

def check_precomputed(name,full_path,gsize=""):
    precomp_dir = "precomputed"
    counts_suffix = ""
    if gsize != "":
        precomp_dir = "{}d_{}".format(gsize,precomp_dir)
        counts_suffix = "{}d".format(gsize)
    precomp_file = "{}/{}.countsO{}".format(precomp_dir,name,counts_suffix)
    if os.path.isfile(precomp_file):
        shutil.copy2(precomp_file,"{}/{}.countsO{}".format(full_path,name,counts_suffix))
        return True
    return False

def convert_orca(file_name):
    G=nx.read_edgelist(file_name,data=False)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    G.remove_nodes_from(list(nx.isolates(G)))
    G = nx.convert_node_labels_to_integers(G,first_label=1)
    nx.write_edgelist(G,file_name+'mod',delimiter='\t',data=False)
    return 1

def count_orca(name):
    current_dir = os.getcwd()
    n = name.split("/")[1]
    indir = name.split("/")[0]
    full_path = "{}/{}".format(current_dir,indir)

    if check_precomputed(n,indir):
        return

    convert_orca(name)
    f=open("{}arg".format(name),'w')
    f.write("{}/{}mod\n".format(full_path,n))
    f.write("{}/{}.countsO\n".format(full_path,n))
    f.write('orbits\n')
    f.write('five\n')
    f.close()
    os.system('Rscript Subgraphcounts5.R {}arg'.format(name))
    shutil.copy2("{}/{}.countsO".format(full_path,n),"precomputed/")
    os.system("rm {}arg".format(name))
    os.system("rm {}/{}mod".format(full_path,n))

def count_undirected(indir,n_threads,graphlet_size):
    queries = get_queries(indir)
    p = Pool(n_threads)
    c=p.map(count_orca, queries)


def prepare_graph(name):
    g = nx.read_edgelist(name,create_using=nx.DiGraph(),data=False,nodetype=int)
    g.remove_edges_from(list(nx.selfloop_edges(g)))
    g.remove_nodes_from(list(nx.isolates(g)))
    g = nx.convert_node_labels_to_integers(g,first_label=1)
    n = g.number_of_edges()
    nx.write_edgelist(g,name+'.mod1',data=False)
    l = ["1" for i in range(n)]
    l.append("")
    open(name+'.mod2',"w").write("\n".join(l))
    os.system('paste -d " " '+name+'.mod1 '+name+'.mod2 > '+name+'.mod3')
    os.system('sort '+name+'.mod3 > '+name+'.mod')
    return list(g.nodes())[-1]

def count_size2dir(name):
    n = name.split("/")[1]
    current_dir = os.getcwd()
    indir = name.split("/")[0]
    if check_precomputed(n,indir,2):
        return

    os.system('python count_2d.py {}.mod {}.countsO2d'.format(name,name))
    shutil.copy2("{}/{}.countsO2d".format(current_dir,name),"{}/2d_precomputed/".format(current_dir))

def count_gtries(name,gsize,n_threads):
    n = name.split("/")[1]
    current_dir = os.getcwd()
    indir = name.split("/")[0]
    if check_precomputed(n,indir,gsize):
        return

    while True:
        os.system('gtscanner/GTScanner -s {} -d -g {}.mod -m gtrie gtscanner/gtries/or_dir{}.gt -or -th {} -countsO {}.countsO{}d'.format(gsize,name,gsize,n_threads,name,gsize))
        if os.stat("{}.countsO{}d".format(name,gsize)).st_size > 0:
            break
        print("Failed to compute 4D orbits for {}, trying again".format(name))
    shutil.copy2("{}/{}.countsO{}d".format(current_dir,name,gsize),"{}/{}d_precomputed/".format(current_dir,gsize))

def count_directed(indir,n_threads,graphlet_size):
    queries = get_queries(indir)
    for name in queries:
        prepare_graph(name)
        count_size2dir(name)
        for gsize in range(3,graphlet_size+1):
            count_gtries(name,gsize,n_threads)

        if os.path.isfile("{}.mod".format(name)):
            os.system("rm {}.mod".format(name))
            os.system("rm {}.mod1".format(name))
            os.system("rm {}.mod2".format(name))
            os.system("rm {}.mod3".format(name))
