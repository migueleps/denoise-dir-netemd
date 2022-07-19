/* -------------------------------------------------
      _       _     ___                            
 __ _| |_ _ _(_)___/ __| __ __ _ _ _  _ _  ___ _ _ 
/ _` |  _| '_| / -_)__ \/ _/ _` | ' \| ' \/ -_) '_|
\__, |\__|_| |_\___|___/\__\__,_|_||_|_||_\___|_|  
|___/                                          
    
gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

Pedro Ribeiro & David Aparício - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
G-Trie Implementation and associated methods

Last Update: 04/06/2014
---------------------------------------------------- */

#include "GTrie.h"
#include "GraphMatrix.h"
#include "GraphTree.h"
#include "Isomorphism.h"
#include "GraphUtils.h"
#include "Conditions.h"
#include "Random.h"
#include "Error.h"
#include "Timer.h"
#include <iostream>
#include <string.h>

#include <assert.h>

/************* Pthreads ************/
#include <vector>
#include <queue>
#include <unordered_map>

const int NO_MORE_WORK   = -1;

const int CANT_BE_CALLED = -2;
const int NOT_CALLED     = -1;
/***********************************/

bool **GTrieNode::adjM;
int **GTrieNode::fastnei;
int *GTrieNode::numnei;
int GTrieNode::numNodes;
bool GTrieNode::isdir;
double *GTrieNode::prob;
/******* Graphlets *******/
int GTrieNode::graphlet_counter = 0;
int GTrie::total_orbits;
long long int **GTrie::orbit_frequency;
//int census_size;
/*************************/

list< list<iPair> >::const_iterator jj, jjend;
list<iPair>::const_iterator kk, kkend;

/************* Pthreads ************/
struct thread_data_gtrie{
   int tid;

   int* mymap;
   bool* used;
   int glk;

   bool answered;
   bool sorry;   // sorry, thread had no work
   bool waiting;
   int  called;

   GTrieNode* root;

   long long* frequencies;                             // frequencies for each gtrie node
   long long int** orbit_frequency;
   unordered_map<GTrieNode*, vector<int> > work;  // gtrie_location -> nodes
   queue<GTrieNode*> workOrder;                   // 
   unordered_map<GTrieNode*, int> limits;         // gtrie_location -> max index limit (symmetry breaking)
};

struct thread_data_gtrie *gtrie_threads;

GraphType   graph_type;

int         gtrie_counter = 0;

/**********************************/

GTrieNode::GTrieNode(int d) {

  depth     = d;

  is_graph  = false;
  nconn     = 0;
  frequency = 0;
  /************* Pthreads ************/
  id = gtrie_counter++;
  /***********************************/

  if (d!=0) {
	/******* Graphlets *******/
	if(Global::compute_orbits){
      full_graph  = new char[(d+1) * (d+1)];
      orbits      = new long long int[d];
      graphlet_id = -1;
    }
    /*************************/
    in   = new bool[d];
    out  = new bool[d];
    conn = new int[d];
  } else {
    in = out = NULL;
    conn = NULL;
  }

  total_in = total_out = total_edges = 0;

  child.clear();
  cond.clear();
  this_node_cond.clear();

  cond_ok = false;
  cond_this_ok = false;
}

GTrieNode::~GTrieNode() {
  if (depth!=0) {
    if (in  != NULL) delete[] in;
    if (out != NULL) delete[] out;
    if (conn != NULL) delete[] conn;
    /******* Graphlets *******/
    if(Global::compute_orbits){
	  if (full_graph != NULL) delete[] full_graph;
      if (orbits != NULL) delete[] orbits;
    }
    /*************************/
  }
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    delete *ii;
}

void GTrieNode::showAsText(FILE *f) {
  int i;
  for (i=0; i<depth; i++) fprintf(f,"   ");
  fprintf(f,"[");
  for (i=0; i<depth; i++) 
    fprintf(f,"%s", out[i]?"1":"0");
  fprintf(f,"][");
  for (i=0; i<depth; i++) 
    fprintf(f,"%s", in[i]?"1":"0");
  fprintf(f,"] ");

  if (cond_ok) fprintf(f,"{}");
  else {    
    list< list<iPair> >:: iterator jj;
    list<iPair>:: iterator kk;
    for (jj=cond.begin(); jj!=cond.end(); jj++) {
      fputc('|', f);
      for (kk=(*jj).begin(); kk!=(*jj).end(); kk++) {
	if (kk!=(*jj).begin()) fputc(' ', f);
	fprintf(f,"%d<%d", kk->first, kk->second);
      }
      fputc('|', f);
    }
  }

  fputc('+', f);

  if (cond_this_ok) fprintf(f,"{}");
  else {
    list< list<int> >:: iterator jjj;
    list<int>:: iterator kkk;
    for (jjj=this_node_cond.begin(); jjj!=this_node_cond.end(); jjj++) {
      fputc('|', f);
      if ((*jjj).size()==0) fprintf(f,"{}");
      else {
	for (kkk=(*jjj).begin(); kkk!=(*jjj).end(); kkk++) {
	  if (kkk!=(*jjj).begin()) fputc(',', f);
	  fprintf(f,"%d", *kkk);	
	}
	fprintf(f,"<%d", depth-1);
      }
      fputc('|', f);
    }
  }
      
  fputc(' ', f);
        
  fprintf(f,"%s\n", is_graph?"isGraph":"");

  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->showAsText(f);
}

float GTrieNode::frequencyGraph(Graph *g) {
  if (g->numNodes()==depth)
    return frequency;
  list<GTrieNode *>::iterator ii;
  int i;
  for(ii=child.begin(); ii!=child.end(); ii++) {
    for (i=0; i<=depth; i++)
      if (g->hasEdge(depth, i) != (*ii)->out[i] ||
	  g->hasEdge(i, depth) != (*ii)->in[i])
	break;
    if (i>depth) break;
  }
  if (ii!=child.end())
    return (*ii)->frequencyGraph(g);
  else
    return -1;
}

void GTrieNode::insert(Graph *g) {
  int i;


  if (g->numNodes() == depth) {
    is_graph = true;
  } else {
    list<GTrieNode *>::iterator ii;
    for(ii=child.begin(); ii!=child.end(); ii++) {
      for (i=0; i<=depth; i++)
	if (g->hasEdge(depth, i) != (*ii)->out[i] ||
	    g->hasEdge(i, depth) != (*ii)->in[i])
	  break;
      if (i>depth) break;
    }
    if (ii==child.end()) {
      GTrieNode *gt = new GTrieNode(depth+1);
      for (i=0; i<=depth; i++) {
      	gt->out[i] = g->hasEdge(depth, i);	
      	gt->in[i]  = g->hasEdge(i, depth);
	if (g->isConnected(depth, i)) 
	  gt->conn[(gt->nconn)++] = i;

	if (gt->in[i])  { (gt->total_in)++;  (gt->total_edges)++;}
	if (gt->out[i]) { (gt->total_out)++; (gt->total_edges)++;}
      }
      gt->insert(g);
      child.push_back(gt);
    } else
      (*ii)->insert(g);
  }
}


void GTrieNode::insertCond(Graph *g, list<iPair> *cond) {
  int i;
  GTrieNode *gt;

  if (g->numNodes() == depth) {
    is_graph = true;
  } else {
    list<GTrieNode *>::iterator ii;
    for(ii=child.begin(); ii!=child.end(); ii++) {
      for (i=0; i<=depth; i++)
	if (g->hasEdge(depth, i) != (*ii)->out[i] ||
	    g->hasEdge(i, depth) != (*ii)->in[i])
	  break;
      if (i>depth) break;
    }
    if (ii==child.end()) {
      gt = new GTrieNode(depth+1);
      for (i=0; i<=depth; i++) {
      	gt->out[i] = g->hasEdge(depth, i);	
      	gt->in[i]  = g->hasEdge(i, depth);
	if (g->isConnected(depth, i)) 
	  gt->conn[(gt->nconn)++] = i;
	if (gt->in[i])  { (gt->total_in)++;  (gt->total_edges)++;}
	if (gt->out[i]) { (gt->total_out)++; (gt->total_edges)++;}
      }
      
    /******* Graphlets *******/
    if(Global::compute_orbits)
      if(depth >= 1)
        for(i = 0; i <= depth; i++)
          for(int j = 0; j <= depth; j++)
            gt->full_graph[i * (depth+1) + j] = (char)(((int)'0') + g->hasEdge(i, j));
    /*************************/
      
      child.push_back(gt);
    } else {
      gt = (*ii);
    }

    gt->insertConditionsFiltered(cond);
    gt->insertCond(g, cond);

  }
}

// Assume a and b are in ascending order
bool GTrieNode::_is_intset_included(list<int> a, list<int> b) {
  list<int>::iterator aa, bb;

  aa = a.begin();
  bb = b.begin();

  while (aa!=a.end() && bb!=b.end()) {
    if (*aa == *bb) {
      aa++;
      bb++;
    } else {
      bb++;
    }      
  }

  if (aa == a.end()) return true;
  else return false;
}


// Assume a and b are in ascending order
bool GTrieNode::_is_pairset_included(list<iPair> a, list<iPair> b) {
  list<iPair>::iterator aa, bb;

  aa = a.begin();
  bb = b.begin();

  while (aa!=a.end() && bb!=b.end()) {
    if ((aa->first == bb->first) && (aa->second == bb->second)) {
      aa++;
      bb++;
    } else {
      bb++;
    }      
  }

  if (aa == a.end()) return true;
  else return false;
}

// Insert Symmetry conditions pertaining to this node
void GTrieNode::insertConditionsFiltered(list<iPair> *conditions) {

  // Already has "empty set" of conditions
  if (cond_ok && cond_this_ok) return;

  list<iPair> aux;
  list<int>   aux_this_node;
  list<iPair>::iterator jj, kk;
  iPair p;


  /*for (jj=conditions->begin(); jj!=conditions->end(); jj++)
    printf("%d < %d  |  ", jj->first, jj->second);
    putchar('\n');*/
  

  // A bit slow but works (search a<b and b<c: remove a<c)
  for (jj=conditions->begin(); jj!=conditions->end(); jj++)
    for (kk=conditions->begin(); kk!=conditions->end(); kk++)
      if (jj->second == kk->first) {
	p.first  = jj->first;
	p.second = kk->second;
	conditions->remove(p); // does not invalidade jj and kk
      }

  /*for (jj=conditions->begin(); jj!=conditions->end(); jj++)
    printf("%d < %d  |  ", jj->first, jj->second);
    putchar('\n');*/

  for (jj=conditions->begin(); jj!=conditions->end(); jj++) {
    if (max(jj->first, jj->second)<=depth-1) // HERE!
      aux.push_back(*jj);
    if (jj->second == depth-1)
      aux_this_node.push_back(jj->first);
 } 

  // Deal with ancestor conditions
  if (!cond_ok) {

    if (aux.size()==0) {
      cond_ok = true;
      cond.clear();
    } else {

      list< list<iPair> >::iterator mm;
      bool is_contained = false;
      // A bit slow, but saves a lot of conditions
      for (mm=cond.begin(); mm!=cond.end();) {
	if (_is_pairset_included(*mm, aux)) {
	  is_contained = true;
	  break;
	}
	else if (_is_pairset_included(aux, *mm)) {
	  //printf("Erasing old condition (ancestors) !!\n");
	  mm = cond.erase(mm);
	  if (mm == cond.end()) break;
	}
	else mm++;
      }

      if (!is_contained) cond.push_back(aux);
      //else printf("avoid insertion (ancestors)!\n");
    }
  }

  return;

 // Deal with this node conditions {
  if (!cond_this_ok) {
    
    if (aux_this_node.size()==0) {
      cond_this_ok = true;
      this_node_cond.clear();
      } else {
    
      list< list<int> >::iterator mm;

      bool is_contained = false;
      // A bit slow, but saves a lot of conditions
      for (mm=this_node_cond.begin(); mm!=this_node_cond.end(); mm++) {
        if (_is_intset_included(*mm, aux_this_node)) {
            is_contained = true;
            break;
        }
        if (_is_intset_included(aux_this_node, *mm)) {
          //printf("Erasing old condition (this_node) !!\n");
          mm = this_node_cond.erase(mm);
       }
      }

      if (!is_contained) this_node_cond.push_back(aux_this_node);
      //else printf("avoid insertion (this_node)!\n");      
    }
  }

   
}

int GTrieNode::countNodes() {
  int aux=1;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countNodes();

  return aux;
}

void GTrieNode::zeroFrequency() {
  frequency = 0;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->zeroFrequency();
}

void GTrieNode::showFrequency() {
  if (is_graph) printf("%lld \n", frequency);
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->showFrequency();
}

int GTrieNode::maxDepth() {
  int aux = 0;
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux=max(aux, 1+(*ii)->maxDepth());
  return aux;
}

void GTrie::census(Graph *g) {
  int i;
  int subgraph_size = maxDepth();
  SharedParallel::num_nodes = g->numNodes();

  _root->zeroFrequency();

  /************* Pthreads ************/
  GTrieNode::numNodes = SharedParallel::num_nodes;
  GTrieNode::fastnei  = g->matrixNeighbours();
  GTrieNode::adjM     = g->adjacencyMatrix();
  GTrieNode::numnei   = g->arrayNumNeighbours(); 

  /******* Graphlets *******/
  if(Global::compute_orbits){
    orbit_frequency = (long long int **)malloc(g->numNodes() * sizeof(long long int *));
    for(int i = 0; i < g->numNodes(); i++) orbit_frequency[i] = (long long int*)malloc(total_orbits * sizeof(long long int));
    for(int i = 0; i < g->numNodes(); i++) for(int j = 0; j < total_orbits; j++) orbit_frequency[i][j] = 0;
  }
  /*************************/

  if (g->type() == DIRECTED) { GTrieNode::isdir = true,  graph_type = DIRECTED;   }
  else                       { GTrieNode::isdir = false, graph_type = UNDIRECTED; }

  Timer::start(1); // parallel census time
  
  gtrie_threads = new thread_data_gtrie[Global::num_threads];
  pthread_t thread[Global::num_threads];

  SharedParallel::waiting_threads = 0;
  SharedParallel::no_more_work    = false;

  for (i = 0; i < Global::num_threads; i++){
    gtrie_threads[i].tid          = i;
    gtrie_threads[i].mymap        = new int[subgraph_size];
    gtrie_threads[i].used         = new bool[SharedParallel::num_nodes];
    gtrie_threads[i].sorry        = false;
    gtrie_threads[i].waiting      = false;
    gtrie_threads[i].called       = CANT_BE_CALLED;
    gtrie_threads[i].root         = _root;
    gtrie_threads[i].frequencies  = new long long[gtrie_counter];
    if(Global::compute_orbits){
      gtrie_threads[i].orbit_frequency = (long long int **)malloc(g->numNodes() * sizeof(long long int *));
      for(int j = 0; j < g->numNodes(); j++) gtrie_threads[i].orbit_frequency[j] = (long long int*)malloc(total_orbits * sizeof(long long int));
      for(int j = 0; j < g->numNodes(); j++) for(int k = 0; k < total_orbits; k++) gtrie_threads[i].orbit_frequency[j][k] = 0;
    }
  }
  
  for (i = 0; i < Global::num_threads; i++){
    pthread_create(&thread[i], NULL, GTrie::startCensus, (void *) &gtrie_threads[i]);
  }

  for (i = 0; i < Global::num_threads; i++)
    pthread_join(thread[i], NULL);

  Timer::stop(1);
  printf("Census time: %.2f\n", Timer::elapsed(1));
    
  GTrieNode *c = *(_root->child.begin());
  list<GTrieNode *>::iterator  ii;
  for(ii=c->child.begin(); ii!=c->child.end(); ii++)
    (*ii)->calcFrequency();
   
  /************* Graphlets ***********/    
  if(Global::compute_orbits){
    for(int j = 0; j < g->numNodes(); j++){
      for(int o = 0; o < total_orbits; o++){
	for (int k = 0; k < Global::num_threads; k++) {
	    if(gtrie_threads[k].orbit_frequency[j][o] != 0)
	      GTrie::orbit_frequency[j][o] += gtrie_threads[k].orbit_frequency[j][o];
	}
      }
    }
   }
  /***********************************/

  for (i = 0; i < Global::num_threads; i++){
    delete [] gtrie_threads[i].mymap;
    delete [] gtrie_threads[i].used;
    delete [] gtrie_threads[i].frequencies;
    if(Global::compute_orbits){
	  for(int j = 0; j < g->numNodes(); j++) free(gtrie_threads[i].orbit_frequency[j]);
      free(gtrie_threads[i].orbit_frequency);
    }
  }

  delete [] gtrie_threads;
  /**********************************/
}

/************* Pthreads ************/
void* GTrie::startCensus(void *data){
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  
  my_data->called = NOT_CALLED;

  Timer::start(my_data->tid+1);

  if(my_data->tid >= SharedParallel::num_nodes){
    sleep(2);
    askForWork(data);
  }
  else{
	if(!Global::compute_orbits)
      censusOnNodes(my_data->tid, SharedParallel::num_nodes-1, data);
    else
      censusOnNodesGraphlets(my_data->tid, SharedParallel::num_nodes-1, data);
  }
  
  return NULL;
}
/**********************************/

/************* Pthreads ************/
void GTrie::censusOnNodes(int begin, int end, void *data) {
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;

  memset(my_data->used, 0, sizeof(*my_data->used) * SharedParallel::num_nodes);
  memset(my_data->frequencies, 0, sizeof(*my_data->frequencies) * gtrie_counter);
  
  GTrieNode *c = *(my_data->root->child.begin());
  list<GTrieNode *>::iterator ii;
  my_data->glk = 1;

  for (int node=begin; node<=end; node+=Global::num_threads){
    if(my_data->called != NOT_CALLED){
      my_data->glk = 0;
      createWorkFromInterval(data, node, end, c);
      break;
    }
    my_data->mymap[0] = node;
    my_data->used[node] = true;
    if (graph_type == DIRECTED)
      for(ii=c->child.begin(); ii!=c->child.end(); ii++)
        (*ii)->goCondDir(data);
    else
      for(ii=c->child.begin(); ii!=c->child.end(); ii++)
        (*ii)->goCondUndir(data);
    my_data->used[node] = false;
  }

  if(my_data->called != -1){
    giveWork(data);
    goBackToWork(data);
  }
  askForWork(data);
}


void GTrie::censusOnNodesGraphlets(int begin, int end, void *data) {
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;

  memset(my_data->used, 0, sizeof(*my_data->used) * SharedParallel::num_nodes);
  memset(my_data->frequencies, 0, sizeof(*my_data->frequencies) * gtrie_counter);
  
  GTrieNode *c = *(my_data->root->child.begin());
  list<GTrieNode *>::iterator ii;
  my_data->glk = 1;

  for (int node=begin; node<=end; node+=Global::num_threads){
    if(my_data->called != NOT_CALLED){
      my_data->glk = 0;
      createWorkFromInterval(data, node, end, c);
      break;
    }
    my_data->mymap[0] = node;
    my_data->used[node] = true;
    if (graph_type == DIRECTED)
      for(ii=c->child.begin(); ii!=c->child.end(); ii++)
        (*ii)->goCondDirGraphlets(data);
    else
      for(ii=c->child.begin(); ii!=c->child.end(); ii++)
        (*ii)->goCondUndirGraphlets(data);
    my_data->used[node] = false;
  }

  if(my_data->called != -1){
    giveWork(data);
    goBackToWorkGraphlets(data);
  }
  askForWork(data);
}
/***********************************/
/***********************************/

/************* Pthreads ************/
void GTrie::askForWork(void* data){
  int to_call;
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data; 
  
  // thread checks if it was itself called; checks if there is no more work
  pthread_mutex_lock(&SharedParallel::waiting_mutex);

  if(my_data->called != NOT_CALLED) { 
    pthread_mutex_lock(&SharedParallel::answer_mutex);
    gtrie_threads[my_data->called].sorry  = true;
    gtrie_threads[my_data->called].answered  = true;
    
    my_data->called = -1;
    pthread_mutex_unlock(&SharedParallel::answer_mutex);
    pthread_cond_signal(&SharedParallel::answer_cond);
  }

  SharedParallel::waiting_threads++;
  my_data->waiting = true;
  
  if(SharedParallel::waiting_threads == Global::num_threads)
    SharedParallel::no_more_work = true;

  pthread_mutex_unlock(&SharedParallel::waiting_mutex);
  //

  // calls a thread for work
  pthread_mutex_lock(&SharedParallel::call_mutex);
  
  if(SharedParallel::no_more_work){
    pthread_mutex_unlock(&SharedParallel::call_mutex);
    pthread_exit(NULL);
  }

  my_data->answered = false;

  pthread_mutex_lock(&SharedParallel::waiting_mutex);
  to_call = getThreadToCall(my_data->tid);

  if(to_call == NO_MORE_WORK){
    pthread_mutex_unlock(&SharedParallel::call_mutex);
    pthread_mutex_unlock(&SharedParallel::waiting_mutex);
    pthread_exit(NULL);
  }

  gtrie_threads[to_call].called = my_data->tid;
  pthread_mutex_unlock(&SharedParallel::waiting_mutex);
  //
  
  // waits for answer
  pthread_mutex_lock(&SharedParallel::answer_mutex);
  
  while(my_data->answered == false)
    pthread_cond_wait(&SharedParallel::answer_cond, &SharedParallel::answer_mutex);
  
  pthread_mutex_unlock(&SharedParallel::answer_mutex);
  pthread_mutex_unlock(&SharedParallel::call_mutex);
  //
  
  // handles work received
  pthread_mutex_lock(&SharedParallel::waiting_mutex);
  SharedParallel::waiting_threads--;
  my_data->waiting = false;
  pthread_mutex_unlock(&SharedParallel::waiting_mutex);
  
  if (my_data->sorry) { 
    my_data->sorry = false; 
    askForWork(data);
  } else{
	/********* Graphlets *********/
    if(Global::compute_orbits) goBackToWorkGraphlets(data);
    else                       goBackToWork(data);
    /******************************/
  }
}
/***********************************/

/************* Pthreads ************/
int GTrie::getThreadToCall(int tid){
  if(SharedParallel::no_more_work)
    return NO_MORE_WORK;
  else{
    int to_call = Random::getInteger(0, Global::num_threads-1);
    
    if(to_call != tid && gtrie_threads[to_call].waiting == false && 
        gtrie_threads[to_call].called == NOT_CALLED)
      return to_call;
    
    return getThreadToCall(tid);
  }
}
/**********************************/

void GTrieNode::goCondUndir(void* data) {  
  int i, j, ci, mylim, glaux;
  int ncand;
  int *p;

  /************* Pthreads ************/
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  /***********************************/

  mylim = INT_MAX;
  if (!cond_ok) {
    i = 1;
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      glaux = -1;
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
        if (kk->second<my_data->glk && my_data->mymap[kk->first] > my_data->mymap[kk->second])
          break;
        else if (kk->second==my_data->glk && my_data->mymap[kk->first]>glaux)
          glaux = my_data->mymap[kk->first];
      if (kk==kkend) {
        i = 0;
        if (glaux < mylim) mylim=glaux;
      }
    }
    if (i) return;
  }
  if (mylim == INT_MAX) mylim = 0;
    
  ncand=0;
  j=ci=INT_MAX;
  for (i=0; i<nconn; i++) {
    glaux = numnei[my_data->mymap[conn[i]]];
    if (glaux<j) {
      ci = my_data->mymap[conn[i]];
      j = glaux;
    }
  }

  glaux = j;
  ncand = ci;
  for (p=&fastnei[ncand][j-1], ci= glaux-1; ci>=0; ci--, p--) { 
    i = *p;
    if (i<mylim) break;
    if (my_data->used[i]) continue;
    /************* Pthreads ************/
    if(my_data->called != NOT_CALLED){
      GTrie::createWorkFromArray(data, ci, &fastnei[ncand][0], this, mylim);
      return;
    }
    /***********************************/
    my_data->mymap[my_data->glk] = i;
    
    bool *b = &adjM[i][0];
    for (j=0; j<my_data->glk; j++)
      if (out[j] != *(b+my_data->mymap[j]))
	break;
    if (j<my_data->glk) continue;
    
    if (is_graph) {
      ++my_data->frequencies[id];
      if (Global::show_occ) {
        for (int k = 0; k<=my_data->glk; k++)
          for (int l = 0; l<=my_data->glk; l++)
            fputc(adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file);
        fputc(':', Global::occ_file);
        for (int k = 0; k<=my_data->glk; k++)
          fprintf(Global::occ_file, " %d", my_data->mymap[k]+1);
        fputc('\n', Global::occ_file);
      }
    }

    my_data->used[i]=true;
    my_data->glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      (*ii)->goCondUndir(data);
    my_data->glk--;
    my_data->used[i]=false;
  }
}

void GTrieNode::goCondDir(void* data) {  
  int i, j, ci, mylim, glaux;
  int ncand;
  int *p;
  /************* Pthreads ************/
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  /***********************************/

  mylim = INT_MAX;
  if (!cond_ok) {
    i = 1;
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      glaux = -1;
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
        if (kk->second<my_data->glk && my_data->mymap[kk->first] > my_data->mymap[kk->second])
          break;
        else if (kk->second==my_data->glk && my_data->mymap[kk->first]>glaux)
          glaux = my_data->mymap[kk->first];
      if (kk==kkend) {
        i = 0;
        if (glaux < mylim) mylim=glaux;
      }
    }
    if (i) return;
  }
  if (mylim == INT_MAX) mylim = 0;
    
  ncand=0;
  j=ci=INT_MAX;
  for (i=0; i<nconn; i++) {
    glaux = numnei[my_data->mymap[conn[i]]];
    if (glaux<j) {
      ci = my_data->mymap[conn[i]];
      j = glaux;
    }
  }

  glaux = j;
  ncand = ci;
  for (p=&fastnei[ncand][j-1], ci= glaux-1; ci>=0; ci--, p--) { 
    i = *p;
    if (i<mylim) break;
    if (my_data->used[i]) continue;
    /************* Pthreads ************/
    if(my_data->called != NOT_CALLED){
      GTrie::createWorkFromArray(data, ci, &fastnei[ncand][0], this, mylim);
      return;
    }
    /***********************************/
    
    my_data->mymap[my_data->glk] = i;

    for (j=0; j<my_data->glk; j++)
      if (in[j] != adjM[my_data->mymap[j]][i])
      break;
    if (j<my_data->glk) continue;
    bool *b = &adjM[i][0];
    for (j=0; j<my_data->glk; j++)
      if (out[j] != *(b+my_data->mymap[j]))
        break;
    if (j<my_data->glk) continue;
    
    if (is_graph) {
      ++my_data->frequencies[id];
      if (Global::show_occ) {
        for (int k = 0; k<=my_data->glk; k++)
          for (int l = 0; l<=my_data->glk; l++)
            fputc(adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file_t[my_data->tid]);
        fputc(':', Global::occ_file_t[my_data->tid]);
        for (int k = 0; k<=my_data->glk; k++)
          fprintf(Global::occ_file_t[my_data->tid], " %d", my_data->mymap[k]);
        fputc('\n', Global::occ_file_t[my_data->tid]);
      }
    }

    my_data->used[i]=true;
    my_data->glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      (*ii)->goCondDir(data);
    my_data->glk--;
    my_data->used[i]=false;
  }
}

void GTrieNode::goCondUndirGraphlets(void* data) {  
  int i, j, ci, mylim, glaux;
  int ncand;
  int *p;

  /************* Pthreads ************/
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  /***********************************/

  mylim = INT_MAX;
  if (!cond_ok) {
    i = 1;
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      glaux = -1;
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
        if (kk->second<my_data->glk && my_data->mymap[kk->first] > my_data->mymap[kk->second])
          break;
        else if (kk->second==my_data->glk && my_data->mymap[kk->first]>glaux)
          glaux = my_data->mymap[kk->first];
      if (kk==kkend) {
        i = 0;
        if (glaux < mylim) mylim=glaux;
      }
    }
    if (i) return;
  }
  if (mylim == INT_MAX) mylim = 0;
    
  ncand=0;
  j=ci=INT_MAX;
  for (i=0; i<nconn; i++) {
    glaux = numnei[my_data->mymap[conn[i]]];
    if (glaux<j) {
      ci = my_data->mymap[conn[i]];
      j = glaux;
    }
  }

  glaux = j;
  ncand = ci;
  for (p=&fastnei[ncand][j-1], ci= glaux-1; ci>=0; ci--, p--) { 
    i = *p;
    if (i<mylim) break;
    if (my_data->used[i]) continue;
    /************* Pthreads ************/
    if(my_data->called != NOT_CALLED){
      GTrie::createWorkFromArray(data, ci, &fastnei[ncand][0], this, mylim);
      return;
    }
    /***********************************/
    my_data->mymap[my_data->glk] = i;
    
    bool *b = &adjM[i][0];
    for (j=0; j<my_data->glk; j++)
      if (out[j] != *(b+my_data->mymap[j]))
	break;
    if (j<my_data->glk) continue;
    
    if (is_graph) {
      ++my_data->frequencies[id];
      /******* Graphlets *******/
      if(Global::compute_orbits){
        for(int o = 0; o < depth; o++) 
          (my_data->orbit_frequency[my_data->mymap[o]][orbits[o]])++;
	  }
      /*************************/
      if (Global::show_occ) {
        for (int k = 0; k<=my_data->glk; k++)
          for (int l = 0; l<=my_data->glk; l++)
            fputc(adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file);
        fputc(':', Global::occ_file);
        for (int k = 0; k<=my_data->glk; k++)
          fprintf(Global::occ_file, " %d", my_data->mymap[k]+1);
        fputc('\n', Global::occ_file);
      }
    }

    my_data->used[i]=true;
    my_data->glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      (*ii)->goCondUndirGraphlets(data);
    my_data->glk--;
    my_data->used[i]=false;
  }
}

void GTrieNode::goCondDirGraphlets(void* data) {  
  int i, j, ci, mylim, glaux;
  int ncand;
  int *p;
  /************* Pthreads ************/
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  /***********************************/

  mylim = INT_MAX;
  if (!cond_ok) {
    i = 1;
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      glaux = -1;
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
        if (kk->second<my_data->glk && my_data->mymap[kk->first] > my_data->mymap[kk->second])
          break;
        else if (kk->second==my_data->glk && my_data->mymap[kk->first]>glaux)
          glaux = my_data->mymap[kk->first];
      if (kk==kkend) {
        i = 0;
        if (glaux < mylim) mylim=glaux;
      }
    }
    if (i) return;
  }
  if (mylim == INT_MAX) mylim = 0;
    
  ncand=0;
  j=ci=INT_MAX;
  for (i=0; i<nconn; i++) {
    glaux = numnei[my_data->mymap[conn[i]]];
    if (glaux<j) {
      ci = my_data->mymap[conn[i]];
      j = glaux;
    }
  }

  glaux = j;
  ncand = ci;
  for (p=&fastnei[ncand][j-1], ci= glaux-1; ci>=0; ci--, p--) { 
    i = *p;
    if (i<mylim) break;
    if (my_data->used[i]) continue;
    /************* Pthreads ************/
    if(my_data->called != NOT_CALLED){
      GTrie::createWorkFromArray(data, ci, &fastnei[ncand][0], this, mylim);
      return;
    }
    /***********************************/
    
    my_data->mymap[my_data->glk] = i;

    for (j=0; j<my_data->glk; j++)
      if (in[j] != adjM[my_data->mymap[j]][i])
      break;
    if (j<my_data->glk) continue;
    bool *b = &adjM[i][0];
    for (j=0; j<my_data->glk; j++)
      if (out[j] != *(b+my_data->mymap[j]))
        break;
    if (j<my_data->glk) continue;
    
    if (is_graph) {
      ++my_data->frequencies[id];
      /******* Graphlets *******/
      if(Global::compute_orbits)
        for(int o = 0; o < depth; o++) (my_data->orbit_frequency[my_data->mymap[o]][orbits[o]])++;
      /*************************/
      if (Global::show_occ) {
        for (int k = 0; k<=my_data->glk; k++)
          for (int l = 0; l<=my_data->glk; l++)
            fputc(adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file_t[my_data->tid]);
        fputc(':', Global::occ_file_t[my_data->tid]);
        for (int k = 0; k<=my_data->glk; k++)
          fprintf(Global::occ_file_t[my_data->tid], " %d", my_data->mymap[k]);
        fputc('\n', Global::occ_file_t[my_data->tid]);
      }
    }

    my_data->used[i]=true;
    my_data->glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      (*ii)->goCondDirGraphlets(data);
    my_data->glk--;
    my_data->used[i]=false;
  }
}


/************* Pthreads ************/
void GTrie::createWorkFromInterval(void* data, int current, int last, GTrieNode* gtrie_location){
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;

  vector<int> nodes;

  for (int k = current; k <= last; k+= Global::num_threads)
    nodes.push_back(k);

  my_data->workOrder.push(gtrie_location);
  my_data->work[gtrie_location] = nodes;
}
/**********************************/

/************* Pthreads ************/
void GTrie::createWorkFromArray(void* data, int index, int* nodes_to_do, GTrieNode* gtrie_location, int mylim){
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  int node;
  
  vector<int> nodes;

  for (; index>=0; index--, nodes_to_do++) { 
    node = *nodes_to_do;
    nodes.push_back(node);
  }

  my_data->workOrder.push(gtrie_location);
  my_data->work[gtrie_location] = nodes;
  my_data->limits[gtrie_location] = mylim;
}
/**********************************/

/************* Pthreads ************/

void GTrie::giveWork(void* data){
  int node;
  bool my_turn = true;
  bool more_work = false;
  
  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;
  
  pthread_mutex_lock(&SharedParallel::answer_mutex);

  queue<GTrieNode*> myWorkOrder;
  queue<GTrieNode*> yourWorkOrder;

  unordered_map<GTrieNode*, vector<int> > mywork;
  unordered_map<GTrieNode*, vector<int> > yourwork;

  my_data->workOrder = correctOrder(my_data->workOrder);
  int my_depth = my_data->workOrder.front()->depth - 1;

  // diagonal task splitting
  while(!my_data->workOrder.empty()){
    GTrieNode* gtrie_node = my_data->workOrder.front();
    
    vector<int> my_nodes;
    vector<int> your_nodes;

    if(my_data->work[gtrie_node].size() > 0) more_work = true;
    const int size = my_data->work[gtrie_node].size();
    
    for(int i = 0; i < size; i++){
      node = my_data->work[gtrie_node][i];
      
      if(my_turn) my_nodes.push_back(node);
      else        your_nodes.push_back(node);
 
      my_turn = !my_turn;
    }
    my_data->workOrder.pop();

    if(!my_nodes.empty()){
      myWorkOrder.push(gtrie_node);
      mywork[gtrie_node] = my_nodes;
    }

    if(!your_nodes.empty()){
      yourWorkOrder.push(gtrie_node);
      yourwork[gtrie_node] = your_nodes;
    }
  }
  
  memcpy(gtrie_threads[my_data->called].mymap, my_data->mymap, sizeof(int) * my_depth);

  my_data->work.clear();
  my_data->work = mywork;
  my_data->workOrder = myWorkOrder;
  gtrie_threads[my_data->called].work = yourwork;
  gtrie_threads[my_data->called].workOrder = yourWorkOrder;
  gtrie_threads[my_data->called].limits.clear();
  gtrie_threads[my_data->called].limits = my_data->limits;
  
  gtrie_threads[my_data->called].answered = true;
  
  if(!more_work)  gtrie_threads[my_data->called].sorry  = true;
  
  my_data->called = -1;

  pthread_mutex_unlock(&SharedParallel::answer_mutex);
  pthread_cond_signal(&SharedParallel::answer_cond);
  
  if(!more_work) askForWork(data);
}
/***********************************/

/************* Pthreads ************/
queue<GTrieNode*> GTrie::correctOrder(queue<GTrieNode*> old_order){
  queue<GTrieNode*> new_order;
  int max_depth = 0;

  while(!old_order.empty()){
    GTrieNode* this_node = old_order.front();
    if(this_node->depth > max_depth) { max_depth = this_node->depth; }
    old_order.pop();
    new_order.push(this_node);
  }

  while(true){
    if(new_order.front()->depth < max_depth){
      GTrieNode* this_node = new_order.front();
      new_order.pop();
      new_order.push(this_node);
    } else
      break;
  }

  return new_order;
}
/***********************************/

/************* Pthreads ************/
void GTrie::goBackToWork(void* data){
  int i, j, node;
  bool did_work = false;

  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;

  list<GTrieNode *>::iterator ii;
  
  while((ssize_t)my_data->workOrder.size() > 0){
    GTrieNode* gtrie_node = my_data->workOrder.front();
    my_data->glk = gtrie_node->depth - 1;

    memset(my_data->used, 0, sizeof(*my_data->used) * SharedParallel::num_nodes);
    for(i=0; i < gtrie_node->depth - 1; i++)
      my_data->used[my_data->mymap[i]] = 1;

    while(!my_data->work[gtrie_node].empty()){
      node = my_data->work[gtrie_node].back();
      
      if (node < my_data->limits[gtrie_node]) { my_data->work[gtrie_node].clear(); break; }

      if (!my_data->used[node]) {
        if(my_data->called != NOT_CALLED && did_work == true){
          break;
        }
        
        my_data->mymap[my_data->glk] = node;
        if (graph_type == DIRECTED){
		      for (j=0; j<my_data->glk; j++)
            if (gtrie_node->in[j] != gtrie_node->adjM[my_data->mymap[j]][node])
              break;
          if (j<my_data->glk) { my_data->work[gtrie_node].pop_back(); continue; }
	      }
        bool *b = &(gtrie_node->adjM[node][0]);
        for (j=0; j<my_data->glk; j++)
          if (gtrie_node->out[j] != *(b+my_data->mymap[j]))
            break;
        if (j<my_data->glk) { my_data->work[gtrie_node].pop_back(); continue; }
        
        if (gtrie_node->is_graph) {
          ++my_data->frequencies[gtrie_node->id];
          if (Global::show_occ) {
            for (int k = 0; k<=my_data->glk; k++)
              for (int l = 0; l<=my_data->glk; l++)
                fputc(gtrie_node->adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file_t[my_data->tid]);
            fputc(':', Global::occ_file_t[my_data->tid]);
            for (int k = 0; k<=my_data->glk; k++)
              fprintf(Global::occ_file_t[my_data->tid], " %d", my_data->mymap[k]);
            fputc('\n', Global::occ_file_t[my_data->tid]);
          }
        }

        my_data->used[node] = true;
        my_data->glk++;
        did_work = true;
        if (graph_type == DIRECTED)
          for(ii=gtrie_node->child.begin(); ii!=gtrie_node->child.end(); ii++)
            (*ii)->goCondDir(data);
        else
          for(ii=gtrie_node->child.begin(); ii!=gtrie_node->child.end(); ii++)
            (*ii)->goCondUndir(data);
         my_data->glk--; 
         my_data->used[node] = false;
      }
      my_data->work[gtrie_node].pop_back();
    }
    if(my_data->called != -1){
      giveWork(data);
      goBackToWork(data);
    }
    my_data->workOrder.pop();
  }
  askForWork(data);
}


void GTrie::goBackToWorkGraphlets(void* data){
  int i, j, node;
  bool did_work = false;

  struct thread_data_gtrie *my_data;
  my_data = (struct thread_data_gtrie *) data;

  list<GTrieNode *>::iterator ii;
  
  while((ssize_t)my_data->workOrder.size() > 0){
    GTrieNode* gtrie_node = my_data->workOrder.front();
    my_data->glk = gtrie_node->depth - 1;

    memset(my_data->used, 0, sizeof(*my_data->used) * SharedParallel::num_nodes);
    for(i=0; i < gtrie_node->depth - 1; i++)
      my_data->used[my_data->mymap[i]] = 1;

    while(!my_data->work[gtrie_node].empty()){
      node = my_data->work[gtrie_node].back();
      
      if (node < my_data->limits[gtrie_node]) { my_data->work[gtrie_node].clear(); break; }

      if (!my_data->used[node]) {
        if(my_data->called != NOT_CALLED && did_work == true){
          break;
        }
        
        my_data->mymap[my_data->glk] = node;
        if (graph_type == DIRECTED){
		      for (j=0; j<my_data->glk; j++)
            if (gtrie_node->in[j] != gtrie_node->adjM[my_data->mymap[j]][node])
              break;
          if (j<my_data->glk) { my_data->work[gtrie_node].pop_back(); continue; }
	      }
        bool *b = &(gtrie_node->adjM[node][0]);
        for (j=0; j<my_data->glk; j++)
          if (gtrie_node->out[j] != *(b+my_data->mymap[j]))
            break;
        if (j<my_data->glk) { my_data->work[gtrie_node].pop_back(); continue; }
        
        if (gtrie_node->is_graph) {
          ++my_data->frequencies[gtrie_node->id];
          /******* Graphlets *******/
          for(int o = 0; o < gtrie_node->depth; o++) (my_data->orbit_frequency[my_data->mymap[o]][gtrie_node->orbits[o]])++;
          /*************************/
          if (Global::show_occ) {
            for (int k = 0; k<=my_data->glk; k++)
              for (int l = 0; l<=my_data->glk; l++)
                fputc(gtrie_node->adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file_t[my_data->tid]);
            fputc(':', Global::occ_file_t[my_data->tid]);
            for (int k = 0; k<=my_data->glk; k++)
              fprintf(Global::occ_file_t[my_data->tid], " %d", my_data->mymap[k]);
            fputc('\n', Global::occ_file_t[my_data->tid]);
          }
        }

        my_data->used[node] = true;
        my_data->glk++;
        did_work = true;
        if (graph_type == DIRECTED)
          for(ii=gtrie_node->child.begin(); ii!=gtrie_node->child.end(); ii++)
            (*ii)->goCondDirGraphlets(data);
        else
          for(ii=gtrie_node->child.begin(); ii!=gtrie_node->child.end(); ii++)
            (*ii)->goCondUndirGraphlets(data);
         my_data->glk--; 
         my_data->used[node] = false;
      }
      my_data->work[gtrie_node].pop_back();
    }
    if(my_data->called != -1){
      giveWork(data);
      goBackToWorkGraphlets(data);
    }
    my_data->workOrder.pop();
  }
  askForWork(data);
}
/**********************************/

/************* Pthreads ************/
void GTrieNode::calcFrequency(){
  int i;

  if (is_graph)
    for(i = 0; i < Global::num_threads; i++)
      frequency += gtrie_threads[i].frequencies[id];

  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->calcFrequency();
}
/**********************************/


int GTrieNode::countGraphsApp() {
  int aux=0;
  if (is_graph && frequency>0) aux++;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countGraphsApp();

  return aux;
}

int GTrieNode::countGraphs() {
  int aux=0;
  if (is_graph) aux=1;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countGraphs();

  return aux;
}

int GTrieNode::countGraphPaths() {
  int aux=0;
  if (is_graph) aux=depth;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countGraphPaths();

  return aux;
}

long long GTrieNode::countOccurrences() {
  long long aux=0;
  if (is_graph && frequency>0) aux+=frequency;
  
  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    aux += (*ii)->countOccurrences();

  return aux;
}

void GTrieNode::writeToFile(FILE *f) {
  int i, bits, aux, nbytes;

  int nchild       = child.size();
  nbytes=0;
  while (nchild>0) {
    nbytes++;
    nchild/=BASE_FORMAT;
  }

  int first_byte   = nbytes<<1;

  // first bit: isgraph ?
  if (is_graph) BIT_SET(first_byte,0);
  else          BIT_CLEAR(first_byte,0);
  fputc(BASE_FIRST+first_byte,f);

  nchild       = child.size();
  while (nchild>0) {
    fputc(BASE_FIRST+(nchild%BASE_FORMAT),f);
    nchild/=BASE_FORMAT;
  }

  // Connections: outgoing
  for (i=0, bits=0, aux=0; i<depth; i++, bits++) {
    if (bits==BASE_BITS) {
      fputc(BASE_FIRST+aux, f);
      aux=0;
      bits=0;
    }
    if (out[i]) BIT_SET(aux, bits);
  }
  fputc(BASE_FIRST+aux, f);


  // Connections: ingoing
  for (i=0, bits=0, aux=0; i<depth; i++, bits++) {
    if (bits==BASE_BITS) {
      fputc(BASE_FIRST+aux, f);
      aux=0;
      bits=0;
    }
    if (in[i]) BIT_SET(aux, bits);
  }
  fputc(BASE_FIRST+aux, f);

  // Previous Conditions
  if (cond_ok) aux = 0;
  else         aux = cond.size();
  fputc(BASE_FIRST+aux, f);

  if (aux>0) {
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {      
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk) {
	fputc(BASE_FIRST+1+(kk->first),f);
	fputc(BASE_FIRST+1+(kk->second),f);
      }
      fputc(BASE_FIRST, f);
    }
  }
  
  /****** Graphlets ******/
  if(Global::compute_orbits){
    if(depth >= 2 && is_graph){
      fprintf(f, "%d", graphlet_id);
      fputc('>', f);

      for(i = 0; i < depth; i++){
        fprintf(f, "%lld", orbits[i]);
        if(i < depth -1) fputc(';', f);
      }
      fputc('#', f);
    }
  }
  /***********************/
  
  fputc('\n', f);

  list<GTrieNode *>::iterator ii;
  for(ii=child.begin(); ii!=child.end(); ii++)
    (*ii)->writeToFile(f);

}

void GTrieNode::readFromFile(FILE *f) {
  int nchilds, i, j, pos, bits, ncond;
  iPair p;
  char aux, buf[MAX_BUF];
  GTrieNode *c;

  if (fgets(buf, MAX_BUF, f)) {
    aux = buf[0]-BASE_FIRST;
    if (BIT_VALUE(aux, 0)) is_graph=true;
    else                   is_graph=false;
    int nbytes = aux>>1;

    pos = 1;

    for (i=0, j=1, nchilds=0; i<nbytes; i++, j*=BASE_FORMAT) {
      aux = buf[pos++]-BASE_FIRST;
      nchilds += int(aux)*j;
    }

    // Connections: outgoing
    aux = buf[pos++]-BASE_FIRST;
    for (i=0, bits=0; i<depth; i++, bits++) {
      if (bits==BASE_BITS) {
	aux = buf[pos++]-BASE_FIRST;
	bits=0;
      }
      if (BIT_VALUE(aux, bits)) out[i] = true;
      else                      out[i] = false;
    }

    // Connections: ingoing
    aux = buf[pos++]-BASE_FIRST;
    for (i=0, bits=0; i<depth; i++, bits++) {
      if (bits==BASE_BITS) {
	aux = buf[pos++]-BASE_FIRST;
	bits=0;
      }
      if (BIT_VALUE(aux, bits)) in[i] = true;
      else                      in[i] = false;
    }

    // Previous Conditions
    aux = buf[pos++]-BASE_FIRST;
    if (aux==0) cond_ok = true;
    else        cond_ok = false;

    if (aux>0) {
      ncond = aux;
      for (i=0; i<ncond; i++) {
	list<iPair> newcond;
	while(1) {
	  aux = buf[pos++]-BASE_FIRST-1;
	  if (aux<0) break;
	  p.first = aux;
	  aux = buf[pos++]-BASE_FIRST-1;
	  p.second = aux;
	  newcond.push_back(p);
	}
	cond.push_back(newcond);
      }
    }
    
    /****** Graphlets ******/
    if(Global::compute_orbits) pos = getGraphletData(pos, buf);
    /***********************/
    
    if (buf[pos]!='\n') {
      fprintf(stderr, "ERROR: [%s] !%d!%c!\n", buf, pos, buf[pos]);
      fprintf(stderr,"[%d](%s) |", depth, is_graph?"X":" ");
      for (i=0, bits=0; i<depth; i++, bits++)
	fprintf(stderr,"%s", out[i]?"1":"0");
      fprintf(stderr,"|\n");
    }

    // Conn and nconn variables (was missing)
    for (i=0; i<depth; i++)
      if (out[i] || in[i]) 
      conn[(nconn)++] = i;
    
    
    for (i=0; i<nchilds; i++) {
      c = new GTrieNode(depth+1);
      c->readFromFile(f);
      child.push_back(c);
    }

  }

}

void GTrieNode::populateGraphTree(GraphTree *tree, char *s, int size) {
  int i, pos=depth-1;
  
  for (i=0;i<depth;i++) {
    s[pos*size+i]=out[i]?'1':'0';
    s[i*size+pos]=in[i]?'1':'0';
  }

  if (is_graph)
    tree->setString(s, frequency);

  list<GTrieNode *>::const_iterator ii, iiend;
  for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
    (*ii)->populateGraphTree(tree, s, size);
}


void GTrieNode::populateMap(mapStringLL *m, char *s, int size) {
  int i, pos=depth-1;

  for (i=0;i<depth;i++) {
    s[pos*size+i]=out[i]?'1':'0';
    s[i*size+pos]=in[i]?'1':'0';
  }

  if (is_graph && frequency>0) (*m)[s]=frequency;

  list<GTrieNode *>::const_iterator ii, iiend;
  for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
    (*ii)->populateMap(m, s, size);
}

// -------------------------------------

GTrie::GTrie() {
  _root = new GTrieNode(0);
  _root->cond_ok = _root->cond_this_ok = true;
}

GTrie::~GTrie() {
  delete _root;
}

void GTrie::insertGraphCond(Graph *g, list<iPair> *cond) {
  _root->insertCond(g, cond);
}

void GTrie::insertGraph(Graph *g) {
  _root->insert(g);
}

void GTrie::showAsText(FILE *f) {
  _root->showAsText(f);
}

float GTrie::frequencyGraphString(int size, const char *s) {
  char larger[size*size+1];
  Isomorphism::canonicalBasedNauty(s, larger, size);

  Graph *g = new GraphMatrix();
  GraphUtils::strToGraph(g, larger, size, DIRECTED); // May change later
  float aux = _root->frequencyGraph(g);
  delete g;
  
  return aux;
}


void GTrie::insertGraphString(int size, const char *s) {
  char larger[size*size+1];
  Isomorphism::canonicalBigger(s, larger, size);


  Graph *g = new GraphMatrix();
  GraphUtils::strToGraph(g, larger, size, DIRECTED); // May change later
  list<iPair> *cond = new list<iPair>;
  Conditions::symmetryConditions(g, cond);

  insertGraphCond(g, cond);

  delete cond;
  delete g;
}

// Populate g-trie with subgraphs of 'size' read from file 's'
void GTrie::readSubgraphs(int size, bool dir, char *s) {
  char buf[MAX_BUF];

  FILE *f = fopen(s, "r");
  if (!f) Error::msg(NULL);
  while (fscanf(f, "%s", buf)==1) {
    insertGraphNautyString(size, buf, dir, 1);
  }
  fclose(f);
  cleanConditions();  
    
  /******* Graphlets *******/
  if(Global::compute_orbits){
    appendOrbits();
    printOrbits();
  }
  /*************************/
}

void GTrie::insertGraphNautyString(int size, const char *s, bool dir, int label) {
  char larger[size*size+1];

  if (label==2)
    Isomorphism::canonicalBigger(s, larger, size);
  else if (label==3)
    Isomorphism::canonicalRandom(s, larger, size);
  else if (label==4)
    Isomorphism::canonicalBasedNautyReverse(s, larger, size);
  else
    Isomorphism::canonicalBasedNauty(s, larger, size);

  Graph *g = new GraphMatrix();

  if (dir)  GraphUtils::strToGraph(g, larger, size, DIRECTED);
  else      GraphUtils::strToGraph(g, larger, size, UNDIRECTED);

  g->makeArrayNeighbours();

  list<iPair> *cond = new list<iPair>;
  Conditions::symmetryConditions(g, cond);  

  insertGraphCond(g, cond);

  delete cond;
  delete g;
}

double GTrie::compressionRate() {
  int nodes  = _root->countNodes()-1;
  int graphs = _root->countGraphPaths();
  return 1 - (double)(nodes) / (double)(graphs);
}

int GTrie::maxDepth() {
  return _root->maxDepth();
}

void GTrie::showFrequency() {
  return _root->showFrequency();
}

double GTrie::countOccurrences() {
  return _root->countOccurrences();
}

int GTrie::countGraphsApp() {
  return _root->countGraphsApp();
}

int GTrie::countGraphs() {
  return _root->countGraphs();
}

void GTrie::writeToFile(char *s) {
  FILE *f;
  
  f=fopen(s,"w");
  if (f!=NULL) {
    fputs("GTRIEFORMAT VERSION 1\n", f);
    /******* Graphlets *******/
    if(Global::compute_orbits){
	  fprintf(f, "%d\n", Isomorphism::getTotalOrbits());
    }
    _root->writeToFile(f);
    /*************************/
    fclose(f);
  } else
    Error::msg("Unable to open g-trie output file \"%s\")", s);    
}

void GTrie::readFromFile(char *s) {
  FILE *f;
  char *dummy, buf[MAX_BUF];
  
  f=fopen(s,"r");
  if (!f) Error::msg(NULL);
  dummy=fgets(buf, MAX_BUF, f);
  /******* Graphlets *******/
  if(Global::compute_orbits){
    dummy=fgets(buf, MAX_BUF, f);

    total_orbits = atoi(dummy);
    printf("\tGraphlet-Trie*\n");
    _root->readFromFile(f);
  }
  /*************************/
  _root->readFromFile(f);
  fclose(f);
}

void GTrie::populateGraphTree(GraphTree *tree, int size) {
  char s[size*size+1];
  s[size*size]=0;
  _root->populateGraphTree(tree, s, size);
}

void GTrie::populateMap(mapStringLL *m, int size) {
  char s[size*size+1];
  s[size*size]=0;
  _root->populateMap(m, s, size);
}


/*void GTrie::censusSample(Graph *g, double *p) {
  int i;
  int subgraph_size = maxDepth();
  int SharedParallel::num_nodes = g->numNodes();

  _root->zeroFrequency();
  
  GTrieNode::my_data->mymap = new int[subgraph_size];
  GTrieNode::my_data->used  = new bool[SharedParallel::num_nodes];
  GTrieNode::numNodes = SharedParallel::num_nodes;
  GTrieNode::fastnei  = g->matrixNeighbours();
  GTrieNode::adjM     = g->adjacencyMatrix();
  GTrieNode::numnei   = g->arrayNumNeighbours(); 
  GTrieNode::prob     = p;

  if (g->type() == DIRECTED) GTrieNode::isdir = true;
  else                       GTrieNode::isdir = false;
  for (i=0; i<SharedParallel::num_nodes; i++)
    GTrieNode::my_data->used[i]=false;


  list<GTrieNode *>::iterator ii, iiend;
  GTrieNode::my_data->glk=0;
  for(ii=_root->child.begin(), iiend = _root->child.end(); ii!=iiend; ii++)
    if (Random::getDouble()<=GTrieNode::prob[0]) {
      (*ii)->goCondSample();
    }
}


void GTrieNode::goCondSample() {
  int i, j, ci, mylim, glaux;
  int ncand;

  if (!cond_ok) {
    list< list<iPair> >::const_iterator jj, jjend;
    list<iPair>::const_iterator kk, kkend;
    for (jj=cond.begin(), jjend=cond.end(); jj!=jjend; ++jj) {
      for (kk=(*jj).begin(), kkend=(*jj).end(); kk!=kkend; ++kk)
	if ( my_data->mymap[kk->first]>my_data->mymap[kk->second])
	  break;
      if (kk==kkend) break;
    }
    if (jj==jjend) return;
  }

  ncand=0;
  j=ci=INT_MAX;
  if (nconn == 0) ncand = numNodes;
  else {
    for (i=0; i<nconn; i++) {
      glaux = numnei[my_data->mymap[conn[i]]];
      if (glaux<j) {
	ci = my_data->mymap[conn[i]];
	j = glaux;
      }
    }
    ncand = j;
  }
  int cand[ncand];

  ncand=0;
  if (nconn == 0) { // We are at a node with no connections to ancestors
    for (i=numNodes-1; i>=0; i--)
      if (!my_data->used[i]) 
	cand[ncand++]=i;
  } else {
    for (i=0; i<j; i++) {
      glaux = fastnei[ci][i];
      if (!my_data->used[glaux])
	cand[ncand++]=glaux;
    }
  }

  if (cond_this_ok) mylim = 0;
  else {
    list< list<int> >::const_iterator jjj;
    list<int>::const_iterator kkk;
    mylim = INT_MAX;
    for (jjj=this_node_cond.begin(); jjj!=this_node_cond.end(); ++jjj) {
      glaux = -1;
      for (kkk=(*jjj).begin(); kkk!=(*jjj).end(); ++kkk)
	if (my_data->mymap[*kkk]>glaux) glaux = my_data->mymap[*kkk];
      if (glaux<mylim) mylim=glaux;
    }
  }

  for (ci=ncand-1; ci>=0; ci--) {
    i = cand[ci];
    if (i<mylim) break;
    my_data->mymap[my_data->glk] = i;

    if (isdir) {  
      for (j=0; j<my_data->glk; j++)
	if (in[j]  != adjM[my_data->mymap[j]][i] ||
	    out[j] != adjM[i][my_data->mymap[j]])
	  break;
    } else {
      for (j=0; j<my_data->glk; j++)
	if (in[j]  != adjM[my_data->mymap[j]][i])
	  break;
    }
    if (j<my_data->glk) continue;

    if (is_graph) {      
      if (Global::show_occ) {
	for (int k = 0; k<=my_data->glk; k++)
	  for (int l = 0; l<=my_data->glk; l++)
	    fputc(adjM[my_data->mymap[k]][my_data->mymap[l]]?'1':'0', Global::occ_file);
	fputc(':', Global::occ_file);
	for (int k = 0; k<=my_data->glk; k++)
	  fprintf(Global::occ_file, " %d", my_data->mymap[k]+1);
	fputc('\n', Global::occ_file);
      }
      frequency[my_data->tid]++;;
    }

    my_data->used[i]=true;
    my_data->glk++;
    list<GTrieNode *>::const_iterator ii, iiend;
    for(ii=child.begin(), iiend = child.end(); ii!=iiend; ++ii)
      if (Random::getDouble() <= prob[my_data->glk]) {
	(*ii)->goCondSample();
      }
    my_data->glk--;
    my_data->used[i]=false;
  }
}
*/

void GTrieNode::clean(int a, int b) {
  list< list<iPair> >::iterator ii;
  list<iPair> ::iterator jj;

  if (cond.size()>0) {

    for (ii=cond.begin(); ii!=cond.end(); ii++) {
      for (jj=ii->begin(); jj!=ii->end();)
	  if (jj->first==a && jj->second==b)
	    jj = ii->erase(jj);
	  else
	    jj++;
      if (ii->size()==0) {
	ii->clear();
	cond_ok = true;
      }
    }    
  }

  list<GTrieNode *>::iterator cc;
  for(cc=child.begin(); cc!=child.end(); cc++)
    (*cc)->clean(a, b);	

}

void GTrieNode::cleanConditions() {
  int a, b;
  list< list<iPair> >::iterator ii;
  list<iPair> ::iterator jj, kk;

  if (cond.size()>0) {
    for (jj=cond.begin()->begin(); jj!=cond.begin()->end();jj++) {
      a=jj->first;
      b=jj->second;
      for (ii=cond.begin(); ii!=cond.end(); ii++) {
	for (kk=ii->begin(); kk!=ii->end(); kk++)
	  if (kk->first==a && kk->second==b) break;
	if (kk==ii->end()) break;
      }
      if (ii==cond.end())  {
	list<GTrieNode *>::iterator ii;
	for(ii=child.begin(); ii!=child.end(); ii++)
	  (*ii)->clean(a, b);	
      }
    }
  }

  list<GTrieNode *>::iterator cc;
  for(cc=child.begin(); cc!=child.end(); cc++)
    (*cc)->cleanConditions();
}


void GTrie::cleanConditions() {
  _root->cleanConditions();
}

/******* Graphlets *******/
vector<GTrieNode *> GTrieNode::computeOrbits(vector<GTrieNode *> nodes){
  int i;

  //printf("Depth: %d\n", depth);
  if(depth >= 2 && is_graph){
    for(i = 0; i < (depth) * (depth); i++){
      //printf("%c ", full_graph[i]);
      //if((i+1) % (depth) == 0) printf("| ");
    }

    graphlet_id = graphlet_counter++;

    int *ret_orbits = Isomorphism::computeOrbits(full_graph, depth);
    for(i = 0; i < depth; i++) orbits[i] = ret_orbits[i];
  }

  list<GTrieNode *>::iterator cc;

  for(cc=child.begin(); cc!=child.end(); cc++) nodes.push_back(*cc);

  return nodes;
}

void GTrie::appendOrbits() {
  vector<GTrieNode *> nodes;
  nodes.push_back(_root);

  buildNodeList(nodes);
}

void GTrie::buildNodeList(vector<GTrieNode *> nodes){
  int size = nodes.size();

  for(int i = 0; i < size; i++) {
    nodes = (*nodes.begin())->computeOrbits(nodes);
    nodes.erase(nodes.begin());
  }

  if(nodes.size() > 0) buildNodeList(nodes);
}


void GTrieNode::printOrbits(){
  int i;

  //printf("Depth: %d\n", depth);
  if(depth >= 2 && is_graph){
    for(i = 0; i < (depth) * (depth); i++){
      printf("%c ", full_graph[i]);
      if((i+1) % (depth) == 0) printf("| ");
    }

    printf("  || G(%d) => Orbits: [", graphlet_id);
    for(i = 0; i < depth; i++)
        printf("%lld ", orbits[i]);
    printf("]\n");
  }

  list<GTrieNode *>::iterator cc;
  for(cc=child.begin(); cc!=child.end(); cc++)
      (*cc)->printOrbits();
}

void GTrie::printOrbits() {
  _root->printOrbits();
}

/****** Graphlets ******/
int GTrieNode::getGraphletData(int pos, const char* s){
  int i = pos;
  int j = 0;
  int orbit_count = 0;

  int num_digits = GTrie::total_orbits > 0 ? (int) log10 ((double) GTrie::total_orbits) + 1 : 1;

  char num[num_digits];

  if(depth >= 2 && is_graph){
    while(s[i] != '>'){
        num[j] = s[i];
        i++; j++;
    }
    memset(num+j, 0, num_digits-j);
    graphlet_id = atoi(num);
    //printf("G%2d", graphlet_id);

    i++; j = 0;
    pos = i;

    //printf(" => [");
    while(i > 0){
      while(s[i] != ';'){
        if(s[i] == '#') {
          i = -1;
          pos++;
          break;
        }
        num[j] = s[i];
        i++; j++;
        pos++;
      }
      memset(num+j, 0, num_digits-j);
      orbits[orbit_count] = atoi(num);
      //printf("%lld ", orbits[orbit_count]);
      if(i == -1) break;
      i++; j = 0;
      pos++;
      orbit_count++;
    }
    //printf("] (Orbits)");
  }

  return pos;
}


/*************************/

/******** Graphlets ******/
void GTrie::freeOrbitFrequencies(int num_nodes){
  for(int i = 0; i < num_nodes; i++) free(orbit_frequency[i]);
  free(orbit_frequency);
}

void GTrie::writeGraphletCountToFile(int numNodes, FILE* out) {
  std::map<int, long long int> distinct;

  for(int i = 0; i < numNodes; i++){
    for(int j = 0; j < total_orbits; j++){
      fprintf(out,"%lld ",orbit_frequency[i][j]);
    }
    fprintf(out,"\n");
  }
}

/*************************/

// -------------------------------------


