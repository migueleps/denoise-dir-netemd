/* -------------------------------------------------

//                                                 
//  88888888888           ad88888ba   88888888888  
//  88                   d8"     "8b  88           
//  88                   Y8,          88           
//  88aaaaa  ,adPPYYba,  `Y8aaaaa,    88aaaaa      
//  88"""""  ""     `Y8    `"""""8b,  88"""""      
//  88       ,adPPPPP88          `8b  88           
//  88       88,    ,88  Y8a     a8P  88           
//  88       `"8bbdP"Y8   "Y88888P"   88888888888  
//                                                 
//

Pedro {Paredes, Ribeiro} - DCC/FCUP (algorithm)
David Aparício           - DCC/FCUP (parallel implementation)
----------------------------------------------------
Base FaSE implementation

---------------------------------------------------- */

#include "Fase.h"
#include "Common.h"
#include "Timer.h"

long long int Fase::MotifCount = 0;
Graph*     Fase::G;
int        Fase::K;
int        Fase::typeLabel;
bool       Fase::directed;
int        Fase::subNum;
int        Fase::graphSize;
double     Fase::samp_prob[MAX_MOTIF_SIZE];
GraphTree* Fase::_sg;

/************* Pthreads ************/
#include <pthread.h>

struct thread_data_fase{
   int tid;

   int* mymap;
   int** extCpy;

   bool answered;
   bool sorry; //sorry, I had no work
   bool waiting;
   int  called;
   char *labelStr;
   int sharedAtSize;
   bool did_work;
   
   trie* pointer;
   
   Graph *G;
   FaseGTrie* myGtrie;
   LSLabeling* myLabeling;
   
   Isomorphism* mynauty;

   long long int myMotifCount;
   unordered_map<trie*, vector<int> > work; // gtrie_location -> nodes
   queue<trie*> workOrder;       // map doesn't keep order
   unordered_map<trie*, int > trieDepth;
};

struct thread_data_fase *fase_threads;

int NO_MORE_WORK = -1;

pthread_mutex_t nauty_mutex     = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t gtrie_mutex     = PTHREAD_MUTEX_INITIALIZER;

const int MAX_LOCKS = 500000;
pthread_mutex_t leaf_mutexes[MAX_LOCKS];

trie* topNode;

FaseGTrie* parGtrie;
LSLabeling* parLabeling;

/***********************************/

/*! List all k-subgraphs in a larger graph
    \param G the graph to be explored
    \param K the size of the subgraphs*/
void Fase::census(Graph *_G, int _K, GraphTree *sg){
  int i, j;
  
  K     = _K;
  G     = _G;
  _sg   = sg;
  
  directed = (G->type() == DIRECTED);
  sg->zeroFrequency();
  graphSize = G->numNodes();
  
  parGtrie    = new FaseGTrie();
  parLabeling = new LSLabeling();
  parGtrie->init(Global::num_threads);
  parLabeling->init(G);

  //printf("Top Pointer: %p\n", parGtrie->top);

  for(i = 0; i < MAX_LOCKS; i++)
    pthread_mutex_init(&leaf_mutexes[i], NULL);
  
  char s[50]; 
  
  /************* Pthreads ************/
  //printf("Using %d thread(s)\n", Global::num_threads);
  
  SharedParallel::num_nodes = G->numNodes();
  
  fase_threads = new thread_data_fase[Global::num_threads];
  pthread_t thread[Global::num_threads];
  
  for (i = 0; i < Global::num_threads; i++){
    fase_threads[i].tid          = i;
    fase_threads[i].mymap        = new int[K];
    fase_threads[i].extCpy       = new int*[K];
    fase_threads[i].myGtrie      = parGtrie;    //new GTrie();      // 
    fase_threads[i].myLabeling   = parLabeling; //new LSLabeling(); // 
    fase_threads[i].myMotifCount = 0;
    fase_threads[i].labelStr     = new char[MAX_LABEL];
    fase_threads[i].waiting      = false;
    fase_threads[i].sorry        = false;
    fase_threads[i].called       = -2;
    fase_threads[i].sharedAtSize = 0;
    fase_threads[i].mynauty      = new Isomorphism(); // CAUSES LEAK
    fase_threads[i].pointer      = parGtrie->top;
  }
  for (i = 0; i < Global::num_threads; i++) //cannot be in the above cycle -> thread calls one that has no struct yet
    pthread_create(&thread[i], NULL, Fase::startCensus, (void *) &fase_threads[i]);
    
  for (i = 0; i < Global::num_threads; i++)
    pthread_join(thread[i], NULL);
  
  parGtrie->populateGraphTree(_sg, fase_threads[0].mynauty);
  
  for (i = 0; i < Global::num_threads; i++){
	fase_threads[i].mynauty->finishNauty();
	delete fase_threads[i].mynauty;
    delete [] fase_threads[i].mymap;
    for (j = 0; j < K; j++)
      delete [] fase_threads[i].extCpy[j];
    delete [] fase_threads[i].extCpy;
    delete [] fase_threads[i].labelStr;
    MotifCount += fase_threads[i].myMotifCount;
  }
  delete [] fase_threads;
  
  /***********************************/
}

/************* Pthreads ************/
void* Fase::startCensus(void *data){
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  my_data->called = -1;
  
  timeval t1, t2;
  gettimeofday(&t1, NULL);
  
  my_data->mynauty->initNauty(Fase::K, Fase::directed);
  
  for (int i = 0; i < Fase::K; i++)
    my_data->extCpy[i] = new int[SharedParallel::num_nodes];
  
  if(my_data->tid >= SharedParallel::num_nodes){ //too many threads for this graph, wait for work
    sleep(2);
    askForWork(data);
  }
  else{
    censusOnNodes(my_data->tid, SharedParallel::num_nodes-1, data);
  }
  
  gettimeofday(&t2, NULL);
  double z = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
  //printf("%f -> Tid(%d): I'm done\n", z, my_data->tid);
  
  return NULL;
}
/************************************/

/************* Pthreads ************/
void Fase::censusOnNodes(int begin, int end, void *data) {
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  int j, node, neighNum, extNum;
  int *neigh;
  
  int *Vext             = new int[Fase::G->numNodes()];
  
  my_data->trieDepth[my_data->pointer] = 1; 
  
  for (node=begin; node<=end; node += Global::num_threads){
    if(my_data->called != -1 && my_data->did_work){
      break;
    }
    my_data->mymap[0]   = node;
    neigh               = Fase::G->arrayNeighbours(node);
    neighNum            = Fase::G->numNeighbours(node);
    extNum = 0;
    for (j = 0; j < neighNum; j++)
      if (neigh[j] > node)
        Vext[extNum++] = neigh[j];
        
    Fase::parallelExtendSubgraph(my_data, Vext, extNum, 1, node);
  }
  
  delete [] Vext;
  
  if(my_data->called != -1 && my_data->did_work){
    //fprintf(fp[my_data->tid], "!!!!!!!STOPPED!!!!!!\n");
    my_data->pointer = (my_data->myGtrie)->jump_th(my_data->pointer);
    createWorkTop(my_data, node, end);
    giveWork(data); // simpleAnswer(my_data);// 
    goBackToWork(my_data);
    return;
  }
  
  askForWork(data);
}
/************************************/

int intRand(const int & min, const int & max) {
    static std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}

/************* Pthreads ************/
void Fase::parallelExtendSubgraph(void *data, int *Vext, int extNum, int subSize, int v){
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  //printf("%p | %p\n", my_data->pointer, (my_data->myGtrie)->jump_th(my_data->pointer));
  
  //printf("subsize: %d\n", subSize);
  
  int i;
  
  if (subSize == K - 1){
    for (i = 0; i < extNum; i++){
      my_data->mymap[subSize] = Vext[i];
      my_data->myMotifCount++;
      
      if (Global::show_occ) {
        for (int kk = 0; kk<=subSize; kk++)
          for (int ll = 0; ll<=subSize; ll++)
            fputc(G->hasEdge(my_data->mymap[kk],my_data->mymap[ll])?'1':'0', Global::occ_file);
        fputc(':', Global::occ_file);
        for (int kk = 0; kk<=subSize; kk++)
          fprintf(Global::occ_file, " %d", my_data->mymap[kk]+1);
        fputc('\n', Global::occ_file);
      }
      
      char* addLabel = (my_data->myLabeling)->Label(my_data->mymap, subSize, Vext[i], typeLabel, my_data->labelStr, directed);
      trie* toinsert = (my_data->myGtrie)->checkInsert_th1(addLabel, my_data->pointer, my_data->tid);
      
      if(toinsert == NULL){
        size_t a = (size_t)my_data->pointer;
        pthread_mutex_lock(&leaf_mutexes[a % MAX_LOCKS]);
        //pthread_mutex_lock(&leaf_mutexes[0]);
        toinsert = (my_data->myGtrie)->checkInsert_th1(addLabel, my_data->pointer, my_data->tid);
        if(toinsert == NULL){
          my_data->pointer = (my_data->myGtrie)->insert_th(my_data->pointer, addLabel, subSize);
        }
        else{
          my_data->pointer = toinsert;
        }
        pthread_mutex_unlock(&leaf_mutexes[a % MAX_LOCKS]);
      }
      else{
        my_data->pointer = toinsert;
      }
      
      /*for (int k = 0; k<=subSize; k++)
        fprintf(fp[my_data->tid], " %d", my_data->mymap[k]);
      fputc('\n', fp[my_data->tid]);*/
      
      if (!((my_data->myGtrie)->getCurrentLeaf_th(my_data->pointer))){ //if it's a new g-trie node, set canonical label
        char s[K * K + 1];
        s[0] = '\0';
        
        pthread_mutex_lock(&nauty_mutex);
        (my_data->mynauty)->canonicalStrNauty(G, my_data->mymap, s);
        pthread_mutex_unlock(&nauty_mutex);
        
        pthread_mutex_lock(&gtrie_mutex);
        my_data->pointer = (my_data->myGtrie)->setCanonicalLabel_th(s, my_data->pointer, my_data->tid);
        pthread_mutex_unlock(&gtrie_mutex);
      }
      (my_data->myGtrie)->incCount(my_data->tid, my_data->pointer);
      my_data->pointer = (my_data->myGtrie)->jump_th(my_data->pointer); //goes back to the parent g-trie node to compute remaining Vext
    }
    return;
  }
  
  int j, o;
  int extCpyNum;

  for (j = 1; j < extNum; j++)
    my_data->extCpy[subSize][extNum - j - 1] = Vext[j]; //populate extCpy[subsize] with Vext: e.g. [1,5,7] => [7, 5]

  for (i = 0; i < extNum; i++){
    if(my_data->called != -1 && my_data->did_work){
      createWork(my_data, Vext, i, extNum, subSize, v);
      my_data->sharedAtSize++;
      return;
    }
  
    extCpyNum    = extNum - i - 1;
    int *neigh   = G->arrayNeighbours(Vext[i]); //  N(u)
    int nneigh   = G->numNeighbours(Vext[i]);   // |N(u)|
    
    for (j = 0; j < nneigh; j++){
      //w = neigh[j]
      if (neigh[j] <= v) // w > Vs[0] requirement
        continue;
      for (o = 0; o < subSize; o++)
        if (neigh[j] == my_data->mymap[o] || G->isConnected(neigh[j], my_data->mymap[o])) // w not in Vs and w in Nexc(u, Vs)
          break;
      if (o == subSize)
        my_data->extCpy[subSize][extCpyNum++] = neigh[j]; //adds valid nodes w to Vext: e.g [7,5] U [8,12]
    }

    char* addLabel = (my_data->myLabeling)->Label(my_data->mymap, subSize, Vext[i], typeLabel, my_data->labelStr, directed);
    
    //printf("-> Current: %p\n", my_data->pointer);
    trie* toinsert = (my_data->myGtrie)->checkInsert_th1(addLabel, my_data->pointer, my_data->tid);
    if(toinsert == NULL){
      size_t a = (size_t)my_data->pointer;
      pthread_mutex_lock(&leaf_mutexes[a % MAX_LOCKS]);
      //pthread_mutex_lock(&leaf_mutexes[0]);
      toinsert = (my_data->myGtrie)->checkInsert_th1(addLabel, my_data->pointer, my_data->tid);
      if(toinsert == NULL){
        my_data->pointer = (my_data->myGtrie)->insert_th(my_data->pointer, addLabel, subSize);
      }
      else{
        my_data->pointer = toinsert;
      }
      pthread_mutex_unlock(&leaf_mutexes[a % MAX_LOCKS]);
    }
    else{
      my_data->pointer = toinsert;
    }

    my_data->mymap[subSize] = Vext[i];
    parallelExtendSubgraph(my_data, my_data->extCpy[subSize], extCpyNum, subSize + 1, v);
    my_data->did_work = true;
    my_data->pointer = (my_data->myGtrie)->jump_th(my_data->pointer);
  }
}
/************************************/

/************* Pthreads ************/
void Fase::askForWork(void* data){
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  int to_call;
  
  #if (DB_LEVEL >= 2)
  printf("Tid(%d): Need work\n", my_data->tid);
  #endif
  
  ////////////checks if it was itself called; checks if there is no more work
  pthread_mutex_lock(&SharedParallel::waiting_mutex);
  
  if(my_data->called != -1) { 
    pthread_mutex_lock(&SharedParallel::answer_mutex);
    #if (DB_LEVEL >= 2)
    printf("Tid(%d): I was called but I have no work\n", my_data->tid);
    #endif
    fase_threads[my_data->called].sorry  = true;
    fase_threads[my_data->called].answered  = true;
    
    my_data->called = -1;
    pthread_mutex_unlock(&SharedParallel::answer_mutex);
    pthread_cond_signal(&SharedParallel::answer_cond);
  }

  SharedParallel::waiting_threads++;
  my_data->waiting = true;
  #if (DB_LEVEL >= 1)
  printf("Tid(%d): Threads waiting: %d\n", my_data->tid, SharedParallel::waiting_threads);
  #endif
  
  if(SharedParallel::waiting_threads == Global::num_threads){
    #if (DB_LEVEL >= 1)
    printf("Tid(%d): There is no more work, we should all stop!\n", my_data->tid);
    //askForWork(data);
    #endif
    SharedParallel::no_more_work = true;
  }
  
  pthread_mutex_unlock(&SharedParallel::waiting_mutex);
  ////////////
  
  ////////////calls a thread for work
  pthread_mutex_lock(&SharedParallel::call_mutex);
  
  if(SharedParallel::no_more_work){
    #if (DB_LEVEL >= 2)
    printf("Tid(%d): No more work!\n", my_data->tid);
    #endif
    pthread_mutex_unlock(&SharedParallel::call_mutex);
    return; //pthread_exit(NULL);
  }
  
  my_data->answered = false;
  
  pthread_mutex_lock(&SharedParallel::waiting_mutex);
  to_call = getThreadToCall(my_data->tid);
  
  if(to_call == NO_MORE_WORK){
    #if (DB_LEVEL >= 2)
    printf("Tid(%d): No more work!\n", my_data->tid);
    #endif
    pthread_mutex_unlock(&SharedParallel::call_mutex);
    pthread_mutex_unlock(&SharedParallel::waiting_mutex);
    return; //pthread_exit(NULL);
  }
  
  fase_threads[to_call].called = my_data->tid;
  pthread_mutex_unlock(&SharedParallel::waiting_mutex);

  #if (DB_LEVEL >= 1)
  printf("Tid(%d): Calling Tid(%d)\n", my_data->tid, to_call); 
  #endif
  
  /////////////////////////////////// waits for answer
  pthread_mutex_lock(&SharedParallel::answer_mutex);
  while(my_data->answered == false){
    pthread_cond_wait(&SharedParallel::answer_cond, &SharedParallel::answer_mutex);
  }  
  pthread_mutex_unlock(&SharedParallel::answer_mutex);
  pthread_mutex_unlock(&SharedParallel::call_mutex);
  ////////////////////////////////////
  
  pthread_mutex_lock(&SharedParallel::waiting_mutex); 
  SharedParallel::waiting_threads--;
  my_data->waiting = false;
  
  #if (DB_LEVEL >= 1)
  printf("Tid(%d): Threads waiting: %d\n", my_data->tid, SharedParallel::waiting_threads);
  #endif
  pthread_mutex_unlock(&SharedParallel::waiting_mutex);
  
  #if (DB_LEVEL >= 1)
  printf("Tid(%d): Sleeping...\n", my_data->tid); 
  #endif
  
  //sleep(2);
  //askForWork(data);
  
  if (my_data->sorry) { 
    #if (DB_LEVEL >= 2)
    printf("Tid(%d): Thread had no work\n", my_data->tid);
    #endif
    my_data->sorry = false; 
    askForWork(data);
  } else { 
    #if (DB_LEVEL >= 1)
    printf("Tid(%d): I received work\n", my_data->tid); 
    #endif
    goBackToWork(data);
    //askForWork(data);
  }
}
/***********************************/

/************* Pthreads ************/
int Fase::getThreadToCall(int tid){
  if(SharedParallel::no_more_work){
    return NO_MORE_WORK;
  }
  else{
    int to_call = tid;
    to_call = Random::getInteger(0, Global::num_threads-1);
    
    if(to_call != tid)
      if(fase_threads[to_call].waiting == false) 
        if(fase_threads[to_call].called == -1)
          return to_call;
    
    return getThreadToCall(tid);
  }
}
/**********************************/

/************* Pthreads ************/
void Fase::createWork(void *data, int* Vext, int current, int end, int subSize, int v){
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  int i;

  trie* gtrie_location = my_data->pointer;
  
  #if (DB_LEVEL >= 3)
  printf("=============================\n");
  
  printf("Tid(%d): In gtrie position: %p\n", my_data->tid, gtrie_location);
  
  printf("Mymap: ");
  for(i = 0; i < subSize; i++){
    printf("%d ", my_data->mymap[i]);
  }
  printf("\n");
  
  printf("Vext: ");
  for(i = current; i < end; i++){
    printf("%d ", Vext[i]);
  }
  printf("\n");
  printf("=============================\n");
  #endif
  
  vector<int> nodes;

  for(i = current; i < end; i++){
    nodes.push_back(Vext[i]);
  }

  my_data->workOrder.push(gtrie_location);
  my_data->work[gtrie_location] = nodes;
  my_data->trieDepth[gtrie_location] = subSize;
  
}
/***********************************/

/************* Pthreads ************/
void Fase::createWorkTop(void* data, int current, int end){
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  int i;

  trie* gtrie_location = my_data->pointer;
  
  #if (DB_LEVEL >= 3)
  printf("=============================\n");
  
  printf("Tid(%d): In gtrie position: %p\n", my_data->tid, gtrie_location);
  
  printf("Top nodes: ");
  for(i = current; i <= end; i+= Global::num_threads){
    printf("%d ", i);
  }
  printf("\n");
  
  printf("=============================\n");
  #endif
  
  vector<int> nodes;

  for (i = current; i <= end; i+= Global::num_threads){
    nodes.push_back(i);
  }

  my_data->workOrder.push(gtrie_location);
  my_data->work[gtrie_location] = nodes;
}
/***********************************/

/************* Pthreads ************/
void Fase::goBackToWork(void* data){
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  my_data->did_work = false;
  int node, i, j, neighNum, extNum;
  
  int *Vext = new int[Fase::G->numNodes()];

  my_data->workOrder = correctOrder(my_data->workOrder, my_data->trieDepth);

  while((ssize_t)my_data->workOrder.size() > 0){
    trie* gtrie_node = my_data->workOrder.front();
    trie* new_node;
    
    //printf("GtrieNode: %p\n", gtrie_node);
    
    std::reverse(my_data->work[gtrie_node].begin(), my_data->work[gtrie_node].end()); //necessary?
    
    while((ssize_t)my_data->work[gtrie_node].size() > 0){
      if(my_data->called != -1 && my_data->did_work){
        break;
      }
      if(gtrie_node == NULL){
        if((ssize_t)my_data->work[gtrie_node].size() == 1)
          my_data->workOrder.pop();
        node =  my_data->work[gtrie_node].back();
        my_data->mymap[0] = node;
        my_data->work[gtrie_node].pop_back();
        
        int * neigh         = Fase::G->arrayNeighbours(node);
        neighNum            = Fase::G->numNeighbours(node);
        extNum = 0;
        for (j = 0; j < neighNum; j++)
          if (neigh[j] > node)
            Vext[extNum++] = neigh[j];
        new_node = parGtrie->top;
        //printf("!!New node: %p -> depth %d\n!!", new_node, my_data->trieDepth[new_node]);
      }
      else{
        my_data->workOrder.pop();
        extNum = 0;
        while((ssize_t)my_data->work[gtrie_node].size() > 0){
          Vext[extNum] = my_data->work[gtrie_node].back();
          my_data->work[gtrie_node].pop_back();
          extNum++;
        }
        new_node = gtrie_node;
      }
      
      my_data->pointer = new_node;
      
      /*printf("Tid(%d): |mymap: ", my_data->tid);
      for(int l=0; l < my_data->trieDepth[new_node]; l++)
        printf("%d ", my_data->mymap[l]);
      printf("|\n");
      
      printf("Tid(%d): |Vext: ", my_data->tid);
      for(int l=0; l < extNum; l++)
        printf("%d ", Vext[l]);
      printf("|\n");*/
      
      Fase::parallelExtendSubgraph(my_data, Vext, extNum, my_data->trieDepth[new_node], my_data->mymap[0]);
    }
    
    if(my_data->called != -1 && my_data->did_work){
      delete[] Vext;
      
      giveWork(data); //simpleAnswer(my_data);//
      goBackToWork(data);
      return;
    }
    //printf("Tid(%d): !Depth %d done! -> gtrie %p\n", my_data->tid, my_data->trieDepth[gtrie_node], new_node);
  }
  
  delete[] Vext;
  askForWork(data);
}
/***********************************/

/************* Pthreads ************/
queue<trie*> Fase::correctOrder(queue<trie*> old_order, unordered_map<trie*, int > trieDepth){
  queue<trie*> new_order;
  int max_depth = 0;

  while(!old_order.empty()){
    trie* this_node = old_order.front();
    if(trieDepth[this_node] > max_depth) { max_depth = trieDepth[this_node]; }
    old_order.pop();
    new_order.push(this_node);
  }

  while(true){
    if(trieDepth[new_order.front()] < max_depth){
      trie* this_node = new_order.front();
      new_order.pop();
      new_order.push(this_node);
    } else
      break;
  }
  return new_order;
}
/*************************************/

void Fase::giveWork(void* data){
  int node;
  bool my_turn   = true;
  bool more_work = false;
  
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  #if (DB_LEVEL >= 1)
  printf("Tid(%d): Called by Tid(%d)\n", my_data->tid, my_data->called);
  #endif

  pthread_mutex_lock(&SharedParallel::answer_mutex);

  queue<trie*> myWorkOrder;
  queue<trie*> yourWorkOrder;

  unordered_map<trie*, vector<int> > mywork;
  unordered_map<trie*, vector<int> > yourwork;
  
  my_data->workOrder = correctOrder(my_data->workOrder, my_data->trieDepth);
  int my_depth = my_data->trieDepth[my_data->workOrder.front()]; //maybe -1
  
  while(!my_data->workOrder.empty()){
    trie* gtrie_node = my_data->workOrder.front();
    
    #if (DB_LEVEL >= 3)
    printf("Tid(%d) Gtrie location: %p(depth:%d) ==> ", my_data->tid, gtrie_node, my_data->trieDepth[gtrie_node]);
    #endif
  
    vector<int> my_nodes;
    vector<int> your_nodes;
    
    if(my_data->work[gtrie_node].size() > 0) more_work = true;
    
    if(my_data->trieDepth[gtrie_node] == 0){
      const int size = my_data->work[gtrie_node].size();
    
      for(int i = 0; i < size; i++){
        node = my_data->work[gtrie_node][i];
        if(my_turn){
          my_nodes.push_back(node);
          #if (DB_LEVEL >= 3)
          printf("I get: %d! ", node);
          #endif
        }
        else{
          your_nodes.push_back(node);
          #if (DB_LEVEL >= 3)
          printf("You get: %d! ", node);
          #endif
        }
        my_turn = !my_turn;
      }
    }
    else{
      const int size = my_data->work[gtrie_node].size();
      
      for(int i = 0; i < size; i++){
        node = my_data->work[gtrie_node][i];
        if(my_turn){
          my_nodes.push_back(node);
          #if (DB_LEVEL >= 3)
          printf("I get: %d! ", node);
          #endif
        }
        else{
          your_nodes.push_back(node);
          #if (DB_LEVEL >= 3)
          printf("You get: %d! ", node);
          #endif
        }
      }
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
    
    #if (DB_LEVEL >= 3)
    printf("\n");
    #endif
  }
  
  memcpy(fase_threads[my_data->called].mymap, my_data->mymap, sizeof(int) * my_depth);

  my_data->work.clear();
  my_data->work = mywork;
  my_data->workOrder = myWorkOrder;
  fase_threads[my_data->called].work = yourwork;
  fase_threads[my_data->called].workOrder = yourWorkOrder;
  //fase_threads[my_data->called].trieDepth.clear();
  fase_threads[my_data->called].trieDepth = my_data->trieDepth;
  fase_threads[my_data->called].answered = true;
  
  if(!more_work) { 
    #if (DEBUG >= 1)
    printf("Tid(%d): No more work to give\n", my_data->tid);
    #endif
    fase_threads[my_data->called].sorry  = true;
  }
  my_data->called = -1;

  pthread_mutex_unlock(&SharedParallel::answer_mutex);
  pthread_cond_signal(&SharedParallel::answer_cond);
  
  if(!more_work) askForWork(data);
}

/************* Pthreads ************/
void Fase::simpleAnswer(void *data){ //debug function
  struct thread_data_fase *my_data;
  my_data = (struct thread_data_fase *) data;
  
  pthread_mutex_lock(&SharedParallel::answer_mutex);
  
  fase_threads[my_data->called].answered = true;

  #if (DB_LEVEL >= 1)
  printf("Tid(%d): Answered Tid(%d)\n", my_data->tid, my_data->called); 
  #endif
  
  my_data->called = -1;

  pthread_mutex_unlock(&SharedParallel::answer_mutex);
  pthread_cond_signal(&SharedParallel::answer_cond);
}
/***********************************/

long long int Fase::getTypes()
{
  return parGtrie->getCanonicalNumber();
}

long long int Fase::getLeafs()
{
  return parGtrie->getClassNumber();
}

long long int Fase::getNodes()
{
  return parGtrie->getNodeNumber();
}

void Fase::listTree(FILE* f)
{
  return parGtrie->listGtrie(f);
}

void Fase::listClasses(FILE* f, double prob)
{
  return parGtrie->listClasses(f, prob);
}

void Fase::destroy()
{
  parGtrie->destroy();
  delete parLabeling;
}
