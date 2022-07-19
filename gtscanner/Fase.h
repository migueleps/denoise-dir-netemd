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

#ifndef _FASE_
#define _FASE_

#include "Common.h"
#include "Graph.h"
#include "LSLabeling.h"
#include "FaseGTrie.h"
#include "Isomorphism.h"
#include "Random.h"
#include "GraphTree.h"
/************* Pthreads ************/
#include "SharedParallel.h"
/***********************************/

class GraphTree; // forward declaration

/*! This class implements the basic FaSE subgraph enumeration algorithm */
class Fase {
 private:
  static int graphSize;
  static int *seqMap;
  static int subNum;
  static char globStr[MAX_LABEL];
  static int **seqExtCpy;
  static char* LSLabel(int w, int subSize);
  static void ExtendSubgraph(int *ext, int extNum, int subSize, int v);
  static GraphTree *_sg;

 public:
  /************* Pthreads ************/
  static int K;
  static Graph *G;
  /***********************************/
  static int typeLabel;
  static long long int MotifCount;
  static bool directed;
  static double samp_prob[MAX_MOTIF_SIZE];
  static long long int getTypes();
  static long long int getLeafs();
  static long long int getNodes();
  static void destroy();
  static void listClasses(FILE* f, double prob);
  static void listTree(FILE* f);
  static void census(Graph *_G, int _K, GraphTree *_sg);
  /************* Pthreads ************/
  static void parallelExtendSubgraph(void* data, int *ext, int extNum, int subSize, int v); 

  static void *startCensus(void *data);
  static void censusOnNodes(int begin, int end, void *data);
  static void askForWork(void *data);
  static int  getThreadToCall(int tid);
  static void simpleAnswer(void *data);
  static void giveWork(void *data);
  static void createWorkTop(void *data, int current, int end);
  static void createWork(void *data, int* Vext, int current, int end, int subSize, int v);
  static void goBackToWork(void *data);
  static std::queue<trie*> correctOrder(std::queue<trie*> old_order, unordered_map<trie*, int > trieDepth);
  /***********************************/

};

#endif
