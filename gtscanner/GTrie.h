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

#ifndef _GTRIE_
#define _GTRIE_

#include "Common.h"
#include "Graph.h"
/************* Pthreads ************/
#include "SharedParallel.h"
/***********************************/

#define BASE_FORMAT      95
#define BASE_FIRST       ' '
#define BASE_BITS        6

class GraphTree; // forward declaration

class GTrieNode {
 private:
  bool _is_intset_included(list<int> a, list<int> b);
  bool _is_pairset_included(list<iPair> a, list<iPair> b);

 public:

  static bool **adjM;
  static int **fastnei;
  static int *numnei;
  static int numNodes;
  static bool isdir;
  static double *prob;
  /******* Graphlets *******/
  static int graphlet_counter;
  /*************************/


  list< list<int> >   this_node_cond; // This node must be bigger than all these nodes
  list< list<iPair> > cond;           // List of symmetry breaking conditions
  list<GTrieNode *> child;            // List of child g-trie nodes

  bool cond_ok;                       // no need to check for conditions
  bool cond_this_ok;                  // no need to check for this node conditions

  /************* Pthreads ************/
  long id;
  /***********************************/
  int depth;          // Depth of g-trie node

  bool is_graph;       // Is this node the end of a subGraph?
  unsigned long long frequency;        // Frequency of this particular subGraph

  bool *in;           // Outward edges
  bool *out;          // Inward edges
  /******* Graphlets *******/
  char *full_graph;   // Full graph representation
  long long int *orbits;
  int graphlet_id;
  /*************************/
  int total_in;       // Number of inward  edges
  int total_out;      // Number of outward edges
  int total_edges;    // Number of inward + outward edges

  int nconn;          // Number of connected nodes
  int *conn;          // Connected nodes


  GTrieNode(int d);   // Create g-trie node with depth 'd'
  ~GTrieNode();

  void insert(Graph *g);
  void insertCond(Graph *g, list<iPair> *cond);
  void showAsText(FILE *f);

  int countNodes();
  int countGraphs();
  int countGraphPaths();
  int maxDepth();

  void zeroFrequency();
  void showFrequency();

  float frequencyGraph(Graph *g);

  /************* Pthreads ************/
  void goCondDir(void* data);
  void goCondUndir(void* data);
  void goCondDirGraphlets(void* data);
  void goCondUndirGraphlets(void* data);
  void calcFrequency();
  /***********************************/
  void goCondSample();

  void insertConditionsFiltered(list<iPair> *cond);

  long long countOccurrences();
  int countGraphsApp();

  void populateGraphTree(GraphTree *tree, char *s, int size);
  void populateMap(mapStringLL *m, char *s, int size);

  void writeToFile(FILE *f);
  void readFromFile(FILE *f);

  void cleanConditions();
  void clean(int a, int b);
  void makeConditionsArray();

  /******* Graphlets *******/
  vector<GTrieNode *> computeOrbits(vector<GTrieNode *> nodes);
  void printOrbits();
  int getGraphletData(int pos, const char* s);
  /*************************/
};

class GTrie {
 private:
  GTrieNode *_root;

 public:
  GTrie();
  ~GTrie();

  /******* Graphlets *******/
  static int total_orbits;
  static long long int** orbit_frequency;
  /*************************/

  void insertGraph(Graph *g);
  void insertGraphCond(Graph *g, list<iPair> *cond);
  void insertGraphString(int size, const char* s);
  void insertGraphNautyString(int size, const char* s, bool dir, int label);

  float frequencyGraphString(int size, const char *s);

  void showAsText(FILE *f);
  double compressionRate();
  int maxDepth();

  void census(Graph *g);
  void censusSample(Graph *g, double *p);

  void showFrequency();

  void unlist();

  int countGraphsApp();
  int countGraphs();
  double countOccurrences();

  void cleanConditions();

  void writeToFile(char *s);
  void readFromFile(char *s);
  void readSubgraphs(int size, bool dir, char *s);

  void populateGraphTree(GraphTree *tree, int size);
  void populateMap(mapStringLL *m, int size);

  /******* Graphlets *******/
  void appendOrbits();
  void buildNodeList(vector<GTrieNode*> nodes);
  void printOrbits();
  void freeOrbitFrequencies(int num_nodes);
  void writeGraphletCountToFile(int numNodes, FILE* out);
  /*************************/

  /************* Pthreads ************/
  static void* startCensus(void *data);
  static void  censusOnNodes(int begin, int end, void *data);
  static void  censusOnNodesGraphlets(int begin, int end, void *data);

  static void askForWork(void *data);
  static int  getThreadToCall(int tid);
 
  static void createWorkFromInterval(void *data, int current, int last, GTrieNode* gtrie_location);
  static void createWorkFromArray(void *data, int index, int* nodes_to_do, GTrieNode* gtrie_location, int mylim);

  static void giveWork(void* data);
  static queue<GTrieNode *> correctOrder(queue<GTrieNode *> old_order);
  static void goBackToWork(void* data);
  static void goBackToWorkGraphlets(void* data);
  /***********************************/
};


#endif

