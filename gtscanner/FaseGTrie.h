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

Pedro {Paredes, Ribeiro} - DCC/FCUP

----------------------------------------------------
G-Tries implementation

---------------------------------------------------- */

#ifndef _FASEGTRIE_
#define _FASEGTRIE_

#include "Common.h"
#include "Isomorphism.h"

class trie;
class childTrie;
class labelTrie;
class GraphTree; // forward declaration

/*! This class creates a G-Trie to store the labelings */
class FaseGTrie
{
 public:

  FILE* f;
  labelTrie** labelList;
  labelTrie* labelRoot;
  int numLabel;
  int maxLabel;
  char globS[MAX_LABEL];
  void init();
  /********** PThreads ********/
  int currentid;
  int threads;
  void init(int threads);
  //trie** current_th;
  trie*  top;
  long long int** canCount;
  /****************************/
  void destroy();
  childTrie* initChild();
  childTrie* searchChild(childTrie* node, char* s);
  childTrie* searchChild1(childTrie* node, char* s);
  char* getLabel(int key);
  void insert(char* s);
  void setCanonicalLabel(char *s);
  /********** PThreads ********/
  int searchLabel_th(labelTrie* node, char* s, int tid);
  trie* jump_th(trie* pointer);
  trie* setCanonicalLabel_th(char *s, trie* pointer, int tid);
  trie* checkInsert_th1(char* s, trie* pointer, int tid);
  trie* insert_th(trie* pointer, char* s, int depth);
  int getCurrentLeaf_th(trie* pointer);
  void incCount(int tid, trie* trieLeaf);
  void populateGraphTree(GraphTree* sg, Isomorphism* mynauty);
  void populateGraphTreeRecursive(GraphTree* sg, Isomorphism* mynauty, labelTrie *node, char* label, int sz);
  //long long int getCurrentId(int tid);
  /****************************/
  void jump();
  void augment();
  long long int getCanonicalNumber();
  long long int getClassNumber();
  long long int getNodeNumber();
  void listClasses(FILE* out, double prob);
  void listGtrie(FILE* out);
  int getCurrentLeaf();

 private:
  long long int count;
  long long int nodes;
  trie* current;
  void printClass(labelTrie *node, char* label, int sz, double prob);
  void printGtrie(childTrie *node, char* label, int sz);
};

#endif
