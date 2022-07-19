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
G-Trie implementation

---------------------------------------------------- */

#include "FaseGTrie.h"
#include "GraphTree.h"

class GTrlabelTrie;
class childTrie;
class trie;

class labelTrie
{
  public:
    labelTrie* childs[2];
    labelTrie* parent;
    int num;
    char let;
};

class trie {
  public:
    trie* parent;
    childTrie* childs;
    int leaf;
    long long int count;
    bool canon;
    trie() {}
    trie(childTrie* childs1, trie* parent1, int leaf1, bool canon1) { 
      parent = parent1;
      childs = childs1;
      leaf = leaf1;
      canon = canon1;
    }
};

class childTrie {
  public:
    childTrie* list[MAX_MOTIF_SIZE + 1];
    trie* ref;
};

childTrie* FaseGTrie::initChild()
{
  childTrie* nn = new childTrie();
  return nn;
}

void FaseGTrie::init(int threads1)
{
  currentid = 0;
  threads = threads1;
  nodes = 0;
  count = 0;
  numLabel = 1;
  maxLabel = 1000;
  labelList = (labelTrie **)malloc(maxLabel * sizeof(labelTrie *));
  canCount = (long long int**)malloc(threads * sizeof(long long int*));
  for(int i = 0; i < threads; i++){
    canCount[i] = (long long int *)malloc(maxLabel * sizeof(long long int));
    for(int j = 0; j < maxLabel; j++){
      canCount[i][j] = 0;
    }
  }
  /*for(int i = 0; i < threads; i++){
    for(int j = 0; j < maxLabel; j++){
      printf("%d ", canCount[i][j]);
    }
  }*/
  labelRoot = new labelTrie();
  labelRoot->parent = NULL;
  //current_th = new trie*[threads];
  top = new trie();
  top->parent = NULL;
  top->childs = new childTrie();
  top->count = 0;
  top->canon = false;
  /*for(int i = 0; i < threads; i++)
    (current_th[i]) = top;*/
  currentid++;
}

void FaseGTrie::destroy()
{
  free(labelList);
  for(int i = 0; i < threads; i++)
    free(canCount[i]);
  free(canCount);
}

void FaseGTrie::augment()
{
  maxLabel *= 2;
  //printf("maxlabel: %d\n", maxLabel);
  labelList = (labelTrie **)realloc(labelList, maxLabel * sizeof(labelTrie *));
  for(int i = 0; i < threads; i++){
    canCount[i] = (long long int *)realloc(canCount[i], maxLabel * sizeof(long long int));
    for(int j = numLabel; j < maxLabel; j++)
      canCount[i][j] = 0; // take out trash values 
  }
}

long long int FaseGTrie::getCanonicalNumber()
{
  return numLabel - 1;
}

long long int FaseGTrie::getNodeNumber()
{
  return nodes;
}

long long int FaseGTrie::getClassNumber()
{
  return count;
}

/************* PThreads ***********/
childTrie* FaseGTrie::searchChild1(childTrie* node, char* s)
{
  while (*s) {
    if (node->list[*s] == NULL) return NULL;
    node = node->list[*s];
    s++;
  }
  return node;
}

trie* FaseGTrie::checkInsert_th1(char* s, trie* pointer, int tid)
{
  childTrie* temp = searchChild1(pointer->childs, s);
  if (temp != NULL)
  {
    //temp->ref->count++;
    //int l = temp->ref->leaf;
    //if(temp->ref != NULL) canCount[tid][temp->ref->leaf]++;
    //pointer = temp->ref;
    return temp->ref;
  }
  else{
    return NULL;
  }
}

childTrie* FaseGTrie::searchChild(childTrie* node, char* s)
{
  while (*s) {
    if (node->list[*s] == NULL) node->list[*s] = initChild();
    node = node->list[*s];
    s++;
  }
  return node;
}

trie* FaseGTrie::insert_th(trie* pointer, char* s, int depth)
{
  //nodes++;
  childTrie* temp = searchChild(pointer->childs, s);
  childTrie* childs = new childTrie();
  
  temp->ref = new trie(childs, pointer, 0, false);
  /*
  temp->ref = new trie();
  temp->ref->childs = childs;
  temp->ref->leaf = 0;
  temp->ref->parent = pointer;
  temp->ref->canon = false;*/
  //temp->ref->count = 1;
  //temp->ref->id = ++currentid;
  //printf("%p -> Count: %d\n", temp->ref, nodes);
  return temp->ref;
}
/***********************************/

int FaseGTrie::searchLabel_th(labelTrie* node, char* s, int tid)
{
  if (*s == 0)
  {
    if (node->num == 0)
    {
      node->num = numLabel;
      //canCount[tid][numLabel] = 1;
      labelList[numLabel++] = node;
      if (numLabel == maxLabel - 2)
        augment();
      //for(int i = 0; i < threads; i++) //some threads may not reach node -> trash values in array
        //for(int j = numLabel; j < maxLabel; j++) // i=0?
          //canCount[i][j] = 0;
    }
    //else
      //canCount[tid][node->num]++;

    return node->num;
  }
  else
  {
    if (node->childs[*s - '0'] == NULL)
    {
      node->childs[*s - '0'] = new labelTrie();
      node->childs[*s - '0']->num = 0;
      node->childs[*s - '0']->let = *s;
      node->childs[*s - '0']->parent = node;
    }
    return searchLabel_th(node->childs[*s - '0'], s + 1, tid);
  }
}

/************* Pthreads **************/
trie* FaseGTrie::setCanonicalLabel_th(char* s, trie* pointer, int tid)
{
  if(pointer->canon == false){
    count++;
    //printf("Tid(%d): Pointer(%p) -> count: %d\n",tid, pointer, count);
    int vl = searchLabel_th(labelRoot, s, tid);
    pointer->leaf = vl;
    pointer->canon = true;
  }
  return pointer;
}

int FaseGTrie::getCurrentLeaf_th(trie* pointer)
{
  return pointer->leaf;
}
/************************************/

int FaseGTrie::getCurrentLeaf()
{
  return current->leaf;
}

void FaseGTrie::jump()
{
  current = current->parent;
}

/********** PThreads ********/
trie* FaseGTrie::jump_th(trie* pointer)
{
  return pointer->parent;
}
/***************************/

char* FaseGTrie::getLabel(int key)
{
  int size = 0;
  labelTrie* dp = labelList[key];
  while (dp->parent != NULL)
  {
    globS[size++] = dp->let;
    dp = dp->parent;
  }
  globS[size] = 0;
  return globS;
}

void FaseGTrie::printGtrie(childTrie *node, char* label, int sz)
{
  if (node->ref != NULL)
  {
    if (sz != 0)
      label[sz++] = '+';
    label[sz] = '\0';
//      printf("Class: %s, Found: %lld, Sample Model: %s\n", label, node->ref->count, getLabel(node->ref->leaf));
    if (!node->ref->leaf)
      printGtrie(node->ref->childs, label, sz);
    else
    {
      label[sz - 1] = '\0';
      fprintf(f, "%s\n", label);
    }
    if (sz != 0)
      label[--sz] = '\0';
  }
  int i;
  for (i = 1; i < MAX_MOTIF_SIZE + 1; i++)
  {
    label[sz] = i + '0' - 1;
    label[sz + 1] = '\0';
    if (node->list[i] != NULL)
      printGtrie(node->list[i], label, sz + 1);
  }
}

void FaseGTrie::listGtrie(FILE* out)
{
  /*f = out;
  fprintf(f, "=\n");
  while (current_th[0]->parent != NULL)
    GTrie::jump();
  char tmpStr[MAX_LABEL];
  tmpStr[0] = '\0';
  GTrie::printGtrie(current_th[0]->childs, tmpStr, 0);
  fprintf(f, "done\n");*/
}


void FaseGTrie::printClass(labelTrie *node, char* label, int sz, double prob)
{
  if (node->num != 0)
  {
    label[sz] = '\0';
    
    long long int totalcount = 0;
    for(int i = 0; i < threads; i++){
      totalcount+= canCount[i][node->num];
    }
    printf("Canonical Class: %s - Occurrences: %lld\n", label, (long long int)(totalcount));
    return;
  }
  if(node->childs[0] != NULL)
  {
    label[sz] = '0';
    printClass(node->childs[0], label, sz + 1, prob);
  }
  if(node->childs[1] != NULL)
  {
    label[sz] = '1';
    printClass(node->childs[1], label, sz + 1, prob);
  }
}
 
void FaseGTrie::listClasses(FILE* out, double prob)
{
  f = out;
  char tmpStr[MAX_LABEL];
  tmpStr[0] = '\0';
  FaseGTrie::printClass(labelRoot, tmpStr, 0, prob);
}

void FaseGTrie::incCount(int tid, trie* trieLeaf){
  //while(trieLeaf->leaf == NULL);
  //int x = trieLeaf->leaf;
  canCount[tid][trieLeaf->leaf]++;
}

void FaseGTrie::populateGraphTree(GraphTree* sg, Isomorphism* mynauty) {
  char tmpStr[MAX_LABEL];
  tmpStr[0] = '\0';
  FaseGTrie::populateGraphTreeRecursive(sg, mynauty, labelRoot, tmpStr, 0);
}

void FaseGTrie::populateGraphTreeRecursive(GraphTree* sg, Isomorphism* mynauty, labelTrie *node, char* label, int sz){
  if (node->num != 0)
  {
    label[sz] = '\0';
    
    long long int totalcount = 0;
    for(int i = 0; i < threads; i++){
      totalcount+= canCount[i][node->num];
    }
	//char canonicalLabel[MAX_LABEL];
    //mynauty->canonicalStrNautyMatrix(label, canonicalLabel);
    
    for(int c = 0; c < (long long int)(totalcount); c++) sg->incrementString(label);
    return;
  }
  if(node->childs[0] != NULL)
  {
    label[sz] = '0';
    populateGraphTreeRecursive(sg, mynauty, node->childs[0], label, sz + 1);
  }
  if(node->childs[1] != NULL)
  {
    label[sz] = '1';
    populateGraphTreeRecursive(sg, mynauty, node->childs[1], label, sz + 1);
  }
}

/*long long int GTrie::getCurrentId(int tid){
  return current_th[tid]->id;
}*/
