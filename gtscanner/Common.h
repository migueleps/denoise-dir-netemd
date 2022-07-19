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
Common definitions

Last Update: 04/06/2014
----------------------------------------------------
*/

#ifndef _COMMON_
#define _COMMON_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <list>
#include <map>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>

#include <vector>
#include <queue>
#include <unordered_map>

#include <random>

#define MAX_BUF 1024 // Maximum string buffer size

#define VERSION      "0.1"
#define PROGRAM_NAME "gtrieScanner"

#define SEPARATOR    "------------------------------------------"

// Limits for subgraph size
#define MIN_MOTIF_SIZE  3
#define MAX_MOTIF_SIZE 50
#define MAX_LABEL      200

#define ERROR_HEADER "Error: "

#define INVALID -1

#define EPSILON 0.00000001

#define BIT_SET(n,i)   ((n)|=(1<<(i)))
#define BIT_CLEAR(n,i) ((n)&=~(1<<(i)))
#define BIT_VALUE(n,i) (((n)>>(i))&1)

#define smallNode char

#define INVALID_FILE   "__NULL__"    // Invalid file name string

#define DEFAULT_RESULTS "results"    // Default name for results file
#define DEFAULT_OCC "occ"            // Default name for occurrences file

typedef enum {NOMETHOD, ESU, GTRIE, SUBGRAPHS, FASE} MethodType;
typedef enum {NOOUTPUT, TEXT, HTML}                  OutputType;
typedef enum {NOFORMAT, SIMPLE, SIMPLE_WEIGHT}       FormatType;

using namespace std; // Could be avoided if wanted

typedef std::pair<int,int> iPair;
typedef std::vector<int> VInt;
typedef std::vector<VInt> VVInt;
typedef vector<smallNode *> VVsmallNode;
typedef map< string, int> mapStringInt;
typedef map< string, long long> mapStringLL;

typedef struct {
  char *s;
  int f_original;
  double avg_random;
  double dev_random;
  double z_score;
} ResultType;


// Class for "global" variables
class Global {
 public:
  static bool show_occ;  // Show occurrences?
  static FILE *occ_file; // FILE handler for dumping occurrences;
  /************* Pthreads ************/
  static int num_threads;   // Number of threads
  static FILE **occ_file_t; // FILE handler for dumping occurrences for each thread
  /***********************************/
  /************ Graphlets ************/
  static bool compute_orbits;
  /***********************************/
};

#endif
