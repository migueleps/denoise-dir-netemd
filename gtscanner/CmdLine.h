/* -------------------------------------------------

  / _ \/__   \    / _\ ___ __ _ _ __  _ __   ___ _ __ 
 / /_\/  / /\/____\ \ / __/ _` | '_ \| '_ \ / _ \ '__|
/ /_\\  / / |_____|\ \ (_| (_| | | | | | | |  __/ |   
\____/  \/        \__/\___\__,_|_| |_|_| |_|\___|_|   

gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

Pedro Ribeiro & David Aparício - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
Command Line Client functions

Last Update: 13/07/2016
---------------------------------------------------- */

#ifndef _CMDLINE_
#define _CMDLINE_

#include "GraphMatrix.h"
#include "GraphUtils.h"
#include "Error.h"
#include "Common.h"
#include "Esu.h"
#include "Fase.h"
#include "GTrie.h"
#include "Isomorphism.h"
#include "Timer.h"
#include "Random.h"
/************* Pthreads ************/
#include "SharedParallel.h"
/***********************************/
/************ Graphlets ************/
#include "GDA.h"
/***********************************/

typedef std::vector< std::map<int, float> > GDD;

class CmdLine {
 private:
  static char graph_file[MAX_BUF];
  static char gtrie_file[MAX_BUF];
  static char subgraphs_file[MAX_BUF];
  static char output_file[MAX_BUF];
  static char occ_file[MAX_BUF];
  /************* Graphlets ************/
  static char orbits_file[MAX_BUF];
  static char orbits_dir[MAX_BUF];
  /************************************/

  static bool dir;
  static bool occurrences;
  static bool create;

  static int motif_size;
  static int random_number;
  static int random_seed;
  static int random_exchanges;
  static int random_tries;

  static double time_original;
  static double time_random;

  static MethodType method;
  static FormatType format;
  static OutputType output;
  static Graph *g;

  static FILE *f_output;
  static FILE *f_occ;
  /************* Pthreads ************/
  static FILE **f_occ_t;
  /***********************************/
  /************ Graphlets ************/
  static FILE *f_orbits;

  static bool compute_orbits;
  static bool compute_gda;
  static GDD gdd;
  /***********************************/

  static GTrie *gt_original;
  static GTrie *gt;

  static GraphTree sg_original;

  static time_t t_start;
  
  static void about();
  static void defaults();
  static void parse_cmdargs(int argc, char **argv);
  static void run_esu(Graph *g, GraphTree *sg);
  static void run_gtrie(Graph *g, GraphTree *sg);
  static void run_fase(Graph *g, int motif_size, GraphTree *sg);
  static void run_subgraphs(Graph *g, GraphTree *sg);

  /************ Graphlets ************/
  static GDA gda;

  static void computeGDD(long long int** orbit_frequencies, int num_nodes, int total_orbits, FILE* out);
  static void compute_gda_dir();
  /***********************************/

  static MethodType str_to_method(char *s);
  static FormatType str_to_format(char *s);
  static OutputType str_to_output(char *s);

  static int compare_results(const void *a, const void *b);

  static void prepare_graph();
  static void prepare_files();
  static void compute_original();
  static void compute_results();
  static void show_results(ResultType *res, int nres);

  static void create_gtrie();

 public:
  static void init(int argc, char **argv);
  static void finish();
  static void decide_action();  
};

#endif
