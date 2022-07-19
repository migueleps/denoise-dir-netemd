/* -------------------------------------------------
David Aparício - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
Common definitions

Last Update: 26/06/2016
----------------------------------------------------
*/

#ifndef _SHAREDPARALLEL_
#define _SHAREDPARALLEL_

#include <pthread.h>

class SharedParallel {
 public:
  static bool no_more_work;
  static int  num_nodes; 
  static int  waiting_threads;

  static pthread_mutex_t answer_mutex;
  static pthread_cond_t  answer_cond;

  static pthread_mutex_t call_mutex;
  static pthread_mutex_t waiting_mutex;
};

#endif
