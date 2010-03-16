/**
 * @file cycle_thread.cc
 * Thread performin a cycle
 */

#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <signal.h>
#include <stdheader.h>
#include <util/debug.h>
#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

#include "cycle_thread.h"


/**
 * Consttuctor
 */
util::CycleThread::CycleThread() {

  // For the while loop in run()
  keeprunning = true;
  // Initialize all the mutex and conditons
  pthread_mutex_init(&mutex_cycle, NULL);
  pthread_mutex_init(&mutex_calc, NULL);
  pthread_cond_init(&cond_cycle, NULL);
  pthread_cond_init(&cond_calc, NULL);
  // Needed to avoid spurious wakeup
  do_a_calc = false;
  do_a_cycle = false;
  // The next cycle will be the first
  first = true;
  DEBUG(15, "CycleThread: Initialized CycleThread variables")
}

/**
 * Called by the derrived class to terminate the thread
 */
void util::CycleThread::terminate_cycle() {
  // Terminate while loop in run()
  keeprunning = false;
  // Make sure it terminates by doing one more cycle
  do_a_calc = true;
  pthread_cond_signal(&cond_calc);
  // wait for thread to finish
  wait();
  pthread_mutex_destroy(&mutex_cycle);
  pthread_mutex_destroy(&mutex_calc);
  pthread_cond_destroy(&cond_cycle);
  pthread_cond_destroy(&cond_calc);
  //pthread_cond_destroy(&cond_wait);
  DEBUG(15, "CycleThread: Destroyed all CycleThread variables")
}

/**
 * Where the cycle (while loop) is placed
 */
void util::CycleThread::run() {

  init_run();
  while (keeprunning) {

    cycle();
    DEBUG(15, "CycleThread: Will send cycle signal")
    pthread_mutex_lock(&mutex_calc);
    DEBUG(15, "CycleThread: Waiting for the signal to calculate")
    do_a_cycle = true;
    pthread_cond_signal(&cond_cycle);
    do_a_calc = false;
    while(!do_a_calc){    // Needed to avoid "Spirious wakeup"
      pthread_cond_wait(&cond_calc, &mutex_calc);
    }
    pthread_mutex_unlock(&mutex_calc);
  }
  do_a_cycle = true;
  pthread_cond_signal(&cond_cycle);
  end_run();
}

/**
 * Initiate a cycle
 */
void util::CycleThread::do_cycle() {
  // Is it the first cycle?
  if (first) {
    pthread_mutex_lock(&mutex_cycle);
    start();
    DEBUG(15, "CycleThread: Just started, waiting for cycle signal");
    do_a_cycle = false;
    while (!do_a_cycle) { // Needed to avoid "Spirious wakeup"
      pthread_cond_wait(&cond_cycle, &mutex_cycle);
    }
    pthread_mutex_unlock(&mutex_cycle);
    first = false;
  } else {
    pthread_mutex_lock(&mutex_cycle);
    do_a_cycle = false;
    do_a_calc = true;
    pthread_cond_signal(&cond_calc);
    DEBUG(15, "CycleThread: Sent calculation signal, waiting for cycle signal")
    while (!do_a_cycle) { // Needed to avoid "Spirious wakeup"
      pthread_cond_wait(&cond_cycle, &mutex_cycle);
    }
    pthread_mutex_unlock(&mutex_cycle);
  }

}