/**
 * @file cycle_thread.cc
 * Thread performin a cycle
 */

#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <signal.h>
#include "../stdheader.h"
#include "../util/debug.h"
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
  // Initialize all the barriers
  pthread_barrier_init(&barrier_start, NULL, 2);
  pthread_barrier_init(&barrier_end, NULL, 2);
  // The next cycle will be the first
  DEBUG(15, "CycleThread: Initialized CycleThread variables")
}

/**
 * Called by the derrived class to terminate the thread
 */
void util::CycleThread::terminate_cycle() {
  // Terminate while loop in run()
  keeprunning = false;
  // Make sure it terminates by doing one more cycle
  do_cycle();
  // wait for thread to finish
  wait();
  pthread_barrier_destroy(&barrier_start);
  pthread_barrier_destroy(&barrier_end);
  DEBUG(15, "CycleThread: Destroyed all CycleThread variables")
}

/**
 * Where the cycle (while loop) is placed
 */
void util::CycleThread::run() {
  init_run();
  while (keeprunning) {
    // wait till told to start
    pthread_barrier_wait(&barrier_start);
    cycle();
    // wait till finished
    pthread_barrier_wait(&barrier_end);
  }
  end_run();
}

/**
 * Initiate a cycle
 */
void util::CycleThread::do_cycle() {
  // start the cycle 
  pthread_barrier_wait(&barrier_start);
  // wait for the cycle to finish
  pthread_barrier_wait(&barrier_end);
}
