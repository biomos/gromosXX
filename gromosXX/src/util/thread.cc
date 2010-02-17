/**
 * @file thread.cc
 * basic thread class (uses pthread)
 */

#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <signal.h>

#include "thread.h"

extern "C" void * util::launchThread(void * t) {
  Thread * thread = (Thread*)(t);
  thread->run();
  return NULL;
}

void util::Thread::start() {
  pthread_create(&id, NULL, launchThread, (void*) this);
}

void util::Thread::kill() {
  pthread_kill(id, SIGTERM);
}

void util::Thread::wait() {
  pthread_join(id, NULL);
};