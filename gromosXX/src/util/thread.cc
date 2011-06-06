/**
 * @file thread.cc
 * basic thread class (uses pthread)
 */

#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <signal.h>
#include "../../config.h"
#include "thread.h"

extern "C" void * util::launchThread(void * t) {
  Thread * thread = (Thread*)(t);
  thread->run();
  return NULL;
}

void util::Thread::start() {
#ifdef OMP
  pthread_create(&id, NULL, launchThread, (void*) this);
#endif
}

void util::Thread::kill() {
#ifdef OMP
  pthread_kill(id, SIGTERM);
#endif
}

void util::Thread::wait() {
#ifdef OMP
  pthread_join(id, NULL);
#endif
};
