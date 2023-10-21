/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
