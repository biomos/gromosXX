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
 * @file thread.h
 * basic thread class (uses pthread)
 */

#ifndef _THREADS_H
#define	_THREADS_H

namespace util {

  /**
   * @class Thread
   * basic thread class
   */
  class Thread {
  public:
    /**
     * start the thread
     */
    void start();

    /**
     * in this function the genuine purpose of every child class
     * can be implemented
     */
    virtual void run() = 0;

    /**
     * kill the thread
     */
    void kill();

    /**
     * interrupt the thread
     */
    void wait();

    /**
     * destructor
     */
    virtual ~Thread() {}

  private:
    /**
     * thread id number
     */
    pthread_t id;
  };

  /**
   * Used to start the thread
   */
  extern "C" void * launchThread(void * t);

}
#endif	/* _THREADS_H */

