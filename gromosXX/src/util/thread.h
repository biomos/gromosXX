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

