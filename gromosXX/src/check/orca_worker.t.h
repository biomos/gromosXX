/**
 * @file orca_worker.t.h
 * Class to perform tests on the Orca QM worker
 */
#ifndef INCLUDED_ORCA_WORKER_T_H
#define	INCLUDED_ORCA_WORKER_T_H

#include "qm_worker.t.h"

class Orca_Worker_Test : public QM_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  Orca_Worker_Test();


  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~Orca_Worker_Test();

};

#endif	/* ORCA_WORKER_T_H */
