/**
 * @file orca_worker_test.h
 * Class to perform tests on the Orca QM worker
 */
#ifndef INCLUDED_ORCA_WORKER_T_H
#define	INCLUDED_ORCA_WORKER_T_H

#include "qm_worker_test.h"

namespace testing {

/**
 * @brief A class that tests the Orca_Worker using an electrostatic embedding
 * 
 */
class Orca_Worker_Electrostatic_Test : public QM_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  Orca_Worker_Electrostatic_Test();


  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~Orca_Worker_Electrostatic_Test() = default;

protected:

  /**
   * Called before the execution of each test
   * 
   */
  void SetUp() override;

  /**
   * Called after the execution of each test
   * 
   */
  void TearDown() override;

};

/**
 * @brief A class that tests the Orca_Worker using an mechanical embedding
 * 
 */
class Orca_Worker_Mechanical_Test : public QM_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  Orca_Worker_Mechanical_Test();


  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~Orca_Worker_Mechanical_Test() = default;

protected:

  /**
   * Called before the execution of each test
   * 
   */
  void SetUp() override;

  /**
   * Called after the execution of each test
   * 
   */
  void TearDown() override;

};

}

#endif	/* ORCA_WORKER_T_H */
