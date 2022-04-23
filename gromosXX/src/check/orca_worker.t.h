/**
 * @file orca_worker.t.h
 * Class to perform tests on the Orca QM worker
 */
#ifndef INCLUDED_ORCA_WORKER_T_H
#define	INCLUDED_ORCA_WORKER_T_H

#include "qm_worker.t.h"

namespace test {

class Orca_Worker_Test : public QM_Worker_Test {

public:

  /**
   * Construct a new qm worker test object
   * 
   */
  Orca_Worker_Test(std::string binary_name
               , std::string test_title
               , std::string topology_file
               , std::string simulation_file
               , std::string configuration_file
               , std::string qmmm_file
               , std::string coordinate_trajectory_file
               , std::string energy_trajectory_file
               , std::string final_configuration_file);


  /**
   * Destroy the qm worker test object
   * 
   */
  virtual ~Orca_Worker_Test() = default;

};

}

#endif	/* ORCA_WORKER_T_H */
