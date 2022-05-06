/**
 * @file qm_worker_test.h
 * Base class to perform tests on QM workers
 */
#ifndef INCLUDED_QM_WORKER_T_H
#define	INCLUDED_QM_WORKER_T_H

#include <gtest/gtest.h>

#include "../stdheader.h"

#include "test_simulation.h"
#include "test_parameters.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_worker.h"
#include "../interaction/qmmm/qm_zone.h"

#include "../configuration/configuration.h"
#include "../topology/topology.h"

#include "../math/periodicity.h"

namespace testing {

class QM_Worker_Test : public ::testing::Test {

public:

  /**
   * @brief Construct a new qm worker test object
   * 
   * @param parameter 
   * @param results 
   */
  QM_Worker_Test(const Parameter& parameter);

  /**
   * @brief Destroy the qm worker test object
   * 
   */
  virtual ~QM_Worker_Test() = default;

protected:

  /**
   * Called before the execution of each test
   * 
   */
  virtual void SetUp() override = 0;

  /**
   * Called after the execution of each test
   * 
   */
  virtual void TearDown() override = 0;

 /**
   * @brief Initializes the parameters
   * 
   */
  virtual void init_parameters() = 0;

  /**
   * @brief Initializes the expected results
   * 
   */
  virtual void init_results() = 0;

  /**
   * @brief Initializes the expected results for energies calculated
   * 
   */
  virtual void init_results_energies() = 0;

  /**
   * @brief Initializes the expected results for forces calculated
   * 
   */
  virtual void init_results_forces() = 0;

  /**
   * @brief Initializes the expected results for forces calculated
   * 
   */
  virtual void init_results_velocities() = 0;

  /**
   * @brief Initializes the expected results for forces calculated
   * 
   */
  virtual void init_results_positions() = 0;

  virtual void init_results_baths() = 0;

  virtual void init_results_bonded_terms() = 0;

  virtual void init_results_nonbonded_terms() = 0;

  virtual void init_results_special_terms() = 0;

  virtual void init_results_mass() = 0;

  virtual void init_results_temperature() = 0;

  virtual void init_results_volume() = 0;

  virtual void init_results_pressure() = 0;

  /**
   * @brief Initializes the expected results for parameters read from the input file
   * 
   */
  virtual void init_results_parameters() = 0;

  /**
   * @brief Initializes the expected results for units and conversion factors
   * 
   */
  virtual void init_results_units() = 0;

  /**
   * @brief Initializes the expected results for element and iac atom mapping
   * 
   */
  virtual void init_results_elements() = 0;

  /**
   * @brief Initializes the expected results for the initial QM zone
   * 
   */
  virtual void init_results_qm_zone_init() = 0;

  /**
   * @brief Checks if Gromos has read in the input parameters correctly
   * 
   */
  virtual void check_parameter_init() = 0;

  /**
   * @brief Checks if Gromos has calculated energies, forces, coordinates correctly
   * 
   */
  virtual void check_simulation_results();

  /**
   * @brief Checks if Gromos has calculated energies correcly
   * 
   */
  virtual void check_simulation_results_energies();

  /**
   * @brief Checks if Gromos has calculated forces correcly
   * 
   */
  virtual void check_simulation_results_forces();

  /**
   * @brief Checks if Gromos has calculated velocities correcly
   * 
   */
  virtual void check_simulation_results_velocities();

  /**
   * @brief Checks if Gromos has calculated positions correcly
   * 
   */
  virtual void check_simulation_results_positions();

  virtual void check_simulation_results_temperature_baths();

  virtual void check_simulation_results_bonded_terms();

  virtual void check_simulation_results_nonbonded_terms();

  virtual void check_simulation_results_special_terms();

  virtual void check_simulation_results_mass();

  virtual void check_simulation_results_temperature();

  virtual void check_simulation_results_volume();

  virtual void check_simulation_results_pressure();

  /**
   * @brief Checks if the units conversion factors have been read in correctly
   * 
   */
  virtual void check_unit_conversion_factors() = 0;

  /**
   * @brief Checks if the qm/mm parameter of the simulation input file have been read in correcly
   * 
   */
  virtual void check_qmmm_parameter();

  /**
   * @brief Checks if the parameters for the QM zone have been read in correcly
   * 
   */
  virtual void check_qm_zone_param();

  /**
   * @brief Checks if the pointers to the QM interaction object, QM worker, and QM zone have been initialized
   * 
   */
  virtual void check_qm_interaction_ptr();

  /**
   * @brief Gathers coordinates of atoms according to out_configuration.cc. Should be implemented elsewhere, probably out_configuration.cc
   * 
   * @tparam b 
   * @param conf Configuration of the simulation
   * @param topo Topology of the simulation
   * @param pos_init Positions as they come out of Gromos
   * @param pos_processed Positions as they would be written out to the .trc file
   */
  template <math::boundary_enum b>
  void preprocess_positions(configuration::Configuration& conf
    , topology::Topology& topo
    , math::VArray& pos_init
    , math::VArray& pos_processed) {

    math::Periodicity<b> periodicity(conf.current().box);
    math::Vec v, v_box, trans, r;

    // put charge groups into the box (on the fly)
    topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
      cg_to = topo.chargegroup_end();

    // solute charge groups
    unsigned int i = 0, count = 0;
    for ( ; i < topo.num_chargegroups(); ++cg_it, ++i) {
      // gather on first atom
      v = pos_init(*cg_it.begin());
      v_box = v;
      periodicity.put_into_positive_box(v_box);
      trans = v_box - v;

      // atoms in a charge group
      topology::Atom_Iterator at_it = cg_it.begin(),
        at_to = cg_it.end();

      for (; at_it != at_to; ++at_it, ++count) {
        r = pos_init(*at_it) + trans;
        pos_processed.push_back(r);
      }
    }

    // solvent charge groups
    unsigned int s = 0;
    unsigned int mol = 0;

    for (; cg_it != cg_to; ++cg_it, ++mol) {
      v = pos_init(**cg_it);
      v_box = v;
      periodicity.put_into_positive_box(v_box);
      trans = v_box - v;

      if (mol >= topo.num_solvent_molecules(s))++s;

      // loop over the atoms
      topology::Atom_Iterator at_it = cg_it.begin(),
        at_to = cg_it.end();

      for (unsigned int atom = 0; at_it != at_to; ++at_it, ++atom) {
        r = pos_processed(*at_it) + trans;
        pos_processed.push_back(r);
      }
    }
  }

  /**
   * Instance of a test simulation based on the input files
   * 
   */
  Test_Simulation test_sim_;

  /**
   * Structure that stores the expected test results
   * 
   */
  Results results_;

  /**
   * Pointer to the qmmm_interaction object of the test simulation
   * 
   */
  interaction::QMMM_Interaction* qmmm_interaction_ptr;

  /**
   * Pointer to the qm_worker object of the test simulation
   * 
   */
  interaction::QM_Worker* qm_worker_ptr;

  /**
   * Pointer to the qm_zone object of the test simulation
   * 
   */
  interaction::QM_Zone* qm_zone_ptr;

  /**
   * Absolute error tolerance for comparing floating point numbers
   * 
   */
  const double epsilon_ = 0.005;

};

}

#endif	/* QM_WORKER_T_H */
