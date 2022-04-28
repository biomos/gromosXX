/**
 * @file orca_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"
#include "../algorithm/algorithm.h"
#include "../algorithm/algorithm/algorithm_sequence.h"

#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../io/configuration/out_configuration.h"

#include "check.h"
#include "test_parameters.h"

#include "qm_worker_test.h"
#include "orca_worker_test.h"
#include "orca_worker_test_electrostatic.h"
#include "orca_worker_test_mechanical.h"

#include "../math/gmath.h"

#include "../interaction/qmmm/orca_worker.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

void Orca_Worker_Test::SetUp() {
  // initialize test simulation and get a pointer to qm/mm relevant objects
  int err = test_sim_.init_simulation();
  ASSERT_EQ(err, 0) << "Initialization of the simulation unsuccessful. Error code: " << err;
  qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  qm_worker_ptr = dynamic_cast<interaction::Orca_Worker*>(qmmm_interaction_ptr->m_worker);
  qm_zone_ptr = qmmm_interaction_ptr->m_qm_zone;
}

void Orca_Worker_Test::TearDown() {
  // set all pointers to nullptr again
  qmmm_interaction_ptr = nullptr;
  qm_worker_ptr = nullptr;
  qm_zone_ptr = nullptr;
}

void Orca_Worker_Test::init_results_binary_name() {
  // binary and worker name
  results_.strs_[Key::binary] = "/home/fpultar/opt/orca-5.0.3/orca"; // modify this
  results_.strs_[Key::qm_worker_name] = "Orca Worker";
}

void Orca_Worker_Test::init_results_units() {
  // units and conversion factors
  results_.doubles_[Key::unit_factor_length] = 0.1;
  results_.doubles_[Key::unit_factor_energy] = 2625.5;
  results_.doubles_[Key::unit_factor_force] = 49641;
  results_.doubles_[Key::unit_factor_charge] = 1.0;
}

void Orca_Worker_Test::init_results_elements() {
  // element symbol and iac atom mapping
  results_.strs_[Key::symbol_h] = "H";
  results_.strs_[Key::symbol_c] = "C";
  results_.strs_[Key::symbol_o] = "O";
  results_.strs_[Key::symbol_p] = "P";
  results_.strs_[Key::symbol_pd] = "Pd";

  results_.ints_[Key::element_1] = 1;
  results_.ints_[Key::element_6] = 6;
  results_.ints_[Key::element_8] = 8;
  results_.ints_[Key::element_15] = 15;
  results_.ints_[Key::element_20] = 20;
  results_.ints_[Key::element_26] = 26;
  results_.ints_[Key::element_29] = 29;
  results_.ints_[Key::element_46] = 46;
}

void Orca_Worker_Test::check_parameter_init() {
  check_qmmm_parameter();
  check_unit_conversion_factors();
  check_element_mapping();
  check_binary_name();
  check_files();
  check_qm_zone_param();
  check_qm_interaction_ptr();
  check_qm_atoms_init();
  check_mm_atoms_init();
}

void Orca_Worker_Test::check_unit_conversion_factors() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, double>& doubles_res = results_.doubles_;
  // units and conversion factors
  EXPECT_EQ(param.qmmm.orca.unit_factor_length, doubles_res[Key::unit_factor_length]);
  EXPECT_EQ(param.qmmm.orca.unit_factor_energy, doubles_res[Key::unit_factor_energy]); // tolerance (float)
  EXPECT_EQ(param.qmmm.orca.unit_factor_force, doubles_res[Key::unit_factor_force]);
  EXPECT_EQ(param.qmmm.orca.unit_factor_charge, doubles_res[Key::unit_factor_charge]);
}

void Orca_Worker_Test::check_element_mapping() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, int>& ints_res = results_.ints_;
  std::unordered_map<Key, std::string>& strs_res = results_.strs_;
    // some tests for elements and iac_elements map
  EXPECT_EQ(param.qmmm.orca.elements.at(1), strs_res[Key::symbol_h]);
  EXPECT_EQ(param.qmmm.orca.elements.at(6), strs_res[Key::symbol_c]);
  EXPECT_EQ(param.qmmm.orca.elements.at(8), strs_res[Key::symbol_o]);
  EXPECT_EQ(param.qmmm.orca.elements.at(15), strs_res[Key::symbol_p]);
  EXPECT_EQ(param.qmmm.orca.elements.at(46), strs_res[Key::symbol_pd]);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(0), ints_res[Key::element_8]); // check if (-1) is subtracted from input file IAC
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(23), ints_res[Key::element_29]);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(24), ints_res[Key::element_29]);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(25), ints_res[Key::element_26]);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(28), ints_res[Key::element_20]);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(69), ints_res[Key::element_15]);
}

void Orca_Worker_Test::check_binary_name() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, std::string>& strs_res = results_.strs_;
  // binary and input file names
  EXPECT_EQ(param.qmmm.orca.binary, strs_res[Key::binary]); // modify this
  // test the qm worker
  EXPECT_EQ(qm_worker_ptr->name(), strs_res[Key::qm_worker_name]);
}

void Orca_Worker_Test::check_files() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, std::string>& strs_res = results_.strs_;
  EXPECT_EQ(param.qmmm.orca.input_file, strs_res[Key::input_file]);
  EXPECT_EQ(param.qmmm.orca.input_coordinate_file, strs_res[Key::input_coordinate_file]);
  EXPECT_EQ(param.qmmm.orca.input_pointcharges_file, strs_res[Key::input_pointcharges_file]);
  EXPECT_EQ(param.qmmm.orca.output_file, strs_res[Key::output_file]);
  EXPECT_EQ(param.qmmm.orca.output_gradient_file, strs_res[Key::output_gradient_file]);
  EXPECT_EQ(param.qmmm.orca.output_mm_gradient_file, strs_res[Key::output_mm_gradient_file]);
}

} // namespace testing