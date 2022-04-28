/**
 * @file orca_worker_test_mechanical.cc
 */

#include "orca_worker_test_mechanical.h"
#include "../interaction/qmmm/orca_worker.h"

namespace testing {

Orca_Worker_Test_Mechanical::Orca_Worker_Test_Mechanical() : Orca_Worker_Test(orca_parameter_mechanical, orca_results_mechanical) {
  // initialize test files and expected results
  init_parameters();
  init_results();
}

void Orca_Worker_Test_Mechanical::SetUp() {
  // initialize test simulation and get a pointer to qm/mm relevant objects
  int err = test_sim_.init_simulation();
  ASSERT_EQ(err, 0) << "Initialization of the simulation unsuccessful. Error code: " << err;
  qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  qm_worker_ptr = dynamic_cast<interaction::Orca_Worker*>(qmmm_interaction_ptr->m_worker);
  qm_zone_ptr = qmmm_interaction_ptr->m_qm_zone;
}

void Orca_Worker_Test_Mechanical::TearDown() {
  // set all pointers to nullptr again
  qmmm_interaction_ptr = nullptr;
  qm_worker_ptr = nullptr;
  qm_zone_ptr = nullptr;
}

void Orca_Worker_Test_Mechanical::init_parameters() {
  test_sim_.parameter().add_input("topo", "md.top");
  test_sim_.parameter().add_input("input", "md_mc.imd");
  test_sim_.parameter().add_input("conf", "md.cnf");
  test_sim_.parameter().add_input("qmmm", "md.qmmm");
  test_sim_.parameter().add_input("trc", "md.trc");
  test_sim_.parameter().add_input("tre", "md.tre");
  test_sim_.parameter().add_input("fin", "md_final.cnf");
}

void Orca_Worker_Test_Mechanical::init_results() {
  // call helper functions
  init_results_parameters();
  init_results_binary_name();
  init_results_units();
  init_results_elements();
  init_results_files();
  init_results_qm_zone_init();
}

void Orca_Worker_Test_Mechanical::init_results_parameters() {
  // parameters from the input file
  results_.ints_[Key::qmmm] = simulation::qmmm_mechanical;
  results_.ints_[Key::qm_software] = simulation::qm_orca;
  results_.ints_[Key::qm_ch] = simulation::qm_ch_dynamic;
  results_.doubles_[Key::qm_cutoff] = 1.2;
  results_.ints_[Key::qm_lj] = simulation::qm_lj_off;
  results_.ints_[Key::qm_constraint] = simulation::qm_constr_off;
  results_.doubles_[Key::qm_mm_scale] = -1.0;
}

void Orca_Worker_Test_Mechanical::init_results_files() {
  // input file names
  results_.strs_[Key::input_file] = "xphos.inp";
  results_.strs_[Key::input_coordinate_file] = "xphos.xyz";
  results_.strs_[Key::input_pointcharges_file] = "xphos.pc";
  results_.strs_[Key::output_file] = "xphos.out";
  results_.strs_[Key::output_gradient_file] = "xphos.engrad";
  results_.strs_[Key::output_mm_gradient_file] = "xphos.pcgrad";
}

void Orca_Worker_Test_Mechanical::init_results_qm_zone_init() {
  // QM zone
  results_.ints_[Key::qm_zone_size_init] = 90;
  results_.ints_[Key::mm_zone_size_init] = 0;
  results_.ints_[Key::qm_zone_charge] = 0;
  results_.ints_[Key::qm_zone_spin_mult] = 1;
  results_.doubles_[Key::qm_energy_init] = 0.0;
  results_.doubles_[Key::qm_atom_charge_init] = 0.0;
  results_.doubles_[Key::qm_atom_force_init] = 0.0;
  results_.doubles_[Key::mm_atom_force_init] = 0.0;
  results_.doubles_[Key::mm_atom_cos_force_init] = 0.0;
}

/**
 * Test if initialization is successful
 * 
 */
TEST_F(Orca_Worker_Test_Mechanical, check_init) {
  // check if input vile has been loaded correctly
  check_parameter_init(); 
}

} // namespace testing