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

#include "../math/gmath.h"

#include "../interaction/qmmm/orca_worker.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

Orca_Worker_Electrostatic_Test::Orca_Worker_Electrostatic_Test() : QM_Worker_Test(Parameter("Orca QM/MM Tests Electrostatic Embedding", "orca_worker_test", "src/check/data/orca")) {
  test_sim_.parameter().add_input("topo", "md.top");
  test_sim_.parameter().add_input("input", "md_el.imd");
  test_sim_.parameter().add_input("conf", "md.cnf");
  test_sim_.parameter().add_input("qmmm", "md.qmmm");
  test_sim_.parameter().add_input("trc", "md.trc");
  test_sim_.parameter().add_input("tre", "md.tre");
  test_sim_.parameter().add_input("fin", "md_final.cnf");
}

void Orca_Worker_Electrostatic_Test::SetUp() {
  int err = test_sim_.init_simulation();
  ASSERT_EQ(err, 0) << "Initialization of the simulation unsuccessful. Error code: " << err;
  qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  qm_worker_ptr = dynamic_cast<interaction::Orca_Worker*>(qmmm_interaction_ptr->m_worker);
  qm_zone_ptr = qmmm_interaction_ptr->m_qm_zone;
}

void Orca_Worker_Electrostatic_Test::TearDown() {
  qmmm_interaction_ptr = nullptr;
  qm_worker_ptr = nullptr;
  qm_zone_ptr = nullptr;
}

/**
 * Test if initialization is successful
 * 
 */
TEST_F(Orca_Worker_Electrostatic_Test, check_init) {
  // references to shorten the code
  const topology::Topology& topo = this->test_sim_.topo();
  const simulation::Simulation& sim = this->test_sim_.sim();
  const configuration::Configuration& conf = this->test_sim_.conf();
  const simulation::Parameter& param = sim.param();
  const auto interaction = this->qmmm_interaction_ptr;
  const interaction::Orca_Worker* worker = dynamic_cast<interaction::Orca_Worker*>(this->qm_worker_ptr);
  const auto qm_zone = this->qm_zone_ptr;

  // check if input vile has been loaded correctly
  EXPECT_EQ(param.qmmm.qmmm, simulation::qmmm_electrostatic);
  EXPECT_EQ(param.qmmm.qm_ch, simulation::qm_ch_constant);
  EXPECT_EQ(param.qmmm.software, simulation::qm_orca);
  EXPECT_EQ(param.qmmm.cutoff, 1.2);
  EXPECT_EQ(param.qmmm.qm_lj, simulation::qm_lj_off);
  EXPECT_EQ(param.qmmm.qm_constraint, simulation::qm_constr_off);
  EXPECT_EQ(param.qmmm.mm_scale, -1.0);

  // units and conversion factors
  EXPECT_EQ(param.qmmm.orca.unit_factor_length, 0.1);
  EXPECT_EQ(param.qmmm.orca.unit_factor_energy, 2625.5); // tolerance (float)
  EXPECT_EQ(param.qmmm.orca.unit_factor_force, 49641);
  EXPECT_EQ(param.qmmm.orca.unit_factor_charge, 1.0);

  // binary and input file names
  EXPECT_EQ(param.qmmm.orca.binary, "/home/fpultar/opt/orca-5.0.3/orca"); // modify this
  EXPECT_EQ(param.qmmm.orca.input_file, "xphos.inp");
  EXPECT_EQ(param.qmmm.orca.input_coordinate_file, "xphos.xyz");
  EXPECT_EQ(param.qmmm.orca.input_pointcharges_file, "xphos.pc");
  EXPECT_EQ(param.qmmm.orca.output_file, "xphos.out");
  EXPECT_EQ(param.qmmm.orca.output_gradient_file, "xphos.engrad");
  EXPECT_EQ(param.qmmm.orca.output_mm_gradient_file, "xphos.pcgrad");

  // some tests for elements and iac_elements map
  EXPECT_EQ(param.qmmm.orca.elements.at(1), "H");
  EXPECT_EQ(param.qmmm.orca.elements.at(6), "C");
  EXPECT_EQ(param.qmmm.orca.elements.at(8), "O");
  EXPECT_EQ(param.qmmm.orca.elements.at(15), "P");
  EXPECT_EQ(param.qmmm.orca.elements.at(46), "Pd");
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(0), 8); // check if (-1) is subtracted from input file IAC
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(23), 29);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(24), 29);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(25), 26);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(28), 20);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(69), 15);
  
  // check if the qm interaction, qm_worker, and qm_zone have been initialized
  EXPECT_TRUE(this->qmmm_interaction_ptr != nullptr) << "qmmm_interaction not initialized";
  EXPECT_TRUE(this->qm_worker_ptr != nullptr) << "qm_worker not initialized";
  EXPECT_TRUE(this->qm_zone_ptr != nullptr) << "qm_zone not initialized";
  EXPECT_EQ("Orca Worker", this->qm_worker_ptr->name()) << "QM worker name does not match expected value";

  // test the qm zone
  EXPECT_EQ(qm_zone->qm.size(), 90) << "Size of the QM zone has incorrect value: " << qm_zone->qm.size();
  EXPECT_EQ(qm_zone->mm.size(), 1194) << "Size of the MM zone has incorrect value: " << qm_zone->mm.size();
  EXPECT_EQ(qm_zone->charge(), 0) << "Charge of the QM zone has incorrect value " << qm_zone->charge();
  EXPECT_EQ(qm_zone->spin_mult(), 1) << "Spin multiplicity of the QM zone has incorrect value: " << qm_zone->spin_mult();
  EXPECT_EQ(qm_zone->QM_energy(), 0) << "Energy of the QM zone is not zero after initialization.";

  // test the atoms of the QM and MM zone
  for (const auto& atom : qm_zone->qm) {
    EXPECT_EQ(atom.qm_charge, 0.0) << "QM charge is not 0.0 at beginning of simulation for atom " << atom.index; // should be set for mechanical embedding
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force(i), 0.0) << "Force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
    }
    switch (atom.index) { // test some selected atoms
      case 0:  EXPECT_EQ(atom.atomic_number, 1); break;
      case 48: EXPECT_EQ(atom.atomic_number, 15); break;
      case 49: EXPECT_EQ(atom.atomic_number, 46); break;
      case 89: EXPECT_EQ(atom.atomic_number, 1); break;
    }
  }
  for (const auto& atom : qm_zone->mm) {
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force[i], 0.0) << "Force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
      EXPECT_EQ(atom.cos_force[i], 0.0) << "COS force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
    }
    switch (atom.index) { // test some selected atoms
      case 105:  EXPECT_EQ(atom.atomic_number, 1); break; 
      case 106:  EXPECT_EQ(atom.atomic_number, 6); break; 
      case 3232: EXPECT_EQ(atom.atomic_number, 8); break;
      case 3233: EXPECT_EQ(atom.atomic_number, 1); break;
    }
  }
}

Orca_Worker_Mechanical_Test::Orca_Worker_Mechanical_Test() : QM_Worker_Test(Parameter("Orca QM/MM Tests Mechanical Embedding", "orca_worker_test", "src/check/data/orca")) {
  test_sim_.parameter().add_input("topo", "md.top");
  test_sim_.parameter().add_input("input", "md_mc.imd");
  test_sim_.parameter().add_input("conf", "md.cnf");
  test_sim_.parameter().add_input("qmmm", "md.qmmm");
  test_sim_.parameter().add_input("trc", "md.trc");
  test_sim_.parameter().add_input("tre", "md.tre");
  test_sim_.parameter().add_input("fin", "md_final.cnf");
}

void Orca_Worker_Mechanical_Test::SetUp() {
  int err = test_sim_.init_simulation();
  ASSERT_EQ(err, 0) << "Initialization of the simulation unsuccessful. Error code: " << err;
  qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  qm_worker_ptr = dynamic_cast<interaction::Orca_Worker*>(qmmm_interaction_ptr->m_worker);
  qm_zone_ptr = qmmm_interaction_ptr->m_qm_zone;
}

void Orca_Worker_Mechanical_Test::TearDown() {
  qmmm_interaction_ptr = nullptr;
  qm_worker_ptr = nullptr;
  qm_zone_ptr = nullptr;
}

/**
 * Test if initialization is successful
 * 
 */
TEST_F(Orca_Worker_Mechanical_Test, check_init) {
  // references to shorten the code
  const topology::Topology& topo = this->test_sim_.topo();
  const simulation::Simulation& sim = this->test_sim_.sim();
  const configuration::Configuration& conf = this->test_sim_.conf();
  const simulation::Parameter& param = sim.param();
  const auto interaction = this->qmmm_interaction_ptr;
  const interaction::Orca_Worker* worker = dynamic_cast<interaction::Orca_Worker*>(this->qm_worker_ptr);
  const auto qm_zone = this->qm_zone_ptr;

  // check if input vile has been loaded correctly
  EXPECT_EQ(param.qmmm.qmmm, simulation::qmmm_mechanical);
  EXPECT_EQ(param.qmmm.qm_ch, simulation::qm_ch_dynamic);
  EXPECT_EQ(param.qmmm.software, simulation::qm_orca);
  EXPECT_EQ(param.qmmm.cutoff, 1.2);
  EXPECT_EQ(param.qmmm.qm_lj, simulation::qm_lj_off);
  EXPECT_EQ(param.qmmm.qm_constraint, simulation::qm_constr_off);
  EXPECT_EQ(param.qmmm.mm_scale, -1.0);

  // units and conversion factors
  EXPECT_EQ(param.qmmm.orca.unit_factor_length, 0.1);
  EXPECT_EQ(param.qmmm.orca.unit_factor_energy, 2625.5); // tolerance (float)
  EXPECT_EQ(param.qmmm.orca.unit_factor_force, 49641);
  EXPECT_EQ(param.qmmm.orca.unit_factor_charge, 1.0);

  // binary and input file names
  EXPECT_EQ(param.qmmm.orca.binary, "/home/fpultar/opt/orca-5.0.3/orca"); // modify this
  EXPECT_EQ(param.qmmm.orca.input_file, "xphos.inp");
  EXPECT_EQ(param.qmmm.orca.input_coordinate_file, "xphos.xyz");
  EXPECT_EQ(param.qmmm.orca.input_pointcharges_file, "xphos.pc");
  EXPECT_EQ(param.qmmm.orca.output_file, "xphos.out");
  EXPECT_EQ(param.qmmm.orca.output_gradient_file, "xphos.engrad");
  EXPECT_EQ(param.qmmm.orca.output_mm_gradient_file, "xphos.pcgrad");

  // some tests for elements and iac_elements map
  EXPECT_EQ(param.qmmm.orca.elements.at(1), "H");
  EXPECT_EQ(param.qmmm.orca.elements.at(6), "C");
  EXPECT_EQ(param.qmmm.orca.elements.at(8), "O");
  EXPECT_EQ(param.qmmm.orca.elements.at(15), "P");
  EXPECT_EQ(param.qmmm.orca.elements.at(46), "Pd");
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(0), 8); // check if (-1) is subtracted from input file IAC
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(23), 29);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(24), 29);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(25), 26);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(28), 20);
  EXPECT_EQ(param.qmmm.orca.iac_elements.at(69), 15);
  
  // check if the qm interaction, qm_worker, and qm_zone have been initialized
  EXPECT_TRUE(this->qmmm_interaction_ptr != nullptr) << "qmmm_interaction not initialized";
  EXPECT_TRUE(this->qm_worker_ptr != nullptr) << "qm_worker not initialized";
  EXPECT_TRUE(this->qm_zone_ptr != nullptr) << "qm_zone not initialized";
  EXPECT_EQ("Orca Worker", this->qm_worker_ptr->name()) << "QM worker name does not match expected value";

  // test the qm zone
  EXPECT_EQ(qm_zone->qm.size(), 90) << "Size of the QM zone has incorrect value: " << qm_zone->qm.size();
  EXPECT_EQ(qm_zone->mm.size(), 0) << "Size of the MM zone has incorrect value: " << qm_zone->mm.size();
  EXPECT_EQ(qm_zone->charge(), 0) << "Charge of the QM zone has incorrect value " << qm_zone->charge();
  EXPECT_EQ(qm_zone->spin_mult(), 1) << "Spin multiplicity of the QM zone has incorrect value: " << qm_zone->spin_mult();
  EXPECT_EQ(qm_zone->QM_energy(), 0) << "Energy of the QM zone is not zero after initialization.";

  // test the atoms of the QM and MM zone
  for (const auto& atom : qm_zone->qm) {
    EXPECT_EQ(atom.qm_charge, 0.0) << "QM charge is not 0.0 at beginning of simulation for atom " << atom.index; // should be set for mechanical embedding
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force(i), 0.0) << "Force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
    }
    switch (atom.index) { // test some selected atoms
      case 0:  EXPECT_EQ(atom.atomic_number, 1); break;
      case 48: EXPECT_EQ(atom.atomic_number, 15); break;
      case 49: EXPECT_EQ(atom.atomic_number, 46); break;
      case 89: EXPECT_EQ(atom.atomic_number, 1); break;
    }
  }
  for (const auto& atom : qm_zone->mm) {
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force[i], 0.0) << "Force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
      EXPECT_EQ(atom.cos_force[i], 0.0) << "COS force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
    }
    switch (atom.index) { // test some selected atoms
      case 105:  EXPECT_EQ(atom.atomic_number, 1); break; 
      case 106:  EXPECT_EQ(atom.atomic_number, 6); break; 
      case 3232: EXPECT_EQ(atom.atomic_number, 8); break;
      case 3233: EXPECT_EQ(atom.atomic_number, 1); break;
    }
  }
}

}