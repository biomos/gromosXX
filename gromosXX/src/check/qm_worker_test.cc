/**
 * @file qm_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"

#include "qm_worker_test.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

QM_Worker_Test::QM_Worker_Test(const Parameter& parameter, const Results& results) : test_sim_(parameter), results_(results) {}

void QM_Worker_Test::check_parameter_init() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, int>& ints_res = results_.ints_;
  std::unordered_map<Key, double>& doubles_res = results_.doubles_;
  std::unordered_map<Key, std::string>& strs_res = results_.strs_;
  EXPECT_EQ(param.qmmm.qmmm, ints_res[Key::qmmm]);
  EXPECT_EQ(param.qmmm.qm_ch, ints_res[Key::qm_ch]);
  EXPECT_EQ(param.qmmm.software, ints_res[Key::qm_software]);
  EXPECT_EQ(param.qmmm.cutoff, doubles_res[Key::qm_cutoff]);
  EXPECT_EQ(param.qmmm.qm_lj, ints_res[Key::qm_lj]);
  EXPECT_EQ(param.qmmm.qm_constraint, ints_res[Key::qm_constraint]);
  EXPECT_EQ(param.qmmm.mm_scale, doubles_res[Key::qm_mm_scale]);

  // units and conversion factors
  EXPECT_EQ(param.qmmm.orca.unit_factor_length, doubles_res[Key::unit_factor_length]);
  EXPECT_EQ(param.qmmm.orca.unit_factor_energy, doubles_res[Key::unit_factor_energy]); // tolerance (float)
  EXPECT_EQ(param.qmmm.orca.unit_factor_force, doubles_res[Key::unit_factor_force]);
  EXPECT_EQ(param.qmmm.orca.unit_factor_charge, doubles_res[Key::unit_factor_charge]);

  // binary and input file names
  EXPECT_EQ(param.qmmm.orca.binary, strs_res[Key::binary]); // modify this
  EXPECT_EQ(param.qmmm.orca.input_file, strs_res[Key::input_file]);
  EXPECT_EQ(param.qmmm.orca.input_coordinate_file, strs_res[Key::input_coordinate_file]);
  EXPECT_EQ(param.qmmm.orca.input_pointcharges_file, strs_res[Key::input_pointcharges_file]);
  EXPECT_EQ(param.qmmm.orca.output_file, strs_res[Key::output_file]);
  EXPECT_EQ(param.qmmm.orca.output_gradient_file, strs_res[Key::output_gradient_file]);
  EXPECT_EQ(param.qmmm.orca.output_mm_gradient_file, strs_res[Key::output_mm_gradient_file]);

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

    // test the qm worker
  EXPECT_EQ(qm_worker_ptr->name(), strs_res[Key::qm_worker_name]);

  // test the qm zone
  EXPECT_EQ(qm_zone_ptr->qm.size(), ints_res[Key::qm_zone_size_init]);
  EXPECT_EQ(qm_zone_ptr->mm.size(), ints_res[Key::mm_zone_size_init]);
  EXPECT_EQ(qm_zone_ptr->charge(), ints_res[Key::qm_zone_charge]);
  EXPECT_EQ(qm_zone_ptr->spin_mult(), ints_res[Key::qm_zone_spin_mult]);
  EXPECT_EQ(qm_zone_ptr->QM_energy(), doubles_res[Key::qm_energy_init]);
  // test the atoms of the QM and MM zone
  for (const auto& atom : qm_zone_ptr->qm) {
    EXPECT_EQ(atom.qm_charge, doubles_res[Key::qm_atom_charge_init]);
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force(i), doubles_res[Key::qm_atom_force_init]);
    }
    switch (atom.index) { // test some selected atoms
      case 0:  EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break;
      case 48: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_15]); break;
      case 49: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_46]); break;
      case 89: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break;
    }
  }
  for (const auto& atom : qm_zone_ptr->mm) {
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force[i], doubles_res[Key::mm_atom_force_init]);
      EXPECT_EQ(atom.cos_force[i], doubles_res[Key::mm_atom_cos_force_init]);
    }
    switch (atom.index) { // test some selected atoms
      case 105:  EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break; 
      case 106:  EXPECT_EQ(atom.atomic_number, ints_res[Key::element_6]); break; 
      case 3232: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_8]); break;
      case 3233: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break;
    }
  }

  // check if the qm interaction, qm_worker, and qm_zone have been initialized
  EXPECT_TRUE(qmmm_interaction_ptr != results_.UNINITIALIZED) << "qmmm_interaction not initialized";
  EXPECT_TRUE(qm_worker_ptr != results_.UNINITIALIZED) << "qm_worker not initialized";
  EXPECT_TRUE(qm_zone_ptr != results_.UNINITIALIZED) << "qm_zone not initialized";
}

}