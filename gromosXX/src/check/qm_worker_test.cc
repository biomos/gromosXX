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

void QM_Worker_Test::check_qmmm_parameter() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, int>& ints_res = results_.ints_;
  std::unordered_map<Key, double>& doubles_res = results_.doubles_;
  EXPECT_EQ(param.qmmm.qmmm, ints_res[Key::qmmm]);
  EXPECT_EQ(param.qmmm.qm_ch, ints_res[Key::qm_ch]);
  EXPECT_EQ(param.qmmm.software, ints_res[Key::qm_software]);
  EXPECT_EQ(param.qmmm.cutoff, doubles_res[Key::qm_cutoff]);
  EXPECT_EQ(param.qmmm.qm_lj, ints_res[Key::qm_lj]);
  EXPECT_EQ(param.qmmm.qm_constraint, ints_res[Key::qm_constraint]);
  EXPECT_EQ(param.qmmm.mm_scale, doubles_res[Key::qm_mm_scale]);
}

void QM_Worker_Test::check_unit_conversion_factors() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, double>& doubles_res = results_.doubles_;
  // units and conversion factors
  EXPECT_EQ(param.qmmm.orca.unit_factor_length, doubles_res[Key::unit_factor_length]);
  EXPECT_EQ(param.qmmm.orca.unit_factor_energy, doubles_res[Key::unit_factor_energy]); // tolerance (float)
  EXPECT_EQ(param.qmmm.orca.unit_factor_force, doubles_res[Key::unit_factor_force]);
  EXPECT_EQ(param.qmmm.orca.unit_factor_charge, doubles_res[Key::unit_factor_charge]);
}

void QM_Worker_Test::check_qm_zone_param() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key, int>& ints_res = results_.ints_;
  std::unordered_map<Key, double>& doubles_res = results_.doubles_;
  // test the qm zone
  EXPECT_EQ(qm_zone_ptr->qm.size(), ints_res[Key::qm_zone_size_init]);
  EXPECT_EQ(qm_zone_ptr->mm.size(), ints_res[Key::mm_zone_size_init]);
  EXPECT_EQ(qm_zone_ptr->charge(), ints_res[Key::qm_zone_charge]);
  EXPECT_EQ(qm_zone_ptr->spin_mult(), ints_res[Key::qm_zone_spin_mult]);
  EXPECT_EQ(qm_zone_ptr->QM_energy(), doubles_res[Key::qm_energy_init]);
}

void QM_Worker_Test::check_qm_interaction_ptr() {
  // check if the qm interaction, qm_worker, and qm_zone have been initialized
  EXPECT_TRUE(qmmm_interaction_ptr != results_.UNINITIALIZED) << "qmmm_interaction not initialized";
  EXPECT_TRUE(qm_worker_ptr != results_.UNINITIALIZED) << "qm_worker not initialized";
  EXPECT_TRUE(qm_zone_ptr != results_.UNINITIALIZED) << "qm_zone not initialized";
}

}