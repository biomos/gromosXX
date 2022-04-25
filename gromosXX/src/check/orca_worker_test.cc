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
#include "error_codes.h"

#include "qm_worker_test.h"
#include "orca_worker_test.h"

#include "../math/gmath.h"

#include "../interaction/qmmm/orca_worker.h"
#include "../interaction/qmmm/qm_atom.h"

namespace testing {

// test files live here
const static std::string orca_test_files = "src/check/data/orca/";

Orca_Worker_Test::Orca_Worker_Test() : QM_Worker_Test(
  "orca_worker_test"
, "Orca QM/MM Tests"
, TOP_SOURCE_DIR "/" + orca_test_files + "md.top"
, TOP_SOURCE_DIR "/" + orca_test_files + "md.imd"
, TOP_SOURCE_DIR "/" + orca_test_files + "md.cnf"
, TOP_SOURCE_DIR "/" + orca_test_files + "md.qmmm"
, TOP_SOURCE_DIR "/" + orca_test_files + "md.trc"
, TOP_SOURCE_DIR "/" + orca_test_files + "md.tre"
, TOP_SOURCE_DIR "/" + orca_test_files + "md_final.cnf") {};

void Orca_Worker_Test::SetUp() {
  int err = this->test_sim_.init_simulation();
  ASSERT_EQ(err, 0) << "Initialization of the simulation unsuccessful. Error code: " << err;
  this->qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  this->qm_worker_ptr = (interaction::Orca_Worker*)this->qmmm_interaction_ptr->m_worker;
  this->qm_zone_ptr = this->qmmm_interaction_ptr->m_qm_zone;
}

void Orca_Worker_Test::TearDown() {}

/**
 * Test if initialization is successful
 * 
 */
TEST_F(Orca_Worker_Test, check_init) {
  // references to shorten the code
  const auto& interaction = this->qmmm_interaction_ptr;
  const auto& worker = (interaction::Orca_Worker*)this->qm_worker_ptr;
  const auto& qm_zone = this->qm_zone_ptr;
  // same basic sanity checks
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
    EXPECT_EQ(atom.qm_charge, 0) << "QM charge is not 0 for atom " << atom.index;
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force[i], 0.0) << "Force component is not initialized to 0.0 for atom " << atom.index << " at index " << i;
    }
  }
  
}

}