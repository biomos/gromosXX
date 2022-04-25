/**
 * @file orca_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"
#include "check.h"

#include "qm_worker_test.h"
#include "orca_worker_test.h"

namespace testing {
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
  this->test_sim_.init_simulation();
  this->qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  this->qm_worker_ptr = this->qmmm_interaction_ptr->m_worker;
  this->qm_zone_ptr = this->qmmm_interaction_ptr->m_qm_zone;
}

void Orca_Worker_Test::TearDown() {

}

TEST_F(Orca_Worker_Test, check_worker_init) {
  // same basic sanity checks
  EXPECT_TRUE(this->qmmm_interaction_ptr != nullptr) << "qmmm_interaction not initialized";
  EXPECT_TRUE(this->qm_worker_ptr != nullptr) << "qm_worker not initialized";
  EXPECT_TRUE(this->qm_zone_ptr != nullptr) << "qm_zone not initialized";
  EXPECT_EQ("Orca Worker", this->qm_worker_ptr->name()) << "QM worker name does not match expected value";
}

TEST_F(Orca_Worker_Test, check_qm_zone) {
  // test the qm zone
  const auto& qm_zone = this->qm_zone_ptr;
  EXPECT_EQ(qm_zone->qm.size(), 90) << "Size of the QM zone has incorrect value: " << qm_zone->qm.size();
  qm_zone->mm.size();
}

}