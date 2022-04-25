/**
 * @file qm_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"

#include "qm_worker_test.h"

#include "../interaction/qmmm/qmmm_interaction.h"

namespace testing {

QM_Worker_Test::QM_Worker_Test(std::string binary_name
               , std::string test_title
               , std::string topology_file
               , std::string simulation_file
               , std::string configuration_file
               , std::string qmmm_file
               , std::string coordinate_trajectory_file
               , std::string energy_trajectory_file
               , std::string final_configuration_file) : 
                   test_sim_(binary_name
                 , test_title
                 , topology_file
                 , simulation_file
                 , configuration_file
                 , qmmm_file
                 , coordinate_trajectory_file
                 , energy_trajectory_file
                 , final_configuration_file) {}

void QM_Worker_Test::SetUp() {

}

void QM_Worker_Test::TearDown() {
  
}

}