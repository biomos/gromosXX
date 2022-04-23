#include "../stdheader.h"

#include "qm_worker.t.h"

#include "../interaction/qmmm/qmmm_interaction.h"

test::QM_Worker_Test::QM_Worker_Test(std::string binary_name
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

void test::QM_Worker_Test::test_friend() {
  this->test_sim_.init_simulation();
  this->qmmm_interaction_ptr = interaction::QMMM_Interaction::pointer();
  this->qm_worker_ptr = this->qmmm_interaction_ptr->m_worker;
  this->qm_zone_ptr = this->qmmm_interaction_ptr->m_qm_zone;
  std::cout << this->qmmm_interaction_ptr << '\n'
            << this->qm_worker_ptr << '\n'
            << this->qm_zone_ptr << std::endl;
  std::cout << this->qm_worker_ptr->name() << std::endl;
}