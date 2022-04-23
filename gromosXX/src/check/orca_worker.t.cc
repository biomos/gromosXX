#include "../stdheader.h"

#include "check.h"

#include "qm_worker.t.h"
#include "orca_worker.t.h"

test::Orca_Worker_Test::Orca_Worker_Test(std::string binary_name
               , std::string test_title
               , std::string topology_file
               , std::string simulation_file
               , std::string configuration_file
               , std::string qmmm_file
               , std::string coordinate_trajectory_file
               , std::string energy_trajectory_file
               , std::string final_configuration_file)
               : test::QM_Worker_Test::QM_Worker_Test(binary_name, test_title, topology_file, simulation_file, configuration_file, qmmm_file, coordinate_trajectory_file, energy_trajectory_file, final_configuration_file) {}

int main(int argc, char** argv) {
  // set up simulation files
  #define ROOT_FOLDER "src/check/data/orca/"
  std::string topo_f, conf_f, input_f, qmmm_f, trc_f, tre_f, fin_f;
  GETFILEPATH(topo_f, "md.top", ROOT_FOLDER);
  GETFILEPATH(conf_f, "md.cnf", ROOT_FOLDER);
  GETFILEPATH(input_f, "md.imd", ROOT_FOLDER);
  GETFILEPATH(qmmm_f, "md.qmmm", ROOT_FOLDER);
  GETFILEPATH(trc_f, "md.trc", ROOT_FOLDER);
  GETFILEPATH(tre_f, "md.tre", ROOT_FOLDER);
  GETFILEPATH(fin_f, "md_final.cnf", ROOT_FOLDER);
  
  // initialize test
  test::QM_Worker_Test orca_worker_test(argv[0], "Orca QM/MM Test", topo_f, input_f, conf_f, qmmm_f, trc_f, tre_f, fin_f);
  orca_worker_test.test_friend();
}