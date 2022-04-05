/**
 * @file check_orca.cc
 * definition of tests and initializations
 */

#include <stdheader.h>

#include <cassert>

#include <check.h>
#include <io/argument.h>
#include <simulation/simulation.h>
#include <check_orca.h>

void addFile(std::string filetype, std::string ending, io::Argument& args) {
  std::string filename;
  GETFILEPATH(filename, "menthol." + ending, "src/check/data/orca/");
  std::vector<std::string> helper;
  helper.push_back(filename);
  args.put(filetype, helper);  
}

int checkSoftware(simulation::Simulation& sim) {
  if (sim.param().qmmm.software == simulation::qm_orca) {
    std::cout << "Software enum matches expected value" << std::endl;
    return 0;
  } 
  else {
    std::cout << "Software enum does NOT match expected value" << std::endl;
    return 1;
  }
}

int checkBinary(simulation::Simulation& sim) {
  if (sim.param().qmmm.orca.binary == "/home/fpultar/opt/orca-5.0.3/orca") {
    std::cout << "QM binary matches expected value" << std::endl;
    return 0;
  }
  else {
    std::cout << "QM binary does NOT match expected value" << std::endl;
    return 1;
  }
}

int checkOrcaFileNames(simulation::Simulation& sim) {
  std::cout << sim.param().qmmm.orca.input_file << std::endl;
  std::cout << sim.param().qmmm.orca.output_file << std::endl;
}