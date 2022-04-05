/**
 * @file check_orca.h
 * definition of tests and initializations
 */

#ifndef CHECK_ORCA_H
#define CHECK_ORCA_H

#include <stdheader.h>
#include <io/argument.h>
#include <simulation/simulation.h>

template <typename T1, typename T2>
int check(T1& first, T2& second, const std::string& message) {
  if (first == second) {
    return 0;
  }
  else {
    std::cerr << message << std::endl;
    return 1;
  }
}

void addFile(std::string filetype, std::string ending, io::Argument& args);

int checkSoftware(simulation::Simulation& sim);

int checkBinary(simulation::Simulation& sim);

int checkOrcaFileNames(simulation::Simulation& sim);

#endif /* CHECK_ORCA_H */