/**
 * @file orca.t.cc
 * tests for the orca qm worker
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>
#include <signal.h>

#include <check.h>
#include <check_orca.h>

#ifdef OMP
  #include <omp.h>
#endif


int main(int argc, char** argv) {
#ifdef OMP
  omp_set_num_threads(1);
#endif
 
  util::Known knowns;
  knowns << "topo" << "conf" << "input" << "verb" << "qmmm";
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);
  usage += "#\n\n";

  io::Argument args;

  if (args.parse(argc, argv, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }
    
  // parse the verbosity flag and set debug levels
  if (util::parse_verbosity(args)){
    std::cerr << "could not parse verbosity argument" << std::endl;
    return 1;
  }

  std::string stopo, sconf, sinput, sqmmm;

  if(args.count("topo") != 1) {
    addFile("topo", "top", args);   
  }

  if(args.count("conf") != 1) {
    addFile("conf", "cnf", args);
  }
  
  if(args.count("input") != 1) {
    addFile("input", "imd", args);
  }

  if(args.count("qmmm") != 1) {
    addFile("qmmm", "qmmm", args);
  }

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  simulation::Simulation sim;
  algorithm::Algorithm_Sequence md;

  if (io::read_input(args, topo, conf, sim, md)) {
    io::messages.display(std::cout);
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }

  // initialises all algorithms (and therefore also the forcefield)
  md.init(topo, conf, sim);

  return checkSoftware(sim);

  return checkBinary(sim);

  return checkOrcaFileNames(sim);
}