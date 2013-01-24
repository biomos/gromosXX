/**
 * @file aladip.t.cc
 * tests using aladip
 */


#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../util/error.h"

#include "../interaction/interaction_types.h"
#include "../io/instream.h"
#include "../util/parse_tcouple.h"
#include "../io/blockinput.h"
#include "../io/topology/in_topology.h"

#include "../algorithm/integration/leap_frog.h"
#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/berendsen_thermostat.h"
#include "../algorithm/pressure/pressure_calculation.h"
#include "../algorithm/pressure/berendsen_barostat.h"

#include "../interaction/forcefield/create_forcefield.h"

#include "../util/create_simulation.h"
#include "../algorithm/create_md_sequence.h"

#include <ctime>

#include "check.h"

#include "check_forcefield.h"
#include "check_state.h"

#ifdef OMP
  #include <omp.h>
#endif

int main(int argc, char* argv[])
{

#ifdef OMP
  omp_set_num_threads(1);
#endif

  int total = 0;

  util::Known knowns;
  knowns << "topo" << "pttopo" << "conf" << "input" << "verb";
    
  std::string usage = argv[0];
  usage += "\n\t[@topo    <topology>]\n";
  usage += "\t[@pttopo  <perturbation topology>]\n";
  usage += "\t[@conf    <starting configuration>]\n";
  usage += "\t[@input   <input>]\n";
  usage += "\t[@verb    <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);
      
  std::string stopo, spttopo, sconf, sinput;
  bool quiet = true;

  if (args.count("verb") != -1) quiet = false;
      
  if(args.count("topo") == 1)
    stopo = args["topo"];
  else
    GETFILEPATH(stopo, "aladip.topo", "src/check/data/");

  if(args.count("pttopo") == 1)
    spttopo = args["pttopo"];
  else
    GETFILEPATH(spttopo, "aladip.pttopo", "src/check/data/");
      
  if(args.count("conf") == 1)
    sconf = args["conf"];
  else
    GETFILEPATH(sconf, "aladip.conf", "src/check/data/");

  if(args.count("input") == 1)
    sinput = args["input"];
  else
    GETFILEPATH(sinput, "aladip_atomic.in", "src/check/data/");

  if (!quiet)
    std::cout << "\n\n"
	      << "topology :      " << stopo << "\n"
	      << "perturbation :  " << spttopo << "\n"
	      << "input :         " << sinput << "\n"
	      << "configuration : " << sconf << "\n"
	      << std::endl;

  util::simulation_struct aladip_sim;
  io::In_Topology in_topo;

  in_topo.quiet = quiet;
      
  if (util::create_simulation(stopo,
			      spttopo,
			      sconf,
			      sinput,
			      aladip_sim,
			      in_topo,
			      "", "", "", "", "", "",
			      quiet
			      )
      != 0){
    std::cerr << "creating simulation failed!" << std::endl;
    return 1;
  }

  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
      
  if (interaction::create_g96_forcefield(*ff, 
					 aladip_sim.topo,
					 aladip_sim.sim,
					 in_topo,
					 std::cout,
					 quiet)
      != 0){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }

  ff->init(aladip_sim.topo, aladip_sim.conf, aladip_sim.sim, std::cout, quiet);
      
  total += check::check_atomic_virial(aladip_sim.topo, aladip_sim.conf, aladip_sim.sim, *ff);
    
  return total;
}
