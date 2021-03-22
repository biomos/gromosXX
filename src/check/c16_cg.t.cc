/**
 * @file c16_cg.t.cc
 * tests using c16cg
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

#include <time.h>

#include "check.h"

#include "check_forcefield.h"
#include "check_state.h"

void hard_coded_values(std::map<std::string, double> & m){
  m["QuarticBond"] = 8.028505;
//  m["PerturbedQuarticBond"] = 1.149568;
  m["Angle"] = 0.675853;
//  m["PerturbedAngle"] = 0.714818;
  m["ImproperDihedral"] = 0.000000;
//  m["PerturbedImproperDihedral"] = 0.000000;
  m["Dihedral"] = 0.000000;
//  m["PerturbedDihedral"] = 13.314602;
  m["NonBonded_cg"] = -1.340336;
//  m["NonBonded"] = -50.196817;
//  m["NonBonded_atomic"] =  -49.912;
 // m["DistanceRestraint"] = 257.189539;
//  m["PerturbedDistanceRestraint"] = 195.899012;
//  m["DihedralRestraint"] = 2127.910749;
//  m["PerturbedDihedralRestraint"] = 279.207857;
}

#ifdef OMP
  #include <omp.h>
#endif

int main(int argc, char* argv[]) {

#ifdef OMP
  omp_set_num_threads(1);
#endif

  int total = 0;

  util::Known knowns;
  knowns << "topo" <<  "conf" << "input" << "verb";

  std::string usage = argv[0];
  usage += "\n\t[@topo    <topology>]\n";
  usage += "\t[@conf    <starting configuration>]\n";
  usage += "\t[@input   <input>]\n";
  usage += "\t[@verb   <[module:][submodule:]level>]\n";

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
    GETFILEPATH(stopo, "cg16.topo", "src/check/data/");
    
  if(args.count("conf") == 1)
    sconf = args["conf"];
  else
    GETFILEPATH(sconf, "cg16.conf", "src/check/data/");

  if(args.count("input") == 1)
    sinput = args["input"];
  else
    GETFILEPATH(sinput, "cg16.in", "src/check/data/");

  if (!quiet)
    std::cout << "\n\n"
	      << "topology :      " << stopo << "\n"
	      << "input :         " << sinput << "\n"
	      << "configuration : " << sconf << "\n"
	      << std::endl;

  // set hard coded values to compare to
  std::map<std::string, double> ref_values;
  hard_coded_values(ref_values);

  util::simulation_struct c16_cg_sim;
  io::In_Topology in_topo;

  in_topo.quiet = quiet;
      
  if (util::create_simulation(stopo,
			      spttopo,
			      sconf,
			      sinput,
			      c16_cg_sim,
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
					 c16_cg_sim.topo,
					 c16_cg_sim.sim,
					 in_topo,
					 std::cout,
					 quiet)
      != 0){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }
  
  ff->init(c16_cg_sim.topo, c16_cg_sim.conf, c16_cg_sim.sim, std::cout,  quiet);
  
  // first check the forcefield
  total += check::check_forcefield(c16_cg_sim.topo, c16_cg_sim.conf,
				   c16_cg_sim.sim, *ff, ref_values);
    
  return total;
}
