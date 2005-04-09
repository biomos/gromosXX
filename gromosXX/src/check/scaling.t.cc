/**
 * @file scaling.t.cc
 * tests scaling using aladip
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
#include <util/error.h>

#include <interaction/interaction_types.h>
#include <io/instream.h>
#include <util/parse_tcouple.h>
#include <io/blockinput.h>
#include <io/topology/in_topology.h>
#include <io/print_block.h>

#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>
#include <algorithm/pressure/pressure_calculation.h>
#include <algorithm/pressure/berendsen_barostat.h>

#include <interaction/forcefield/create_forcefield.h>

#include <util/create_simulation.h>
#include <algorithm/create_md_sequence.h>

#include <ctime>

#include <config.h>

#include "check.h"

#include "check_forcefield.h"
#include "check_state.h"

int main(int argc, char* argv[])
{

  int total = 0;
  
  char *knowns[] = 
    {
      "verb"
    };
  
  int nknowns = 1;
    
  std::string usage = argv[0];
  usage += "\t[@verb   <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, nknowns, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);
      
  std::string stopo, spttopo, slambdadeppttopo, sconf, sinput, slambdadepinput;
  bool quiet = true;

  if (args.count("verb") != -1) quiet = false;

  GETFILEPATH(stopo, "aladip.topo", "src/check/data/");
  GETFILEPATH(spttopo, "scaling.pttopo", "src/check/data/");
  GETFILEPATH(slambdadeppttopo, "scaling.lambdadep.pttopo", "src/check/data/");
  GETFILEPATH(sconf, "aladip.conf", "src/check/data/");
  GETFILEPATH(sinput, "scaling.in", "src/check/data/");
  GETFILEPATH(slambdadepinput, "scaling.lambdadep.in", "src/check/data/");

  if (!quiet)
    std::cout << "\n\n"
	      << "topology :               " << stopo << "\n"
	      << "perturbation :           " << spttopo << "\n"
	      << "lambdadep perturbation : " << slambdadeppttopo << "\n"
	      << "input :                  " << sinput << "\n"
	      << "lambdadep input :        " << slambdadepinput << "\n"
	      << "configuration :          " << sconf << "\n"
	      << std::endl;

  util::simulation_struct aladip_sim;
  util::simulation_struct aladip_lambdadep_sim;
      
  io::In_Topology in_topo;
  io::In_Topology in_lambdadep_topo;

  in_topo.quiet = quiet;
  in_lambdadep_topo.quiet = quiet;
      
  if (util::create_simulation(stopo,
			      spttopo,
			      sconf,
			      sinput,
			      aladip_sim,
			      in_topo,
			      "",
			      quiet
			      )
      != 0){
    std::cerr << "creating simulation failed!" << std::endl;
    return 1;
  }
      
  // std::cout << "aladip sim" << std::endl;
  // io::messages.display(std::cout);
  // io::messages.clear();

  if (util::create_simulation(stopo,
			      slambdadeppttopo,
			      sconf,
			      slambdadepinput,
			      aladip_lambdadep_sim,
			      in_lambdadep_topo,
			      "",
			      quiet
			      )
      != 0){
    std::cerr << "creating lambda dependent simulation failed!" << std::endl;
    return 1;
  }

  // std::cout << "aladip lambdadep sim" << std::endl;
  // io::messages.display(std::cout);
  // io::messages.clear();

  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
  interaction::Forcefield *lambdadep_ff = new interaction::Forcefield;
      
  if (interaction::create_g96_forcefield(*ff, 
					 aladip_sim.topo,
					 aladip_sim.sim,
					 aladip_sim.conf,
					 in_topo,
					 std::cout,
					 quiet)
      != 0){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }

  ff->init(aladip_sim.topo, aladip_sim.conf, aladip_sim.sim, std::cout, quiet);

  // std::cout << "aladip forcefield" << std::endl;
  // io::messages.display(std::cout);
  // io::messages.clear();

  if (interaction::create_g96_forcefield(*lambdadep_ff, 
					 aladip_lambdadep_sim.topo,
					 aladip_lambdadep_sim.sim,
					 aladip_lambdadep_sim.conf,
					 in_lambdadep_topo,
					 std::cout,
					 quiet)
      != 0){
    std::cerr << "creating lambda dependent forcefield failed!" << std::endl;
    return 1;
  }

  lambdadep_ff->init(aladip_lambdadep_sim.topo, aladip_lambdadep_sim.conf, aladip_lambdadep_sim.sim, std::cout, quiet);
  
  // std::cout << "aladip lambdadep forcefield" << std::endl;
  // io::messages.display(std::cout);
  // io::messages.clear();

  int res;
  total = 0;


  CHECKING("different lambda dependence (nonbonded)", res);

  ff->apply(aladip_sim.topo, 
	    aladip_sim.conf,
	    aladip_sim.sim);

  aladip_sim.conf.current().perturbed_energy_derivatives.calculate_totals();
  if (!quiet){
    io::print_ENERGY(std::cout,
		     aladip_sim.conf.current().perturbed_energy_derivatives,
		     aladip_sim.topo.energy_groups());
  }
      
  // save the un - scaled energies
  std::vector<std::vector<double> > lj_energy, crf_energy;
  lj_energy.resize(aladip_sim.conf.current().perturbed_energy_derivatives.lj_energy.size());
  crf_energy.resize(aladip_sim.conf.current().perturbed_energy_derivatives.crf_energy.size());

  for(unsigned int i=0; i<lj_energy.size(); ++i){
    lj_energy[i] = aladip_sim.conf.current().perturbed_energy_derivatives.lj_energy[i];
    crf_energy[i] = aladip_sim.conf.current().perturbed_energy_derivatives.crf_energy[i];
  }
      
  // set lambda correctly
  const double lp = aladip_sim.topo.lambda();
  double dlp;

  // std::cout << "\nscaled only : " << aladip_lambdadep_sim.sim.param().perturbation.scaled_only << std::endl;

  // for all energy groups which are scaled
  std::map<std::pair<int, int>, std::pair<int, double> >::const_iterator
    it = aladip_lambdadep_sim.topo.energy_group_lambdadep().begin(),
    to = aladip_lambdadep_sim.topo.energy_group_lambdadep().end();
      
  for( ; it != to; ++it){
    // account for i,j and j,i pairs being present (only due it once...)
    if (it->first.first > it->first.second) continue;

    const double alpha = it->second.second;

    if (!quiet){
      std::cout << "\nenergy group " << it->first.first << " - " 
		<< it->first.second << std::endl;
      std::cout << "\talpha = " << alpha << std::endl;
    }
	
    if (alpha != 0.0){
      const double l = (alpha - 1 + sqrt((1-alpha)*(1-alpha) + 4 * alpha * lp)) 
	/ (2 * alpha);
	  
      aladip_lambdadep_sim.topo.lambda(l);
      dlp = (2 * l - 1.0) * alpha + 1;

      if (!quiet){
	std::cout << "\tdl'/dl = " << dlp << std::endl;
	std::cout << "\t=> l = " << l << " for l' = " << lp << std::endl;
      }
    }
    else{
      aladip_lambdadep_sim.topo.lambda(lp);
      dlp = (2 * lp - 1.0) * alpha + 1;	  
      if (!quiet)
	std::cout << "\tdl'/dl = " << dlp << std::endl;
    }

    lambdadep_ff->apply(aladip_lambdadep_sim.topo,
			aladip_lambdadep_sim.conf,
			aladip_lambdadep_sim.sim);

    aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives.calculate_totals();
    if (!quiet)
      io::print_ENERGY(std::cout,
		       aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives,
		       aladip_lambdadep_sim.topo.energy_groups());
	
    // the total (nonbonded) energy lambda derivative (including the dl'/dl term)
    double e =
      aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives.
      lj_energy[it->first.first][it->first.second] +
      aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives.
      crf_energy[it->first.first][it->first.second];

    // the non - scaled dE/dl
    double normal_e =
      lj_energy[it->first.first][it->first.second] +
      crf_energy[it->first.first][it->first.second];

    // add the j,i part (if atoms in the second energy group are perturbed)
    if (it->first.first != it->first.second){
      e +=
	aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives.
	lj_energy[it->first.second][it->first.first] +
	aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives.
	crf_energy[it->first.second][it->first.first];
	  
      normal_e +=
	lj_energy[it->first.second][it->first.first] +
	crf_energy[it->first.second][it->first.first];
    }

    // get rid of the dl'/dl derivative
    e /= dlp;
    CHECK_APPROX_EQUAL(normal_e, e, 0.0000001, res);
  }

  RESULT(res, total);

  return total;
}
