/**
 * @file scaling.t.cc
 * tests scaling using aladip
 */


#include <util/stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

#include <simulation/parameter.h>
#include <interaction/interaction_types.h>
#include <io/instream.h>
#include <util/parse_tcouple.h>
#include <io/blockinput.h>
#include <io/topology/in_topology.h>

#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>
#include <algorithm/pressure/pressure_calculation.h>
#include <algorithm/pressure/berendsen_barostat.h>

#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <util/create_simulation.h>
#include <algorithm/create_md_sequence.h>

#include <time.h>

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

    try{
      io::Argument args(argc, argv, nknowns, knowns, usage, true);

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
					     quiet)
	  != 0){
	std::cerr << "creating forcefield failed!" << std::endl;
	return 1;
      }

      // std::cout << "aladip forcefield" << std::endl;
      // io::messages.display(std::cout);
      // io::messages.clear();

      if (interaction::create_g96_forcefield(*lambdadep_ff, 
					     aladip_lambdadep_sim.topo,
					     aladip_lambdadep_sim.sim,
					     aladip_lambdadep_sim.conf,
					     in_lambdadep_topo,
					     quiet)
	  != 0){
	std::cerr << "creating lambda dependent forcefield failed!" << std::endl;
	return 1;
      }

      // std::cout << "aladip lambdadep forcefield" << std::endl;
      // io::messages.display(std::cout);
      // io::messages.clear();

      int res, total = 0;


      CHECKING("different lambda dependence (nonbonded)", res);

      ff->apply(aladip_sim.topo, 
		aladip_sim.conf,
		aladip_sim.sim);

      aladip_sim.conf.current().perturbed_energy_derivatives[0].
	calculate_totals();

      // set lambda correctly
      const double lp = aladip_sim.topo.lambda();
      double nonbonded_der = 0;

      for(size_t s=0; 
	  s < aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives.size();
	  ++s){
	
	const double alpha = 
	  aladip_lambdadep_sim.topo.perturbed_energy_derivative_alpha()[s];

	// std::cerr << "\nalpha = " << alpha << std::endl;
	
	if (alpha != 0.0){
	  const double l = (alpha - 1 + sqrt((1-alpha)*(1-alpha) + 4 * alpha * lp)) 
	    / (2 * alpha);
	  
	  // std::cerr << "setting lambda to " << l << std::endl;
	  
	  aladip_lambdadep_sim.topo.lambda(l);
	}
	else{
	  // std::cerr << "normal lambda dependency: l = " << lp << std::endl;
	  aladip_lambdadep_sim.topo.lambda(lp);
	}
	
	lambdadep_ff->apply(aladip_lambdadep_sim.topo,
			    aladip_lambdadep_sim.conf,
			    aladip_lambdadep_sim.sim);

	aladip_lambdadep_sim.conf.current().perturbed_energy_derivatives[s].
	  calculate_totals();
	
	// std::cerr << "nonbonded_total[" << s << "] = "
	// << aladip_lambdadep_sim.conf.current().
	// perturbed_energy_derivatives[s].nonbonded_total << std::endl;

	nonbonded_der += aladip_lambdadep_sim.conf.current().
	  perturbed_energy_derivatives[s].nonbonded_total;

      }

      const double normal_nonbonded_der = 
	aladip_sim.conf.current().perturbed_energy_derivatives[0].nonbonded_total;
      
      // std::cerr << "normal nonbonded der = " << normal_nonbonded_der << std::endl;

      CHECK_APPROX_EQUAL(normal_nonbonded_der, nonbonded_der, 0.0000001, res);
      RESULT(res, total);

    }
    catch (std::string s){
      std::cout << s << std::endl;
      return 1;
    }

    return total;
}
