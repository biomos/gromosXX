/**
 * @file rasn.t.cc
 * tests using rasn
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

int main(int argc, char* argv[])
{
  
    char *knowns[] = 
      {
        "topo", "conf", "input", "pttopo", "verb"
      };
    
    int nknowns = 5;
    
    std::string usage = argv[0];
    usage += "\n\t[@topo    <topology>]\n";
    usage += "\t[@pttopo <perturbation topology>]\n";
    usage += "\t[@conf    <starting configuration>]\n";
    usage += "\t[@input   <input>]\n";
    usage += "\t[@verb   <[module:][submodule:]level>]\n";

    try{
      io::Argument args(argc, argv, nknowns, knowns, usage, true);

      // parse the verbosity flag and set debug levels
      util::parse_verbosity(args);
      
      std::string stopo, spttopo, sconf, sinput;
      bool quiet = true;

      if (args.count("verb") != -1) quiet = false;
      
      if(args.count("topo") == 1)
	stopo = args["topo"];
      else
	GETFILEPATH(stopo, "rasn.topo", "src/check/data/");

      if(args.count("pttopo") == 1)
	spttopo = args["pttopo"];
      else
	GETFILEPATH(spttopo, "rasn.pttopo", "src/check/data/");
      
      if(args.count("conf") == 1)
	sconf = args["conf"];
      else
	GETFILEPATH(sconf, "rasn.conf", "src/check/data/");

      if(args.count("input") == 1)
	sinput = args["input"];
      else
	GETFILEPATH(sinput, "rasn.in", "src/check/data/");

      if (!quiet)
	std::cout << "\n\n"
		  << "topology :      " << stopo << "\n"
		  << "perturbation :  " << spttopo << "\n"
		  << "input :         " << sinput << "\n"
		  << "configuration : " << sconf << "\n"
		  << std::endl;

      util::simulation_struct rasn_sim;
      io::In_Topology in_topo;

      in_topo.quiet = quiet;
      
      if (util::create_simulation(stopo,
				  spttopo,
				  sconf,
				  sinput,
				  rasn_sim,
				  in_topo,
				  quiet
				  )
	  != 0){
	std::cerr << "creating simulation failed!" << std::endl;
	return 1;
      }

      if (false){
	// create a forcefield
	interaction::Forcefield *ff = new interaction::Forcefield;
	
	if (interaction::create_g96_forcefield(*ff, 
					       rasn_sim.topo,
					       rasn_sim.sim,
					       rasn_sim.conf,
					       in_topo,
					       quiet)
	  != 0){
	  std::cerr << "creating forcefield failed!" << std::endl;
	  return 1;
	}
      }
      // vtune: run the thing...
      else{
	
	algorithm::create_md_sequence(rasn_sim.md,
				      rasn_sim.topo,
				      rasn_sim.conf,
				      rasn_sim.sim,
				      in_topo);

	double end_time = rasn_sim.sim.param().step.t0 + 
	  rasn_sim.sim.time_step_size() * 
	  rasn_sim.sim.param().step.number_of_steps;
    
	int error;

	while(rasn_sim.sim.time() < end_time){
      
	  if ((error = rasn_sim.md.run(rasn_sim.topo, 
				       rasn_sim.conf,
				       rasn_sim.sim))){

	    std::cout << "\nError during MD run!\n" << std::endl;
	    // try to save the final structures...
	    break;
	  }

	  // update the energies
	  rasn_sim.conf.old().energies.calculate_totals();
	  rasn_sim.conf.current().energy_averages.
	    update(rasn_sim.conf.old().energies,
		   rasn_sim.conf.old().energy_averages,
		   rasn_sim.sim.time_step_size());

      // perturbed energy derivatives
	  if (rasn_sim.sim.param().perturbation.perturbation){
	    rasn_sim.conf.old().perturbed_energy_derivatives.calculate_totals();

	    rasn_sim.conf.current().perturbed_energy_derivative_averages.
	      update(rasn_sim.conf.old().perturbed_energy_derivatives,
		     rasn_sim.conf.old().perturbed_energy_derivative_averages,
		     rasn_sim.sim.time_step_size(),
		     rasn_sim.sim.param().perturbation.dlamt);
	  }
	  
	  rasn_sim.sim.time() +=  rasn_sim.sim.time_step_size();
	  ++ rasn_sim.sim.steps();

	}
	
      }
    }
    catch (std::string s){
      std::cout << s << std::endl;
      return 1;
    }

    return 0;
}
