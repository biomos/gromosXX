/**
 * @file create_md_sequence.cc
 */

#include <stdheader.h>
#include <fstream>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>
#include <interaction/interaction_types.h>

#include <io/argument.h>
#include <io/blockinput.h>
#include <io/instream.h>
#include <io/topology/in_topology.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>
#include <algorithm/pressure/pressure_calculation.h>
#include <algorithm/pressure/berendsen_barostat.h>

#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <math/periodicity.h>

#include <algorithm/constraints/create_constraints.h>
#include <algorithm/constraints/remove_com_motion.h>

#include <algorithm/integration/slow_growth.h>
#include <algorithm/integration/steepest_descent.h>

#include <io/print_block.h>

#include <algorithm/create_md_sequence.h>


#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

int algorithm::create_md_sequence(algorithm::Algorithm_Sequence &md_seq,
				  topology::Topology &topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  io::In_Topology &it)
{

  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
  interaction::create_g96_forcefield(*ff, topo, sim, conf, it);
  
  //==================================================
  // construct the md algorithm
  //==================================================

  // center of mass removal
  algorithm::Remove_COM_Motion * rcom =
    new algorithm::Remove_COM_Motion;
  
  if (sim.param().centreofmass.remove_trans ||
      sim.param().centreofmass.remove_rot)
    md_seq.push_back(rcom);
  else{
    rcom->apply(topo, conf, sim);
    delete rcom;
  }
  
  // add the forcefield
  md_seq.push_back(ff);

  // energy minimisation or MD?
  if (sim.param().minimise.ntem == 1){
    algorithm::Steepest_Descent * sd = new algorithm::Steepest_Descent;
    md_seq.push_back(sd);
  }
  else{
    md_seq.push_back(new algorithm::Leap_Frog_Velocity);

    // temperature scaling? -> has to be done before temperature calculation!!!
    if (sim.param().multibath.couple){
      algorithm::Berendsen_Thermostat * tcoup =
	new algorithm::Berendsen_Thermostat;
      md_seq.push_back(tcoup);
    }
    
    md_seq.push_back(new algorithm::Leap_Frog_Position);
  }
  
  // CONSTRAINTS
  create_constraints(md_seq, topo, conf, sim, it);

  // temperature calculation (always!)
  {
    algorithm::Temperature_Calculation * tcalc =
      new algorithm::Temperature_Calculation;
    // calculate initial temperature
    tcalc->apply(topo, conf, sim);

    // do we scale the initial temperatures?
    if (sim.param().replica.scale){
      std::cout << "\tscale initial velocities (replica exchange)\n";
      
      algorithm::Berendsen_Thermostat tcoup;
      tcoup.calc_scaling(topo, conf, sim, true);
      tcoup.scale(topo, conf, sim);

      tcalc->apply(topo, conf, sim);
    }

    io::print_MULTIBATH_COUPLING(std::cout, sim.multibath());

    io::print_DEGREESOFFREEDOM(std::cout, sim.multibath());
    
    io::print_MULTIBATH(std::cout, sim.multibath(),
			conf.old().energies,
			"INITIAL TEMPERATURES");


    DEBUG(7, tcalc->name);
    md_seq.push_back(tcalc);

  }
  
  // pressure calculation?
  io::print_PCOUPLE(std::cout, sim.param().pcouple.calculate,
		    sim.param().pcouple.scale,
		    sim.param().pcouple.pres0,
		    sim.param().pcouple.compressibility,
		    sim.param().pcouple.tau,
		    sim.param().pcouple.virial);
  
  if (sim.param().pcouple.calculate){
    algorithm::Pressure_Calculation * pcalc =
      new algorithm::Pressure_Calculation;
    md_seq.push_back(pcalc);
  }

  //  pressure scaling
  if (sim.param().pcouple.scale != math::pcouple_off){
    algorithm::Berendsen_Barostat * pcoup =
      new algorithm::Berendsen_Barostat;
    md_seq.push_back(pcoup);
  }

  // slow growth
  if (sim.param().perturbation.perturbation){
    if (sim.param().perturbation.dlamt){
      algorithm::Slow_Growth *sg =
	new algorithm::Slow_Growth;
      md_seq.push_back(sg);
    }
    
    std::cout << "PERTURBATION\n"
	      << "\tlambda         : " << sim.param().perturbation.lambda << "\n"
	      << "\texponent       : " << sim.param().perturbation.lambda_exponent << "\n"
	      << "\tdlambda        : " << sim.param().perturbation.dlamt << "\n"
	      << "\tscaling        : ";

    if (sim.param().perturbation.scaling){
      if (sim.param().perturbation.scaled_only)
	std::cout << "perturbing only scaled interactions\n";
      else
	std::cout << "on\n";
    }
    else
      std::cout << "off\n";
    std::cout << "END\n";
  }
  else
    std::cout << "PERTURATION OFF\n";
  
  return 0;

}

