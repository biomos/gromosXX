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
#include <algorithm/constraints/position_constraints.h>
#include <algorithm/integration/energy_calculation.h>
#include <algorithm/integration/leap_frog.h>

#include <io/blockinput.h>
#include <io/instream.h>
#include <io/configuration/inframe.h>
#include <io/configuration/in_configuration.h>

#include <algorithm/integration/analyze.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/nosehoover_thermostat.h>
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
				  simulation::Simulation & sim,
				  io::In_Topology &it,
				  std::ostream & os)
{

  // analyze trajectory:
  // overwrite current configuration with trajectory data
  if (sim.param().analyze.analyze){
    algorithm::Analyze_Step * as = 
      new algorithm::Analyze_Step(sim.param().analyze.trajectory);
    md_seq.push_back(as);
  }

  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
  interaction::create_g96_forcefield(*ff, topo, sim, it, os);
  
  //==================================================
  // construct the md algorithm
  //==================================================

  // center of mass removal
  if (sim.param().centreofmass.remove_trans ||
      sim.param().centreofmass.remove_rot){

    algorithm::Remove_COM_Motion * rcom =
      new algorithm::Remove_COM_Motion(os);
      
    md_seq.push_back(rcom);
  }
  
  // add the forcefield
  md_seq.push_back(ff);

  // position constraints?
  if (sim.param().posrest.posrest == 3){
    algorithm::Position_Constraints * pc = new algorithm::Position_Constraints;
    md_seq.push_back(pc);
  }

  // energy minimisation or MD?
  if (sim.param().minimise.ntem == 1){
    algorithm::Steepest_Descent * sd = new algorithm::Steepest_Descent;
    md_seq.push_back(sd);
  }
  else{

    if(sim.param().integrate.method == simulation::integrate_leap_frog){
      md_seq.push_back(new algorithm::Leap_Frog_Velocity);
    }
    else{
      std::cout << "\tno integration (velocities) selected!\n";
    }

    // temperature scaling? -> has to be done before temperature calculation!!!
    if (sim.param().multibath.couple){

      if (sim.param().multibath.nosehoover == 0){
	algorithm::Berendsen_Thermostat * tcoup =
	  new algorithm::Berendsen_Thermostat;
	md_seq.push_back(tcoup);
      }
      else if (sim.param().multibath.nosehoover >= 1){
	algorithm::NoseHoover_Thermostat *tcoup =
	  new algorithm::NoseHoover_Thermostat;
	md_seq.push_back(tcoup);
      }
    }
    
    if(sim.param().integrate.method == simulation::integrate_leap_frog){
      md_seq.push_back(new algorithm::Leap_Frog_Position);
    }
    else{
      std::cout << "\tno integration (position) selected!\n";
    }
  }
  
  // CONSTRAINTS
  create_constraints(md_seq, topo, sim, it);

  // temperature calculation (always!)
  {
    algorithm::Temperature_Calculation * tcalc =
      new algorithm::Temperature_Calculation;

    DEBUG(7, tcalc->name);
    md_seq.push_back(tcalc);

  }
  
  // pressure calculation?
  io::print_PCOUPLE(os, sim.param().pcouple.calculate,
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
    
    os << "PERTURBATION\n"
       << "\tlambda         : " << sim.param().perturbation.lambda << "\n"
       << "\texponent       : " << sim.param().perturbation.lambda_exponent << "\n"
       << "\tdlambda        : " << sim.param().perturbation.dlamt << "\n"
       << "\tscaling        : ";

    if (sim.param().perturbation.scaling){
      if (sim.param().perturbation.scaled_only)
	os << "perturbing only scaled interactions\n";
      else
	os << "on\n";
    }
    else
      os << "off\n";
    os << "END\n";
  }
  else
    os << "PERTURATION OFF\n";
  
  // total energy calculation and energy average update
  {
    algorithm::Energy_Calculation * ec = 
      new algorithm::Energy_Calculation();
    md_seq.push_back(ec);
  }

  return 0;

}

