/**
 * @file script.cc
 * example for a simple md script
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction.h>
#include <interaction/interaction_types.h>
#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

#include <time.h>

#include <io/instream.h>
#include <util/parse_tcouple.h>
#include <io/blockinput.h>
#include <io/topology/in_topology.h>
#include <io/configuration/out_configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>
#include <algorithm/pressure/pressure_calculation.h>
#include <algorithm/pressure/berendsen_barostat.h>

#include <math/periodicity.h>
#include <math/volume.h>

#include <algorithm/constraints/shake.h>
#include <algorithm/constraints/perturbed_shake.h>
#include <algorithm/constraints/lincs.h>
#include <algorithm/constraints/flexible_constraint.h>
#include <algorithm/constraints/perturbed_flexible_constraint.h>

#include <interaction/forcefield/create_forcefield.h>

#include <util/create_simulation.h>
#include <algorithm/create_md_sequence.h>

#include <time.h>

#include <config.h>

#include <sstream>
#include <iomanip>

int main(int argc, char *argv[]){

  char *knowns[] = 
    {
      "topo", "pttopo", "conf", "input", "fin", "verb", "limit", "fcon", "stepsize"
    };
    
  int nknowns = 9;
    
  std::string usage = argv[0];
  usage += "\n\t@topo        <topology>\n";
  usage += "\t@pttopo      <perturbation topology>\n";
  usage += "\t@conf        <configuration>\n";
  usage += "\t@input       <input file>\n";
  usage += "\t@fin         <final configuration>\n";
  usage += "\t@fcon        <flexible constraints trajectory>\n";
  usage += "\t[@limit      <convergence limit (default 0.001)>]\n";
  usage += "\t[@stepsize   <step size (default 0.1)>]\n";
  usage += "\t[@verb   <verbosity>]\n";

  io::Argument args;
  if (args.parse(argc, argv, nknowns, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }
    
  double limit = 0.001;
  if (args.count("limit") > 0){
    std::istringstream is(args["limit"]);
    is >> limit;
  }

  double stepsize = 0.1;
  if(args.count("stepsize") > 0){
    std::istringstream is(args["stepsize"]);
    is >> stepsize;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);

  std::cout << "\n\n"
	    << "topology :      " << args["topo"] << "\n"
	    << "perturbation :  " << args["pttopo"] << "\n"
	    << "input :         " << args["input"] << "\n"
	    << "configuration : " << args["conf"] << "\n"
	    << std::endl;


  util::simulation_struct sys;
  io::In_Topology in_topo;

  if (util::create_simulation(args["topo"],
			      args["pttopo"],
			      args["conf"],
			      args["input"],
			      sys,
			      in_topo
			      )
      != 0){
    std::cerr << "creating simulation failed!" << std::endl;
    return 1;
  }
  
  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
	
  if (interaction::create_g96_forcefield(*ff, 
					 sys.topo,
					 sys.sim,
					 sys.conf,
					 in_topo
					 )
      != 0){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }

  ff->init(sys.topo, sys.conf, sys.sim);

  // need leap frog and flex shake
  algorithm::Leap_Frog_Velocity leap_frog_v;
  algorithm::Leap_Frog_Position leap_frog_p;
      
  algorithm::Perturbed_Flexible_Constraint  
    perturbed_flexible_constraints(sys.sim.param().constraint.solute.shake_tolerance,
				   1000,
				   ff);

  in_topo.read_harmonic_bonds(perturbed_flexible_constraints.parameter());
    
  // create output files...
  io::Out_Configuration traj("GromosXX\n");
  traj.title("GromosXX\n" + sys.sim.param().title);
  traj.init(args, sys.sim.param());

  std::ofstream fcon(args["fcon"].c_str());
  traj._print_title("GromosXX\n" + sys.sim.param().title,
		    "flexible constraints",
		    fcon);
  
  io::messages.display();

  // ==================================================
  bool convergence = false;
  int error;
  int iteration = 0;
  
  while (!convergence){
    
    ++iteration;
    std::cout << "=== iteration " << std::setw(3) << iteration
	      << " ==================================================" << std::endl;

    std::cout << "calculating free force..." << std::endl;
    
    ff->apply(sys.topo, sys.conf, sys.sim);
    
    std::cout << "unconstrained positions..." << std::endl;
    leap_frog_v.apply(sys.topo, sys.conf, sys.sim);
    leap_frog_p.apply(sys.topo, sys.conf, sys.sim);

    std::cout << "calculating constraint lengths..." << std::endl;
    perturbed_flexible_constraints.calc_distance(sys.topo, sys.conf, sys.sim);

    std::cout << "resetting positions..." << std::endl;
    sys.conf.current().pos = sys.conf.old().pos;
    
    std::cout << "flexible shaking positions..." << std::endl;
    perturbed_flexible_constraints.solute(sys.topo, sys.conf, sys.sim, error);
    if (error) return 1;
    
    perturbed_flexible_constraints._store_lengths(sys.conf);

    sys.conf.current().vel = sys.conf.old().vel;

    double max_dev = 0.0, av_dev = 0.0;
    for(unsigned int i=0; i<sys.topo.num_solute_atoms(); ++i){
      double dev = math::abs2(sys.conf.current().pos(i) - sys.conf.old().pos(i));
      av_dev += dev;
      if (dev > max_dev) max_dev = dev;
    }
    av_dev /= sys.topo.num_solute_atoms();
    av_dev = sqrt(av_dev);
    max_dev = sqrt(max_dev);

    std::cout << "average deviation : " << av_dev << std::endl;
    std::cout << "maximum deviation : " << max_dev << std::endl;

    sys.conf.exchange_state();

    for(int i=0; i<sys.topo.num_solute_atoms(); ++i){

      sys.conf.current().pos(i) = 
	sys.conf.current().pos(i) +
	stepsize * 
	(sys.conf.old().pos(i) - sys.conf.current().pos(i));
      
    }

    for(int k=0; k<sys.conf.special().flexible_constraint.flexible_vel.size(); ++k)
      sys.conf.special().flexible_constraint.flexible_vel[k] = 0.0;
  
    sys.conf.current().vel = sys.conf.old().vel;
  
    traj._print_timestep(sys.sim, fcon);
    traj._print_flexv(sys.conf, sys.topo, fcon);

    sys.sim.time() += sys.sim.time_step_size();
    ++sys.sim.steps();
  
    if (max_dev < limit)
      convergence = true;
  }

  std::cout << "writing final configuration" << std::endl;
  traj.write(sys.conf, sys.topo, sys.sim, io::final);


  // ==================================================

  return 0;
}

