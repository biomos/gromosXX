
#include <util/stdheader.h>

#include <topology/core/core.h>
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

#include <io/read_input.h>
#include <io/read_special.h>

#include <io/configuration/out_configuration.h>

int main(int argc, char *argv[]){
  try{
    
    char *knowns[] = 
      {
        "topo", "conf", "input", "verb", "alg", "pttopo",
        "trj", "fin", "trv", "trf", "tre", "trg", "print", "trp",
	"posres"
      };
    
    int nknowns = 15;
    
    std::string usage = argv[0];
    usage += "\n\t@topo    <topology>\n";
    usage += "\t[@pttopo <perturbation topology>]\n";
    usage += "\t@conf    <starting configuration>\n";
    usage += "\t@input   <input>\n";
    usage += "\t@trj     <trajectory>\n";
    usage += "\t@fin     <final structure>\n";
    usage += "\t[@trv    <velocity trajectory>]\n";
    usage += "\t[@trf    <force trajectory>]\n";
    usage += "\t[@tre    <energy trajectory>]\n";
    usage += "\t[@trg    <free energy trajectory>]\n";
    usage += "\t[@posres <position restraints data>]\n";
    usage += "\t[@alg    <RK|LF>]\n";
    usage += "\t[@print  <pairlist/force>]\n";
    usage += "\t[@trp    <print file>]\n";
    usage += "\t[@verb   <[module:][submodule:]level>]\n";

    io::Argument args(argc, argv, nknowns, knowns, usage);

    // parse the verbosity flag and set debug levels
    util::parse_verbosity(args);

    // create the simulation classes
    topology::Topology topo;
    configuration::Configuration conf;
    algorithm::Algorithm_Sequence md;
    simulation::Simulation sim;

    io::Out_Configuration traj("GromosXX");
    
    io::read_input(args, topo, conf, sim,  md);
    io::read_special(args, topo, conf, sim);

    if (args.count("fin") > 0)
      traj.final_configuration(args["fin"]);
    else throw std::string("argument fin for final configuration required!");
    if (args.count("trj") > 0)
      traj.trajectory(args["trj"], sim.param().write.position);
    else if (sim.param().write.position)
      throw std::string("write trajectory but no trj argument");
    if (args.count("trv") > 0)
      traj.velocity_trajectory(args["trv"], sim.param().write.velocity);
    else if (sim.param().write.velocity)
      throw std::string("write velocity trajectory but no trv argument");
    if (args.count("trf") > 0)
      traj.force_trajectory(args["trf"], 1);
    //else if (sim.param().write.force)
    //  throw std::string("write force trajectory but no trf argument");
    if (args.count("tre") > 0)
      traj.energy_trajectory(args["tre"], sim.param().write.energy);
    else if (sim.param().write.energy)
      throw std::string("write energy trajectory but no tre argument");
    if (args.count("trg") > 0)
      traj.free_energy_trajectory(args["trg"], sim.param().write.free_energy);
    else if (sim.param().write.free_energy)
      throw std::string("write free energy trajectory but no trg argument");

    std::cout << "\nMESSAGES FROM INITIALIZATION\n";
    io::messages.display(std::cout);
    io::messages.clear();

    std::cout.precision(5);
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    std::cout << "\nenter the next level of molecular "
	      << "dynamics simulations\n" << std::endl;

    double end_time = sim.param().step.t0 + 
      sim.time_step_size() * sim.param().step.number_of_steps;
    

    std::cout << "MD loop\n\tstart t=" << sim.time() 
	      << "\tend t=" << end_time << "\n" << std::endl;

    int error;
    
    while(sim.time() < end_time){
      // std::cout << "\tmd step " << sim.time() << std::endl;
      
      traj.write(conf, topo, sim, io::reduced);

      if ((error = md.run(topo, conf, sim))){

	if (error == E_MINIMUM_REACHED){
	  error = 0; // clear error condition
	  break;
	}
	
	std::cout << "\nError during MD run!\n\n";
	// try to save the final structures...
	break;
      }

      // update the energies
      conf.old().energies.calculate_totals();
      conf.current().energy_averages.update(conf.old().energies,
					    conf.old().energy_averages,
					    sim.time_step_size());
      // perturbed energy derivatives
      if (sim.param().perturbation.perturbation){
	conf.old().perturbed_energy_derivatives.calculate_totals();
	conf.current().perturbed_energy_derivative_averages.update
	  (conf.old().perturbed_energy_derivatives,
	   conf.old().perturbed_energy_derivative_averages,
	   sim.time_step_size());
      }

      traj.print(topo, conf, sim);

      sim.time() += sim.time_step_size();
      ++sim.steps();

    }
    

    std::cout << "writing final configuration" << std::endl;
    
    traj.write(conf, topo, sim, io::final);
    traj.print_final(topo, conf, sim);
    
    std::cout << "\nMESSAGES FROM SIMULATION\n";
    io::messages.display(std::cout);

    if (error)
      std::cout << "\nErrors encountered during run - check above!\n" << std::endl;
    else
    std::cout << "\nGromosXX finished successfully\n" << std::endl;
    
  }
  catch (std::string s){
    io::messages.display();
    std::cerr << "there was something wrong:\n" << s << std::endl;
    return 1;
  }
  
    return 0;
}

