/**
 * @file md.cc
 * the main md program
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

#include <io/read_input.h>
#include <io/read_special.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>

#ifdef OMP
#include <omp.h>
#endif

#pragma hdrstop

#include "BUILD_NUMBER"

int main(int argc, char *argv[]){

  const double start = util::now();
  // double init_start = util::now();

  try{
    
    char *knowns[] = 
      {
        "topo", "conf", "input", "verb", "pttopo",
        "trj", "fin", "trv", "trf", "tre", "trg", "print", "trp",
	"bae", "bag", "posres", "version"
      };
    
    int nknowns = 17;
    
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
    usage += "\t[@bae    <block averaged energy trajectory>]\n";
    usage += "\t[@bag    <block averaged free energy trajectory>]\n";    
    usage += "\t[@posres <position restraints data>]\n";
    usage += "\t[@print  <pairlist/force>]\n";
    usage += "\t[@trp    <print file>]\n";
    usage += "\t[@verb   <[module:][submodule:]level>]\n";
    usage += "\t[@version]\n";

    io::Argument args(argc, argv, nknowns, knowns, usage);
    
    if (args.count("version") >= 0){
      
#ifdef NDEBUG
#ifndef BZDEBUG
      std::cout << "\033[1;32m";
#else
      std::cout << "\033[1;31m";
#endif
#else
      std::cout << "\033[1;31m";
#endif
      std::cout << "\n\nGromosXX 0.1.3 development\033[22;0m\n\n"
		<< "1st June 2004\n";
    }
    else
      std::cout << "\n\nGromosXX 0.1.2 development\n\n"
		<< "20. January 2004\n";
    
    std::cout << "build date    " << BUILD_DATE << "\n"
	      << "build number  " << BUILD_NUMBER << "\n\n";
    
#ifdef NDEBUG
    std::cout << "standard library debugging disabled.\n";
#else
    std::cout << "standard library debugging enabled.\n";
#endif
#ifdef BZDEBUG
    std::cout << "Blitz debugging enabled.\n";
#else
    std::cout << "Blitz debugging disabled.\n";
#endif
    
    std::cout << "\nGruppe fuer Informatikgestuetzte Chemie\n"
	      << "Professor W. F. van Gunsteren\n"
	      << "Swiss Federal Institute of Technology\n"
	      << "Zuerich\n\n"
	      << "Bugreports to http://www.igc.ethz.ch:5555\n\n";
    
    if (args.count("version") >= 0){
      return 0;
    }
    
    // parse the verbosity flag and set debug levels
    util::parse_verbosity(args);

    // create the simulation classes
    topology::Topology topo;
    configuration::Configuration conf;
    algorithm::Algorithm_Sequence md;
    simulation::Simulation sim;

    io::Out_Configuration traj("GromosXX\n");

    io::read_input(args, topo, conf, sim,  md);

    io::read_special(args, topo, conf, sim);

    traj.title("GromosXX\n" + sim.param().title);

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
    if (args.count("bae") > 0)
      traj.block_averaged_energy(args["bae"], sim.param().write.block_average);
    else if (sim.param().write.block_average && sim.param().write.energy)
      throw std::string("write block averaged energy but no bae argument");
	if (sim.param().perturbation.perturbation){
			if (args.count("bag") > 0)
				traj.block_averaged_free_energy(args["bag"], sim.param().write.block_average);
			else if (sim.param().write.block_average && sim.param().write.free_energy)
				throw std::string("write block averaged free energy but no bag argument");
	}
    std::cout << "\nMESSAGES FROM INITIALIZATION\n";
    if (io::messages.display(std::cout) >= io::message::error){
      // exit
      std::cout << "\nErrors during initialization!\n" << std::endl;
      return 1;
    }
    
    io::messages.clear();

    std::cout.precision(5);
    std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
    std::cout << "\nenter the next level of molecular "
	      << "dynamics simulations\n" << std::endl;


    // some omp stuff
#ifdef OMP
    int nthreads, tid;
#pragma omp parallel private(nthreads, tid)
    {
      tid = omp_get_thread_num();
      if (tid == 0){
	nthreads = omp_get_num_threads();
	std::cout << "OpenMP code enabled\n"
		  << "\tshared memory parallelization\n"
		  << "\twww.openmp.org\n\n"
		  << "\tusing "
		  << omp_get_num_threads() << " threads\n"
		  << "\tthis can be adjusted by setting the\n"
		  << "\tOMP_NUM_THREADS environment variable\n"
		  << std::endl;
      }
   
    }
#endif
    
    double end_time = sim.param().step.t0 + 
      sim.time_step_size() * sim.param().step.number_of_steps;
    
    
    
    std::cout << "==================================================\n"
	      << " MAIN MD LOOP\n"
	      << "==================================================\n"
	      << std::endl;

    int error;

    const double init_time = util::now() - start;
    
    while(sim.time() < end_time){
      // std::cout << "\tmd step " << sim.time() << std::endl;
      
      traj.write(conf, topo, sim, io::reduced);

      if ((error = md.run(topo, conf, sim))){

	if (error == E_MINIMUM_REACHED){
	  conf.old().energies.calculate_totals();
	  traj.print_timestep(sim, traj.output());
	  io::print_ENERGY(traj.output(), conf.old().energies, 
			   topo.energy_groups(),
			   "MINIMUM ENERGY", "EMIN_");
	  
	  error = 0; // clear error condition
	  break;
	}

	std::cout << "\nError during MD run!\n" << std::endl;
	// try to save the final structures...
	break;
      }

      // update the energies
      if (conf.old().energies.calculate_totals()){
	std::cout << "\nError during MD run!\n" << std::endl;
	// try to save the final structures...
	break;	
      }

      // perturbed energy derivatives
      if (sim.param().perturbation.perturbation){
	conf.old().perturbed_energy_derivatives.calculate_totals();
      }
      
      conf.current().averages.apply(topo, conf, sim);
      
      traj.print(topo, conf, sim);

      sim.time() += sim.time_step_size();
      ++sim.steps();

    }
    
    std::cout << "writing final configuration" << std::endl;
    
    traj.write(conf, topo, sim, io::final);
    traj.print_final(topo, conf, sim);
    
    std::cout << "\nMESSAGES FROM SIMULATION\n";
    io::messages.display(std::cout);

    std::cout << "\n\n";
    
    md.print_timing(std::cout);

    std::cout << "Overall time used:\t" << util::now() - start << "\n"
	      << "(initialization took " << init_time << ")\n\n";

    const time_t time_now = time_t(util::now());
    std::cout << ctime(&time_now) << "\n\n";
    
    if (error)
      std::cout << "\nErrors encountered during run - check above!\n" << std::endl;
    else
    std::cout << "\nGromosXX finished successfully\n" << std::endl;
    
  }
  catch (std::string s){
    io::messages.display();
    std::cerr << "error encountered while running md:\n" << s << std::endl;
    return 1;
  }
  
    return 0;
}

