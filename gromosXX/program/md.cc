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
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>

#ifdef OMP
#include <omp.h>
#endif

#include <sys/utsname.h>

#pragma hdrstop

#include "BUILD_NUMBER"

void print_title(bool color = false);

int main(int argc, char *argv[]){

  const double start = util::now();

  util::Known knowns;
  knowns << "topo" << "conf" << "input" << "verb" << "pttopo"
	 << "trj" << "fin" << "trv" << "trf" << "tramd" << "tre" << "trg"
	 << "bae" << "bag" << "posres" <<"distrest" << "dihrest" << "jval"
	 << "anatrj" << "print" 
	 << "version";
  
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);
  usage += "#\n\n";

  io::Argument args;

  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }
    
  if (args.count("version") >= 0){
    print_title(true);
    return 0;
  }
  else print_title();
    
  // parse the verbosity flag and set debug levels
  if (util::parse_verbosity(args)){
    std::cerr << "could not parse verbosity argument" << std::endl;
    return 1;
  }

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  simulation::Simulation sim;
  algorithm::Algorithm_Sequence md;

  io::Out_Configuration traj("GromosXX\n");

  if (io::read_input(args, topo, conf, sim,  md)){
    io::messages.display(std::cout);
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }

  traj.title("GromosXX\n" + sim.param().title);

  // create output files...
  traj.init(args, sim.param());

  // initialises all algorithms (and therefore also the forcefield)
  md.init(topo, conf, sim);

  std::cout << "\nMESSAGES FROM INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }
    
  io::messages.clear();

  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    
  std::cout << "\nenter the next level of molecular "
	    << "dynamics simulations\n" << std::endl;

    
  int percent = 0;
  double end_time = sim.time() + 
    sim.time_step_size() * (sim.param().step.number_of_steps - 1);
    
  std::cout << "==================================================\n"
	    << " MAIN MD LOOP\n"
	    << "==================================================\n"
	    << std::endl;

  int error;

  const double init_time = util::now() - start;
    
  while(sim.time() < end_time + math::epsilon){
      
    traj.write(conf, topo, sim, io::reduced);

    // run a step
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
      else { 
	// try to print energies anyway
	// if (error == E_NAN){
	io::print_ENERGY(traj.output(), conf.current().energies,
			 topo.energy_groups(),
			 "OLDERROR", "OLDERR_");
	
	io::print_ENERGY(traj.output(), conf.old().energies, 
			 topo.energy_groups(),
			 "ERROR", "ERR_");
      }

      std::cout << "\nError during MD run!\n" << std::endl;
      std::cout << "\tat step " << sim.steps() << " (time " << sim.time() << ")\n" << std::endl;
      // try to save the final structures...
      break;
    }

    traj.print(topo, conf, sim);

    sim.time() += sim.time_step_size();
    ++sim.steps();
    
    
    if ((sim.param().step.number_of_steps / 10 > 0) &&
	(sim.steps() % (sim.param().step.number_of_steps / 10) == 0)){
      ++percent;
      const double spent = util::now() - start;
      const int hh = int(spent / 3600);
      const int mm = int((spent - hh * 3600) / 60);
      const int ss = int(spent - hh * 3600 - mm * 60);

      std::cerr << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
      std::cout << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
      std::cerr << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;
      std::cout << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;

      const double eta_spent = spent / sim.steps() * sim.param().step.number_of_steps - spent;
      const int eta_hh = int(eta_spent / 3600);
      const int eta_mm = int((eta_spent - eta_hh * 3600) / 60);
      const int eta_ss = int(eta_spent - eta_hh * 3600 - eta_mm * 60);
      
      std::cerr << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
      std::cout << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
    }
  }
    
  std::cout << "writing final configuration" << std::endl;
    
  traj.write(conf, topo, sim, io::final);
  traj.print_final(topo, conf, sim);
    
  std::cout << "\nMESSAGES FROM SIMULATION\n";
  io::message::severity_enum err_msg = io::messages.display(std::cout);

  std::cout << "\n\n";
    
  md.print_timing(std::cout);

  std::cout << "Overall time used:\t" << util::now() - start << "\n"
	    << "(initialization took " << init_time << ")\n\n";

  const time_t time_now = time_t(util::now());
  std::cout << ctime(&time_now) << "\n\n";
    
  if (error){
    std::cout << "\nErrors encountered during run - check above!\n" << std::endl;
    return 1;
  }
  else if(err_msg > io::message::notice){
    std::cout << "\nGromosXX finished. "
	      << "Check the messages for possible problems during the run."
	      << std::endl;
    return 0;
  }
  else{
    
    std::cout << "\nGromosXX finished successfully\n" << std::endl;
  }
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////////////////////////////////

void print_title(bool color)
{
  if (color){

#ifdef NDEBUG
#ifndef BZDEBUG
    std::cout << "\033[1;32m";
#else
    std::cout << "\033[1;31m";
#endif
#else
    std::cout << "\033[1;31m";
#endif
    std::cout << "\n\nGromosXX 0.3.0 development\033[22;0m\n\n"
	      << "12 February 2008\n";
  }
  else
    std::cout << "\n\nGromosXX 0.3.0 development\n\n"
	      << "12 February 2008\n";
  
  std::cout << "build date    " << BUILD_DATE << "\n"
	    << "build number  " << BUILD_NUMBER << "\n\n";
  
#ifdef NDEBUG
  std::cout << "standard library debugging disabled.\n";
#else
  std::cout << "standard library debugging enabled.\n";
#endif

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
  
  std::cout << "\nGruppe fuer Informatikgestuetzte Chemie\n"
	    << "Professor W. F. van Gunsteren\n"
	    << "Swiss Federal Institute of Technology\n"
	    << "Zuerich\n\n"
	    << "Bugreports to https://gromos/svn/trac/gromosXXc++\n\n";

  struct utsname sysinf;
  if (uname(&sysinf) != -1){
    std::cout << "running on"
	      << "\n\t" << sysinf.nodename
	      << "\n\t" << sysinf.sysname
	      << " " << sysinf.release
	      << " " << sysinf.version
	      << " " << sysinf.machine
	      << "\n\n";
  }
}
