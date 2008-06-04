/**
 * @file repex.cc
 * the main md program for replica exchange simulations
 */

#ifdef XXMPI
#include <mpi.h>
#endif

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
#include <unistd.h>

#include <io/configuration/out_configuration.h>

#include <gsl/gsl_rng.h>
#include <util/replica_exchange.h>

#ifdef OMP
#include <omp.h>
#endif

#include <sys/utsname.h>

#pragma hdrstop

#include "BUILD_NUMBER"

void print_title(io::Argument &args, bool color = false);

int main(int argc, char *argv[]){

#ifdef REPEX

  util::Known knowns;
  knowns << "topo" << "cg_topo" << "conf" << "cg_conf"
     << "input" << "cg_input" << "verb" << "pttopo" << "cg_pttopo"
	 << "trj" << "cg_trj" << "fin" << "cg_fin" << "trv" << "trf" << "trs" << "tre" << "re" << "trg"
	 << "bae" << "bag" << "posres" <<"distrest" << "jval" << "friction"
	 << "rep" << "master" << "slave" << "control"
	 << "version" << "gzip";
  
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);
  usage += "#\n\n";

  io::Argument args;

  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }
    
  if (args.count("version") >= 0){
    print_title(args, true);
    return 0;
  }
  else print_title(args);
    
  // parse the verbosity flag and set debug levels
  if (util::parse_verbosity(args)){
    std::cerr << "could not parse verbosity argument" << std::endl;
    return 1;
  }

  // #ifdef XXMPI

  // MPI_Init(&argc, &argv);
  // MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  
  int error = 0;

  if (args.count("master") >= 0){

    // the master
    util::Replica_Exchange_Master rep_master;

    std::cout << "starting master thread" << std::endl;
    error = rep_master.run(args);
    
  }
  else if (args.count("slave") >= 0){
    
    // and the slave
    std::cout << "repex: starting slave" << std::endl;
    util::Replica_Exchange_Slave rep_slave;

    error = rep_slave.run(args);

  }
  else if (args.count("control") >= 0){
    std::cout << "repex: starting control!" << std::endl;
    util::Replica_Exchange_Control rep_control;
    error = rep_control.run(args);
  }
  else{
    std::cout << "repex: either @master @slave or @interactive required" << std::endl;
  }

  if (error){
    std::cerr << "errors during repex run! (error code " << error << ")\n"
	      << std::endl;

    io::messages.display();
    
    return 1;
  }
  
  return 0;

#else
  std::cout << argv[0] << " needs REPEX to run\n\tuse --enable-repex at configure\n" << std::endl;
  return 1;
#endif
  
}

////////////////////////////////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////////////////////////////////

void print_title(io::Argument & args, bool color)
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

#ifdef XXMPI

  std::cout << "MPI code enabled\n"
	    << "\tdistributed memory parallelization\n"
	    << "\twww.mpi-group.org\n\n"
	    << std::endl;
#endif

  std::cout << "Replica Exchange Method\n\t";

  if (args.count("master") >= 0)
    std::cout << "master process";
  else
    std::cout << "slave process";
    
  
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
