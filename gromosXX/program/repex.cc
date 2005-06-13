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
#include <unistd.h>

#include <io/configuration/out_configuration.h>

#include <gsl/gsl_rng.h>
#include <util/replica_exchange.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#pragma hdrstop

#include "BUILD_NUMBER"

void print_title(io::Argument &args, bool color = false);

int main(int argc, char *argv[]){

  const double start = util::now();

  util::Known knowns;
  knowns << "topo" << "conf" << "input" << "verb" << "pttopo"
	 << "trj" << "fin" << "trv" << "trf" << "tre" << "trg"
	 << "bae" << "bag" << "posres" <<"distrest" << "jval"
	 << "rep" << "master" << "slave"
	 << "version";
  
  
  std::string usage;
  util::get_usage(knowns, usage, argv[0]);

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

#ifdef XXMPI

  MPI_Init(&argc, &argv);

  if (args.count("master") >= 0){

    // the master
    util::Replica_Exchange_Master rep_master;

    std::cout << "starting master thread" << std::endl;
    rep_master.run(args);
    
  }
  else{
    
    // and the slaves
    util::Replica_Exchange_Slave rep_slave;

    rep_slave.run(args);

  }
  
#else
  
  std::cout << "MPI required! (use configure --enable-mpi)" << std::endl;
  return 1;

#endif
  
  return 0;
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
    std::cout << "\n\nGromosXX 0.2.1 development\033[22;0m\n\n"
	      << "26th October 2004\n";
  }
  else
    std::cout << "\n\nGromosXX 0.2.1 development\n\n"
	      << "26th October 2004\n";
  
  std::cout << "build date    " << BUILD_DATE << "\n"
	    << "build number  " << BUILD_NUMBER << "\n\n";
  
#ifdef NDEBUG
  std::cout << "standard library debugging disabled.\n";
#else
  std::cout << "standard library debugging enabled.\n";
#endif

  // some omp stuff
#ifdef XXMPI

  std::cout << "MPI code enabled\n"
	    << "\tdistributed memory parallelization\n"
	    << "\twww.mpi-group.org\n\n"
	    << std::endl;

  std::cout << "Replica Exchange Method\n\t";

  if (args.count("master") >= 0)
    std::cout << "master process";
  else
    std::cout << "slave process";
    
#endif
  
  std::cout << "\nGruppe fuer Informatikgestuetzte Chemie\n"
	    << "Professor W. F. van Gunsteren\n"
	    << "Swiss Federal Institute of Technology\n"
	    << "Zuerich\n\n"
	    << "Bugreports to http://www.igc.ethz.ch:5555\n\n";

}
