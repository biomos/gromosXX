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
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>

#include <util/replica_exchange.h>

#ifdef REM
#include <omp.h>
#endif

#pragma hdrstop

#include "BUILD_NUMBER"

void print_title(bool color = false);

int main(int argc, char *argv[]){

  const double start = util::now();

  char *knowns[] = 
    {
      "topo", "conf", "input", "verb", "pttopo",
      "trj", "fin", "trv", "trf", "tre", "trg", "print", "trp",
      "bae", "bag", "posres", "jval", "rep", "version"
    };
    
  int nknowns = 19;
    
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
  usage += "\t[@jval   <jvalue restraints data>]\n";
  usage += "\t[@rep    <replica exchange final data>]\n";
  usage += "\t[@print  <pairlist/force>]\n";
  usage += "\t[@trp    <print file>]\n";
  usage += "\t[@verb   <[module:][submodule:]level>]\n";
  usage += "\t[@version]\n";

  io::Argument args;

  if (args.parse(argc, argv, nknowns, knowns)){
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

#ifdef REM
  int nthreads, tid;
#pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();

    if (tid == 0)
      nthreads = omp_get_num_threads();
  }
  
  std::vector<util::Replica_Exchange> rep(nthreads);
  
  util::replica_master = &rep[0];
  std::cerr << "starting REM on " << nthreads << " threads" <<  std::endl;

#pragma omp parallel private(tid)
  {

    tid = omp_get_thread_num();

    assert(rep.size() > tid);
    std::cerr << "running " << tid+1 << " of " << nthreads << std::endl;
    rep[tid].run(args, tid, nthreads);

  }

#else

  std::cout << "OpenMP required! (use configure --enable-rem)" << std::endl;
  return 1;

#endif
  
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
#ifdef REM
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
	    << "Bugreports to http://www.igc.ethz.ch:5555\n\n";

}
