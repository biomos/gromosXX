
#include "../stdheader.h"
#include "../simulation/simulation.h"

#include "../io/argument.h"
#include "../io/usage.h"

#include <config.h>

void print_title(bool color = false);

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  
  Known knowns;
  knowns << "topo" << "conf" << "input" << "verb" << "pttopo"
	 << "trj" << "fin" << "trv" << "trf" << "tre" << "trg"
	 << "bae" << "bag" << "posres" <<"distrest" << "jval"
	 << "anatrj" << "print"
	 << "version";

  std::string usage;
  get_usage(knowns, usage, argv[0]);

  Argument args;

  if (args.parse(argc, argv, knowns)){
    std::cerr << usage << std::endl;
    return 1;
  }
    
  if (args.count("version") >= 0){
    int rank = MPI::COMM_WORLD.Get_rank();
    if (rank == 0)
      print_title(true);
    MPI_Finalize();
    return 0;
  }
  else print_title();

  MPI_Finalize();
  
  return 0;
}


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
    std::cout << "\n\nddmd 0.0.1 development\033[22;0m\n\n"
	      << "12th June 2005\n";
  }
  else
    std::cout << "\n\nddmd 0.0.1 development\n\n"
	      << "12th June 2005\n";
  
#ifdef NDEBUG
  std::cout << "standard library debugging disabled.\n";
#else
  std::cout << "standard library debugging enabled.\n";
#endif

  int nthreads = 0, tid = 0;
  #pragma omp parallel private(tid)
  {
    tid = omp_get_thread_num();
    if (tid == 0){
      nthreads = omp_get_num_threads();
      std::cout << "OpenMP code enabled\n"
		<< "\tshared memory parallelization\n"
		<< "\twww.openmp.org\n\n"
		<< "\tusing "
		<< nthreads << " threads\n"
		<< "\tthis can be adjusted by setting the\n"
		<< "\tOMP_NUM_THREADS environment variable\n"
		<< std::endl;

      std::cout << "MPI code enabled\n"
		<< "\tdistributed memory parallelization\n"
		<< "\tusing "
		<< MPI::COMM_WORLD.Get_size()
		<< " nodes\n\n";
    }
  }

  std::cout << "\nGruppe fuer Informatikgestuetzte Chemie\n"
	    << "Professor W. F. van Gunsteren\n"
	    << "Swiss Federal Institute of Technology\n"
	    << "Zuerich\n\n"
	    << "Bugreports to http://www.igc.ethz.ch:5555\n\n";
  
}
