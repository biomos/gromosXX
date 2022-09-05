/**
 * @file usage.cc
 * get usage string
 */

#include "../stdheader.h"
#include "usage.h"
#include <sys/utsname.h>

void util::get_usage(util::Known const &knowns, std::string &usage, std::string name)
{
  usage = "#\n# " + name + "\n\n";

  if (knowns.count("topo")){
    usage += "\t# topology data\n";
    usage += "\t@topo        filename\n\n";
  }
  
  if (knowns.count("pttopo")){
    usage += "\t# perturbation topology data\n";
    usage += "\t# @pttopo    filename\n\n";
  }
  
  if (knowns.count("conf")){
    usage += "\t# coordinates\n";
    usage += "\t@conf        filename\n\n";
  }
  
  if (knowns.count("input")){
    usage += "\t# input parameter\n";
    usage += "\t@input       filename\n\n";
  }

  if (knowns.count("fin")){
    usage += "\t# output final coordinates\n";
    usage += "\t@fin         filename\n\n";
  }
  
  if (knowns.count("trc")){
    usage += "\t# output coordinates trajectory\n";
    usage += "\t@trc         filename\n\n";
  }
  
  if (knowns.count("trv")){
    usage += "\t# output velocity trajectory\n";
    usage += "\t# @trv       filename\n\n";
  }
  
  if (knowns.count("trf")){
    usage += "\t# output force trajectory\n";
    usage += "\t# @trf       filename\n\n";
  }

  if (knowns.count("trs")){
    usage += "\t# output special trajectory\n";
    usage += "\t# @trs       filename\n\n";
  }  
  
  if (knowns.count("tre")){
    usage += "\t# output energy trajectory\n";
    usage += "\t# @tre       filename\n\n";
  }

  if (knowns.count("re")){
    usage += "\t# output replica energy trajectory (per switch)\n";
    usage += "\t# @re        filename\n\n";
  }
  
  if (knowns.count("bae")){
    usage += "\t# output block averaged energy trajectory\n";
    usage += "\t# @bae       filename\n\n";
  }
  
  if (knowns.count("trg")){
    usage += "\t# output free-energy trajectory\n";
    usage += "\t# @trg       filename\n\n";
  }
  
  if (knowns.count("bag")){
    usage += "\t# output block averaged free-energy trajectory\n";
    usage += "\t# @bag       filename\n\n";    
  }
  
  if (knowns.count("posresspec")){
    usage += "\t# position restraints specification\n";
    usage += "\t# @posresspec    filename\n\n";
  }
  
  if (knowns.count("refpos")){
    usage += "\t# position restraints reference positions\n";
    usage += "\t# @refpos    filename\n\n";
  }
  
  if (knowns.count("distrest")){
    usage += "\t# distance restraints specification\n";
    usage += "\t# @distrest  filename\n\n";
  }

  if (knowns.count("angrest")){
    usage += "\t# angle restraints specification\n";
    usage += "\t# @angrest  filename\n\n";
  }

  if (knowns.count("dihrest")){
    usage += "\t# dihedral restraints specification\n";
    usage += "\t# @dihrest  filename\n\n";
  }
  
  if (knowns.count("jval")){
    usage += "\t# J-value restraints specification\n";
    usage += "\t# @jval      filename\n\n";
  }

  if (knowns.count("xray")){
    usage += "\t# X-ray restraints specification\n";
    usage += "\t# @xray      filename\n\n";
  }
  
  if (knowns.count("sym")){
    usage += "\t# symmetry restraints specification\n";
    usage += "\t# @sym       filename\n\n";
  }

  if (knowns.count("lud")){
    usage += "\t# local elevation umbrella database file\n";
    usage += "\t# @lud       filename\n\n";
  }

  if (knowns.count("led")){
    usage += "\t# local elevation coordinate definition file\n";
    usage += "\t# @led       filename\n\n";
  }
  
  if (knowns.count("bsleus")) {
    usage += "\t# B&S-LEUS definition file\n";
    usage += "\t# @bsleus  filename\n\n";
  }

  if (knowns.count("order")){
    usage += "\t# order parameter restraints specification\n";
    usage += "\t# @order     filename\n\n";
  }
  
  if (knowns.count("rdc")){
    usage += "\t# RDC restraints specification\n";
    usage += "\t# @rdc       filename\n\n";
  }

  if (knowns.count("friction")){
    usage += "\t# atomic friction coefficients\n";
    usage += "\t# @friction   filename\n\n";
  }

  if (knowns.count("qmmm")){
    usage += "\t# QM/MM specification file\n";
    usage += "\t# @qmmm       filename\n\n";
  }
  
  if (knowns.count("repout")){
    usage += "\t# output file for replica exchange\n";
    usage += "\t# @repout    filename\n\n";
  }

  if (knowns.count("repdat")){
    usage += "\t# data file for replica exchange\n";
    usage += "\t# @repdat     filename\n\n";
  }
  
  if (knowns.count("print")){
    usage += "\t# print additional information\n";
    usage += "\t# @print     <pairlist/force>\n\n";
  }
  
  if (knowns.count("trp")){
    usage += "\t# write additional information to file\n";
    usage += "\t# @trp       filename\n\n";
  }
  
  if (knowns.count("anatrj")){
    usage += "\t# re-analyze trajectory\n";
    usage += "\t# @anatrj    filename\n\n";
  }
  
  if (knowns.count("verb")){
    usage += "\t# control verbosity (in debug builds)\n";
    usage += "\t# @verb      <[module:][submodule:]level>\n\n";
  }
  
  if (knowns.count("version")){
    usage += "\t# print version information\n";
    usage += "\t# @version\n\n";
  }
  
  if (knowns.count("develop")){
    usage += "\t# run untested code under development\n";
    usage += "\t# @develop\n\n";
  } 
  
  // usage += "#\n\n";

}

#ifdef OMP
#include <omp.h>
#endif
#ifdef XXMPI
#include <mpi.h>
#endif

void util::print_title(bool mpi, std::ostream & os, bool repex)
{
  os << GROMOSXX << "\n" 
     << "========\n";

  os << "version    :     " << MD_VERSION << "\n"
     << "build date :     " << MD_DATE << "\n\n";

#ifdef NDEBUG
  os << "Debugging is disabled.\n\n";
#else
  os << "This is a debug code. It will run much slower than the optimized code.\n"
          << "Use --disable-debug in order to get a faster code.\n\n";
#endif

  // some omp stuff
#ifdef OMP
  int nthreads = 0, tid = 0;
#pragma omp parallel private(nthreads, tid)
  {
    tid = omp_get_thread_num();
    if (tid == 0){
      nthreads = omp_get_num_threads();
      os << "OpenMP code enabled\n"
		<< "\tusing "
		<< omp_get_num_threads() << " threads.\n"
		<< "\tThis can be adjusted by setting the\n "
		<< "\tOMP_NUM_THREADS environment variable.\n"
		<< std::endl;
    }
  }
#endif
#ifdef XXMPI
  if (mpi) {
    os << "MPI code enabled\n";
    int size = MPI::COMM_WORLD.Get_size();
    if (size > 1)
      os << "\trunning on " << size << " nodes.\n";
    else
      os << "\trunning on " << size << " node.\n";
    int rank = MPI::COMM_WORLD.Get_rank();
    os << "\tthis is node " << rank << ".\n\n";
  }
#endif

  if (repex) {
    os << "using the replica exchange method.\n\t";
  }

  os << "\nGruppe fuer Informatikgestuetzte Chemie\n"
	    << "Professor W. F. van Gunsteren\n"
	    << "ETH Swiss Federal Institute of Technology\n"
	    << "Zurich, Switzerland\n\n"
	    << "Bugreports to http://www.gromos.net\n\n";

  struct utsname sysinf;
  if (uname(&sysinf) != -1){
    os << "Running on " << sysinf.nodename << " (" << sysinf.sysname
	      << " " << sysinf.release
	      << " " << sysinf.version
	      << " " << sysinf.machine
	      << ").\n\n";
  }
}

