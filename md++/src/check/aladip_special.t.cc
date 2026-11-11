/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file aladip.t.cc
 * tests using aladip
 */


#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../interaction/interaction.h"
#include "../interaction/forcefield/forcefield.h"

#include "../io/argument.h"
#include "../util/parse_verbosity.h"
#include "../util/usage.h"
#include "../util/error.h"

#include "../interaction/interaction_types.h"
#include "../io/instream.h"
#include "../util/parse_tcouple.h"
#include "../io/blockinput.h"
#include "../io/topology/in_topology.h"
#include "../io/message.h"

#include "../algorithm/integration/leap_frog.h"
#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/berendsen_thermostat.h"
#include "../algorithm/pressure/pressure_calculation.h"
#include "../algorithm/pressure/berendsen_barostat.h"

#include "../interaction/forcefield/create_forcefield.h"

#include "../util/create_simulation.h"
#include "../algorithm/create_md_sequence.h"

#include <time.h>

#include "check.h"
#include "check_forcefield.h"
#include "check_state.h"

#ifdef XXMPI
  #include <mpi.h>
#endif
#ifdef OMP
  #include <omp.h>
#endif

void hard_coded_values(std::map<std::string, double> & m){
  // Cartesian reference: d(5,9), angle(2,6,10), dihedral(2,4,6,8).
  // Harmonic radians; perturbed terms include 4*lambda*(1-lambda), lambda=0.125.
  m["ColvarRestraint"] = 6.17781013380949;
  m["PerturbedColvarRestraint"] = 6.63043981085383;
  m["XrayRestraint"] = 5.9411e+03;
  m["Local Elevation"] = 3.5284e+01;
  m["OrderParameterRestraint"] = 3.316416e-02;
}


int main(int argc, char* argv[]) {

  #ifdef XXMPI
  MPI_Init(&argc,&argv);
  #endif
  #ifdef OMP
    //omp_set_num_threads(1);
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      if (tid == 0){
        std::cout << "OpenMP code enabled; using " << omp_get_num_threads() << " threads." << std::endl;
      }
    }
  #endif

  int total = 0;

  util::Known knowns;
  knowns << "topo" << "pttopo" << "conf" << "input" << "colvarres" << "xray" << "led" << "lud" << "order" << "verb";
    
  std::string usage = argv[0];
  usage += "\n\t[@topo      <topology>]\n";
  usage += "\t[@pttopo   <perturbation topology>]\n";
  usage += "\t[@conf      <starting configuration>]\n";
  usage += "\t[@input     <input>]\n";
  usage += "\t[@colvarres <CV restraint specification>]\n";
  usage += "\t[@xray      <xray>]\n";
  usage += "\t[@led       <le definition file>]\n";
  usage += "\t[@lud       <le umbrella file>]\n";
  usage += "\t[@order     <order parameter specification file]\n";
  usage += "\t[@verb     <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);
      
  std::string stopo, spttopo, sconf, sinput, scolvarres, sxray, sled, slud, sorder;
  bool quiet = true;

  if (args.count("verb") != -1) quiet = false;
      
  if(args.count("topo") == 1)
    stopo = args["topo"];
  else
    GETFILEPATH(stopo, "aladip.topo", "src/check/data/");

  if(args.count("pttopo") == 1)
    spttopo = args["pttopo"];
  else
    GETFILEPATH(spttopo, "aladip.pttopo", "src/check/data/");
      
  if(args.count("conf") == 1)
    sconf = args["conf"];
  else
    GETFILEPATH(sconf, "aladip.conf", "src/check/data/");

  if(args.count("input") == 1)
    sinput = args["input"];
  else
    GETFILEPATH(sinput, "aladip_special.in", "src/check/data/");

  if(args.count("colvarres") == 1)
    scolvarres = args["colvarres"];
  else
    GETFILEPATH(scolvarres, "aladip.colvarres", "src/check/data/");

  if(args.count("led") ==1)
    sled = args["led"];
  else
    GETFILEPATH(sled, "aladip.led", "src/check/data/");
  //sled="";
  
  if(args.count("lud") ==1)
    slud = args["lud"];
  else
    GETFILEPATH(slud, "aladip.lud", "src/check/data/");
  //slud="";
  
  if(args.count("order") ==1)
    sorder = args["order"];
  else
    GETFILEPATH(sorder, "aladip.order", "src/check/data/");
  //sorder="";
    
  #ifdef HAVE_CLIPPER
    if(args.count("xray") ==1)
      sxray = args["xray"];
    else
      GETFILEPATH(sxray, "aladip.xrs", "src/check/data/");
  #else
    sxray = "";
  #endif
  
  if (!quiet)
    std::cout << "\n\n"
	      << "topology :      " << stopo << "\n"
	      << "perturbation :  " << spttopo << "\n"
	      << "input :         " << sinput << "\n"
              << "colvarres :      " << scolvarres << "\n"
        << "xray :          " << sxray << "\n"
	      << "led :           " << sled << "\n"
	      << "lud :           " << slud << "\n"
        << "order :         " << sorder << "\n"
	      << "configuration : " << sconf << "\n"
	      << std::endl;

  // set hard coded values to compare to
  std::map<std::string, double> ref_values;
  hard_coded_values(ref_values);

  util::simulation_struct aladip_sim;
  io::In_Topology in_topo;

  in_topo.quiet = quiet;

  if (util::create_simulation(stopo,
			      spttopo,
			      sconf,
			      sinput,
			      aladip_sim,
			      in_topo,
			      "",
			      "",
			      "",
            sxray,
            "", // qmmm
			      sled,
			      slud,
            sorder,  
			      quiet, scolvarres
			      ) != 0 ){
    std::cerr << "creating simulation failed!" << std::endl;
    return 1;
  }
  io::messages.display(std::cout);
  io::messages.clear();

  #ifndef HAVE_CLIPPER
    aladip_sim.sim.param().xrayrest.xrayrest = simulation::xrayrest_off;
  #endif
      
  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;

  if (interaction::create_g96_forcefield(*ff, 
					 aladip_sim.topo,
					 aladip_sim.sim,
					 in_topo,
					 std::cout,
					 quiet) != 0 ){
    std::cerr << "creating forcefield failed!" << std::endl;
    return 1;
  }
  io::messages.display(std::cout);
  io::messages.clear();

  bool ordinary_cv = false, perturbed_cv = false;
  for (interaction::Forcefield::iterator it = ff->begin(); it != ff->end(); ++it) {
    ordinary_cv |= ((*it)->name == "ColvarRestraint");
    perturbed_cv |= ((*it)->name == "PerturbedColvarRestraint");
  }
  if (!ordinary_cv || !perturbed_cv ||
      ff->init(aladip_sim.topo, aladip_sim.conf, aladip_sim.sim, std::cout, quiet) ||
      io::messages.contains(io::message::error) ||
      io::messages.contains(io::message::critical)) {
    io::messages.display(std::cerr);
    std::cerr << "CV forcefield registration/initialization failed!" << std::endl;
    delete ff;
    return 1;
  }

  // first check the forcefield
  total += check::check_forcefield(aladip_sim.topo, aladip_sim.conf, 
				   aladip_sim.sim, *ff, ref_values);
           
  io::messages.display(std::cout);
  io::messages.clear();

  delete ff;
  return total;
}
