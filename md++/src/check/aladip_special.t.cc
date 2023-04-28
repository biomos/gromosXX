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

#ifdef OMP
  #include <omp.h>
#endif

void hard_coded_values(std::map<std::string, double> & m){
  m["DistanceRestraint"] = 257.189539;
  m["PerturbedDistanceRestraint"] = 195.899012;
  m["AngleRestraint"] = 54.5267;
  m["PerturbedAngleRestraint"] = 4.60917;
  m["DihedralRestraint"] = 2127.910749;
  m["PerturbedDihedralRestraint"] = 399.820792;
  m["XrayRestraint"] = 5.9411e+03;
  m["Local Elevation"] = 3.5284e+01;
  m["OrderParameterRestraint"] = 3.316416e-02;
  m["TFRDCRestraint"] = 0.276657772809125;
}


int main(int argc, char* argv[]) {

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
  knowns << "topo" << "pttopo" << "conf" << "input" << "distanceres" << "angrest"
         << "dihrest" << "xray" << "led" << "lud" << "order" << "tfrdcres" << "verb";
    
  std::string usage = argv[0];
  usage += "\n\t[@topo      <topology>]\n";
  usage += "\t[@pttopo   <perturbation topology>]\n";
  usage += "\t[@conf      <starting configuration>]\n";
  usage += "\t[@input     <input>]\n";
  usage += "\t[@distanceres  <distanceres>]\n";
  usage += "\t[@angrest   <angrest>]\n";
  usage += "\t[@dihrest   <dihrest>]\n";
  usage += "\t[@xray      <xray>]\n";
  usage += "\t[@led       <le definition file>]\n";
  usage += "\t[@lud       <le umbrella file>]\n";
  usage += "\t[@order     <order parameter specification file]\n";
  usage += "\t[@tfrdcres     <tensor-free RDC restraints specification file>]\n";
  usage += "\t[@verb     <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);
      
  std::string stopo, spttopo, sconf, sinput, sdistanceres, sangrest, sdihrest, 
    sxray, sled, slud, sorder, stfrdcres;
  bool quiet = false;

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

  if(args.count("distanceres") == 1)
    sdistanceres = args["distanceres"];
  else
    GETFILEPATH(sdistanceres, "aladip.distrest", "src/check/data/");

  if(args.count("angrest") ==1)
    sdihrest = args["angrest"];
  else
    GETFILEPATH(sangrest, "aladip.angrest", "src/check/data/");

  if(args.count("dihrest") ==1)
    sdihrest = args["dihrest"];
  else
    GETFILEPATH(sdihrest, "aladip.dihrest", "src/check/data/");

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

  if(args.count("tfrdcres") ==1)
    stfrdcres = args["tfrdcres"];
  else
    GETFILEPATH(stfrdcres, "aladip.tfr", "src/check/data/");
  //stfrdcres="";
    
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
	      << "distanceres :   " << sdistanceres << "\n"
	      << "angrest :       " << sangrest << "\n"
	      << "dihrest :       " << sdihrest << "\n"
	      << "xray :          " << sxray << "\n"
	      << "led :           " << sled << "\n"
	      << "lud :           " << slud << "\n"
	      << "order :         " << sorder << "\n"
	      << "tfrdcres :      " << stfrdcres << "\n"
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
			      sdistanceres,
			      sangrest,
			      sdihrest,
            sxray,
            "", // qmmm
			      sled,
			      slud,
            sorder,
			      stfrdcres,
			      quiet
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

  ff->init(aladip_sim.topo, aladip_sim.conf, aladip_sim.sim, std::cout,  quiet);

  // first check the forcefield
  total += check::check_forcefield(aladip_sim.topo, aladip_sim.conf, 
				   aladip_sim.sim, *ff, ref_values);
           
  io::messages.display(std::cout);
  io::messages.clear();

  return total;
}
