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
 * @file aladip_unperturbed.t.cc
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
  m["QuarticBond"] = 19.356071;
  //m["PerturbedQuarticBond"] = 1.149568;
  m["Angle"] = 12.670413;
 // m["PerturbedAngle"] = 0.714818;
  m["ImproperDihedral"] = 1.485471;
 // m["PerturbedImproperDihedral"] = 2.642780;
  m["Dihedral"] = 6.292624;
//  m["PerturbedDihedral"] = 13.314602;
 // m["NonBonded_cg"] = -7.352312;
  m["NonBonded"] = -50.892127;
  m["NonBonded_newRF"] = -67.613934;
  m["NonBonded_atomic"] =  -50.791941;
 // m["DistanceRestraint"] = 257.189539;
//  m["PerturbedDistanceRestraint"] = 195.899012;
 // m["DihedralRestraint"] = 2127.910749;
//  m["PerturbedDihedralRestraint"] = 279.207857;
  //m["LeapFrogPositions"] = ;
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
  knowns << "topo" << "conf" << "input" << "verb";

  std::string usage = argv[0];
  usage += "\n\t[@topo    <topology>]\n";
  usage += "\t[@conf    <starting configuration>]\n";
  usage += "\t[@input   <input>]\n";
  usage += "\t[@verb   <[module:][submodule:]level>]\n";

  io::Argument args;
  if (args.parse(argc, argv, knowns, true)){
    std::cerr << usage << std::endl;
    return 1;
  }

  // parse the verbosity flag and set debug levels
  util::parse_verbosity(args);
      
  std::string stopo, spttopo, sconf, sinput;
  bool quiet = true;

  if (args.count("verb") != -1) quiet = false;
      
  if(args.count("topo") == 1)
    stopo = args["topo"];
  else
    GETFILEPATH(stopo, "aladip.topo", "src/check/data/");
      
  if(args.count("conf") == 1)
    sconf = args["conf"];
  else
    GETFILEPATH(sconf, "aladip.conf", "src/check/data/");

  if(args.count("input") == 1)
    sinput = args["input"];
  else
    GETFILEPATH(sinput, "aladip_unperturbed.in", "src/check/data/");

  if (!quiet)
    std::cout << "\n\n"
	      << "topology :      " << stopo << "\n"
	      << "input :         " << sinput << "\n"
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
			      "", "", "", "", "", "", "", "",
			      quiet
			      ) != 0 ){
    std::cerr << "creating simulation failed!" << std::endl;
    return 1;
  }

  io::messages.display(std::cout);
  io::messages.clear();
      
  if (true){
    // create a forcefield
    interaction::Forcefield *ff = new interaction::Forcefield;
	
    if (interaction::create_g96_forcefield(*ff, 
					   aladip_sim.topo,
					   aladip_sim.sim,
					   in_topo,
					   std::cout,
					   quiet)	!= 0 ){
      std::cerr << "creating forcefield failed!" << std::endl;
      return 1;
    }

    io::messages.display(std::cout);
    io::messages.clear();

    ff->init(aladip_sim.topo, aladip_sim.conf, aladip_sim.sim, std::cout,  quiet);

    // first check the forcefield
    total += check::check_forcefield(aladip_sim.topo, aladip_sim.conf, 
				     aladip_sim.sim, *ff, ref_values);

    // run it with atomic cutoff
    aladip_sim.sim.param().pairlist.atomic_cutoff = true;
    total += check::check_atomic_cutoff(aladip_sim.topo, aladip_sim.conf, 
				aladip_sim.sim, *ff, ref_values);
    aladip_sim.sim.param().pairlist.atomic_cutoff = false;

    // check virial, ...
    total += check::check_state(aladip_sim.topo, aladip_sim.conf, 
				aladip_sim.sim, *ff);
  }
  // vtune: run the thing...
  else{
	
    algorithm::create_md_sequence(aladip_sim.md,
				  aladip_sim.topo,
				  aladip_sim.sim,
				  in_topo);

    aladip_sim.md.init(aladip_sim.topo, 
          aladip_sim.conf,
          aladip_sim.sim);

    double end_time = aladip_sim.sim.param().step.t0 + 
        aladip_sim.sim.time_step_size() * 
        aladip_sim.sim.param().step.number_of_steps;
    
    int error = 0;

    while(aladip_sim.sim.time() < end_time){
      
      if ((error = aladip_sim.md.run(aladip_sim.topo, 
				     aladip_sim.conf,
				     aladip_sim.sim))){

	      std::cout << "\nError during MD run!\n" << std::endl;
	      // try to save the final structures...
	      break;
      }

      // update the energies
      aladip_sim.conf.old().energies.calculate_totals();

      // perturbed energy derivatives
      if (aladip_sim.sim.param().perturbation.perturbation){
	      aladip_sim.conf.old().perturbed_energy_derivatives.calculate_totals();
      }

      aladip_sim.conf.current().averages.apply(aladip_sim.topo,
	      aladip_sim.conf,
	      aladip_sim.sim);
	  
      aladip_sim.sim.time() +=  aladip_sim.sim.time_step_size();
      ++ aladip_sim.sim.steps();

    }

    aladip_sim.md.print_timing(std::cout);
  }
  
  io::messages.display(std::cout);
  io::messages.clear();

  return total;
}
