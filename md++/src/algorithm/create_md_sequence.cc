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
 * @file create_md_sequence.cc
 */

#include "../stdheader.h"
#include <fstream>

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../interaction/interaction.h"
#include "../interaction/interaction_types.h"

#include "../io/argument.h"
#include "../io/blockinput.h"
#include "../io/instream.h"
#include "../io/topology/in_topology.h"

#include "../interaction/forcefield/forcefield.h"
#include "../interaction/forcefield/create_forcefield.h"

#include "../algorithm/algorithm/algorithm_sequence.h"
#include "../algorithm/constraints/position_constraints.h"
#include "../algorithm/integration/energy_calculation.h"
#include "../algorithm/integration/leap_frog.h"
#include "../algorithm/integration/scaled_leap_frog.h"
#include "../algorithm/integration/monte_carlo.h"
#include "../algorithm/integration/stochastic.h"
#include "../algorithm/integration/lattice_shift.h"
#include "../algorithm/integration/multigradient.h"
#include "../algorithm/integration/eds.h"
#include "../algorithm/integration/gamd.h"

#include "../io/blockinput.h"
#include "../io/instream.h"
#include "../io/configuration/inframe.h"
#include "../io/configuration/in_configuration.h"

#include "../algorithm/integration/analyze.h"
#include "../algorithm/temperature/temperature_calculation.h"
#include "../algorithm/temperature/nosehoover_thermostat.h"
#include "../algorithm/temperature/berendsen_thermostat.h"
#include "../algorithm/pressure/pressure_calculation.h"
#include "../algorithm/pressure/berendsen_barostat.h"
#include "../algorithm/virtualatoms/prepare_virtualatoms.h"
#include "../algorithm/virtualatoms/propagate_forces.h"

#include "../math/periodicity.h"

#include "gpu/cuda/manager/cuda_manager.h"

#include "../algorithm/constraints/create_constraints.h"
#include "../algorithm/constraints/remove_com_motion.h"


#include "../algorithm/integration/slow_growth.h"
#include "../algorithm/integration/steepest_descent.h"
#include "../algorithm/integration/conjugate_gradient.h"


#include "../io/print_block.h"

#include "../algorithm/create_md_sequence.h"



#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

int algorithm::create_md_sequence(algorithm::Algorithm_Sequence &md_seq,
				  topology::Topology &topo,
				  simulation::Simulation & sim,
				  io::In_Topology &it,
				  std::ostream & os,
				  bool quiet)
{


  // create a forcefield
  interaction::Forcefield *ff = new interaction::Forcefield;
  interaction::create_g96_forcefield(*ff, topo, sim, it, os, quiet);
  
  //==================================================
  // construct the md algorithm
  //==================================================

  if (sim.param().multigradient.multigradient) {
    md_seq.push_back(new algorithm::Multi_Gradient());
  }

  // center of mass motion printing / removal
  if (sim.param().centreofmass.skip_step ||
      sim.param().print.centreofmass) {
      if (sim.param().analyze.analyze) {
          io::messages.add("COM removal is ignored with anatrj",
                "create_md_sequence", io::message::warning);
      } else {
        algorithm::IAlgorithm * rcom = algorithm::make_algorithm<algorithm::Remove_COM_Motion>(sim, os);
        md_seq.push_back(rcom);
      }
  }
  
  // add the lattice shift tracking 
  if (sim.param().boundary.boundary != math::vacuum) {
    algorithm::IAlgorithm * lst = algorithm::make_algorithm<algorithm::Lattice_Shift_Tracker>(sim);
    md_seq.push_back(lst);
  }
  
  // prepare virtual atoms
  if (sim.param().virtualatoms.virtualatoms) {
    algorithm::IAlgorithm * pva = algorithm::make_algorithm<algorithm::Prepare_VirtualAtoms>(sim);
    md_seq.push_back(pva);
  }

  // add the forcefield
  md_seq.push_back(ff);

  // propagate forces for virtual atoms
  if (sim.param().virtualatoms.virtualatoms) {
      algorithm::IAlgorithm * pf = algorithm::make_algorithm<algorithm::Propagate_Forces>(sim);
      md_seq.push_back(pf);
  }
  
  //add EDS
  if (sim.param().eds.eds) {
    algorithm::IAlgorithm * eds = algorithm::make_algorithm<algorithm::EDS>(sim);
    md_seq.push_back(eds);
  } 

  //ORIOL_GAMD add GAMD
  if (sim.param().gamd.gamd){
    algorithm::IAlgorithm * gamd = algorithm::make_algorithm<algorithm::GAMD>(sim);
    md_seq.push_back(gamd);
  } 
              
  // position constraints?
  if (sim.param().posrest.posrest == 3 && !sim.param().analyze.no_constraints) {
    algorithm::IAlgorithm * pc = algorithm::make_algorithm<algorithm::Position_Constraints>(sim);
    md_seq.push_back(pc);
  }

  // monte-carlo steps?
  if (sim.param().montecarlo.mc){
    algorithm::IAlgorithm * mc = algorithm::make_algorithm<algorithm::Monte_Carlo>(sim, ff);
    md_seq.push_back(mc);
  }

  // energy minimisation or MD?
  algorithm::Stochastic_Dynamics_Vel1 * sd_vel = nullptr;
  if (sim.param().minimise.ntem == simulation::emin_steepest_descent) {
    algorithm::IAlgorithm * sd = algorithm::make_algorithm<algorithm::Steepest_Descent>(sim);
    md_seq.push_back(sd);
  }
  else if (sim.param().minimise.ntem == simulation::emin_conjugate_gradient_fr
            || sim.param().minimise.ntem == simulation::emin_conjugate_gradient_pr) {
    algorithm::IAlgorithm * cg = algorithm::make_algorithm<algorithm::Conjugate_Gradient>(sim, md_seq);
    md_seq.push_back(cg);
  }
  else if (sim.param().analyze.analyze) {
    algorithm::IAlgorithm * as = algorithm::make_algorithm<algorithm::Analyze_Step>(sim, sim.param().analyze.trajectory);
    md_seq.push_back(as);
  } else {    
    // SD ?
    if (sim.param().stochastic.sd) {
      sd_vel = algorithm::make_algorithm<algorithm::Stochastic_Dynamics_Vel1>(sim, sim.param());
      assert(sd_vel != nullptr && "make_algorithm did not return Stochastic_Dynamics_Vel1");
      md_seq.push_back(sd_vel);
    }
    // MD ?
    else if (sim.param().integrate.method == simulation::integrate_leap_frog) {
      algorithm::IAlgorithm * lfv = nullptr;
      if(sim.param().addecouple.adgr > 0)
          lfv = algorithm::make_algorithm<algorithm::Scaled_Leap_Frog_Velocity>(sim);
      else
          lfv = algorithm::make_algorithm<algorithm::Leap_Frog_Velocity>(sim);
      md_seq.push_back(lfv);
    }
    // ??
    else{
      if (!quiet)
        std::cout << "\tno integration (velocities) selected!\n";
    }

    // temperature scaling? -> has to be done before temperature calculation!!!
    if (sim.param().multibath.couple){
      
      if (sim.param().stochastic.sd){
	io::messages.add("temperature coupling only of atoms with gamma=0 not implemented,"
			 " couples ALL ATOMS!",
			 "create_md_sequence",
			 io::message::warning);
      }

      if (sim.param().multibath.algorithm == 0) {
        algorithm::IAlgorithm * tcoup = algorithm::make_algorithm<algorithm::Berendsen_Thermostat>(sim);
        md_seq.push_back(tcoup);
      }
      else if (sim.param().multibath.algorithm >= 1) {
        algorithm::IAlgorithm * tcoup = algorithm::make_algorithm<algorithm::NoseHoover_Thermostat>(sim);
        md_seq.push_back(tcoup);
      }
    }
    
    if (sim.param().stochastic.sd) {
      //calculating the new position without the contribution form the random velocity
      algorithm::IAlgorithm * sdp1 = algorithm::make_algorithm<algorithm::Stochastic_Dynamics_Pos1>(sim);
      md_seq.push_back(sdp1);
    }
    else if (sim.param().integrate.method == simulation::integrate_leap_frog) {
      algorithm::IAlgorithm * lfp = algorithm::make_algorithm<algorithm::Leap_Frog_Position>(sim);
       md_seq.push_back(lfp);
    }
    else{
      if (!quiet)
	std::cout << "\tno integration (position) selected!\n";
    }
  }
  
  // CONSTRAINTS
  create_constraints(md_seq, topo, sim, it);
  if (sim.param().stochastic.sd){
      //getting the Velocities from the restraint conformation
      algorithm::IAlgorithm * sdv2 = algorithm::make_algorithm<algorithm::Stochastic_Dynamics_Vel2>(sim);
      md_seq.push_back(sdv2);
    }

  // temperature calculation (always!)
  {
    algorithm::IAlgorithm * tcalc = algorithm::make_algorithm<algorithm::Temperature_Calculation>(sim);
    DEBUG(7, tcalc->name);
    md_seq.push_back(tcalc);

  }
  
  // pressure calculation?
  if (!quiet){
    io::print_PCOUPLE(os, sim.param().pcouple.calculate,
		      sim.param().pcouple.scale,
		      sim.param().pcouple.pres0,
		      sim.param().pcouple.compressibility,
		      sim.param().pcouple.tau,
		      sim.param().pcouple.virial,
                      sim.param().pcouple.x_semi,
                      sim.param().pcouple.y_semi,
                      sim.param().pcouple.z_semi);
  }
  
  if (sim.param().pcouple.calculate){
    algorithm::Pressure_Calculation * pcalc = algorithm::make_algorithm<algorithm::Pressure_Calculation>(sim);
    md_seq.push_back(pcalc);
    
    // coarse grain factor for pressure correction
    double tot_factor = 0.0;
    for (unsigned int i = 0; i < topo.num_atoms(); ++i)
      tot_factor += topo.cg_factor(i);
    topo.tot_cg_factor() = double(tot_factor) / double(topo.num_atoms());
    DEBUG(10, "atoms = " << topo.num_atoms());
    DEBUG(10, "coarse grain correction factor = " << topo.tot_cg_factor());
  }

  //  pressure scaling
  if (sim.param().pcouple.scale != math::pcouple_off) {
    algorithm::Berendsen_Barostat * pcoup = algorithm::make_algorithm<algorithm::Berendsen_Barostat>(sim);
    md_seq.push_back(pcoup);
  }

  // slow growth
  if (sim.param().perturbation.perturbation) {
    if (sim.param().perturbation.dlamt) {
      algorithm::Slow_Growth * sg = algorithm::make_algorithm<algorithm::Slow_Growth>(sim);
      md_seq.push_back(sg);
    }
  }
  
  // total energy calculation and energy average update
  {
    algorithm::Energy_Calculation * ec = algorithm::make_algorithm<algorithm::Energy_Calculation>(sim);
    md_seq.push_back(ec);
  }
  if (sim.param().stochastic.sd) {
      // getting the new positions including the contribution from the random velocity
      algorithm::IAlgorithm * sdp2 = algorithm::make_algorithm<
                                    algorithm::Stochastic_Dynamics_Pos2>(sim, sd_vel->rng(), 
                                                                  &sd_vel->random_vectors());
      md_seq.push_back(sdp2);
      //do the constraints again without changing the velocities
      create_constraints(md_seq, topo, sim, it);
  }


  return 0;

}

