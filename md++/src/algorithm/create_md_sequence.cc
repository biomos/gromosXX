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

#include "../interaction/forcefield/forcefield.h"
#include "../interaction/forcefield/create_forcefield.h"

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

  // initialize the CudaManager
  // if user asks for CUDA acceleration, we call CUDA-capable variants wherever possible

  // center of mass motion printing / removal
  if (sim.param().centreofmass.skip_step ||
      sim.param().print.centreofmass){
      if (sim.param().analyze.analyze) {
          io::messages.add("COM removal is ignored with anatrj",
                "create_md_sequence", io::message::warning);
      } else {
        // if constexpr (false) {
        //   algorithm::Remove_COM_Motion<util::gpuBackend> * rcom_gpu =
        //     new algorithm::Remove_COM_Motion<util::gpuBackend>(os);
        //   md_seq.push_back(rcom_gpu);
        // } else {
        //   algorithm::Remove_COM_Motion<util::cpuBackend> * rcom =
        //     new algorithm::Remove_COM_Motion<util::cpuBackend>(os);
        //   md_seq.push_back(rcom);
        // }
        // We rather use a factory to create the proper variant
        // backend_factory<algorithm::Remove_COM_Motion>()
        // algorithm::make_algorithm<algorithm::Remove_COM_Motion>(os);
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
  if (sim.param().virtualatoms.virtualatoms)
    md_seq.push_back(new algorithm::Prepare_VirtualAtoms());

  // add the forcefield
  md_seq.push_back(ff);

  // propagate forces for virtual atoms
  if (sim.param().virtualatoms.virtualatoms)
      md_seq.push_back(new algorithm::Propagate_Forces());
  
  //add EDS
  if (sim.param().eds.eds) {
    md_seq.push_back(new algorithm::EDS());
  } 

  //ORIOL_GAMD add GAMD
  if (sim.param().gamd.gamd){
    md_seq.push_back(new algorithm::GAMD());
  } 
              
  // position constraints?
  if (sim.param().posrest.posrest == 3 && !sim.param().analyze.no_constraints) {
    algorithm::Position_Constraints * pc = new algorithm::Position_Constraints;
    md_seq.push_back(pc);
  }

  // monte-carlo steps?
  if (sim.param().montecarlo.mc){
    algorithm::Monte_Carlo * mc = new algorithm::Monte_Carlo(ff);
    md_seq.push_back(mc);
  }

  // energy minimisation or MD?
  algorithm::Stochastic_Dynamics_Vel1 * sd_vel = NULL;
  if (sim.param().minimise.ntem == 1){
    algorithm::Steepest_Descent * sd = new algorithm::Steepest_Descent;
    md_seq.push_back(sd);
  }
  else if (sim.param().minimise.ntem == 2 || sim.param().minimise.ntem == 3){
    algorithm::Conjugate_Gradient * cg = new algorithm::Conjugate_Gradient(md_seq);
    md_seq.push_back(cg);
  }
  else if (sim.param().analyze.analyze){
    algorithm::Analyze_Step * as = 
      new algorithm::Analyze_Step(sim.param().analyze.trajectory);
    md_seq.push_back(as);
  } else {    
    // SD ?
    if (sim.param().stochastic.sd){
      sd_vel = new algorithm::Stochastic_Dynamics_Vel1(sim.param());
      md_seq.push_back(sd_vel);
    }
    // MD ?
    else if(sim.param().integrate.method == simulation::integrate_leap_frog){
      if(sim.param().addecouple.adgr>0)
       md_seq.push_back(new algorithm::Scaled_Leap_Frog_Velocity);
      else
       md_seq.push_back(new algorithm::Leap_Frog_Velocity);
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

      if (sim.param().multibath.algorithm == 0){
	algorithm::Berendsen_Thermostat * tcoup =
	  new algorithm::Berendsen_Thermostat;
	md_seq.push_back(tcoup);
      }
      else if (sim.param().multibath.algorithm >= 1){
	algorithm::NoseHoover_Thermostat *tcoup =
	  new algorithm::NoseHoover_Thermostat;
	md_seq.push_back(tcoup);
      }
    }
    
    if (sim.param().stochastic.sd){
      //calculating the new position without the contribution form the random velocity
      md_seq.push_back(new algorithm::Stochastic_Dynamics_Pos1);
    }
    else if (sim.param().integrate.method == simulation::integrate_leap_frog){
       md_seq.push_back(new algorithm::Leap_Frog_Position);
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
      md_seq.push_back(new algorithm::Stochastic_Dynamics_Vel2);
    }

  // temperature calculation (always!)
  {
    algorithm::Temperature_Calculation * tcalc =
      new algorithm::Temperature_Calculation;

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
    algorithm::Pressure_Calculation * pcalc =
      new algorithm::Pressure_Calculation;
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
  if (sim.param().pcouple.scale != math::pcouple_off){
    algorithm::Berendsen_Barostat * pcoup =
      new algorithm::Berendsen_Barostat;
    md_seq.push_back(pcoup);
  }

  // slow growth
  if (sim.param().perturbation.perturbation){
    if (sim.param().perturbation.dlamt){
      algorithm::Slow_Growth *sg =
	new algorithm::Slow_Growth;
      md_seq.push_back(sg);
    }
  }
  
  // total energy calculation and energy average update
  {
    algorithm::Energy_Calculation * ec = 
      new algorithm::Energy_Calculation();
    md_seq.push_back(ec);
  }
  if (sim.param().stochastic.sd){
      // getting the new positions including the contribution from the random velocity
      md_seq.push_back(new algorithm::Stochastic_Dynamics_Pos2(sd_vel->rng(), 
              &sd_vel->random_vectors()));
      //do the constraints again without changing the velocities
      create_constraints(md_seq, topo, sim, it);
      
  }


  return 0;

}

