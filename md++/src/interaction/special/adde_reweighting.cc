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
 * @file adde_reweighting.cc
 * template methods of Distance_Restraint_Interaction
 */

#include "../../stdheader.h"
#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"
#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"
#include "../../interaction/special/adde_reweighting.h"
#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate reweighting for adiabatic decoupling
 */

int interaction::Adde_Reweighting::calculate_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{ 
  unsigned int eg=sim.param().addecouple.adc_index()[0].eg;
  double vhl=0, evhl = 0.0;
  //int not_adde;
  conf.special().adde.vhh=conf.current().energies.lj_energy[eg][eg]+
          conf.current().energies.crf_energy[eg][eg];
  //for(unsigned int i=0; i<sim.param().multibath.multibath.size(); ++i){
  //  if(i!=sim.param().addecouple.adc_index()[0].tg)
  //    not_adde=i;
  //}
  for(unsigned int i=0; i<topo.energy_groups().size(); ++i){
    if(i==eg)
      continue;
    vhl+= conf.current().energies.lj_energy[i][eg]+
          conf.current().energies.crf_energy[i][eg]+
          conf.current().energies.lj_energy[eg][i]+
          conf.current().energies.crf_energy[eg][i];
  }
  if (sim.steps()==0)
    conf.special().adde.evhl=vhl;
  //double betal=1/
  //        (math::k_Boltzmann*sim.param().multibath.multibath[not_adde].temperature);
 
  evhl=(1-std::exp(-sim.time_step_size()/sim.param().addecouple.tmf))*
          vhl+
          std::exp(-sim.time_step_size()/sim.param().addecouple.tmf)*
          conf.special().adde.evhl;
  //evhl=(1-std::exp(-sim.time_step_size()/sim.param().addecouple.tmf))*
  //        std::exp(betal*(vhl-conf.special().adde.vhl0))+
  //        std::exp(-sim.time_step_size()/sim.param().addecouple.tmf)*
  //        conf.special().adde.evhl;
  
  conf.special().adde.evhl = evhl;  

  
  return 0;
}


