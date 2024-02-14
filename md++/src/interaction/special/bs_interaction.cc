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

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../util/bs_umbrella.h"
#include "bs_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

int interaction::BS_Interaction::init(topology::Topology& topo, 
                                      configuration::Configuration& conf, 
                                      simulation::Simulation& sim, 
                                      std::ostream& os, 
                                      bool quiet)
{
  DEBUG(4, "Initialze the B&S-LEUS interaction.");
  if(!quiet){
    os << "BSLEUS\n";
    os << conf.special().bs_umbrella.str();
    os << "END\n";
  }
  return 0;
}

int interaction::BS_Interaction::calculate_interactions(
                                        topology::Topology& topo, 
                                        configuration::Configuration& conf, 
                                        simulation::Simulation& sim)
{
  m_timer.start(sim);
  conf.special().bs_umbrella.apply(conf, sim);
  m_timer.stop();
  return 0;
}