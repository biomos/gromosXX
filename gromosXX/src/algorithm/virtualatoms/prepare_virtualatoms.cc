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
 * @file prepare_virtualatoms.cc
 * methods of the Prepare_VirtualAtoms
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../util/virtual_atom.h"
#include "../../interaction/interaction.h"
#include "../../interaction/forcefield/forcefield.h"
#include "../../interaction/forcefield/create_forcefield.h"
#include "../../interaction/nonbonded/interaction/nonbonded_interaction.h"
#include "../io/topology/in_topology.h"

#include "prepare_virtualatoms.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE virtualatoms


int algorithm::Prepare_VirtualAtoms::apply(topology::Topology &topo, 
		  configuration::Configuration &conf,
		   simulation::Simulation &sim){
    // loop over virtual atoms and update positions in configuration
    DEBUG(7, "Preparing Virtual Atoms");
    std::map<unsigned int, util::Virtual_Atom>::iterator it;
    for ( it = topo.virtual_atoms_group().atoms().begin(); it != topo.virtual_atoms_group().atoms().end(); it++ )
    {
        int atom_num = it->first;
        util::Virtual_Atom atom = it->second;
        conf.current().pos[atom_num] = atom.pos(conf, topo);
        conf.current().vel[atom_num] = 0.0;
        DEBUG(10 ,"Setting atom " << atom_num << " to position x: " << conf.current().pos[atom_num][0] << " y: " << conf.current().pos[atom_num][1]
                  << " z: " << conf.current().pos[atom_num][2]);                     
    }
    return 0;
}