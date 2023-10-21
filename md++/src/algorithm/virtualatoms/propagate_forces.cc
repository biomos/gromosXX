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
 * @file propagate_forces.cc
 * methods of the Propagate_Forces class
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../util/virtual_atom.h"
#include "../../interaction/interaction.h"

#include "propagate_forces.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE virtualatoms

int algorithm::Propagate_Forces::apply(topology::Topology &topo, 
		  configuration::Configuration &conf,
		   simulation::Simulation &sim){
    DEBUG(7, "Propagate Forces of virtual atoms: apply");
    // loop over virtual atoms and propagate the assigned forces 
    std::map<unsigned int, util::Virtual_Atom>::iterator it;
    for ( it = topo.virtual_atoms_group().atoms().begin(); it != topo.virtual_atoms_group().atoms().end(); it++ )
    {
        int atom_num = it->first;
        util::Virtual_Atom atom = it->second;

        DEBUG(10, "Atom " << atom_num << " force: x = " <<  conf.current().force[atom_num][0] << " y= " <<
                   conf.current().force[atom_num][1] << " z= " << conf.current().force[atom_num][2]);

        // propagate the forces assigned to the virtual atoms and set them to 0
        atom.force(conf, topo, conf.current().force[atom_num]);
        conf.current().force[atom_num] = 0.0;
        conf.old().force[atom_num] = 0.0;            
    }
    return 0;
}
