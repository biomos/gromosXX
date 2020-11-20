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