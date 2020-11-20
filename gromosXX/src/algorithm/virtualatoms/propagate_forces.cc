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
