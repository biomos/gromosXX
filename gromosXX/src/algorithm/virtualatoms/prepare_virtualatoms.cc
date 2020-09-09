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
#include "../../interaction/nonbonded/interaction/nonbonded_interaction.h"

#include "prepare_virtualatoms.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE virtualatoms


int algorithm::Prepare_VirtualAtoms::apply(topology::Topology &topo, 
		  configuration::Configuration &conf,
		   simulation::Simulation &sim){
    // loop over virtual atoms and update positions in configuration
    std::map<unsigned int, util::Virtual_Atom>::iterator it;
    unsigned int atom_counter = 0;
    int num_atoms = topo.num_atoms();
    for ( it = topo.virtual_atoms_group().atoms().begin(); it != topo.virtual_atoms_group().atoms().end(); it++ )
    {
        // check if atom is in the  list of atoms with nonbonded interactions
        int atom_num = it->first;
        util::Virtual_Atom atom = it->second;
        if (atom.has_nonbonded()){
            conf.current().pos[num_atoms + atom_counter] = atom.pos(conf, topo); 
            atom_counter++;                       
        }

    }
    return 0;
}

int algorithm::Prepare_VirtualAtoms::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim, 
        std::ostream &os = std::cout, bool quiet = false)
    {

    std::map<unsigned int, util::Virtual_Atom>::iterator it;
    std::vector<util::Virtual_Atom> nonbonded_virtual;
    bool add_atom = false;
    int num_atoms = topo.num_atoms();
    // loop over all virtual atoms to check which ones have nonbonded interactions
    for ( it = topo.virtual_atoms_group().atoms().begin(); it != topo.virtual_atoms_group().atoms().end(); it++ )
    {
        int atom_num = it->first;
        util::Virtual_Atom atom = it->second;
        //check charge
        if (atom.charge() != 0){
            add_atom = true;
        }
        else{
            // check lj
            int iac = atom.iac();
            // "NonBonded"
            std::vector<std::vector<interaction::lj_parameter_struct> > lj_list;
            interaction::Nonbonded_Interaction * ni = dynamic_cast<interaction::Nonbonded_Interaction *> (m_ff.interaction("NonBonded"));
            if (sim.param().force.interaction_function == simulation::cgrain_func ||
             sim.param().force.interaction_function == simulation::cggromos_func){
                  lj_list = ni->parameter().cg_parameter();
             }
             else{
                 lj_list = ni->parameter().lj_parameter();
             }
            //loop over all pairs
            for (int i = 0; i < lj_list.size(); i++){
                for (int j = 0; j < lj_list.size(); j++){
                    if (i == iac || j == iac){
                        interaction::lj_parameter_struct lj_struct = lj_list[i][j];
                        if (lj_struct.c6 != 0.0 || lj_struct.c12 != 0.0 || lj_struct.cs6 != 0.0 || lj_struct.cs12 != 0.0){
                            add_atom = true;
                        }
                    } 
                }
            }
        } // end check lj
        if (add_atom){
            nonbonded_virtual.push_back(atom);
            atom.set_nonbonded();
            add_atom = false;

        }
    } // end loop over atoms
    // now resize the arrays
    topo.resize(num_atoms + nonbonded_virtual.size());
    conf.current().resize(num_atoms + nonbonded_virtual.size());
    conf.old().resize(num_atoms + nonbonded_virtual.size());
    // add the IAC and charge to the new arrays
    for ( unsigned int i = 0; i < nonbonded_virtual.size(); i++){
        topo.charge()[num_atoms + i] =  nonbonded_virtual[i].charge();
        topo.iac_array()[num_atoms + i] =  nonbonded_virtual[i].iac();
    }
    //os << "nonbonded virtuals " << nonbonded_virtual.size() << "\n";
    //os << "topo size " << num_atoms << "\n";
    //os << "topo size " << topo.charge().size() << "\n";
    return 0;

}