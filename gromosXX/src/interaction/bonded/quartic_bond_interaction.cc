/**
 * @file quartic_bond_interaction.cc
 * template methods of Quartic_bond_interaction.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// interactions
#include "../../interaction/interaction_types.h"
#include "quartic_bond_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate quartic bond forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_quartic_bond_interactions(topology::Topology &topo,
						configuration::Configuration &conf,
						simulation::Simulation &sim)
{

  std::vector<interaction::bond_type_struct> const & bondtypes=topo.bond_types_quart();
  math::Periodicity<B> periodicity(conf.current().box);

  // loop over the bonds
  std::vector<topology::two_body_term_struct>::iterator b_it =
    topo.solute().bonds().begin(),
    b_to = topo.solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double e = 0.0;

  for( ; b_it != b_to; ++b_it){
    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist2 = abs2(v);
    
    assert(unsigned(b_it->type) < bondtypes.size());
    const double r02 = bondtypes[b_it->type].r0 *
      bondtypes[b_it->type].r0;

    DEBUG(7, "bond " << b_it->i << "-" << b_it->j
	  << " type " << b_it->type);
    DEBUG(10, "K " << bondtypes[b_it->type].K
	  << " r02 " << r02);
    DEBUG(10, "DF " << (-bondtypes[b_it->type].K *
			(dist2 - r02)) << "\n" << math::v2s(v));

    f = v * (-bondtypes[b_it->type].K *
	     (dist2 - r02));
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    // if (V == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int c=0; c<3; ++c)
	  conf.current().virial_tensor(a, c) += 
	    v(a) * f(c);
      
      DEBUG(7, "\tatomic virial done");
      // }

    e = 0.25 * bondtypes[b_it->type].K *
      (dist2 -r02) * (dist2 - r02);

    DEBUG(10, "energy: " << e);
    DEBUG(10, "bond energy size: " << unsigned(conf.current().energies.bond_energy.size()));
    DEBUG(10, "energy group size: " << unsigned(topo.atom_energy_group().size()));

    assert(conf.current().energies.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().energies.
      bond_energy[topo.atom_energy_group()[b_it->i]] += e;

    // ORIOL_GAMD
    if(sim.param().gamd.gamd){
      unsigned int gamd_group = topo.gamd_accel_group(b_it->i);
      std::vector<unsigned int> key = {gamd_group, gamd_group};
      unsigned int igroup = topo.gamd_interaction_group(key);
      DEBUG(10, "\tGAMD interaction group is " << igroup);
      conf.special().gamd.total_force[igroup](b_it->i) += f;
      conf.special().gamd.total_force[igroup](b_it->j) -= f;
      conf.current().energies.gamd_potential_total[igroup] += e;
      // virial
      for(int a=0; a<3; ++a){
        for(int bb=0; bb < 3; ++bb){
          conf.special().gamd.virial_tensor[igroup](a, bb) +=  v(a) * f(bb);
        }
      }

    } // end gamd
  }
 
  return 0;
}

int interaction::Quartic_Bond_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start(sim);

  SPLIT_VIRIAL_BOUNDARY(_calculate_quartic_bond_interactions, topo, conf, sim);
  
  m_timer.stop();
  
  return 0;
}

