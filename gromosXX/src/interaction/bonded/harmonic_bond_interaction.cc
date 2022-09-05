/**
 * @file harmonic_bond_interaction.cc
 * template methods of harmonic_bond_interaction.
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
#include "harmonic_bond_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate harmonic bond forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_harmonic_bond_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  std::vector<interaction::bond_type_struct> const & bondtypes=topo.bond_types_harm();
  // loop over the bonds
  std::vector<topology::two_body_term_struct>::const_iterator b_it =
    topo.solute().bonds().begin(),
    b_to = topo.solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy = 0.0, diff = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; b_it != b_to; ++b_it){
    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist = sqrt(abs2(v));
    
    assert(dist != 0.0);
    assert(unsigned(b_it->type) < bondtypes.size());
    
    DEBUG(7, "bond " << b_it->i << "-" << b_it->j << " type " << b_it->type);
    DEBUG(10, "K " << bondtypes[b_it->type].K << " r0 " 
	  << bondtypes[b_it->type].r0);
    DEBUG(10, "pos i " << math::v2s(pos(b_it->i)));
    DEBUG(10, "pos j " << math::v2s(pos(b_it->j)));
    DEBUG(10, "dist " << dist);
    DEBUG(10, "DF " << (-bondtypes[b_it->type].K * 
			(dist - bondtypes[b_it->type].r0) / dist) 
	  << "\nr(ij) " << math::v2s(v));

    diff = dist - bondtypes[b_it->type].r0;

    f = v * (-bondtypes[b_it->type].K *
	     (diff) / dist);
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    // if (V == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    v(a) * f(bb);

      DEBUG(7, "\tatomic virial done");
      // }

    energy = 0.5 * bondtypes[b_it->type].K * diff * diff;
    conf.current().energies.bond_energy[topo.atom_energy_group()
					[b_it->i]] += energy;

    DEBUG(9, "\tenergy = " << energy);
    
  }

  return 0;
  
}


int interaction::Harmonic_Bond_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  m_timer.start();

  SPLIT_VIRIAL_BOUNDARY(_calculate_harmonic_bond_interactions,
			topo, conf, sim);

  m_timer.stop();

  return 0;
  
}
