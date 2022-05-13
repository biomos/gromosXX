/**
 * @file position_restraint_interaction.cc
 * template methods of Position_Restraint_Interaction
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

#include "../../interaction/special/position_restraint_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate position restraint interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_position_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim)
{
  // loop over the position restraints
  std::vector<topology::position_restraint_struct>::const_iterator 
    it = topo.position_restraints().begin(),
    to = topo.position_restraints().end();

  const math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  const math::VArray &ref = conf.special().reference_positions;
  const math::SArray &bfactor = conf.special().bfactors;
  math::Vec v, f;

  double energy = 0.0;

  math::Periodicity<B> periodicity(conf.current().box);

  for( ; it != to; ++it){

    periodicity.nearest_image(pos(it->seq), ref(it->seq), v);

    // double dist = sqrt(abs2(v));
    
    double bf = 1.0;
    if (sim.param().posrest.posrest == simulation::posrest_bfactor)
      bf = bfactor(it->seq);

    f = (- sim.param().posrest.force_constant / bf) * v;

    force(it->seq) += f;
  
    // should there be a contribution of this special ia to the virial?
    // no there should be NO contribution
    /*
    if (V == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    v(a) * f(bb);

      DEBUG(7, "\tatomic virial done");
    }
     */

    energy = 0.5 * sim.param().posrest.force_constant / bf * abs2(v);

    conf.current().energies.posrest_energy[topo.atom_energy_group()
					  [it->seq]] += energy;
    
  }

  return 0;
}

int interaction::Position_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  SPLIT_VIRIAL_BOUNDARY(_calculate_position_restraint_interactions,
			topo, conf, sim);
  
  return 0;
}
