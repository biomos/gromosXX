/**
 * @file position_restraint_interaction.tcc
 * template methods of Position_Restraint_Interaction
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

#include <util/debug.h>

/**
 * calculate position restraint interactions
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_position_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  // loop over the bonds
  std::vector<topology::position_restraint_struct>::const_iterator 
    it = topo.position_restraints().begin(),
    to = topo.position_restraints().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy;

  math::Periodicity<b> periodicity(conf.current().box);

  for( ; it != to; ++it){

    periodicity.nearest_image(pos(it->seq), it->pos, v);

    double dist = sqrt(dot(v, v));
    
    assert(dist != 0.0);
    
    f = (- sim.param().posrest.force_constant / it->bfactor) * v;

    force(it->seq) += f;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    v(a) * f(bb);

      DEBUG(7, "\tatomic virial done");
    }

    energy = 0.5 * sim.param().posrest.force_constant / it->bfactor * dist;

    conf.current().energies.posrest_energy[topo.atom_energy_group()
					   [it->seq]] += energy;
    
  }

  return 0;
  
}


template<typename t_interaction_spec>
int interaction::Position_Restraint_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_position_restraint_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim);
      break;
    case math::triclinic :
      return _calculate_position_restraint_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim);
      break;
    case math::rectangular :
      return _calculate_position_restraint_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
}
