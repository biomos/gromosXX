/**
 * @file harmonic_bond_interaction.tcc
 * template methods of harmonic_bond_interaction.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include <util/debug.h>

/**
 * calculate harmonic bond forces and energies.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_harmonic_bond_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::vector<interaction::bond_type_struct> const & param)
{
  // loop over the bonds
  std::vector<topology::two_body_term_struct>::const_iterator b_it =
    topo.solute().bonds().begin(),
    b_to = topo.solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double energy, diff;

  math::Periodicity<b> periodicity(conf.current().box);

  for( ; b_it != b_to; ++b_it){
    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist = sqrt(dot(v, v));
    
    assert(dist != 0.0);
    assert(unsigned(b_it->type) < param.size());
    
    DEBUG(7, "bond " << b_it->i << "-" << b_it->j << " type " << b_it->type);
    DEBUG(10, "K " << param[b_it->type].K << " r0 " 
	  << param[b_it->type].r0);

    DEBUG(10, "DF " << (-param[b_it->type].K * 
			(dist - param[b_it->type].r0) / dist) 
	  << "\n" << v);

    diff = dist - param[b_it->type].r0;

    f = v * (-param[b_it->type].K *
	     (diff) / dist);
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    v(a) * f(bb);

      DEBUG(7, "\tatomic virial done");
    }

    energy = 0.5 * param[b_it->type].K * diff * diff;
    conf.current().energies.bond_energy[topo.atom_energy_group()
					[b_it->i]] += energy;
    
  }

  return 0;
  
}


template<typename t_interaction_spec>
int interaction::Harmonic_Bond_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{

  const double start = util::now();

  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_harmonic_bond_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::triclinic :
      return _calculate_harmonic_bond_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::rectangular :
      return _calculate_harmonic_bond_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    default:
      throw std::string("Wrong boundary type");
  }

  m_timing += util::now() - start;
  
}
