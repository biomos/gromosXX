/**
 * @file perturbed_harmonic_bond_interaction.cc
 * template methods of Perturbed_Harmonic_Bond_Interaction
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include <util/debug.h>

/**
 * calculate quartic bond forces and energies and lambda derivatives.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_perturbed_hbond_interactions
(  topology::Topology & topo,
   configuration::Configuration & conf,
   simulation::Simulation & sim,
   interaction::Harmonic_Bond_Interaction<t_interaction_spec> const & m_interaction)
{

  DEBUG(7, "perturbed harmonic bond interaction");
  DEBUG(9, "using the bond interaction: " << m_interaction.name);
  DEBUG(9, std::setprecision(5));
  
  // loop over the bonds
  std::vector<topology::perturbed_two_body_term_struct>::const_iterator b_it =
    topo.perturbed_solute().bonds().begin(),
    b_to = topo.perturbed_solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double e, diff, e_lambda;

  math::Periodicity<b> periodicity(conf.current().box);

  for( ; b_it != b_to; ++b_it){

    DEBUG(7, "bond " << b_it->i << "-" << b_it->j
	  << " A-type " << b_it->A_type
	  << " B-type " << b_it->B_type
	  << " lambda " << topo.lambda());

    assert(pos.size() > int(b_it->i) && pos.size() > int(b_it->j));

    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist = sqrt(dot(v, v));

    DEBUG(7, "dist: " << dist);

    assert(dist != 0.0);
    assert(unsigned(b_it->A_type) < m_interaction.parameter().size());
    assert(unsigned(b_it->B_type) < m_interaction.parameter().size());

    const double K = (1-topo.lambda()) *
      m_interaction.parameter()[b_it->A_type].K +
      topo.lambda() *
      m_interaction.parameter()[b_it->B_type].K;
    
    DEBUG(7, "K: " << K);

    const double r0 = ((1-topo.lambda()) *
		       m_interaction.parameter()[b_it->A_type].r0 +
		       topo.lambda() *
		       m_interaction.parameter()[b_it->B_type].r0);
    diff = dist - r0;
    
    DEBUG(9, "r0: " << r0);
    
    DEBUG(9, "DF " << K * (diff) << "\n" << v);
    
    f = v * (-K) * diff/dist;

    DEBUG(9, "FORCE: " << f);
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int bb=0; bb<3; ++bb)
	  conf.current().virial_tensor(a, bb) += 
	    v(a) * f(bb);

      DEBUG(7, "\tatomic virial done");
    }
  
    e = 0.5 * K * diff * diff;

    DEBUG(9, "energy: " << e);

    const double K_diff = m_interaction.parameter()[b_it->B_type].K -
      m_interaction.parameter()[b_it->A_type].K;
    DEBUG(9, "K_diff: " << K_diff);
    
    const double b_diff = m_interaction.parameter()[b_it->B_type].r0 -
      m_interaction.parameter()[b_it->A_type].r0;
    DEBUG(9, "b_diff: " << b_diff);
    
    e_lambda = 0.5 * K_diff*diff*diff - K*diff*b_diff;

    DEBUG(9, "e_lambda: " << e_lambda);
    
    assert(conf.current().energies.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().energies.
      bond_energy[topo.atom_energy_group()
		  [b_it->i]] += e;
    
    assert(conf.current().perturbed_energy_derivatives.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().perturbed_energy_derivatives.
      bond_energy[topo.atom_energy_group()
		  [b_it->i]] += e_lambda;
    
  }
  
  return 0;
  
}


template<typename t_interaction_spec>
int interaction::Perturbed_Harmonic_Bond_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  const double start = util::now();
  
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_perturbed_hbond_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    case math::triclinic :
      return _calculate_perturbed_hbond_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    case math::rectangular :
      return _calculate_perturbed_hbond_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_interaction);
      break;
    default:
      throw std::string("Wrong boundary type");
  }

  m_timing += util::now() - start;
  
}
