/**
 * @file quartic_bond_interaction.tcc
 * template methods of Quartic_bond_interaction.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

/**
 * calculate quartic bond forces and energies.
 */
template<math::boundary_enum b, typename t_interaction_spec>
static int _calculate_quartic_bond_interactions(topology::Topology &topo,
			    configuration::Configuration &conf,
			    simulation::Simulation &sim,
			    std::vector<interaction::bond_type_struct> const & param)
{
  
  math::Periodicity<b> periodicity(conf.current().box);

  // loop over the bonds
  std::vector<topology::two_body_term_struct>::iterator b_it =
    topo.solute().bonds().begin(),
    b_to = topo.solute().bonds().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec v, f;

  double e;

  for( ; b_it != b_to; ++b_it){
    periodicity.nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist2 = dot(v, v);
    
    assert(unsigned(b_it->type) < param.size());
    const double r02 = param[b_it->type].r0 *
      param[b_it->type].r0;

    DEBUG(7, "bond " << b_it->i << "-" << b_it->j
	  << " type " << b_it->type);
    DEBUG(10, "K " << param[b_it->type].K
	  << " r02 " << r02);
    DEBUG(10, "DF " << (-param[b_it->type].K *
			(dist2 - r02)) << "\n" << v);

    f = v * (-param[b_it->type].K *
	     (dist2 - r02));
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int c=0; c<3; ++c)
	  conf.current().virial_tensor(a, c) += 
	    v(a) * f(c);

      DEBUG(7, "\tatomic virial done");
    }

    e = 0.25 * param[b_it->type].K *
      (dist2 -r02) * (dist2 - r02);

    DEBUG(10, "energy: " << e);
    DEBUG(10, "bond energy size: " << conf.current().energies.bond_energy.size());
    DEBUG(10, "energy group size: " << topo.atom_energy_group().size());

    assert(conf.current().energies.bond_energy.size() >
	   topo.atom_energy_group()[b_it->i]);
    
    conf.current().energies.
      bond_energy[topo.atom_energy_group()
		  [b_it->i]] += e;
  }
 
  return 0;
}

template<typename t_interaction_spec>
int interaction::Quartic_Bond_Interaction<t_interaction_spec>
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  const double start = util::now();
  
  switch(conf.boundary_type){
    case math::vacuum :
      return _calculate_quartic_bond_interactions<math::vacuum, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::triclinic :
      return _calculate_quartic_bond_interactions<math::triclinic, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    case math::rectangular :
      return _calculate_quartic_bond_interactions<math::rectangular, t_interaction_spec>
	(topo, conf, sim, m_parameter);
      break;
    default:
      throw std::string("Wrong boundary type");
  }
  
  m_timing += util::now() - start;
  
}

