/**
 * @file quartic_bond_interaction.tcc
 * template methods of Quartic_bond_interaction.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Quartic_bond_interaction<t_simulation, t_interaction_spec>
::Quartic_bond_interaction()
  : Interaction<t_simulation, t_interaction_spec>("QuarticBond")
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Quartic_bond_interaction<t_simulation, t_interaction_spec>
::~Quartic_bond_interaction()
{
}

/**
 * calculate quartic bond forces and energies.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Quartic_bond_interaction<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  // loop over the bonds
  std::vector<simulation::Bond>::iterator b_it =
    sim.topology().solute().bonds().begin(),
    b_to = sim.topology().solute().bonds().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec v, f;

  double e;

  for( ; b_it != b_to; ++b_it){
    sim.system().periodicity()
      .nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist2 = dot(v, v);
    
    assert(unsigned(b_it->type) < m_bond_parameter.size());
    const double r02 = m_bond_parameter[b_it->type].r0 *
      m_bond_parameter[b_it->type].r0;

    DEBUG(7, "bond " << b_it->i << "-" << b_it->j
	  << " type " << b_it->type);
    DEBUG(10, "K " << m_bond_parameter[b_it->type].K
	  << " r02 " << r02);
    DEBUG(10, "DF " << (-m_bond_parameter[b_it->type].K *
			(dist2 - r02)) << "\n" << v);

    f = v * (-m_bond_parameter[b_it->type].K *
	     (dist2 - r02));
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    if (t_interaction_spec::do_virial == atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  sim.system().virial()(a, b) += 
	    v(a) * f(b);

      DEBUG(7, "\tatomic virial done");
    }

    e = 0.25 * m_bond_parameter[b_it->type].K *
      (dist2 -r02) * (dist2 - r02);

    assert(sim.system().energies().bond_energy.size() >
	   sim.topology().atom_energy_group()[b_it->i]);
    
    sim.system().energies().
      bond_energy[sim.topology().atom_energy_group()
		  [b_it->i]] += e;
  }
    
}

/**
 * add bond type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Quartic_bond_interaction<t_simulation, t_interaction_spec>
::add(bond_type_struct s)
{
  m_bond_parameter.push_back(s);
}

/**
 * add bond type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::Quartic_bond_interaction<t_simulation, t_interaction_spec>
::add(double K, double r0)
{
  bond_type_struct s;
  s.K = K;
  s.r0 = r0;
  add(s);
}

/**
 * access bond parameter.
 */
template<typename t_simulation, typename t_interaction_spec>
inline std::vector<interaction::bond_type_struct> const &
interaction::Quartic_bond_interaction<t_simulation, t_interaction_spec>
::parameter()const
{
  return m_bond_parameter;
}
