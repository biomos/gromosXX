/**
 * @file harmonic_bond_interaction.tcc
 * template methods of harmonic_bond_interaction.
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
inline interaction::harmonic_bond_interaction<t_simulation, t_interaction_spec>
::harmonic_bond_interaction()
  : Interaction<t_simulation, t_interaction_spec>("HarmonicBond")
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::harmonic_bond_interaction<t_simulation, t_interaction_spec>
::~harmonic_bond_interaction()
{
}

/**
 * calculate harmonic bond forces and energies.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::harmonic_bond_interaction<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  // loop over the bonds
  std::vector<simulation::Bond>::iterator b_it =
    sim.topology().solute().bonds().begin(),
    b_to = sim.topology().solute().bonds().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec v, f;

  double energy, diff;

  for( ; b_it != b_to; ++b_it){
    sim.system().periodicity()
      .nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist = sqrt(dot(v, v));
    
    assert(dist != 0.0);
    assert(unsigned(b_it->type) < m_bond_parameter.size());
    
    DEBUG(7, "bond " << b_it->i << "-" << b_it->j << " type " << b_it->type);
    DEBUG(10, "K " << m_bond_parameter[b_it->type].K << " r0 " 
	  << m_bond_parameter[b_it->type].r0);

    DEBUG(10, "DF " << (-m_bond_parameter[b_it->type].K * 
			(dist - m_bond_parameter[b_it->type].r0) / dist) 
	  << "\n" << v);

    diff = dist - m_bond_parameter[b_it->type].r0;

    f = v * (-m_bond_parameter[b_it->type].K *
	     (diff) / dist);
    
    force(b_it->i) += f;
    force(b_it->j) -= f;

    energy = 0.5 * m_bond_parameter[b_it->type].K * diff * diff;
    sim.system().energies().bond_energy[sim.topology().atom_energy_group()
					[b_it->i]] += energy;
    
  }
    
}

/**
 * add bond type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::harmonic_bond_interaction<t_simulation, t_interaction_spec>
::add(bond_type_struct s)
{
  m_bond_parameter.push_back(s);
}

/**
 * add bond type.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::harmonic_bond_interaction<t_simulation, t_interaction_spec>
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
interaction::harmonic_bond_interaction<t_simulation, t_interaction_spec>
::parameter()const
{
  return m_bond_parameter;
}
