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
template<typename t_simulation>
inline interaction::harmonic_bond_interaction<t_simulation>
::harmonic_bond_interaction()
  : Interaction<t_simulation>("HarmonicBond")
{
}

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::harmonic_bond_interaction<t_simulation>
::~harmonic_bond_interaction()
{
}

/**
 * calculate harmonic bond forces and energies.
 */
template<typename t_simulation>
inline void interaction::harmonic_bond_interaction<t_simulation>
::calculate_interactions(t_simulation &sim)
{
  // loop over the bonds
  simulation::Bond::iterator b_it =
    sim.topology().solute().bonds().begin();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec v, f;

  for( ; !b_it.eol(); ++b_it){
    sim.system().periodicity()
      .nearest_image(pos(b_it.i()), pos(b_it.j()), v);

    double dist = sqrt(dot(v, v));
    
    assert(dist != 0.0);
    assert(unsigned(b_it.type()) < m_bond_parameter.size());
    
    DEBUG(7, "bond " << b_it.i() << "-" << b_it.j() << " type " << b_it.type());
    DEBUG(10, "K " << m_bond_parameter[b_it.type()].K << " r0 " << m_bond_parameter[b_it.type()].r0);

    DEBUG(10, "DF " << (-m_bond_parameter[b_it.type()].K * (dist - m_bond_parameter[b_it.type()].r0) / dist) << "\n" << v);

    f = v * (-m_bond_parameter[b_it.type()].K *
	     (dist - m_bond_parameter[b_it.type()].r0) / dist);
    
    force(b_it.i()) += f;
    force(b_it.j()) -= f;
  }
    
}

/**
 * add bond type.
 */
template<typename t_simulation>
inline void interaction::harmonic_bond_interaction<t_simulation>
::add(bond_type_struct s)
{
  m_bond_parameter.push_back(s);
}

/**
 * add bond type.
 */
template<typename t_simulation>
inline void interaction::harmonic_bond_interaction<t_simulation>
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
template<typename t_simulation>
inline std::vector<interaction::bond_type_struct> const &
interaction::harmonic_bond_interaction<t_simulation>
::parameter()const
{
  return m_bond_parameter;
}
