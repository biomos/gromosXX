/**
 * @file perturbed_quartic_bond_interaction.tcc
 * template methods of Perturbed_Quartic_Bond_Interaction
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
inline interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
::Perturbed_Quartic_Bond_Interaction(
    interaction::Quartic_bond_interaction<t_simulation>
    & bond_interaction)
  : Interaction<t_simulation>("Perturbed QuarticBond"),
    m_bond_interaction(bond_interaction)
{
}

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
::~Perturbed_Quartic_Bond_Interaction()
{
}

/**
 * calculate quartic bond forces and energies and lambda derivatives.
 */
template<typename t_simulation>
inline void interaction::Perturbed_Quartic_Bond_Interaction<t_simulation>
::calculate_interactions(t_simulation &sim)
{
  // this is repeated code from Quartic_Bond_Interaction !!!

  DEBUG(5, "perturbed quartic bond interaction");
  DEBUG(7, "using the bond interaction: " << m_bond_interaction.name);

  // loop over the bonds
  std::vector<simulation::Perturbed_Bond>::iterator b_it =
    sim.topology().perturbed_solute().bonds().begin(),
    b_to = sim.topology().perturbed_solute().bonds().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec v, f;

  double e, e_lambda;

  for( ; b_it != b_to; ++b_it){

    DEBUG(7, "bond " << b_it->i << "-" << b_it->j
	  << " A-type " << b_it->type
	  << " B-type " << b_it->B_type
	  << " lambda " << sim.topology().lambda());

    assert(pos.size() > b_it->i && pos.size() > b_it->j);
    sim.system().periodicity()
      .nearest_image(pos(b_it->i), pos(b_it->j), v);

    double dist2 = dot(v, v);

    DEBUG(7, "dist2: " << dist2);

    assert(unsigned(b_it->type) < m_bond_interaction.parameter().size());
    assert(unsigned(b_it->B_type) < m_bond_interaction.parameter().size());

    const double K = (1-sim.topology().lambda()) *
      m_bond_interaction.parameter()[b_it->type].K +
      sim.topology().lambda() *
      m_bond_interaction.parameter()[b_it->B_type].K;
    
    DEBUG(7, "K: " << K);

    const double r02 = (1-sim.topology().lambda()) *
      m_bond_interaction.parameter()[b_it->type].r0 *
      m_bond_interaction.parameter()[b_it->type].r0 +
      sim.topology().lambda() *
      m_bond_interaction.parameter()[b_it->B_type].r0 *
      m_bond_interaction.parameter()[b_it->B_type].r0;

    DEBUG(7, "r02: " << r02);
    
    DEBUG(7, "DF " << K * (dist2 - r02) << "\n" << v);
    
    f = v * K * (dist2 - r02);
    
    force(b_it->i) += f;
    force(b_it->j) -= f;
  
    e = 0.25 * K * (dist2 -r02) * (dist2 - r02);

    DEBUG(7, "energy: " << e);

    const double K_diff = m_bond_interaction.parameter()[b_it->B_type].K -
      m_bond_interaction.parameter()[b_it->type].K;
    DEBUG(7, "K_diff: " << K_diff);
    
    const double b_diff = m_bond_interaction.parameter()[b_it->B_type].r0 -
      m_bond_interaction.parameter()[b_it->type].r0;
    DEBUG(7, "b_diff: " << b_diff);
    
    const double b_mix = m_bond_interaction.parameter()[b_it->type].r0 +
      sim.topology().lambda() * b_diff;
    DEBUG(7, "b_mix: " << b_mix);
    
    e_lambda = 0.25 * ( -4 * (m_bond_interaction.parameter()[b_it->type].K +
			      sim.topology().lambda() * K_diff) *
			b_diff * b_mix *
			(dist2 - b_mix * b_mix) +
			K_diff * 
			(dist2 - b_mix * b_mix) * (dist2 - b_mix * b_mix));
    DEBUG(7, "e_lambda: " << e_lambda);
    
    assert(sim.system().energies().bond_energy.size() >
	   sim.topology().atom_energy_group()[b_it->i]);
    
    sim.system().energies().
      bond_energy[sim.topology().atom_energy_group()
		  [b_it->i]] += e;
    
    assert(sim.system().lambda_energies().bond_energy.size() >
	   sim.topology().atom_energy_group()[b_it->i]);
    
    sim.system().lambda_energies().
      bond_energy[sim.topology().atom_energy_group()
		  [b_it->i]] += e_lambda;

  }
    
}

