/**
 * @file perturbed_angle_interaction.tcc
 * template methods of Perturbed_Angle_Interaction
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
inline interaction::Perturbed_Angle_Interaction<t_simulation, t_interaction_spec>
::Perturbed_Angle_Interaction(
    interaction::angle_interaction<t_simulation, t_interaction_spec>
    & angle_interaction)
  : Interaction<t_simulation, t_interaction_spec>("Perturbed Angle"),
    m_angle_interaction(angle_interaction)
{
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_interaction_spec>
inline interaction::Perturbed_Angle_Interaction<t_simulation, t_interaction_spec>
::~Perturbed_Angle_Interaction()
{
}

/**
 * calculate angle forces and energies and lambda derivatives.
 */
template<typename t_simulation, typename t_interaction_spec>
inline void interaction::
Perturbed_Angle_Interaction<t_simulation, t_interaction_spec>
::calculate_interactions(t_simulation &sim)
{
  // this is repeated code from Angle_Interaction !!!

  DEBUG(5, "perturbed angle interaction");
  DEBUG(7, "using the angle interaction: " << m_angle_interaction.name);
  DEBUG(7, setprecision(5));
  
  // loop over the angles
  std::vector<simulation::Perturbed_Angle>::iterator a_it =
    sim.topology().perturbed_solute().angles().begin(),
    a_to = sim.topology().perturbed_solute().angles().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec rij, rkj, fi, fj, fk;

  double energy, e_lambda;

  for( ; a_it != a_to; ++a_it){

    DEBUG(7, "angle " << a_it->i << "-" << a_it->j << "-" << a_it->k
	  << " A-type " << a_it->type
	  << " B-type " << a_it->B_type
	  << " lambda " << sim.topology().lambda());

    assert(pos.size() > a_it->i && pos.size() > a_it->j && 
	   pos.size() > a_it->k);
    sim.system().periodicity().nearest_image(pos(a_it->i), pos(a_it->j), rij);
    sim.system().periodicity().nearest_image(pos(a_it->k), pos(a_it->j), rkj);

    double dij = sqrt(dot(rij, rij));
    double dkj = sqrt(dot(rkj, rkj));
    
    assert(dij != 0.0);
    assert(dkj != 0.0);

    double ip = dot(rij, rkj);
    double cost = ip / (dij * dkj);
    assert(unsigned(a_it->type)   < m_angle_interaction.parameter().size());
    assert(unsigned(a_it->B_type) < m_angle_interaction.parameter().size());
    
    double K    = (1-sim.topology().lambda()) *
      m_angle_interaction.parameter()[a_it->type].K +
      sim.topology().lambda() *
      m_angle_interaction.parameter()[a_it->B_type].K;
    double cos0 =  (1-sim.topology().lambda()) *
      m_angle_interaction.parameter()[a_it->type].cos0 +
      sim.topology().lambda() *
      m_angle_interaction.parameter()[a_it->B_type].cos0;

    const double K_diff = m_angle_interaction.parameter()[a_it->B_type].K - 
      m_angle_interaction.parameter()[a_it->type].K;
    const double cos_diff=m_angle_interaction.parameter()[a_it->B_type].cos0- 
      m_angle_interaction.parameter()[a_it->type].cos0;
    
    DEBUG(10, "K=" << K << " cos0=" << cos0 << " dij=" << dij << " dkj=" << dkj)
;

    double ki = -K * (cost - cos0) / dij;
    double kk = -K * (cost - cos0) / dkj;
    
    DEBUG(10, "cost=" << cost << " ki=" << ki << " kk=" << kk);

    fi = ki*(rkj/dkj - rij/dij * cost);
    fk = kk*(rij/dij - rkj/dkj * cost);
    fj = -1.0 * fi - fk;
    
    force(a_it->i) += fi;
    force(a_it->j) += fj;
    force(a_it->k) += fk;

    energy = 0.5 * K * (cost - cos0) * (cost - cos0);

    e_lambda = 0.5 * ( -2.0 * K * cos_diff * (cost - cos0) +
		       K_diff * (cost - cos0) * (cost - cos0) );

    DEBUG(7, "energy: " << energy);

    DEBUG(7, "K_diff: " << K_diff);

    DEBUG(7, "cos_diff: " << cos_diff);
    
    DEBUG(7, "e_lambda: " << e_lambda);
    
    assert(sim.system().energies().angle_energy.size() >
	   sim.topology().atom_energy_group()[a_it->i]);
    
    sim.system().energies().
      angle_energy[sim.topology().atom_energy_group()
		  [a_it->i]] += energy;
    
    assert(sim.system().lambda_energies().angle_energy.size() >
	   sim.topology().atom_energy_group()[a_it->i]);
    
    sim.system().lambda_energies().
      angle_energy[sim.topology().atom_energy_group()
		  [a_it->i]] += e_lambda;

  }
    
}

