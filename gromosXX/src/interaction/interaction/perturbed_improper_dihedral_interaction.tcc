/**
 * @file perturbed_improper_dihedral_interaction.tcc
 * template methods of Perturbed_Improper_Dihedral_Interaction
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
inline interaction::Perturbed_Improper_Dihedral_Interaction<t_simulation>
::Perturbed_Improper_Dihedral_Interaction(
    interaction::Improper_dihedral_interaction<t_simulation>
    & improper_dihedral_interaction)
  : Interaction<t_simulation>("Perturbed Improper Dihedral"),
    m_improper_dihedral_interaction(improper_dihedral_interaction)
{
}

/**
 * Destructor.
 */
template<typename t_simulation>
inline interaction::Perturbed_Improper_Dihedral_Interaction<t_simulation>
::~Perturbed_Improper_Dihedral_Interaction()
{
}

/**
 * calculate angle forces and energies and lambda derivatives.
 */
template<typename t_simulation>
inline void interaction::Perturbed_Improper_Dihedral_Interaction<t_simulation>
::calculate_interactions(t_simulation &sim)
{
  // this is repeated code from Improper_Dihedral_Interaction !!!

  DEBUG(5, "perturbed improper dihedral interaction");
  DEBUG(7, "using the improper dihedral interaction: " 
	<< m_improper_dihedral_interaction.name);
  DEBUG(7, setprecision(5));
  
  // loop over the angles
  std::vector<simulation::Perturbed_Improper_Dihedral>::iterator i_it =
    sim.topology().perturbed_solute().improper_dihedrals().begin(),
    i_to = sim.topology().perturbed_solute().improper_dihedrals().end();

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();
  math::Vec rij, rkj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, q;
  double energy, e_lambda;

  for( ; i_it != i_to; ++i_it){

    DEBUG(7, "improper dihedral " << i_it->i << "-" << i_it->j << "-" 
	  << i_it->k << "-" << i_it->l
	  << " A-type " << i_it->type
	  << " B-type " << i_it->B_type
	  << " lambda " << sim.topology().lambda());

    assert(pos.size() > i_it->i && pos.size() > i_it->j && 
	   pos.size() > i_it->k && pos.size() > i_it->l);
    sim.system().periodicity().
      nearest_image(pos(i_it->k), pos(i_it->j), rkj);
    sim.system().periodicity().
      nearest_image(pos(i_it->i), pos(i_it->j), rij);
    sim.system().periodicity().
      nearest_image(pos(i_it->k), pos(i_it->l), rkl);

    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = dot(rkj, rkj);
    dmj2 = dot(rmj, rmj);
    dnk2 = dot(rnk, rnk);
    dkj  = sqrt(dkj2);
    dmj  = sqrt(dmj2);
    dnk  = sqrt(dnk2);
    
    DEBUG(15,"dkj="<<dkj<<" dmj="<<dmj<<" dnk="<<dnk);
   

    assert(dmj != 0.0);
    assert(dnk != 0.0);

    ip = dot(rmj, rnk);
   
    q  = acos(ip / (dmj*dnk));

    DEBUG(10, "zeta="<<q);
    
    ip = dot(rij, rnk);
    if(ip < 0) q *= -1.0;
    
    assert(unsigned(i_it->type) < m_improper_dihedral_interaction.parameter().size());

    double K    = (1-sim.topology().lambda()) *
      m_improper_dihedral_interaction.parameter()[i_it->type].K +
      sim.topology().lambda() *
      m_improper_dihedral_interaction.parameter()[i_it->B_type].K;
    double q0 =  (1-sim.topology().lambda()) *
      m_improper_dihedral_interaction.parameter()[i_it->type].q0 +
      sim.topology().lambda() *
      m_improper_dihedral_interaction.parameter()[i_it->B_type].q0;

    const double K_diff = 
      m_improper_dihedral_interaction.parameter()[i_it->B_type].K - 
      m_improper_dihedral_interaction.parameter()[i_it->type].K;
    const double q_diff =
      m_improper_dihedral_interaction.parameter()[i_it->B_type].q0- 
      m_improper_dihedral_interaction.parameter()[i_it->type].q0;
    
    DEBUG(10, "K=" << K << " q0=" << q0 );

    const double ki = -K * (q - q0) * dkj / dmj2;
    const double kl = K * (q - q0) * dkj / dnk2;
    const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
    const double kj2 = dot(rkl, rkj) / dkj2;
    
    fi = ki * rmj;
    fl = kl * rnk;
    fj = kj1 * fi - kj2 * fl;
    fk = -1.0*(fi + fj + fl);
    
    force(i_it->i) += fi;
    force(i_it->j) += fj;
    force(i_it->k) += fk;
    force(i_it->l) += fl;
    energy = 0.5 * K * (q-q0) * (q-q0);
    
    e_lambda = 0.5 * ( -2.0 * K * q_diff * (q-q0) +
		       K_diff * (q-q0) * (q-q0));
    
    assert(sim.system().energies().improper_energy.size() >
	   sim.topology().atom_energy_group()[i_it->i]);
    sim.system().energies().
      improper_energy[sim.topology().atom_energy_group()
		      [i_it->i]] += energy;

    assert(sim.system().lambda_energies().improper_energy.size() >
	   sim.topology().atom_energy_group()[i_it->i]);
    
    sim.system().lambda_energies().
      improper_energy[sim.topology().atom_energy_group()
		  [i_it->i]] += e_lambda;

   
  }
}

    


