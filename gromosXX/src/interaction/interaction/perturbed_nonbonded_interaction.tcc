/**
 * @file perturbed_nonbonded_interaction.tcc
 * template methods of Nonbonded_Interaction.
 */

#undef MODULE
#undef SUBMODULE

#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

/**
 * Constructor.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
						    t_pairlist, t_innerloop>
::Perturbed_Nonbonded_Interaction(t_simulation &sim, interaction::Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop> &nonbonded_interaction)
  : Interaction<t_simulation>("Perturbed NonBonded"),
    m_nonbonded_interaction(nonbonded_interaction)
{
  // this should maybe be done somewhere else, but it seems to work
  m_nonbonded_interaction.alpha_lj(sim.topology().alpha_lj());
  m_nonbonded_interaction.alpha_crf(sim.topology().alpha_crf());
}

/**
 * Destructor.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_pairlist, 
						    t_innerloop>
::~Perturbed_Nonbonded_Interaction()
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::destructor");
}

/**
 * calculate nonbonded forces and energies.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
						   t_pairlist, t_innerloop>
::calculate_interactions(t_simulation &sim)
{
  DEBUG(4, "Perturbed_Nonbonded_Interaction::calculate_interactions");

  // calculate forces / energies
  DEBUG(7, "\tshort range");

  do_perturbed_interactions
    (sim, m_nonbonded_interaction.pairlist().perturbed_begin(),
     m_nonbonded_interaction.pairlist().perturbed_end());
  
  // add 1,4 - interactions
  do_perturbed_14_interactions(sim);
  
  // possibly do the RF contributions due to excluded atoms
  if(sim.nonbonded().RF_exclusion()){
    do_perturbed_RF_excluded_interactions(sim);
  }

  // do the perturbed pairs
  do_perturbed_pair_interactions(sim);

}

/**
 * helper function to calculate perturbed forces and energies, 
 * stores them in the arrays pointed to by parameters
 * to make it usable for longrange calculations.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
						 t_pairlist, t_innerloop>
::do_perturbed_interactions(t_simulation &sim,
			    typename t_pairlist::iterator it,
			    typename t_pairlist::iterator to)
{  
  DEBUG(7, "\tcalculate perturbed interactions");  

  DEBUG(8, "Alpha LJ: " << m_nonbonded_interaction.alpha_lj() 
	<< " Alpha CRF: " << m_nonbonded_interaction.alpha_crf());

  for( ; it != to; ++it){    
    DEBUG(8, "perturbed pair: " << it.i() << " - " << *it);
    
    m_nonbonded_interaction.perturbed_interaction_inner_loop(sim, it.i(), *it);

  }

  // and long-range energy lambda-derivatives
  DEBUG(7, "add long-range lambda-derivatives");

  for(size_t i = 0; 
      i < m_nonbonded_interaction.pairlist().filter()
	.lambda_energies().lj_energy.size(); ++i){
    for(size_t j = 0; j < m_nonbonded_interaction.pairlist()
	  .filter().lambda_energies().lj_energy.size(); ++j){

      assert(sim.system().lambda_energies().lj_energy.size() > i);
      assert(sim.system().lambda_energies().lj_energy[i].size() > j);
      assert(sim.system().lambda_energies().lj_energy.size() > j);
      assert(sim.system().lambda_energies().lj_energy[j].size() > i);
      
      sim.system().lambda_energies().lj_energy[i][j] += 
	m_nonbonded_interaction.pairlist().filter().lambda_energies()
	.lj_energy[i][j];
      sim.system().lambda_energies().crf_energy[i][j] += 
	m_nonbonded_interaction.pairlist().filter().lambda_energies()
	.crf_energy[i][j];
    }
  }

  DEBUG(7, "end of function perturbed nonbonded interaction");
  
}

/**
 * helper function to calculate the forces and energies from the
 * 1,4 interactions.
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, 
					   t_pairlist, t_innerloop>
::do_perturbed_14_interactions(t_simulation &sim)
{
  DEBUG(7, "\tcalculate perturbed 1,4-interactions");

  std::set<int>::const_iterator it, to;
  std::map<size_t, simulation::Perturbed_Atom>::const_iterator 
    mit=sim.topology().perturbed_solute().atoms().begin(), 
    mto=sim.topology().perturbed_solute().atoms().end();
  
  for(; mit!=mto; ++mit){
    it = mit->second.one_four_pair().begin();
    to = mit->second.one_four_pair().end();
    
    for( ; it != to; ++it){

      m_nonbonded_interaction.perturbed_one_four_interaction_inner_loop
	(sim, mit->second.sequence_number(), *it);

    } // loop over 1,4 pairs
  } // loop over solute atoms
}  

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::Perturbed_Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::do_perturbed_RF_excluded_interactions(t_simulation &sim)
{
  DEBUG(7, "\tcalculate perturbed excluded RF interactions");

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();

  math::Vec r, f_rf, A_f_rf, B_f_rf;
  double e_rf, A_e_rf, B_e_rf, de_rf, A_de_rf, B_de_rf;
  const double l=sim.topology().lambda();
  
  std::set<int>::const_iterator it, to;
  std::map<size_t, simulation::Perturbed_Atom>::const_iterator
    mit=sim.topology().perturbed_solute().atoms().begin(),
    mto=sim.topology().perturbed_solute().atoms().end();
  DEBUG(7, "\tSize of perturbed atoms " << sim.topology().perturbed_solute().atoms().size());
  
  for(; mit!=mto; ++mit){
    // self term has already been calculated for state A, correct for that and 
    // calculate it for this lambda
    // only a distance independent part
    r=0.0;
    const int i=mit->second.sequence_number();
    const double q_i_a = mit->second.A_charge();
    const double q_i_b = mit->second.B_charge();

    // now calculate everything
    const double B_l = sim.topology().lambda();
    const double B_ln = pow(B_l, sim.topology().nlam());
    const double B_lnm = pow(B_l, sim.topology().nlam()-1);
    
    const double A_l = 1.0 - sim.topology().lambda();
    const double A_ln = pow(A_l, sim.topology().nlam());
    const double A_lnm = pow(A_l, sim.topology().nlam()-1);

    m_nonbonded_interaction.rf_soft_interaction(r, q_i_a*q_i_a, 
						B_l, A_f_rf, A_e_rf, A_de_rf);
    m_nonbonded_interaction.rf_soft_interaction(r, q_i_b*q_i_b, 
						A_l, B_f_rf, B_e_rf, B_de_rf);

    DEBUG(7, "Self term for atom " << i << " A: " << A_e_rf << " B: " << B_e_rf);
    
    // (1-l)^n * A + l^n *B - A 
    e_rf  = B_ln * B_e_rf  + A_ln * A_e_rf - A_e_rf;
    de_rf = B_ln * B_de_rf + A_ln * A_de_rf 
      + sim.topology().nlam() * B_lnm * B_e_rf 
      - sim.topology().nlam() * A_lnm * A_e_rf;
  
    sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(i)] += 0.5 * e_rf;
    sim.system().lambda_energies().crf_energy
      [sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(i)] += 0.5 * de_rf;
    
    // now loop over the exclusions
    it = mit->second.exclusion().begin();
    to = mit->second.exclusion().end();
    
    for( ;it != to; ++it){
      sim.system().periodicity().nearest_image(pos(i), 
					       pos(*it), r);
      DEBUG(7, "r2 i(" << i << "-" << *it << ") " << dot(r,r));
      
      double q_ij_a;
      double q_ij_b;
      if(sim.topology().perturbed_atom()[*it]){
	q_ij_a = q_i_a * 
	  sim.topology().perturbed_solute().atoms()[*it].A_charge();
	q_ij_b = q_i_b *
	  sim.topology().perturbed_solute().atoms()[*it].B_charge();
      }
      else{
	q_ij_a = q_i_a * sim.topology().charge()(*it);
	q_ij_b = q_i_b * sim.topology().charge()(*it);
      }
      m_nonbonded_interaction.rf_soft_interaction(r, q_ij_a, A_l, 
						  A_f_rf, A_e_rf, A_de_rf);
      m_nonbonded_interaction.rf_soft_interaction(r, q_ij_b, B_l, 
						  B_f_rf, B_e_rf, B_de_rf);
      DEBUG(7, "excluded atoms " << i << " & " << *it << " A: " << A_e_rf << " B: " << B_e_rf);
      e_rf  = B_ln * B_e_rf  + A_ln * A_e_rf;
      de_rf = B_ln * B_de_rf + A_ln * A_de_rf 
	+ sim.topology().nlam() * B_lnm * B_e_rf 
	- sim.topology().nlam() * A_lnm * A_e_rf;
      f_rf = B_ln * B_f_rf + A_ln * A_f_rf;
     
      // and add everything to the correct arrays 
      sim.system().energies().crf_energy 
	[sim.topology().atom_energy_group(i)]
	[sim.topology().atom_energy_group(*it)] += e_rf;
      sim.system().lambda_energies().crf_energy 
	[sim.topology().atom_energy_group(i)]
	[sim.topology().atom_energy_group(*it)] += de_rf;
      force(i) += f_rf;
      force(*it) -=f_rf;
    }
  }
}

/**
 * helper function to calculate the forces and energies from the
 * RF contribution of excluded atoms and self term
 */
template<typename t_simulation, typename t_pairlist, typename t_innerloop>
inline void interaction::
Perturbed_Nonbonded_Interaction<t_simulation, t_pairlist, t_innerloop>
::do_perturbed_pair_interactions(t_simulation &sim)
{

  std::vector<simulation::Perturbed_Atompair>::const_iterator
    it = sim.topology().perturbed_solute().atompairs().begin(),
    to = sim.topology().perturbed_solute().atompairs().end();
  
  math::Vec r, f, A_f, B_f;
  double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
    B_e_lj, B_e_crf, B_de_lj, B_de_crf;
  double e_lj, e_crf, de_lj, de_crf;
  lj_parameter_struct const *A_lj;
  lj_parameter_struct const *B_lj;
  double A_q, B_q;
  const double l = sim.topology().lambda();
  DEBUG(7, "lambda: " << l);
  bool is_perturbed;
  
  for(; it != to; ++it){

    DEBUG(7, "\tperturbed-pair\t" << it->i << "\t" << it->j);
    
    sim.system().periodicity().nearest_image(sim.system().pos()(it->i), 
					     sim.system().pos()(it->j), r);

    // is i perturbed?
    if (sim.topology().perturbed_atom()[it->i] == true){
      assert(sim.topology().perturbed_solute().atoms().count(it->i) == 1);

      is_perturbed = true;
      
      // is j perturbed as well?
      if(sim.topology().perturbed_atom()[it->j] == true){
	assert(sim.topology().perturbed_solute().atoms().count(it->j) == 1);
	
	A_lj =  &m_nonbonded_interaction.lj_parameter(
		sim.topology().perturbed_solute().atoms()[it->i].A_IAC(),
		sim.topology().perturbed_solute().atoms()[it->j].A_IAC());
	B_lj = &m_nonbonded_interaction.lj_parameter(
	        sim.topology().perturbed_solute().atoms()[it->i].B_IAC(),
		sim.topology().perturbed_solute().atoms()[it->j].B_IAC());

	A_q = sim.topology().perturbed_solute().atoms()[it->i].A_charge() * 
	  sim.topology().perturbed_solute().atoms()[it->j].A_charge();
	B_q = sim.topology().perturbed_solute().atoms()[it->i].B_charge() *
	  sim.topology().perturbed_solute().atoms()[it->j].B_charge();
      
      }
      else{ // only i perturbed
	A_lj = &m_nonbonded_interaction.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[it->i].A_IAC(),
	       sim.topology().iac(it->j));
	B_lj = &m_nonbonded_interaction.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[it->i].B_IAC(),
	       sim.topology().iac(it->j));
	A_q = sim.topology().perturbed_solute().atoms()[it->i].A_charge() * 
	  sim.topology().charge()(it->j);
	B_q = sim.topology().perturbed_solute().atoms()[it->i].B_charge() *
	  sim.topology().charge()(it->j);
      }
    }
    else{ // i unperturbed
      // is j perturbed
      if(sim.topology().perturbed_atom()[it->j] == true){
	assert(sim.topology().perturbed_solute().atoms().count(it->j) == 1);
	
	is_perturbed = true;

	A_lj =  &m_nonbonded_interaction.lj_parameter(
		sim.topology().iac(it->i),
		sim.topology().perturbed_solute().atoms()[it->j].A_IAC());
	B_lj = &m_nonbonded_interaction.lj_parameter(
	        sim.topology().iac(it->i),
		sim.topology().perturbed_solute().atoms()[it->j].B_IAC());

	A_q = sim.topology().charge()(it->i) *
	  sim.topology().perturbed_solute().atoms()[it->j].A_charge();
	B_q = sim.topology().charge()(it->j) *
	  sim.topology().perturbed_solute().atoms()[it->j].B_charge();
      
      }
      else{
	// both unperturbed

	is_perturbed = false;
	
	A_lj = &m_nonbonded_interaction.lj_parameter(
				     sim.topology().iac(it->i),
				     sim.topology().iac(it->j));
	B_lj = A_lj;
      
	A_q = sim.topology().charge()(it->i) *
	  sim.topology().charge()(it->j);
	B_q = A_q;
      }
      
      
    }

    // interaction in state A
    // ======================
    switch(it->A_interaction){
      // --------------------------
      case 0: // excluded
	A_e_lj = A_e_crf = A_de_lj = A_de_crf = 0.0;
	A_f = 0.0;
	DEBUG(7, "excluded in A");
	
	break;
	// --------------------------
      case 1: // normal interaction
	DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 
	      << " c12=" << A_lj->c12);
	DEBUG(7, "\tcharges state A i*j = " << A_q);
    
	if (is_perturbed){
	  DEBUG(7, "perturbed interaction");
	  m_nonbonded_interaction.lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
							  A_q, sim.topology().lambda(),
							  A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
	}
	else{
	  DEBUG(7, "non-perturbed interaction");
	  m_nonbonded_interaction.lj_crf_interaction(r, A_lj->c6, A_lj->c12,
						     A_q, A_f, A_e_lj, 
						     A_e_crf);
	  A_de_lj = A_de_crf = 0.0;
	}
 
	DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
	      << A_f << " e_lj: " 
	      << A_e_lj << " e_crf: " << A_e_crf 
	      << " de_lj: " << A_de_lj 
	      << " de_crf: " << A_de_crf);

	break;
      case 2: // 1,4 interaction
	DEBUG(7, "\tlj-parameter state A cs6=" << A_lj->cs6 
	      << " cs12=" << A_lj->cs12);
	DEBUG(7, "\tcharges state A i*j = " << A_q);
    
	if (is_perturbed){
	  DEBUG(7, "perturbed 1,4 interaction");
	  m_nonbonded_interaction.lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
							  A_q, sim.topology().lambda(),
							  A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
	}
	else{
	  DEBUG(7, "non-perturbed 1,4 interaction");
	  m_nonbonded_interaction.lj_crf_interaction(r, A_lj->cs6, A_lj->cs12,
						     A_q, A_f, A_e_lj, 
						     A_e_crf);
	  A_de_lj = A_de_crf = 0.0;
	}
 
	DEBUG(7, "\tcalculated interaction state A:\n\t\tf: "
	      << A_f << " e_lj: " 
	      << A_e_lj << " e_crf: " << A_e_crf 
	      << " de_lj: " << A_de_lj 
	      << " de_crf: " << A_de_crf);

	break;
	// --------------------------
    }

    // interaction in state B
    // ======================
    switch(it->B_interaction){
      // --------------------------
      case 0: // excluded
	DEBUG(7, "B excluded");
	B_e_lj = B_e_crf = B_de_lj = B_de_crf = 0.0;
	B_f = 0.0;
	
	break;
	// --------------------------
      case 1: // normal interaction
	DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 
	      << " c12=" << B_lj->c12);
	DEBUG(7, "\tcharges state B i*j = " << B_q);
    
	if (is_perturbed){
	  DEBUG(7, "perturbed interaction");
	  m_nonbonded_interaction.
	    lj_crf_soft_interaction(r, B_lj->c6, B_lj->c12,
				    B_q, 1.0 - sim.topology().lambda(),
				    B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
	}
	else{
	  DEBUG(7, "non-perturbed interaction");
	  m_nonbonded_interaction.lj_crf_interaction(r, B_lj->c6, B_lj->c12,
						     B_q, B_f, B_e_lj, 
						     B_e_crf);
	  B_de_lj = B_de_crf = 0.0;
	  
	}
 
	DEBUG(7, "\tcalculated interaction state B:\n\t\tf: "
	      << B_f << " e_lj: " 
	      << B_e_lj << " e_crf: " << B_e_crf 
	      << " de_lj: " << B_de_lj 
	      << " de_crf: " << B_de_crf);

	break;
      case 2: // 1,4 interaction
	DEBUG(7, "\tlj-parameter state B cs6=" << B_lj->cs6 
	      << " cs12=" << B_lj->cs12);
	DEBUG(7, "\tcharges state B i*j = " << B_q);
    
	if (is_perturbed){
	  DEBUG(7, "perturbed 1,4 interaction");
	  m_nonbonded_interaction.
	    lj_crf_soft_interaction(r, B_lj->cs6, B_lj->cs12,
				    B_q, 1.0 - sim.topology().lambda(),
				    B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
	}
	else{
	  DEBUG(7, "non-perturbed 1,4 interaction");
	  m_nonbonded_interaction.lj_crf_interaction(r, B_lj->cs6, B_lj->cs12,
						     B_q, B_f, B_e_lj, 
						     B_e_crf);
	  B_de_lj = B_de_crf = 0.0;
	}
	
	DEBUG(7, "\tcalculated interaction state B:\n\t\tf: "
	      << B_f << " e_lj: " 
	      << B_e_lj << " e_crf: " << B_e_crf 
	      << " de_lj: " << B_de_lj 
	      << " de_crf: " << B_de_crf);

	break;
	// --------------------------
    }
    
    // now combine everything
    const double B_l = sim.topology().lambda();
    const double B_ln = pow(B_l, sim.topology().nlam());
    const double B_lnm = pow(B_l, sim.topology().nlam()-1);
    
    const double A_l = 1.0 - sim.topology().lambda();
    const double A_ln = pow(A_l, sim.topology().nlam());
    const double A_lnm = pow(A_l, sim.topology().nlam()-1);
    
    DEBUG(10, "B_l: " << B_l <<
	  " B_ln: " << B_ln <<
	  " A_l: " << A_l <<
	  " A_ln: " << A_ln);

    f      = B_ln * B_f      + A_ln * A_f;
    e_lj   = B_ln * B_e_lj   + A_ln * A_e_lj;
    e_crf  = B_ln * B_e_crf  + A_ln * A_e_crf;
    de_lj  = B_ln * B_de_lj  + A_ln * A_de_lj  
      + sim.topology().nlam() * B_lnm * B_e_lj  
      - sim.topology().nlam() * A_lnm * A_e_lj;

    de_crf = B_ln * B_de_crf + A_ln * A_de_crf 
      + sim.topology().nlam() * B_lnm * B_e_crf 
      - sim.topology().nlam() * A_lnm * A_e_crf;
        
    sim.system().force()(it->i) += f;
    sim.system().force()(it->j) -= f;

    DEBUG(7, "A_lnm: " << A_lnm << " B_lnm: " << B_lnm);
    DEBUG(7, "\tcalculated interaction:\n\t\tf: " << f << " e_lj: " 
	  << e_lj << " e_crf: " << e_crf << " de_lj: " << de_lj 
	  << " de_crf: " << de_crf);
    
    // energy
    //assert(m_storage.energies().lj_energy.size() > 
    //   sim.topology().atom_energy_group(i));
    //assert(m_storage.energies().lj_energy.size() >
    //   sim.topology().atom_energy_group(j));
    //assert(m_storage.energies().crf_energy.size() > 
    //   sim.topology().atom_energy_group(i));
    //assert(m_storage.energies().crf_energy.size() >
    //   sim.topology().atom_energy_group(j));

    sim.system().energies().lj_energy[sim.topology().atom_energy_group(it->i)]
      [sim.topology().atom_energy_group(it->j)] += e_lj;
    
    sim.system().energies().crf_energy[sim.topology().atom_energy_group(it->i)]
      [sim.topology().atom_energy_group(it->j)] += e_crf;
    
    DEBUG(7, "\tenergy i and j " << sim.topology().atom_energy_group(it->i)
	  << " " << sim.topology().atom_energy_group(it->j));

    //assert(m_storage.lambda_energies().lj_energy.size() > 
    //   sim.topology().atom_energy_group(i));
    //assert(m_storage.lambda_energies().lj_energy.size() >
    //   sim.topology().atom_energy_group(j));
    //assert(m_storage.lambda_energies().crf_energy.size() > 
    //   sim.topology().atom_energy_group(i));
    //assert(m_storage.lambda_energies().crf_energy.size() >
    //   sim.topology().atom_energy_group(j));
        
    sim.system().lambda_energies().
      lj_energy[sim.topology().atom_energy_group(it->i)]
      [sim.topology().atom_energy_group(it->j)]
      += de_lj;

    sim.system().lambda_energies().
      crf_energy[sim.topology().atom_energy_group(it->i)]
      [sim.topology().atom_energy_group(it->j)]
      += de_crf;

  }
  
}

