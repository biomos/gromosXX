/**
 * @File perturbed_nonbonded_inner_loop.tcc
 * template methods of Perturbed_Nonbonded_Inner_Loop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

template<typename t_simulation, typename t_storage, bool do_scaling>
interaction::Perturbed_Nonbonded_Inner_Loop<t_simulation, t_storage, do_scaling>
::Perturbed_Nonbonded_Inner_Loop(Nonbonded_Base &base, t_storage &storage)
  : Nonbonded_Inner_Loop<t_simulation, t_storage>(base, storage)
{
}

template<typename t_simulation, typename t_storage, bool do_scaling>
void interaction::Perturbed_Nonbonded_Inner_Loop<t_simulation, t_storage, 
						 do_scaling>
::perturbed_interaction_inner_loop(t_simulation &sim, 
				   size_t const i, size_t const j)
{
    DEBUG(7, "\treally! perturbed-pair\t" << i << "\t" << j);

    math::Vec r, f, A_f, B_f;
    double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
           B_e_lj, B_e_crf, B_de_lj, B_de_crf;
    double e_lj, e_crf, de_lj, de_crf;
    
    
    sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					     sim.system().pos()(j), r);

    lj_parameter_struct const *A_lj;
    lj_parameter_struct const *B_lj;
    double A_q, B_q;
    double alpha_lj, alpha_crf;
    
    const double l = sim.topology().lambda();
    
    DEBUG(7, "lambda: " << l);
    
    //assert(sim.topology().perturbed_atom().size() > j);
    
    assert(sim.topology().perturbed_solute().atoms().count(i) == 1);

    
    if(sim.topology().perturbed_atom()[j] ==true){
      // both perturbed
      assert(sim.topology().perturbed_solute().atoms().count(j) == 1);
      A_lj =  &m_base.lj_parameter(
	        sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	        sim.topology().perturbed_solute().atoms()[j].A_IAC());
      B_lj = &m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().perturbed_solute().atoms()[j].B_IAC());
      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().perturbed_solute().atoms()[j].A_charge();
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().perturbed_solute().atoms()[j].B_charge();
      
      alpha_lj = (sim.topology().perturbed_solute().atoms()[i].LJ_softcore() +
		  sim.topology().perturbed_solute().atoms()[j].LJ_softcore()) /
	2.0;
      alpha_crf = (sim.topology().perturbed_solute().atoms()[i].CRF_softcore() +
		   sim.topology().perturbed_solute().atoms()[j].CRF_softcore()) /
	2.0;
      
    }
    else{
      A_lj = &m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	       sim.topology().iac(j));
      B_lj = &m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().iac(j));
      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().charge()(j);
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().charge()(j);

      alpha_lj = sim.topology().perturbed_solute().atoms()[i].LJ_softcore();
      alpha_crf = sim.topology().perturbed_solute().atoms()[i].CRF_softcore();

    }
    
    DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 << " c12=" << A_lj->c12);
    DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 << " c12=" << B_lj->c12);
    DEBUG(7, "\tcharges state A i*j = " << A_q);
    DEBUG(7, "\tcharges state B i*j = " << B_q);
    
    m_base.lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
				   A_q, sim.topology().lambda(),
				   alpha_lj, alpha_crf,
				   A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
    
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " << A_f << " e_lj: " 
	  << A_e_lj << " e_crf: " << A_e_crf << " de_lj: " << A_de_lj 
	  << " de_crf: " << A_de_crf);
    
    m_base.lj_crf_soft_interaction(r, B_lj->c6, B_lj->c12,
				   B_q, 1.0 - sim.topology().lambda(),
				   alpha_lj, alpha_crf,
				   B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
    
    DEBUG(7, "\tcalculated interaction state B:\n\t\tf: " << B_f << " e_lj: " 
	  << B_e_lj << " e_crf: " << B_e_crf << " de_lj: " << B_de_lj 
	  << " de_crf: " << B_de_crf);
    if (do_scaling){
      
      // check whether we need to do scaling
      // based on energy groups
      DEBUG(7, "scaled interaction: " << i << " " << j);

      std::pair<int, int> 
	energy_group_pair(sim.topology().atom_energy_group(i),
			  sim.topology().atom_energy_group(j));
      if (sim.topology().energy_group_scaling().count(energy_group_pair)){
	
	// YES, we do scale the interactions!
	std::pair<double, double> scaling_pair =
	  sim.topology().energy_group_scaling()[energy_group_pair];
	A_f      *= scaling_pair.first;
	A_e_lj   *= scaling_pair.first;
	A_e_crf  *= scaling_pair.first;
	A_de_lj  *= scaling_pair.first;
	A_de_crf *= scaling_pair.first;
	B_f      *= scaling_pair.second;
	B_e_lj   *= scaling_pair.second;
	B_e_crf  *= scaling_pair.second;
	B_de_lj  *= scaling_pair.second;
	B_de_crf *= scaling_pair.second;
      }
    }
    
    // now combine everything
    const double B_l = sim.topology().lambda();
    const double B_ln = pow(B_l, sim.topology().nlam());
    const double B_lnm = pow(B_l, sim.topology().nlam()-1);
    
    const double A_l = 1.0 - sim.topology().lambda();
    const double A_ln = pow(A_l, sim.topology().nlam());
    const double A_lnm = pow(A_l, sim.topology().nlam()-1);
    
    DEBUG(7, "B_l: " << B_l <<
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
    
    //assert(m_storage.force().size() > i);
    //assert(m_storage.force().size() > j);
    
    m_storage.force()(i) += f;
    m_storage.force()(j) -= f;

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

    m_storage.energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;
    
    m_storage.energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;
    
    DEBUG(7, "\tenergy i and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));

    assert(unsigned(m_storage.lambda_energies().lj_energy.size()) > 
	   sim.topology().atom_energy_group(i));
    assert(unsigned(m_storage.lambda_energies().lj_energy.size()) >
	   sim.topology().atom_energy_group(j));
    assert(unsigned(m_storage.lambda_energies().crf_energy.size()) > 
	   sim.topology().atom_energy_group(i));
    assert(unsigned(m_storage.lambda_energies().crf_energy.size()) >
	   sim.topology().atom_energy_group(j));

    DEBUG(11, "lj dE/dLAMBDA: " << 
	  m_storage.lambda_energies().lj_energy
	  [sim.topology().atom_energy_group(i)]    
	  [sim.topology().atom_energy_group(j)]);
    DEBUG(11, "crf dE/dLAMBDA: " << 
	  m_storage.lambda_energies().crf_energy
	  [sim.topology().atom_energy_group(i)]    
	  [sim.topology().atom_energy_group(j)]);
        
    m_storage.lambda_energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_lj;

    m_storage.lambda_energies().
      crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_crf;

}


template<typename t_simulation, typename t_storage, bool do_scaling>
void interaction::
Perturbed_Nonbonded_Inner_Loop<t_simulation, t_storage, do_scaling>
::perturbed_one_four_interaction_inner_loop(t_simulation &sim,
				  size_t const i, size_t const j)
{
    DEBUG(7, "\tpair\t" << i << "\t" << j);

    math::Vec r, f, A_f, B_f;
    double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
           B_e_lj, B_e_crf, B_de_lj, B_de_crf;
    double e_lj, e_crf, de_lj, de_crf;
    

    //assert(sim.system().pos().size() > i &&
    //   sim.system().pos().size() > j);

    sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					     sim.system().pos()(j), r);

    lj_parameter_struct const *A_lj;
    lj_parameter_struct const *B_lj;
    
    double A_q, B_q;
    double alpha_lj, alpha_crf;
    
    const double l = sim.topology().lambda();

    //assert(sim.topology().perturbed_solute().atoms().size() > i &&
    //	   sim.topology().perturbed_solute().atoms().size() > j);
    
    if(sim.topology().perturbed_atom()[j] ==true){
      A_lj =  &m_base.lj_parameter(
	        sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	        sim.topology().perturbed_solute().atoms()[j].A_IAC());
      B_lj = &m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().perturbed_solute().atoms()[j].B_IAC());

      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().perturbed_solute().atoms()[j].A_charge();
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().perturbed_solute().atoms()[j].B_charge();
    
      alpha_lj = (sim.topology().perturbed_solute().atoms()[i].LJ_softcore() +
		  sim.topology().perturbed_solute().atoms()[j].LJ_softcore()) /
	2.0;
      alpha_crf = (sim.topology().perturbed_solute().atoms()[i].CRF_softcore() +
		   sim.topology().perturbed_solute().atoms()[j].CRF_softcore()) /
	2.0;

    }
    else{
      A_lj = &m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	       sim.topology().iac(j));
      B_lj = &m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().iac(j));
      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().charge()(j);
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().charge()(j);

      alpha_lj = sim.topology().perturbed_solute().atoms()[i].LJ_softcore();
      alpha_crf = sim.topology().perturbed_solute().atoms()[i].CRF_softcore();

   }
    
    DEBUG(7, "\tlj-parameter state A c6=" << A_lj->cs6 << " c12=" << A_lj->cs12);
    DEBUG(7, "\tlj-parameter state B c6=" << B_lj->cs6 << " c12=" << B_lj->cs12);
    DEBUG(7, "\tcharges state A i*j = " << A_q);
    DEBUG(7, "\tcharges state B i*j = " << B_q);
    
    m_base.lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
				   A_q, sim.topology().lambda(),
				   alpha_lj, alpha_crf,
				   A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
    
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " 
	  << A_f << " e_lj: " << A_e_lj << " e_crf: " << A_e_crf 
	  << " de_lj: " << A_de_lj  << " de_crf: " << A_de_crf);
    
    m_base.lj_crf_soft_interaction(r, B_lj->cs6, B_lj->cs12,
				   B_q, 1.0 - sim.topology().lambda(),
				   alpha_lj, alpha_crf,
				   B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
    
    DEBUG(7, "\tcalculated interaction state B:\n\t\tf: " 
	  << B_f << " e_lj: " << B_e_lj << " e_crf: " << B_e_crf 
	  << " de_lj: " << B_de_lj << " de_crf: " << B_de_crf);

    if (do_scaling){
      
      // check whether we need to do scaling
      // based on energy groups

      std::pair<int, int> 
	energy_group_pair(sim.topology().atom_energy_group(i),
			  sim.topology().atom_energy_group(j));
      if (sim.topology().energy_group_scaling().count(energy_group_pair)){
	
	// YES, we do scale the interactions!

	std::pair<double, double> scaling_pair =
	  sim.topology().energy_group_scaling()[energy_group_pair];
	A_f      *= scaling_pair.first;
	A_e_lj   *= scaling_pair.first;
	A_e_crf  *= scaling_pair.first;
	A_de_lj  *= scaling_pair.first;
	A_de_crf *= scaling_pair.first;
	B_f      *= scaling_pair.second;
	B_e_lj   *= scaling_pair.second;
	B_e_crf  *= scaling_pair.second;
	B_de_lj  *= scaling_pair.second;
	B_de_crf *= scaling_pair.second;
      }
    }
    
    // now combine everything
    const double B_l = sim.topology().lambda();
    const double B_ln = pow(B_l, sim.topology().nlam());
    const double B_lnm = pow(B_l, sim.topology().nlam()-1);

    const double A_l = 1.0 - sim.topology().lambda();
    const double A_ln = pow(A_l, sim.topology().nlam());
    const double A_lnm = pow(A_l, sim.topology().nlam()-1);
    
    f      = B_ln * B_f      + A_ln * A_f;
    e_lj   = B_ln * B_e_lj   + A_ln * A_e_lj;
    e_crf  = B_ln * B_e_crf  + A_ln * A_e_crf;
    de_lj  = B_ln * B_de_lj  + A_ln * A_de_lj  
      + sim.topology().nlam() * B_lnm * B_e_lj  
      - sim.topology().nlam() * A_lnm * A_e_lj;
    de_crf = B_ln * B_de_crf + A_ln * A_de_crf 
      + sim.topology().nlam() * B_lnm * B_e_crf 
      - sim.topology().nlam() * A_lnm * A_e_crf;

    DEBUG(11, "B_ln " << B_ln << " B_de_lj " << B_de_lj
	  << " A_ln " << A_ln << " A_de_lj " << A_de_lj
	  << " nlam " << sim.topology().nlam()
	  << " B_l " << B_l << " B_e_lj " << B_e_lj
	  << " A_l " << A_l << " A_e_lj " << A_e_lj);
    DEBUG(11, "B_ln " << B_ln << " B_de_crf " << B_de_crf
	  << " A_ln " << A_ln << " A_de_crf " << A_de_crf
	  << " nlam " << sim.topology().nlam()
	  << " B_l " << B_l << " B_e_crf " << B_e_crf
	  << " A_l " << A_l << " A_e_crf " << A_e_crf);
    

    //assert(m_storage.force().size() > i &&
    //	   m_storage.force().size() > j);
    
    m_storage.force()(i) += f;
    m_storage.force()(j) -= f;

    DEBUG(7, "\tforces stored");
    
    // energy
    //assert(m_storage.energies().lj_energy.size() > 
    //	   sim.topology().atom_energy_group(i));
    //assert(m_storage.energies().lj_energy.size() >
    //	   sim.topology().atom_energy_group(j));

    //    assert(m_storage.energies().crf_energy.size() > 
    //	   sim.topology().atom_energy_group(i));
    //assert(m_storage.energies().crf_energy.size() >
    //	   sim.topology().atom_energy_group(j));

    m_storage.energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;
    
    m_storage.energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;
    
    DEBUG(7, "\ti and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));

    assert(unsigned(m_storage.lambda_energies().lj_energy.size()) > 
    	   sim.topology().atom_energy_group(i));
    assert(unsigned(m_storage.lambda_energies().lj_energy.size()) >
    	   sim.topology().atom_energy_group(j));
    assert(unsigned(m_storage.lambda_energies().crf_energy.size()) > 
    	   sim.topology().atom_energy_group(i));
    assert(unsigned(m_storage.lambda_energies().crf_energy.size()) >
    	   sim.topology().atom_energy_group(j));

    DEBUG(11, "lj dE/dLAMBDA: " << 
	  m_storage.lambda_energies().lj_energy
	  [sim.topology().atom_energy_group(i)]    
	  [sim.topology().atom_energy_group(j)]
	  << "\t\tde_lj: " << de_lj);
    DEBUG(11, "crf dE/dLAMBDA: " << 
	  m_storage.lambda_energies().crf_energy
	  [sim.topology().atom_energy_group(i)]    
	  [sim.topology().atom_energy_group(j)]
	  << "\t\tde_crf: " << de_crf);
    
    m_storage.lambda_energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_lj;
    m_storage.lambda_energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_crf;

}

template<typename t_simulation, typename t_storage, bool do_scaling>
inline void 
interaction::Perturbed_Nonbonded_Inner_Loop<
  t_simulation, t_storage, do_scaling>
::perturbed_RF_excluded_interaction_inner_loop(t_simulation &sim,
					       std::map<size_t, simulation::Perturbed_Atom>::const_iterator const & mit)
{

  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();

  math::Vec r, f_rf, A_f_rf, B_f_rf, f_old_A;
  double e_rf, A_e_rf, B_e_rf, de_rf, A_de_rf, B_de_rf, e_crf_old_A;
  const double l=sim.topology().lambda();
  
  std::set<int>::const_iterator it, to;
  
  // self term has already been calculated for state A, 
  // correct for that and 
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

  m_base.rf_interaction(r, q_i_a*q_i_a, f_old_A, e_crf_old_A);

  m_base.rf_soft_interaction(r, q_i_a*q_i_a, 
			     B_l, mit->second.CRF_softcore(),
			     A_f_rf, A_e_rf, A_de_rf);
  
  m_base.rf_soft_interaction(r, q_i_b*q_i_b, 
			     A_l, mit->second.CRF_softcore(),
			     B_f_rf, B_e_rf, B_de_rf);
  
  if (do_scaling){

    // check whether we need to do scaling
    // based on energy groups
    DEBUG(7, "scaled interaction: (self term) " << i << " " << i);

    std::pair<int, int> 
      energy_group_pair(sim.topology().atom_energy_group(i),
			sim.topology().atom_energy_group(i));
    
    if (sim.topology().energy_group_scaling().count(energy_group_pair)){
      
	// YES, we do scale the interactions!

      std::pair<double, double> scaling_pair =
	sim.topology().energy_group_scaling()[energy_group_pair];
      A_f_rf  *= scaling_pair.first;
      A_e_rf  *= scaling_pair.first;
      A_de_rf *= scaling_pair.first;
      B_f_rf  *= scaling_pair.second;
      B_e_rf  *= scaling_pair.second;
      B_de_rf *= scaling_pair.second;
    }
  }

  DEBUG(7, "Self term for atom " << i << " A: "
	<< A_e_rf << " B: " << B_e_rf);
  
  // (1-l)^n * A + l^n *B - A 
  e_rf  = B_ln * B_e_rf  + A_ln * A_e_rf - e_crf_old_A;
  de_rf = B_ln * B_de_rf + A_ln * A_de_rf 
    + sim.topology().nlam() * B_lnm * B_e_rf 
    - sim.topology().nlam() * A_lnm * A_e_rf;
  
  sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
    [sim.topology().atom_energy_group(i)] += 0.5 * e_rf;
  sim.system().lambda_energies().crf_energy
    [sim.topology().atom_energy_group(i)]
    [sim.topology().atom_energy_group(i)] += 0.5 * de_rf;
  
  // now loop over the exclusions
  // those are fortunately not in the normal exclusions!
  it = mit->second.exclusion().begin();
  to = mit->second.exclusion().end();
  
  for( ;it != to; ++it){
    sim.system().periodicity().nearest_image(pos(i), 
					     pos(*it), r);
    DEBUG(8, "r2 i(" << i << "-" << *it << ") " << dot(r,r));
    
    double q_ij_a;
    double q_ij_b;
    double alpha_crf;
    
    if(sim.topology().perturbed_atom()[*it]){
      // j perturbed
      q_ij_a = q_i_a * 
	sim.topology().perturbed_solute().atoms()[*it].A_charge();
      q_ij_b = q_i_b *
	sim.topology().perturbed_solute().atoms()[*it].B_charge();
      
      alpha_crf = (mit->second.CRF_softcore() +
		   sim.topology().perturbed_solute().
		   atoms()[*it].CRF_softcore()) * 0.5;
    }
    else{
	// only i perturbed
      q_ij_a = q_i_a * sim.topology().charge()(*it);
      q_ij_b = q_i_b * sim.topology().charge()(*it);
      
      alpha_crf = mit->second.CRF_softcore();
      
    }
    DEBUG(8, "q_i_a " << q_i_a << " q_i_b " << q_i_b
	  << " q_ij_a " << q_ij_a << " q_ij_b " << q_ij_b
	  << " A_l " << A_l << " B_l " << B_l);
    
    m_base.rf_soft_interaction(r, q_ij_a, B_l,
			       alpha_crf,
			       A_f_rf, A_e_rf, A_de_rf);

    m_base.rf_soft_interaction(r, q_ij_b, A_l, 
			       alpha_crf,
			       B_f_rf, B_e_rf, B_de_rf);

    DEBUG(8, "alpha_crf : " << alpha_crf);
    DEBUG(7, "excluded atoms " << i << " & " << *it << " A: " 
	  << A_e_rf << " B: " << B_e_rf);
    
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

template<typename t_simulation, typename t_storage, bool do_scaling>
inline void 
interaction::Perturbed_Nonbonded_Inner_Loop<
  t_simulation, t_storage, do_scaling>
::perturbed_pair_interaction_inner_loop(t_simulation &sim,
					std::vector<simulation::Perturbed_Atompair>::const_iterator const &it)
{

  // NO SCALING for PERTURBED PAIRS ??

  math::Vec r, f, A_f, B_f;
  double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
    B_e_lj, B_e_crf, B_de_lj, B_de_crf;
  double e_lj, e_crf, de_lj, de_crf;
  lj_parameter_struct const *A_lj;
  lj_parameter_struct const *B_lj;
  double A_q, B_q;
  double alpha_lj, alpha_crf;
  
  const double l = sim.topology().lambda();
  DEBUG(7, "lambda: " << l);
  bool is_perturbed;


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
	
      A_lj =  &m_base.
	lj_parameter(sim.topology().perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     sim.topology().perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_base.
	lj_parameter(sim.topology().perturbed_solute().atoms()
		     [it->i].B_IAC(),
		     sim.topology().perturbed_solute().atoms()
		     [it->j].B_IAC());

      A_q = sim.topology().perturbed_solute().atoms()[it->i].A_charge() * 
	sim.topology().perturbed_solute().atoms()[it->j].A_charge();

      B_q = sim.topology().perturbed_solute().atoms()[it->i].B_charge() *
	sim.topology().perturbed_solute().atoms()[it->j].B_charge();

      alpha_lj = (sim.topology().perturbed_solute().atoms()
		  [it->i].LJ_softcore() +
		  sim.topology().perturbed_solute().atoms()
		  [it->j].LJ_softcore()) / 2.0;

      alpha_crf = (sim.topology().perturbed_solute().atoms()
		   [it->i].CRF_softcore() +
		   sim.topology().perturbed_solute().atoms()
		   [it->j].CRF_softcore() ) / 2.0;
      
    }
    else{ // only i perturbed

      A_lj = &m_base.
	lj_parameter(sim.topology().perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     sim.topology().iac(it->j));

      B_lj = &m_base.
	lj_parameter(sim.topology().perturbed_solute().atoms()
		     [it->i].B_IAC(),
		     sim.topology().iac(it->j));

      A_q = sim.topology().perturbed_solute().atoms()
	[it->i].A_charge() * 
	sim.topology().charge()(it->j);
      
      B_q = sim.topology().perturbed_solute().atoms()
	[it->i].B_charge() *
	sim.topology().charge()(it->j);
      
      alpha_lj = sim.topology().perturbed_solute().atoms()
	[it->i].LJ_softcore();

      alpha_crf = sim.topology().perturbed_solute().atoms()
	[it->i].CRF_softcore();
      
    }
  }
  else{ // i unperturbed
    // is j perturbed
    if(sim.topology().perturbed_atom()[it->j] == true){
      assert(sim.topology().perturbed_solute().atoms().count(it->j)
	       == 1);
	
      is_perturbed = true;

      A_lj =  &m_base.
	lj_parameter(sim.topology().iac(it->i),
		     sim.topology().perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_base.
	lj_parameter(sim.topology().iac(it->i),
		     sim.topology().perturbed_solute().atoms()
		     [it->j].B_IAC());

      A_q = sim.topology().charge()(it->i) *
	sim.topology().perturbed_solute().atoms()[it->j].A_charge();

      B_q = sim.topology().charge()(it->j) *
	sim.topology().perturbed_solute().atoms()[it->j].B_charge();

      alpha_lj = sim.topology().perturbed_solute().atoms()
	[it->j].LJ_softcore();

      alpha_crf = sim.topology().perturbed_solute().atoms()
	[it->j].CRF_softcore();
      
    }
    else{
      // both unperturbed
      
      is_perturbed = false;
	
      A_lj = &m_base.
	lj_parameter(sim.topology().iac(it->i),
		     sim.topology().iac(it->j));

      B_lj = A_lj;
      
      A_q = sim.topology().charge()(it->i) *
	sim.topology().charge()(it->j);
      B_q = A_q;
    }
          
  } // all parameters done

  // interaction in state A
  // ======================
  switch(it->A_interaction){
    // --------------------------
    case 0: // excluded
      A_e_lj = A_e_crf = A_de_lj = A_de_crf = 0.0;
      A_f = 0.0;
      DEBUG(7, "excluded in A");
      if(sim.nonbonded().RF_exclusion()){
	if (is_perturbed)
	  m_base.rf_soft_interaction(r, A_q, sim.topology().lambda(),
			       alpha_crf,
			       A_f, A_e_crf, A_de_crf);
	else
	  m_base.rf_interaction(r, A_q, A_f, A_e_crf);
      }

      break;
      // --------------------------
    case 1: // normal interaction
      DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 
	    << " c12=" << A_lj->c12);
      DEBUG(7, "\tcharges state A i*j = " << A_q);
      
      if (is_perturbed){
	DEBUG(7, "perturbed interaction");
	m_base.
	  lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
				  A_q, sim.topology().lambda(),
				  alpha_lj, alpha_crf,
				  A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
      }
      else{
	DEBUG(7, "non-perturbed interaction");
	m_base.
	  lj_crf_interaction(r, A_lj->c6, A_lj->c12,
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
	m_base.
	  lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
				  A_q, sim.topology().lambda(),
				  alpha_lj, alpha_crf,
				  A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
      }
      else{
	DEBUG(7, "non-perturbed 1,4 interaction");
	m_base.
	  lj_crf_interaction(r, A_lj->cs6, A_lj->cs12,
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
      if(sim.nonbonded().RF_exclusion()){
	if (is_perturbed)
	  m_base.rf_soft_interaction(r, B_q, sim.topology().lambda(),
			       alpha_crf,
			       B_f, B_e_crf, B_de_crf);
	else
	  m_base.rf_interaction(r, B_q, B_f, B_e_crf);
      }
 
      break;
      // --------------------------
    case 1: // normal interaction
      DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 
	    << " c12=" << B_lj->c12);
      DEBUG(7, "\tcharges state B i*j = " << B_q);
    
      if (is_perturbed){
	DEBUG(7, "perturbed interaction");
	m_base.
	  lj_crf_soft_interaction(r, B_lj->c6, B_lj->c12,
				  B_q, 1.0 - sim.topology().lambda(),
				  alpha_lj, alpha_crf,
				  B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
      }
      else{
	DEBUG(7, "non-perturbed interaction");
	m_base.lj_crf_interaction(r, B_lj->c6, B_lj->c12,
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
	m_base.
	  lj_crf_soft_interaction(r, B_lj->cs6, B_lj->cs12,
				  B_q, 1.0 - sim.topology().lambda(),
				  alpha_lj, alpha_crf,
				  B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
      }
      else{
	DEBUG(7, "non-perturbed 1,4 interaction");
	m_base.lj_crf_interaction(r, B_lj->cs6, B_lj->cs12,
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
  
  sim.system().energies().lj_energy
    [sim.topology().atom_energy_group(it->i)]
    [sim.topology().atom_energy_group(it->j)] += e_lj;
  
  sim.system().energies().crf_energy
    [sim.topology().atom_energy_group(it->i)]
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

