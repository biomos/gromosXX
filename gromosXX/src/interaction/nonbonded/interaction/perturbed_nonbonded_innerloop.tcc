/**
 * @file perturbed_nonbonded_innerloop.tcc
 * template methods of Perturbed_Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

template<typename t_interaction_spec,
	 typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec, t_perturbation_details>
::perturbed_lj_crf_innerloop
( topology::Topology & topo, configuration::Configuration & conf,
  size_t const i, size_t const j, Storage &storage,
  Periodicity_type const & periodicity,
  int pc)
{
  DEBUG(8, "\tperturbed pair\t" << i << "\t" << j << " (inner loop)");

  math::Vec r, f;
  double f1, f6, f12;
  
  double e_lj, e_crf, de_lj, de_crf;
  
  int energy_derivative_index = -1;

  if (t_interaction_spec::do_bekker){
    r = conf.current().pos(i) + periodicity.shift(pc).pos
      - conf.current().pos(j);
  }
  else{
    periodicity.nearest_image(conf.current().pos(i), 
			      conf.current().pos(j), r);
  }
  
  lj_parameter_struct const *A_lj;
  lj_parameter_struct const *B_lj;
  double A_q, B_q;
  double alpha_lj=0, alpha_crf=0;
  
  // const double l = topo.lambda();
  
  if(j < topo.num_solute_atoms() && 
     topo.is_perturbed(j) ==true){

    A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
				 topo.perturbed_solute().atoms()[j].A_IAC());
    B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
				 topo.perturbed_solute().atoms()[j].B_IAC());

    A_q = topo.perturbed_solute().atoms()[i].A_charge() * 
      topo.perturbed_solute().atoms()[j].A_charge();
    B_q = topo.perturbed_solute().atoms()[i].B_charge() *
      topo.perturbed_solute().atoms()[j].B_charge();
    
    alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
		topo.perturbed_solute().atoms()[j].LJ_softcore()) /
      2.0;
    alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
		 topo.perturbed_solute().atoms()[j].CRF_softcore()) /
      2.0;
    
  }
  else{
    A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
				 topo.iac(j));
    B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
				 topo.iac(j));
    DEBUG(10, "\tiac-i (A) : " << topo.perturbed_solute().atoms()[i].A_IAC() 
	  << " iac-i (B) : " << topo.perturbed_solute().atoms()[i].B_IAC()
	  << " iac-j : " << topo.iac(j));
    
    A_q = topo.perturbed_solute().atoms()[i].A_charge() * 
      topo.charge()(j);
    B_q = topo.perturbed_solute().atoms()[i].B_charge() *
      topo.charge()(j);
    
    alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
    alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();
    
  }
  
  DEBUG(8, "\tlj-parameter state A c6=" << A_lj->c6 << " c12=" << A_lj->c12);
  DEBUG(8, "\tlj-parameter state B c6=" << B_lj->c6 << " c12=" << B_lj->c12);
  DEBUG(8, "\tcharges state A i*j = " << A_q);
  DEBUG(8, "\tcharges state B i*j = " << B_q);
    

  if (t_perturbation_details::do_scaling){
    // SCALING ON

    // check whether we need to do scaling
    // based on energy groups
    
    std::pair<int, int> 
      energy_group_pair(topo.atom_energy_group(i),
			topo.atom_energy_group(j));
    bool reset_lambda = false;

    // check whether we have changing lambda dependencies
    if (topo.energy_group_lambdadep().count(energy_group_pair)){
      
      energy_derivative_index = topo.energy_group_lambdadep()[energy_group_pair].first;

      DEBUG(8, "energy_derivative_index=" << energy_derivative_index);
      assert(energy_derivative_index >= 0);

      assert(energy_derivative_index < int(topo.lambda_prime().size()));
      assert(energy_derivative_index < int(topo.lambda_prime_derivative().size()));

      // set lambdas
      DEBUG(8, "lambda dep l=" << topo.lambda() 
	    << " alpha=" << topo.energy_group_lambdadep()[energy_group_pair].second
	    << " lp=" << topo.lambda_prime()[energy_derivative_index]);
      

      set_lambda(topo.lambda_prime()[energy_derivative_index], topo.lambda_exp());
      reset_lambda = true;

    }

    if (topo.energy_group_scaling().count(energy_group_pair)){
    
    // YES, we do scale the interactions!
      lj_crf_scaled_interaction(r, A_lj->c6, A_lj->c12,
				B_lj->c6, B_lj->c12,
				A_q, B_q,
				alpha_lj, alpha_crf,
				topo.energy_group_scaling()[energy_group_pair].first,
				topo.energy_group_scaling()[energy_group_pair].second,
				f1, f6, f12,
				e_lj, e_crf, de_lj, de_crf);
    }
    else{
      // no scaling
      lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
			      B_lj->c6, B_lj->c12,
			      A_q, B_q,
			      alpha_lj, alpha_crf,
			      f1, f6, f12,
			      e_lj, e_crf, de_lj, de_crf);
    }
    
    if (reset_lambda)
      set_lambda(topo.lambda(), topo.lambda_exp());

  } 
  // END OF SCALING ON ---
  //
  else{

    lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
			    B_lj->c6, B_lj->c12,
			    A_q, B_q,
			    alpha_lj, alpha_crf,
			    f1, f6, f12,
			    e_lj, e_crf, de_lj, de_crf);
  }
  //--------------------------------------------------
  // interactions have been calculated
  //--------------------------------------------------
  
  DEBUG(8, "\tcalculated interaction state A:\n\t\tf: " 
	<< f1 << " " << f6 << " " << f12 << " e_lj: " << e_lj 
	<< " e_crf: " << e_crf 
	<< " de_lj: " << de_lj << " de_crf: " << de_crf);
  
  // now combine everything
  f      = (f1 + f6 + f12) * r;
  
  storage.force(i) += f;
  storage.force(j) -= f;
  
  DEBUG(8, "\tforces stored");
  
  if (t_interaction_spec::do_virial == math::molecular_virial){
    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
	storage.virial_tensor(a, b) += 
	  (r(a) - conf.special().rel_mol_com_pos(i)(a) + 
	   conf.special().rel_mol_com_pos(j)(a)) * f(b);
    
    DEBUG(8, "\tvirial done");
  }
  
  if (t_interaction_spec::do_virial == math::atomic_virial){
    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
	storage.virial_tensor(a, b) += 
	  r(a) * f(b);
    
    DEBUG(8, "\tatomic virial done");
  }
  
  // energy
  assert(storage.energies.lj_energy.size() > 
	 topo.atom_energy_group(i));
  assert(storage.energies.lj_energy.size() >
	 topo.atom_energy_group(j));
  
  storage.energies.lj_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_lj;
  
  storage.energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(j)] += e_crf;
  
  DEBUG(8, "\tenergy group: i and j " << topo.atom_energy_group(i)
	<< " " << topo.atom_energy_group(j)
	<< " pert der index = " << energy_derivative_index);
  
  assert(storage.perturbed_energy_derivatives.
	 lj_energy.size() > max(topo.atom_energy_group(i),
				topo.atom_energy_group(j)));
  
  assert(storage.perturbed_energy_derivatives.
	 lj_energy[topo.atom_energy_group(i)].size() 
	 > max(topo.atom_energy_group(i),
	       topo.atom_energy_group(j)));

  if (t_perturbation_details::do_scaling &&
      energy_derivative_index != -1){

    // lambda dependent energy derivatives, add
    // additional d lambda prime / d lambda factor

    storage.perturbed_energy_derivatives.lj_energy
      [topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += 
      de_lj * topo.lambda_prime_derivative()[energy_derivative_index];
    
    storage.perturbed_energy_derivatives.crf_energy
      [topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += 
      de_crf * topo.lambda_prime_derivative()[energy_derivative_index];
  }
  else{
    // standard...

    storage.perturbed_energy_derivatives.lj_energy
      [topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_lj;
    
    storage.perturbed_energy_derivatives.crf_energy
      [topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_crf;
    
  }
  
  DEBUG(7, "\tperturbed lj_crf_innerloop " << i << " - " << j << " done!");
  
}


template<typename t_interaction_spec,
	 typename t_perturbation_details>
void interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec, t_perturbation_details>
::perturbed_one_four_interaction_innerloop
( topology::Topology & topo,
  configuration::Configuration & conf,
  size_t const i, size_t const j,
  Periodicity_type const & periodicity)
{
    DEBUG(7, "\tone four pair\t" << i << "\t" << j);

    math::Vec r, f;

    double e_lj, e_crf, de_lj, de_crf;
    double f1, f6, f12;
    
    periodicity.nearest_image(conf.current().pos(i), 
			      conf.current().pos(j), r);

    lj_parameter_struct const * A_lj;
    lj_parameter_struct const * B_lj;
    double A_q, B_q;
    double alpha_lj=0, alpha_crf=0;
    
    // const double l = topo.lambda();
    
    if(topo.is_perturbed(j) ==true){
      A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
				   topo.perturbed_solute().atoms()[j].A_IAC());
      B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
				   topo.perturbed_solute().atoms()[j].B_IAC());
      A_q = topo.perturbed_solute().atoms()[i].A_charge() * 
	    topo.perturbed_solute().atoms()[j].A_charge();
      B_q = topo.perturbed_solute().atoms()[i].B_charge() *
	    topo.perturbed_solute().atoms()[j].B_charge();

      alpha_lj = (topo.perturbed_solute().atoms()[i].LJ_softcore() +
		  topo.perturbed_solute().atoms()[j].LJ_softcore()) /
	2.0;
      alpha_crf = (topo.perturbed_solute().atoms()[i].CRF_softcore() +
		   topo.perturbed_solute().atoms()[j].CRF_softcore()) /
	2.0;
      
    }
    else{
      A_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].A_IAC(),
				   topo.iac(j));
      B_lj = &m_param->lj_parameter(topo.perturbed_solute().atoms()[i].B_IAC(),
				   topo.iac(j));
      A_q = topo.perturbed_solute().atoms()[i].A_charge() * 
	    topo.charge()(j);
      B_q = topo.perturbed_solute().atoms()[i].B_charge() *
	    topo.charge()(j);

      alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
      alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

   }
    
    DEBUG(7, "\tlj-parameter state A c6=" << A_lj->cs6 
	  << " c12=" << A_lj->cs12);
    DEBUG(7, "\tlj-parameter state B c6=" << B_lj->cs6 
	  << " c12=" << B_lj->cs12);
    DEBUG(7, "\tcharges state A i*j = " << A_q);
    DEBUG(7, "\tcharges state B i*j = " << B_q);
    
    lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
			    B_lj->cs6, B_lj->cs12,
			    A_q, B_q,
			    alpha_lj, alpha_crf,
			    f1, f6, f12, 
			    e_lj, e_crf, de_lj, de_crf);
    
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " 
	  << f1 << " " << f6 << " " << f12 << " e_lj: " << e_lj 
	  << " e_crf: " << e_crf 
	  << " de_lj: " << de_lj  << " de_crf: " << de_crf);
    
    // now combine everything
    f      = (f1 + f6 + f12) * r;

    conf.current().force(i) += f;
    conf.current().force(j) -= f;

    DEBUG(7, "\tforces stored");
    
    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  conf.current().virial_tensor(a, b) += 
	    r(a) * f(b);

      DEBUG(7, "\tatomic virial done");
    }

    // energy
    assert(conf.current().energies.lj_energy.size() > 
	   topo.atom_energy_group(i));
    assert(conf.current().energies.lj_energy.size() >
	   topo.atom_energy_group(j));

    conf.current().energies.lj_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += e_lj;
    
    conf.current().energies.crf_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += e_crf;
    
    DEBUG(7, "\ti and j " << topo.atom_energy_group(i)
	  << " " << topo.atom_energy_group(j));
    DEBUG(20,"de_lj tot (before) " 
	  << conf.current().perturbed_energy_derivatives.lj_energy[topo.atom_energy_group(i)]
	  [topo.atom_energy_group(j)]);
    
    conf.current().perturbed_energy_derivatives.lj_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_lj;
    conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_crf;

}

template<typename t_interaction_spec,
	 typename t_perturbation_details>
inline void 
interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec, t_perturbation_details>
::perturbed_RF_excluded_interaction_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
 std::map<size_t, topology::Perturbed_Atom>::const_iterator const & mit,
 Periodicity_type const & periodicity)
{

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;

  math::Vec r, f_rf;
  double e_rf, de_rf;
  // const double l=topo.lambda();
  
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
  // rf_interaction(r, q_i_a*q_i_a, f_old_A, e_crf_old_A);

  rf_soft_interaction(r, q_i_a*q_i_a, q_i_b * q_i_b,
		      B_lambda(), mit->second.CRF_softcore(),
		      f_rf, e_rf, de_rf, true);
  
  DEBUG(7, "Self term for atom " << i << " = "
	<< e_rf);
    
  conf.current().energies.crf_energy[topo.atom_energy_group(i)]
    [topo.atom_energy_group(i)] += 0.5 * e_rf;
  conf.current().perturbed_energy_derivatives.crf_energy
    [topo.atom_energy_group(i)]
    [topo.atom_energy_group(i)] += 0.5 * de_rf;
  
  // now loop over the exclusions
  // those are fortunately not in the normal exclusions!
  it = mit->second.exclusion().begin();
  to = mit->second.exclusion().end();
  
  for( ;it != to; ++it){
    periodicity.nearest_image(pos(i), 
			      pos(*it), r);

    DEBUG(8, "r2 i(" << i << "-" << *it << ") " << dot(r,r));
    
    double q_ij_a;
    double q_ij_b;
    double alpha_crf=0;
    
    if(unsigned(*it) < topo.num_solute_atoms() && topo.is_perturbed(*it)){
      // j perturbed
      q_ij_a = q_i_a * 
	topo.perturbed_solute().atoms()[*it].A_charge();
      q_ij_b = q_i_b *
	topo.perturbed_solute().atoms()[*it].B_charge();
      
      alpha_crf = (mit->second.CRF_softcore() +
		   topo.perturbed_solute().
		   atoms()[*it].CRF_softcore()) * 0.5;
    }
    else{
	// only i perturbed
      q_ij_a = q_i_a * topo.charge()(*it);
      q_ij_b = q_i_b * topo.charge()(*it);
      
      alpha_crf = mit->second.CRF_softcore();
      
    }
    DEBUG(8, "q_i_a " << q_i_a << " q_i_b " << q_i_b
	  << " q_ij_a " << q_ij_a << " q_ij_b " << q_ij_b
	  << " A_l " << A_lambda() << " B_l " << B_lambda());
    
    rf_soft_interaction(r, q_ij_a, q_ij_b, B_lambda(),
			alpha_crf, f_rf, e_rf, de_rf);

    DEBUG(8, "alpha_crf : " << alpha_crf);
    DEBUG(7, "excluded atoms " << i << " & " << *it << ": " << e_rf);
    DEBUG(7, "\tde_rf: " << de_rf);
    
    // and add everything to the correct arrays 
    conf.current().energies.crf_energy 
      [topo.atom_energy_group(i)]
      [topo.atom_energy_group(*it)] += e_rf;
    conf.current().perturbed_energy_derivatives.crf_energy 
      [topo.atom_energy_group(i)]
      [topo.atom_energy_group(*it)] += de_rf;
    force(i) += f_rf;
    force(*it) -=f_rf;

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  conf.current().virial_tensor(a, b) += 
	    r(a) * f_rf(b);

      DEBUG(7, "\tatomic virial done");
    }

  }
}

