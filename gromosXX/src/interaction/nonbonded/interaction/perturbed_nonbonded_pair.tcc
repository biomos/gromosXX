/**
 * @file perturbed_nonbonded_pair.tcc
 * template methods of Perturbed_Nonbonded_Pair
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

template<typename t_interaction_spec, typename perturbation_details>
interaction::Perturbed_Nonbonded_Pair<t_interaction_spec, perturbation_details>
::Perturbed_Nonbonded_Pair(Nonbonded_Parameter &nbp,
			   Nonbonded_Term & nbt,
			   Perturbed_Nonbonded_Term & pnbt)
  : m_param(&nbp),
    m_nonbonded_term(&nbt),
    m_perturbed_nonbonded_term(&pnbt)
{
}

/**
 * calculate the interactions for the
 * PERTURBED PAIRS
 * (different interaction types in A and in B)
 */
template<typename t_interaction_spec, typename perturbation_details>
inline void interaction::Perturbed_Nonbonded_Pair<
  t_interaction_spec, perturbation_details>
::perturbed_pair_outerloop(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Storage & storage)
{
  DEBUG(8, "perturbed pairs");
  
  Periodicity_type periodicity(conf.current().box);
  
  std::vector<topology::perturbed_two_body_term_struct>::const_iterator
    it = topo.perturbed_solute().atompairs().begin(),
    to = topo.perturbed_solute().atompairs().end();
    
  for(; it != to; ++it){
    perturbed_pair_interaction_innerloop(topo, conf, sim, it, periodicity);
  }
  
}



template<typename t_interaction_spec, typename perturbation_details>
inline void 
interaction::Perturbed_Nonbonded_Pair<
  t_interaction_spec, perturbation_details>
::perturbed_pair_interaction_innerloop
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  std::vector<topology::perturbed_two_body_term_struct>
  ::const_iterator const &it,
  Periodicity_type const & periodicity)
{

  // NO RANGE FILTER FOR PERTURBED PAIRS ??
  // NO SCALING for PERTURBED PAIRS ??
  // NO MOLECULAR VIRIAL CONTRIBUTION ??

  math::Vec r, f, A_f, B_f;
  double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
    B_e_lj, B_e_crf, B_de_lj, B_de_crf;
  double e_lj, e_crf, de_lj, de_crf;
  lj_parameter_struct const *A_lj;
  lj_parameter_struct const *B_lj;
  double A_q, B_q;
  double alpha_lj=0, alpha_crf=0;
  
  const double l = topo.lambda();
  DEBUG(7, "lambda: " << l);
  bool is_perturbed;


  DEBUG(7, "\tperturbed-pair\t" << it->i << "\t" << it->j);
  
  periodicity.nearest_image(conf.current().pos(it->i), 
			    conf.current().pos(it->j), r);

  // is i perturbed?
  if (topo.is_perturbed(it->i)){
    assert(topo.perturbed_solute().atoms().count(it->i) == 1);
    
    is_perturbed = true;
      
    // is j perturbed as well?
    if(topo.is_perturbed(it->j)){
      assert(topo.perturbed_solute().atoms().count(it->j) == 1);
	
      A_lj =  &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     topo.perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].B_IAC(),
		     topo.perturbed_solute().atoms()
		     [it->j].B_IAC());

      A_q = topo.perturbed_solute().atoms()[it->i].A_charge() * 
	topo.perturbed_solute().atoms()[it->j].A_charge();

      B_q = topo.perturbed_solute().atoms()[it->i].B_charge() *
	topo.perturbed_solute().atoms()[it->j].B_charge();

      alpha_lj = (topo.perturbed_solute().atoms()
		  [it->i].LJ_softcore() +
		  topo.perturbed_solute().atoms()
		  [it->j].LJ_softcore()) / 2.0;

      alpha_crf = (topo.perturbed_solute().atoms()
		   [it->i].CRF_softcore() +
		   topo.perturbed_solute().atoms()
		   [it->j].CRF_softcore() ) / 2.0;
      
    }
    else{ // only i perturbed

      A_lj = &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     topo.iac(it->j));

      B_lj = &m_param->
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].B_IAC(),
		     topo.iac(it->j));

      A_q = topo.perturbed_solute().atoms()
	[it->i].A_charge() * 
	topo.charge()(it->j);
      
      B_q = topo.perturbed_solute().atoms()
	[it->i].B_charge() *
	topo.charge()(it->j);
      
      alpha_lj = topo.perturbed_solute().atoms()
	[it->i].LJ_softcore();

      alpha_crf = topo.perturbed_solute().atoms()
	[it->i].CRF_softcore();
      
    }
  }
  else{ // i unperturbed
    // is j perturbed
    if(topo.is_perturbed(it->j)){
      assert(topo.perturbed_solute().atoms().count(it->j)
	       == 1);
	
      is_perturbed = true;

      A_lj =  &m_param->
	lj_parameter(topo.iac(it->i),
		     topo.perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_param->
	lj_parameter(topo.iac(it->i),
		     topo.perturbed_solute().atoms()
		     [it->j].B_IAC());

      A_q = topo.charge()(it->i) *
	topo.perturbed_solute().atoms()[it->j].A_charge();

      B_q = topo.charge()(it->j) *
	topo.perturbed_solute().atoms()[it->j].B_charge();

      alpha_lj = topo.perturbed_solute().atoms()
	[it->j].LJ_softcore();

      alpha_crf = topo.perturbed_solute().atoms()
	[it->j].CRF_softcore();
      
    }
    else{
      // both unperturbed
      
      is_perturbed = false;
	
      A_lj = &m_param->
	lj_parameter(topo.iac(it->i),
		     topo.iac(it->j));

      B_lj = A_lj;
      
      A_q = topo.charge()(it->i) *
	topo.charge()(it->j);
      B_q = A_q;
    }
          
  } // all parameters done

  // interaction in state A
  // ======================
  switch(it->A_type){
    // --------------------------
    case 0: // excluded
      A_e_lj = A_e_crf = A_de_lj = A_de_crf = 0.0;
      A_f = 0.0;
      DEBUG(7, "excluded in A");
      if(sim.param().longrange.rf_excluded){
	if (is_perturbed){
	  m_perturbed_nonbonded_term->
	    rf_soft_interaction(r, A_q, 0, topo.lambda(),
			      alpha_crf,
			      A_f, A_e_crf, A_de_crf);
	}
	else
	  m_nonbonded_term->
	    rf_interaction(r, A_q, A_f, A_e_crf);
      }

      break;
      // --------------------------
    case 1: // normal interaction
      DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 
	    << " c12=" << A_lj->c12);
      DEBUG(7, "\tcharges state A i*j = " << A_q);
      
      if (is_perturbed){
	DEBUG(7, "perturbed interaction");
	double A_f1, A_f6, A_f12;
	
	m_perturbed_nonbonded_term->
	  lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
				0, 0,
				A_q, 0,
				alpha_lj, alpha_crf,
				A_f1, A_f6, A_f12,
				A_e_lj, A_e_crf, A_de_lj, A_de_crf);

	A_f = (A_f1 + A_f6 + A_f12) * r;
	
      }
      else{
	  double A_f1;
	DEBUG(7, "non-perturbed interaction");
	m_nonbonded_term->
	  lj_crf_interaction(r, A_lj->c6, A_lj->c12,
			     A_q, A_f1, A_e_lj, 
			     A_e_crf);

	A_de_lj = A_de_crf = 0.0;
	A_f = A_f1 * r;

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
	double A_f1, A_f6, A_f12;
	
	m_perturbed_nonbonded_term->
	  lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
				  0, 0,
				  A_q, 0,
				  alpha_lj, alpha_crf,
				  A_f1, A_f6, A_f12, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
	A_f = (A_f1 + A_f6 + A_f12) * r;
      }
      else{
	  double A_f1;
	DEBUG(7, "non-perturbed 1,4 interaction");
	m_nonbonded_term->
	  lj_crf_interaction(r, A_lj->cs6, A_lj->cs12,
			     A_q, A_f1, A_e_lj, 
			     A_e_crf);
	A_de_lj = A_de_crf = 0.0;
	A_f = A_f1 * r;
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
  switch(it->B_type){
    // --------------------------
    case 0: // excluded
      DEBUG(7, "B excluded");
      B_e_lj = B_e_crf = B_de_lj = B_de_crf = 0.0;
      B_f = 0.0;
      if(sim.param().longrange.rf_excluded){
	if (is_perturbed)
	  m_perturbed_nonbonded_term->
	    rf_soft_interaction(r, 0, B_q, topo.lambda(),
				alpha_crf,
				B_f, B_e_crf, B_de_crf);
	else
	  m_nonbonded_term->
	    rf_interaction(r, B_q, B_f, B_e_crf);
      }
 
      break;
      // --------------------------
    case 1: // normal interaction
      DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 
	    << " c12=" << B_lj->c12);
      DEBUG(7, "\tcharges state B i*j = " << B_q);
    
      if (is_perturbed){
	DEBUG(7, "perturbed interaction");
	double B_f1, B_f6, B_f12;
	
	m_perturbed_nonbonded_term->
	  lj_crf_soft_interaction(r, 0, 0, B_lj->c6, B_lj->c12,
				  0, B_q,
				  alpha_lj, alpha_crf,
				  B_f1, B_f6, B_f12,
				  B_e_lj, B_e_crf, B_de_lj, B_de_crf);
	B_f = (B_f1 + B_f6 + B_f12) * r;
      }
      else{
	  double B_f1;
	DEBUG(7, "non-perturbed interaction");
	m_nonbonded_term->
	  lj_crf_interaction(r, B_lj->c6, B_lj->c12,
			     B_q, B_f1, B_e_lj, 
			     B_e_crf);
	B_de_lj = B_de_crf = 0.0;
	B_f = B_f1 * r;
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
	double B_f1, B_f6, B_f12;

	m_perturbed_nonbonded_term->
	  lj_crf_soft_interaction(r, 0, 0, B_lj->cs6, B_lj->cs12,
				  0, B_q, alpha_lj, alpha_crf,
				  B_f1, B_f6, B_f12,
				  B_e_lj, B_e_crf, B_de_lj, B_de_crf);
	B_f = (B_f1 + B_f6 + B_f12) * r;
      }
      else{
	DEBUG(7, "non-perturbed 1,4 interaction");
	double B_f1;
	m_nonbonded_term->
	  lj_crf_interaction(r, B_lj->cs6, B_lj->cs12,
			     B_q, B_f1, B_e_lj, 
			     B_e_crf);
	B_de_lj = B_de_crf = 0.0;
	B_f = B_f1 * r;
      }
      
      DEBUG(7, "\tcalculated interaction state B:\n\t\tf: "
	    << B_f << " e_lj: " 
	    << B_e_lj << " e_crf: " << B_e_crf 
	    << " de_lj: " << B_de_lj 
	    << " de_crf: " << B_de_crf);
      
      break;
      // --------------------------
  }
  
  if (perturbation_details::do_scaling){

    // check whether we need to do scaling
    // based on energy groups
    DEBUG(7, "scaled interaction: (perturbed pair) " << it->i << " - " << it->j);

    std::pair<int, int> 
      energy_group_pair(topo.atom_energy_group(it->i),
			topo.atom_energy_group(it->j));
    
    if (topo.energy_group_scaling().count(energy_group_pair)){
      
	// YES, we do scale the interactions!

      std::pair<double, double> scaling_pair =
	topo.energy_group_scaling()[energy_group_pair];
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
  DEBUG(10, "B_l: " << m_perturbed_nonbonded_term->B_lambda() <<
	" B_ln: " << m_perturbed_nonbonded_term->B_lambda_n() <<
	" A_l: " << m_perturbed_nonbonded_term->A_lambda() <<
	" A_ln: " << m_perturbed_nonbonded_term->A_lambda_n());
  
  f      = m_perturbed_nonbonded_term->B_lambda_n() * B_f      + m_perturbed_nonbonded_term->A_lambda_n() * A_f;
  e_lj   = m_perturbed_nonbonded_term->B_lambda_n() * B_e_lj   + m_perturbed_nonbonded_term->A_lambda_n() * A_e_lj;
  e_crf  = m_perturbed_nonbonded_term->B_lambda_n() * B_e_crf  + m_perturbed_nonbonded_term->A_lambda_n() * A_e_crf;
  de_lj  = m_perturbed_nonbonded_term->B_lambda_n() * B_de_lj  + m_perturbed_nonbonded_term->A_lambda_n() * A_de_lj  
    + topo.lambda_exp() * m_perturbed_nonbonded_term->B_lambda_n_1() * B_e_lj  
    - topo.lambda_exp() * m_perturbed_nonbonded_term->A_lambda_n_1() * A_e_lj;
  
  de_crf = m_perturbed_nonbonded_term->B_lambda_n() * B_de_crf + m_perturbed_nonbonded_term->A_lambda_n() * A_de_crf 
    + topo.lambda_exp() * m_perturbed_nonbonded_term->B_lambda_n_1() * B_e_crf 
    - topo.lambda_exp() * m_perturbed_nonbonded_term->A_lambda_n_1() * A_e_crf;
  
  conf.current().force(it->i) += f;
  conf.current().force(it->j) -= f;
  
  if (t_interaction_spec::do_virial == math::atomic_virial){
    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
	conf.current().virial_tensor(a, b) += 
	  r(a) * f(b);
    
    DEBUG(7, "\tatomic virial done");
  }


  DEBUG(7, "A_lnm: " << m_perturbed_nonbonded_term->A_lambda_n_1() << " B_lnm: " << m_perturbed_nonbonded_term->B_lambda_n_1());
  DEBUG(7, "\tcalculated interaction:\n\t\tf: " << f << " e_lj: " 
	<< e_lj << " e_crf: " << e_crf << " de_lj: " << de_lj 
	<< " de_crf: " << de_crf);
  
  // energy
  //assert(m_storage.energies.lj_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.energies.lj_energy.size() >
  //   topo.atom_energy_group(j));
  //assert(m_storage.energies.crf_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.energies.crf_energy.size() >
  //   topo.atom_energy_group(j));
  
  conf.current().energies.lj_energy
    [topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)] += e_lj;
  
  conf.current().energies.crf_energy
    [topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)] += e_crf;
  
  DEBUG(7, "\tenergy i and j " << topo.atom_energy_group(it->i)
	<< " " << topo.atom_energy_group(it->j));
  
  //assert(m_storage.perturbed_energy_derivatives.lj_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.perturbed_energy_derivatives.lj_energy.size() >
  //   topo.atom_energy_group(j));
  //assert(m_storage.perturbed_energy_derivatives.crf_energy.size() > 
  //   topo.atom_energy_group(i));
  //assert(m_storage.perturbed_energy_derivatives.crf_energy.size() >
  //   topo.atom_energy_group(j));
  
  conf.current().perturbed_energy_derivatives.
    lj_energy[topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)]
    += de_lj;
  
  conf.current().perturbed_energy_derivatives.
    crf_energy[topo.atom_energy_group(it->i)]
    [topo.atom_energy_group(it->j)]
    += de_crf;
  
}

