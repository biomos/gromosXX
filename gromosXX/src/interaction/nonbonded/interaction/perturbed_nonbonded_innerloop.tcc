/**
 * @file perturbed_nonbonded_innerloop.tcc
 * template methods of Perturbed_Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include <util/debug.h>

template<typename t_interaction_spec>
interaction::Perturbed_Nonbonded_Innerloop<t_interaction_spec>
::Perturbed_Nonbonded_Innerloop(Nonbonded_Base &base)
  : m_base(base)
{
}

template<typename t_interaction_spec>
template<typename t_storage>
void interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec>
::perturbed_interaction_innerloop
( topology::Topology & topo, configuration::Configuration & conf,
  size_t const i, size_t const j, t_storage &storage,
  math::Periodicity<t_interaction_spec::boundary_type> const & periodicity,
  int pc)
{
    DEBUG(7, "\tperturbed pair\t" << i << "\t" << j);

    math::Vec r, f, A_f, B_f;
    double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
           B_e_lj, B_e_crf, B_de_lj, B_de_crf;
    double e_lj, e_crf, de_lj, de_crf;
    
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
    
    if(topo.is_perturbed(j) ==true){
      A_lj =  &m_base.lj_parameter(
	        topo.perturbed_solute().atoms()[i].A_IAC(),
	        topo.perturbed_solute().atoms()[j].A_IAC());
      B_lj = &m_base.lj_parameter(
	       topo.perturbed_solute().atoms()[i].B_IAC(),
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
      A_lj = & m_base.lj_parameter(
	       topo.perturbed_solute().atoms()[i].A_IAC(),
	       topo.iac(j));
      B_lj = & m_base.lj_parameter(
	       topo.perturbed_solute().atoms()[i].B_IAC(),
	       topo.iac(j));
      A_q = topo.perturbed_solute().atoms()[i].A_charge() * 
	    topo.charge()(j);
      B_q = topo.perturbed_solute().atoms()[i].B_charge() *
	    topo.charge()(j);

      alpha_lj = topo.perturbed_solute().atoms()[i].LJ_softcore();
      alpha_crf = topo.perturbed_solute().atoms()[i].CRF_softcore();

   }
    
    DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 << " c12=" << A_lj->c12);
    DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 << " c12=" << B_lj->c12);
    DEBUG(7, "\tcharges state A i*j = " << A_q);
    DEBUG(7, "\tcharges state B i*j = " << B_q);

    m_base.lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
				   A_q, topo.lambda(),
				   alpha_lj, alpha_crf,
				   A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
    
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " 
	  << A_f << " e_lj: " << A_e_lj << " e_crf: " << A_e_crf 
	  << " de_lj: " << A_de_lj << " de_crf: " << A_de_crf);
    
    m_base.lj_crf_soft_interaction(r, B_lj->c6, B_lj->c12,
				   B_q, 1.0 - topo.lambda(),
				   alpha_lj, alpha_crf,
				   B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
    
    DEBUG(7, "\tcalculated interaction state B:\n\t\tf: " 
	  << B_f << " e_lj: " << B_e_lj << " e_crf: " << B_e_crf 
	  << " de_lj: " << B_de_lj << " de_crf: " << B_de_crf);
    
    if (t_interaction_spec::do_scaling){
      
      // check whether we need to do scaling
      // based on energy groups

      std::pair<int, int> 
	energy_group_pair(topo.atom_energy_group(i),
			  topo.atom_energy_group(j));
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
    f      = m_base.B_lambda_n() * B_f      + m_base.A_lambda_n() * A_f;
    e_lj   = m_base.B_lambda_n() * B_e_lj   + m_base.A_lambda_n() * A_e_lj;
    e_crf  = m_base.B_lambda_n() * B_e_crf  + m_base.A_lambda_n() * A_e_crf;
    de_lj  = m_base.B_lambda_n() * B_de_lj  + m_base.A_lambda_n() * A_de_lj  
      + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_lj  
      - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_lj;
    de_crf = m_base.B_lambda_n() * B_de_crf + m_base.A_lambda_n() * A_de_crf 
      + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_crf 
      - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_crf;
    

    storage.force(i) += f;
    storage.force(j) -= f;

    DEBUG(7, "\tforces stored");

    if (t_interaction_spec::do_virial == math::molecular_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  storage.virial_tensor(a, b) += 
	    (r(a) - conf.special().rel_mol_com_pos(i)(a) + 
	     conf.special().rel_mol_com_pos(j)(a)) * f(b);
      
      DEBUG(7, "\tvirial done");
    }

    if (t_interaction_spec::do_virial == math::atomic_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  storage.virial_tensor(a, b) += 
	    r(a) * f(b);

      DEBUG(7, "\tatomic virial done");
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
    
    DEBUG(7, "\ti and j " << topo.atom_energy_group(i)
	  << " " << topo.atom_energy_group(j));
    
    storage.perturbed_energy_derivatives.lj_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_lj;
    storage.perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_crf;
    
}


template<typename t_interaction_spec>
void interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec>
::perturbed_one_four_interaction_innerloop
( topology::Topology & topo,
  configuration::Configuration & conf,
  size_t const i, size_t const j,
  math::Periodicity<t_interaction_spec::boundary_type> const & periodicity)
{
    DEBUG(7, "\tpair\t" << i << "\t" << j);

    math::Vec r, f, A_f, B_f;
    double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
           B_e_lj, B_e_crf, B_de_lj, B_de_crf;
    double e_lj, e_crf, de_lj, de_crf;
    
    
    periodicity.nearest_image(conf.current().pos(i), 
			      conf.current().pos(j), r);

    lj_parameter_struct const * A_lj;
    lj_parameter_struct const * B_lj;
    double A_q, B_q;
    double alpha_lj=0, alpha_crf=0;
    
    // const double l = topo.lambda();
    
    if(topo.is_perturbed(j) ==true){
      A_lj =  & m_base.lj_parameter(
	        topo.perturbed_solute().atoms()[i].A_IAC(),
	        topo.perturbed_solute().atoms()[j].A_IAC());
      B_lj = & m_base.lj_parameter(
	       topo.perturbed_solute().atoms()[i].B_IAC(),
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
      A_lj = & m_base.lj_parameter(
	       topo.perturbed_solute().atoms()[i].A_IAC(),
	       topo.iac(j));
      B_lj = & m_base.lj_parameter(
	       topo.perturbed_solute().atoms()[i].B_IAC(),
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
    
    m_base.lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
				   A_q, topo.lambda(),
				   alpha_lj, alpha_crf,
				   A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
    
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " 
	  << A_f << " e_lj: " << A_e_lj << " e_crf: " << A_e_crf 
	  << " de_lj: " << A_de_lj  << " de_crf: " << A_de_crf);
    
    m_base.lj_crf_soft_interaction(r, B_lj->cs6, B_lj->cs12,
				   B_q, 1.0 - topo.lambda(),
				   alpha_lj, alpha_crf,
				   B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
    
    DEBUG(7, "\tcalculated interaction state B:\n\t\tf: " 
	  << B_f << " e_lj: " << B_e_lj << " e_crf: " << B_e_crf 
	  << " de_lj: " << B_de_lj << " de_crf: " << B_de_crf);
    
    if (t_interaction_spec::do_scaling){
      
      // check whether we need to do scaling
      // based on energy groups

      std::pair<int, int> 
	energy_group_pair(topo.atom_energy_group(i),
			  topo.atom_energy_group(j));
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
    f      = m_base.B_lambda_n() * B_f      + m_base.A_lambda_n() * A_f;
    e_lj   = m_base.B_lambda_n() * B_e_lj   + m_base.A_lambda_n() * A_e_lj;
    e_crf  = m_base.B_lambda_n() * B_e_crf  + m_base.A_lambda_n() * A_e_crf;
    de_lj  = m_base.B_lambda_n() * B_de_lj  + m_base.A_lambda_n() * A_de_lj  
      + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_lj  
      - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_lj;
    de_crf = m_base.B_lambda_n() * B_de_crf + m_base.A_lambda_n() * A_de_crf 
      + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_crf 
      - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_crf;

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
    
    conf.current().perturbed_energy_derivatives.lj_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_lj;
    conf.current().perturbed_energy_derivatives.crf_energy[topo.atom_energy_group(i)]
      [topo.atom_energy_group(j)] += de_crf;

}

template<typename t_interaction_spec>
inline void 
interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec>
::perturbed_RF_excluded_interaction_innerloop
(topology::Topology & topo, configuration::Configuration & conf,
 std::map<size_t, topology::Perturbed_Atom>::const_iterator const & mit,
 math::Periodicity<t_interaction_spec::boundary_type> const & periodicity)
{

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;

  math::Vec r, f_rf, A_f_rf, B_f_rf, f_old_A;
  double e_rf, A_e_rf, B_e_rf, de_rf, A_de_rf, B_de_rf, e_crf_old_A;
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
  m_base.rf_interaction(r, q_i_a*q_i_a, f_old_A, e_crf_old_A);

  m_base.rf_soft_interaction(r, q_i_a*q_i_a, 
			     m_base.B_lambda(), mit->second.CRF_softcore(),
			     A_f_rf, A_e_rf, A_de_rf);
  
  m_base.rf_soft_interaction(r, q_i_b*q_i_b, 
			     m_base.A_lambda(), mit->second.CRF_softcore(),
			     B_f_rf, B_e_rf, B_de_rf);
  
  if (t_interaction_spec::do_scaling){

    // check whether we need to do scaling
    // based on energy groups
    DEBUG(7, "scaled interaction: (self term) " << i << " " << i);

    std::pair<int, int> 
      energy_group_pair(topo.atom_energy_group(i),
			topo.atom_energy_group(i));
    
    if (topo.energy_group_scaling().count(energy_group_pair)){
      
	// YES, we do scale the interactions!

      std::pair<double, double> scaling_pair =
	topo.energy_group_scaling()[energy_group_pair];
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
  e_rf  = m_base.B_lambda_n() * B_e_rf  + m_base.A_lambda_n() * A_e_rf - e_crf_old_A;
  de_rf = m_base.B_lambda_n() * B_de_rf + m_base.A_lambda_n() * A_de_rf 
    + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_rf 
    - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_rf;
  
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
    
    if(topo.is_perturbed(*it)){
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
	  << " A_l " << m_base.A_lambda() << " B_l " << m_base.B_lambda());
    
    m_base.rf_soft_interaction(r, q_ij_a, m_base.B_lambda(),
			       alpha_crf,
			       A_f_rf, A_e_rf, A_de_rf);

    m_base.rf_soft_interaction(r, q_ij_b, m_base.A_lambda(), 
			       alpha_crf,
			       B_f_rf, B_e_rf, B_de_rf);

    DEBUG(8, "alpha_crf : " << alpha_crf);
    DEBUG(7, "excluded atoms " << i << " & " << *it << " A: " 
	  << A_e_rf << " B: " << B_e_rf);
    
    e_rf  = m_base.B_lambda_n() * B_e_rf  + m_base.A_lambda_n() * A_e_rf;
    de_rf = m_base.B_lambda_n() * B_de_rf + m_base.A_lambda_n() * A_de_rf 
      + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_rf 
      - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_rf;
    f_rf = m_base.B_lambda_n() * B_f_rf + m_base.A_lambda_n() * A_f_rf;
    
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

template<typename t_interaction_spec>
inline void 
interaction::Perturbed_Nonbonded_Innerloop<
  t_interaction_spec>
::perturbed_pair_interaction_innerloop
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  std::vector<topology::perturbed_two_body_term_struct>
  ::const_iterator const &it,
  math::Periodicity<t_interaction_spec::boundary_type> const & periodicity)
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
	
      A_lj =  &m_base.
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     topo.perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_base.
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

      A_lj = &m_base.
	lj_parameter(topo.perturbed_solute().atoms()
		     [it->i].A_IAC(),
		     topo.iac(it->j));

      B_lj = &m_base.
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

      A_lj =  &m_base.
	lj_parameter(topo.iac(it->i),
		     topo.perturbed_solute().atoms()
		     [it->j].A_IAC());

      B_lj = &m_base.
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
	
      A_lj = &m_base.
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
	if (is_perturbed)
	  m_base.rf_soft_interaction(r, A_q, topo.lambda(),
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
				  A_q, topo.lambda(),
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
				  A_q, topo.lambda(),
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
  switch(it->B_type){
    // --------------------------
    case 0: // excluded
      DEBUG(7, "B excluded");
      B_e_lj = B_e_crf = B_de_lj = B_de_crf = 0.0;
      B_f = 0.0;
      if(sim.param().longrange.rf_excluded){
	if (is_perturbed)
	  m_base.rf_soft_interaction(r, B_q, topo.lambda(),
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
				  B_q, 1.0 - topo.lambda(),
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
				  B_q, 1.0 - topo.lambda(),
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
  
  if (t_interaction_spec::do_scaling){

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
  DEBUG(10, "B_l: " << m_base.B_lambda() <<
	" B_ln: " << m_base.B_lambda_n() <<
	" A_l: " << m_base.A_lambda() <<
	" A_ln: " << m_base.A_lambda_n());
  
  f      = m_base.B_lambda_n() * B_f      + m_base.A_lambda_n() * A_f;
  e_lj   = m_base.B_lambda_n() * B_e_lj   + m_base.A_lambda_n() * A_e_lj;
  e_crf  = m_base.B_lambda_n() * B_e_crf  + m_base.A_lambda_n() * A_e_crf;
  de_lj  = m_base.B_lambda_n() * B_de_lj  + m_base.A_lambda_n() * A_de_lj  
    + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_lj  
    - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_lj;
  
  de_crf = m_base.B_lambda_n() * B_de_crf + m_base.A_lambda_n() * A_de_crf 
    + topo.lambda_exp() * m_base.B_lambda_n_1() * B_e_crf 
    - topo.lambda_exp() * m_base.A_lambda_n_1() * A_e_crf;
  
  conf.current().force(it->i) += f;
  conf.current().force(it->j) -= f;
  
  if (t_interaction_spec::do_virial == math::atomic_virial){
    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
	conf.current().virial_tensor(a, b) += 
	  r(a) * f(b);
    
    DEBUG(7, "\tatomic virial done");
  }


  DEBUG(7, "A_lnm: " << m_base.A_lambda_n_1() << " B_lnm: " << m_base.B_lambda_n_1());
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
