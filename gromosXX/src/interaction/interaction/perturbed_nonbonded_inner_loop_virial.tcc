/**
 * @file perturbed_nonbonded_inner_loop_virial.tcc
 * template methods of Perturbed_Nonbonded_Inner_Loop_Virial
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

template<typename t_simulation, typename t_storage>
interaction::Perturbed_Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::Perturbed_Nonbonded_Inner_Loop_Virial(Nonbonded_Base &base, t_storage &storage)
  : interaction::Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>(base,storage)
{
}

template<typename t_simulation, typename t_storage>
void interaction::Perturbed_Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::perturbed_interaction_inner_loop(t_simulation &sim, 
				   size_t const i, size_t const j)
{
    DEBUG(7, "\tpair\t" << i << "\t" << j);

    math::Vec r, f, A_f, B_f;
    double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
           B_e_lj, B_e_crf, B_de_lj, B_de_crf;
    double e_lj, e_crf, de_lj, de_crf;
    
    
    sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					     sim.system().pos()(j), r);

    lj_parameter_struct const *A_lj;
    lj_parameter_struct const *B_lj;
    double A_q, B_q;
    const double l = sim.topology().lambda();
    
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
      
    }
    else{
      A_lj = & m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	       sim.topology().iac(j));
      B_lj = & m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().iac(j));
      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().charge()(j);
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().charge()(j);
   }
    
    DEBUG(7, "\tlj-parameter state A c6=" << A_lj->c6 << " c12=" << A_lj->c12);
    DEBUG(7, "\tlj-parameter state B c6=" << B_lj->c6 << " c12=" << B_lj->c12);
    DEBUG(7, "\tcharges state A i*j = " << A_q);
    DEBUG(7, "\tcharges state B i*j = " << B_q);
    
    m_base.lj_crf_soft_interaction(r, A_lj->c6, A_lj->c12,
				   A_q, sim.topology().lambda(),
				   A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
 
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " << A_f << " e_lj: " 
	  << A_e_lj << " e_crf: " << A_e_crf << " de_lj: " << A_de_lj 
	  << " de_crf: " << A_de_crf);
    


    m_base.lj_crf_soft_interaction(r, B_lj->c6, B_lj->c12,
				   B_q, 1.0 - sim.topology().lambda(),
				   B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
    DEBUG(7, "\tcalculated interaction state B:\n\t\tf: " << B_f << " e_lj: " 
	  << B_e_lj << " e_crf: " << B_e_crf << " de_lj: " << B_de_lj 
	  << " de_crf: " << B_de_crf);

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
    

    m_storage.force()(i) += f;
    m_storage.force()(j) -= f;

    DEBUG(7, "\tforces stored");

    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
	m_storage.virial()(a, b) += 
	  (r(a) - sim.system().rel_mol_com_pos()(i)(a) + 
	          sim.system().rel_mol_com_pos()(j)(a)) * f(b);
    
    // energy
    assert(m_storage.energies().lj_energy.size() > 
	   sim.topology().atom_energy_group(i));
    assert(m_storage.energies().lj_energy.size() >
	   sim.topology().atom_energy_group(j));

    m_storage.energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;
    
    m_storage.energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;
    
    DEBUG(7, "\ti and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));
    
    m_storage.lambda_energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_lj;
    m_storage.lambda_energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_crf;
    
}


template<typename t_simulation, typename t_storage>
void interaction::Perturbed_Nonbonded_Inner_Loop_Virial<t_simulation, 
							t_storage>
::perturbed_one_four_interaction_inner_loop(t_simulation &sim,
				  size_t const i, size_t const j)
{
    DEBUG(7, "\tpair\t" << i << "\t" << j);

    math::Vec r, f, A_f, B_f;
    double A_e_lj, A_e_crf, A_de_lj, A_de_crf, 
           B_e_lj, B_e_crf, B_de_lj, B_de_crf;
    double e_lj, e_crf, de_lj, de_crf;
    
    
    sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					     sim.system().pos()(j), r);

    lj_parameter_struct const * A_lj;
    lj_parameter_struct const * B_lj;
    double A_q, B_q;
    const double l = sim.topology().lambda();
    
    if(sim.topology().perturbed_atom()[j] ==true){
      A_lj =  & m_base.lj_parameter(
	        sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	        sim.topology().perturbed_solute().atoms()[j].A_IAC());
      B_lj = & m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().perturbed_solute().atoms()[j].B_IAC());
      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().perturbed_solute().atoms()[j].A_charge();
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().perturbed_solute().atoms()[j].B_charge();
      
    }
    else{
      A_lj = & m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].A_IAC(),
	       sim.topology().iac(j));
      B_lj = & m_base.lj_parameter(
	       sim.topology().perturbed_solute().atoms()[i].B_IAC(),
	       sim.topology().iac(j));
      A_q = sim.topology().perturbed_solute().atoms()[i].A_charge() * 
	    sim.topology().charge()(j);
      B_q = sim.topology().perturbed_solute().atoms()[i].B_charge() *
	    sim.topology().charge()(j);
   }
    
    DEBUG(7, "\tlj-parameter state A c6=" << A_lj->cs6 << " c12=" << A_lj->cs12);
    DEBUG(7, "\tlj-parameter state B c6=" << B_lj->cs6 << " c12=" << B_lj->cs12);
    DEBUG(7, "\tcharges state A i*j = " << A_q);
    DEBUG(7, "\tcharges state B i*j = " << B_q);
    
    m_base.lj_crf_soft_interaction(r, A_lj->cs6, A_lj->cs12,
				   A_q, sim.topology().lambda(),
				   A_f, A_e_lj, A_e_crf, A_de_lj, A_de_crf);
 
    DEBUG(7, "\tcalculated interaction state A:\n\t\tf: " << A_f << " e_lj: " 
	  << A_e_lj << " e_crf: " << A_e_crf << " de_lj: " << A_de_lj 
	  << " de_crf: " << A_de_crf);
    


    m_base.lj_crf_soft_interaction(r, B_lj->cs6, B_lj->cs12,
				   B_q, 1.0 - sim.topology().lambda(),
				   B_f, B_e_lj, B_e_crf, B_de_lj, B_de_crf);
    DEBUG(7, "\tcalculated interaction state B:\n\t\tf: " << B_f << " e_lj: " 
	  << B_e_lj << " e_crf: " << B_e_crf << " de_lj: " << B_de_lj 
	  << " de_crf: " << B_de_crf);

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

    m_storage.force()(i) += f;
    m_storage.force()(j) -= f;

    DEBUG(7, "\tforces stored");
    
    // energy
    assert(m_storage.energies().lj_energy.size() > 
	   sim.topology().atom_energy_group(i));
    assert(m_storage.energies().lj_energy.size() >
	   sim.topology().atom_energy_group(j));

    m_storage.energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;
    
    m_storage.energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;
    
    DEBUG(7, "\ti and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));
    
    m_storage.lambda_energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_lj;
    m_storage.lambda_energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += de_crf;


    /*    DEBUG(10, "\t1,4-pair\t" << i << "\t" << j);

    math::Vec r, f;
    double e_lj, e_crf;
    
    sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					     sim.system().pos()(j), r);

    const lj_parameter_struct &lj = 
      m_base.lj_parameter(sim.topology().iac(i),
		   sim.topology().iac(j));

    DEBUG(11, "\tlj-parameter cs6=" << lj.cs6 << " cs12=" << lj.cs12);

    m_base.lj_crf_interaction(r, lj.cs6, lj.cs12,
			    sim.topology().charge()(i) * 
			    sim.topology().charge()(j),
			    f, e_lj, e_crf);

    m_storage.force()(i) += f;
    m_storage.force()(j) -= f;

    // energy
    m_storage.energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;

    m_storage.energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;

    DEBUG(11, "\ti and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));
    */
}


