/**
 * @file nonbonded_inner_loop_virial.tcc
 * template methods of Nonbonded_Inner_Loop_Virial
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../debug.h"

template<typename t_simulation, typename t_storage>
interaction::Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::Nonbonded_Inner_Loop_Virial(Nonbonded_Base &base, t_storage &storage)
  : interaction::Nonbonded_Inner_Loop<t_simulation, t_storage>(base,storage)
{
}

template<typename t_simulation, typename t_storage>
void interaction::Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::interaction_inner_loop(t_simulation const &sim, size_t const i, size_t const j)
{
    DEBUG(7, "\tpair\t" << i << "\t" << j);

    math::Vec r, f;
    double e_lj, e_crf;
    
    sim.system().periodicity().nearest_image(sim.system().pos()(i), 
					     sim.system().pos()(j), r);

    const lj_parameter_struct &lj = 
      m_base.lj_parameter(sim.topology().iac(i),
			  sim.topology().iac(j));

    DEBUG(7, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

    m_base.lj_crf_interaction(r, lj.c6, lj.c12,
			      sim.topology().charge()(i) * 
			      sim.topology().charge()(j),
			      f, e_lj, e_crf);

    DEBUG(7, "\tcalculated interaction f: " << f << " e_lj: " << e_lj << " e_crf: " << e_crf);
    
    m_storage.force()(i) += f;
    m_storage.force()(j) -= f;
    DEBUG(7, "\tforces stored");

    for(int a=0; a<3; ++a)
      for(int b=0; b<3; ++b)
	m_storage.virial()(a, b) += 
	  (r(a) - sim.system().rel_mol_com_pos()(i)(a) + 
	          sim.system().rel_mol_com_pos()(j)(a)) * f(b);

    DEBUG(7, "\tvirial done");

    
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
}


template<typename t_simulation, typename t_storage>
void interaction::Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::one_four_interaction_inner_loop(t_simulation &sim,
				  size_t const i, size_t const j)
{
    DEBUG(10, "\t1,4-pair\t" << i << "\t" << j);

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

}

template<typename t_simulation, typename t_storage>
void interaction::Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::perturbed_interaction_inner_loop(t_simulation &sim, 
				   size_t const i, size_t const j)
{
}


template<typename t_simulation, typename t_storage>
void interaction::Nonbonded_Inner_Loop_Virial<t_simulation, t_storage>
::perturbed_one_four_interaction_inner_loop(t_simulation &sim,
				  size_t const i, size_t const j)
{
}


