/**
 * @file nonbonded_innerloop.tcc
 * template methods of Nonbonded_Innerloop
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE interaction

#include "../../../debug.h"

template<typename t_simulation, typename t_nonbonded_spec>
interaction::Nonbonded_Innerloop<t_simulation, t_nonbonded_spec>
::Nonbonded_Innerloop(Nonbonded_Base &base)
  : m_base(base)
{
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_storage>
inline void 
interaction::Nonbonded_Innerloop<t_simulation, t_nonbonded_spec>
::interaction_innerloop(t_simulation const &sim, size_t const i, size_t const j,
			t_storage &storage)
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
    
    storage.force()(i) += f;
    storage.force()(j) -= f;
    DEBUG(7, "\tforces stored");

    if (t_nonbonded_spec::do_virial == molecular_virial){
      for(int a=0; a<3; ++a)
	for(int b=0; b<3; ++b)
	  storage.virial()(a, b) += 
	    (r(a) - sim.system().rel_mol_com_pos()(i)(a) + 
	     sim.system().rel_mol_com_pos()(j)(a)) * f(b);

      DEBUG(7, "\tvirial done");
    }
    
    // energy
    assert(storage.energies().lj_energy.size() > 
	   sim.topology().atom_energy_group(i));
    assert(storage.energies().lj_energy.size() >
	   sim.topology().atom_energy_group(j));

    storage.energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;

    storage.energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;

    DEBUG(7, "\ti and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));
}


template<typename t_simulation, typename t_nonbonded_spec>
void interaction::Nonbonded_Innerloop<t_simulation, t_nonbonded_spec>
::one_four_interaction_innerloop(t_simulation &sim,
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

    sim.system().force()(i) += f;
    sim.system().force()(j) -= f;

    // energy
    sim.system().energies().lj_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_lj;

    sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(j)] += e_crf;

    DEBUG(11, "\ti and j " << sim.topology().atom_energy_group(i)
	  << " " << sim.topology().atom_energy_group(j));

}

template<typename t_simulation, typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_simulation, t_nonbonded_spec>
::RF_excluded_interaction_innerloop(t_simulation &sim,
				    size_t const i)
{
  math::Vec r, f;
  double e_crf;
  
  math::VArray &pos   = sim.system().pos();
  math::VArray &force = sim.system().force();

  std::set<int>::const_iterator it, to;
  it = sim.topology().exclusion(i).begin();
  to = sim.topology().exclusion(i).end();
  
  DEBUG(11, "\tself-term " << i );
  r=0;
  
  // this will only contribute in the energy, the force should be zero.
  m_base.rf_interaction(r,sim.topology().charge()(i) * sim.topology().charge()(i),
			f, e_crf);
  sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
    [sim.topology().atom_energy_group(i)] += 0.5 * e_crf;
  DEBUG(11, "\tcontribution " << 0.5*e_crf);
  
  for( ; it != to; ++it){
    
    DEBUG(11, "\texcluded pair " << i << " - " << *it);
    
    sim.system().periodicity().nearest_image(pos(i), pos(*it), r);
    
    
    m_base.rf_interaction(r, sim.topology().charge()(i) * 
			  sim.topology().charge()(*it),
			  f, e_crf);
    
    force(i) += f;
    force(*it) -= f;
    
    // energy
    sim.system().energies().crf_energy[sim.topology().atom_energy_group(i)]
      [sim.topology().atom_energy_group(*it)] += e_crf;
    DEBUG(11, "\tcontribution " << e_crf);
    
  } // loop over excluded pairs
  
}

template<typename t_simulation, typename t_nonbonded_spec>
inline void 
interaction::Nonbonded_Innerloop<t_simulation, t_nonbonded_spec>
::RF_solvent_interaction_innerloop
(t_simulation &sim, simulation::chargegroup_iterator const & cg_it)
{
  math::Vec r;
  double e_crf;
  
  math::VArray &pos   = sim.system().pos();

  // loop over the atoms
  simulation::Atom_Iterator at_it = cg_it.begin(),
    at_to = cg_it.end();
  
  for ( ; at_it != at_to; ++at_it){
    DEBUG(11, "\tsolvent self term " << *at_it);
    // no solvent self term. The distance dependent part and the forces
    // are zero. The distance independent part should add up to zero 
    // for the energies and is left out.
    
    for(simulation::Atom_Iterator at2_it=at_it+1; at2_it!=at_to; ++at2_it){
      
      DEBUG(11, "\tsolvent " << *at_it << " - " << *at2_it);
      sim.system().periodicity().nearest_image(pos(*at_it), 
					       pos(*at2_it), r);
      
      // for solvent, we don't calculate internal forces (rigid molecules)
      // and the distance independent parts should go to zero
      e_crf = -sim.topology().charge()(*at_it) * 
	sim.topology().charge()(*at2_it) * 
	m_base.coulomb_constant() * 
	m_base.crf_2cut3i() * dot(r,r);
      
      // energy
      sim.system().energies().crf_energy
	[sim.topology().atom_energy_group(*at_it) ]
	[sim.topology().atom_energy_group(*at2_it)] += e_crf;
    } // loop over at2_it
  } // loop over at_it
  
}
