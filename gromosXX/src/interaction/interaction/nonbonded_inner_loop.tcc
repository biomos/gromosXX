/**
 * @file nonbonded_inner_loop.tcc
 * template methods of Nonbonded_Inner_Loop
 */

template<typename t_simulation, typename t_pairlist, typename t_storage>
interaction::Nonbonded_Inner_Loop<t_simulation, t_pairlist, t_storage>
::Nonbonded_Inner_Loop(t_storage &store, Nonbonded_Base &base)
  : m_base(base),
    m_storage(storage)
{
}

template<typename t_simulation, typename t_pairlist, typename t_storage>
void interaction::Nonbonded_Inner_Loop<t_simulation, t_pairlist, t_storage>
::do_interaction(t_simulation &sim, typename t_pairlist::iterator &it)
{
    DEBUG(10, "\tpair\t" << it.i() << "\t" << *it);

    math::Vec r, f;
    double e_lj, e_crf;
    
    sim.system().periodicity().nearest_image(sim.system().pos(it.i()), 
					     sim.system().pos(*it), r);

    const lj_parameter_struct &lj = 
      lj_parameter(sim.topology().iac(it.i()),
		   sim.topology().iac(*it));

    DEBUG(11, "\tlj-parameter c6=" << lj.c6 << " c12=" << lj.c12);

    base.lj_crf_interaction(r, lj.c6, lj.c12,
			    sim.topology().charge()(it.i()) * 
			    sim.topology().charge()(*it),
			    f, e_lj, e_crf);

    m_store.force()(it.i()) += f;
    m_store.force()(*it) -= f;

    // energy
    m_store.energy().lj_energy[sim.topology().atom_energy_group(it.i())]
      [sim.topology().atom_energy_group(*it)] += e_lj;

    m_store.energy().crf_energy[sim.topology().atom_energy_group(it.i())]
      [sim.topology().atom_energy_group(*it)] += e_crf;

    DEBUG(11, "\ti and j " << sim.topology().atom_energy_group(it.i())
	  << " " << sim.topology().atom_energy_group(*it));
}

// this is copy paste code = BAD !!!
template<typename t_simulation, typename t_pairlist, typename t_storage>
void interaction::Nonbonded_Inner_Loop<t_simulation, t_pairlist, t_storage>
::do_one_four_interaction(t_simulation &sim, 
			 typename t_pairlist::iterator &it)
{
    DEBUG(10, "\t1,4-pair\t" << it.i() << "\t" << *it);

    math::Vec r, f;
    double e_lj, e_crf;
    
    sim.system().periodicity().nearest_image(sim.system().pos(it.i()), 
					     sim.system().pos(*it), r);

    const lj_parameter_struct &lj = 
      lj_parameter(sim.topology().iac(it.i()),
		   sim.topology().iac(*it));

    DEBUG(11, "\tlj-parameter cs6=" << lj.sc6 << " cs12=" << lj.cs12);

    base.lj_crf_interaction(r, lj.cs6, lj.cs12,
			    sim.topology().charge()(it.i()) * 
			    sim.topology().charge()(*it),
			    f, e_lj, e_crf);

    m_store.force()(it.i()) += f;
    m_store.force()(*it) -= f;

    // energy
    m_store.energy().lj_energy[sim.topology().atom_energy_group(it.i())]
      [sim.topology().atom_energy_group(*it)] += e_lj;

    m_store.energy().crf_energy[sim.topology().atom_energy_group(it.i())]
      [sim.topology().atom_energy_group(*it)] += e_crf;

    DEBUG(11, "\ti and j " << sim.topology().atom_energy_group(it.i())
	  << " " << sim.topology().atom_energy_group(*it));

}


