/**
 * @file standard_pairlist_algorithm2.tcc
 * create an atomic pairlist with a
 * chargegroup or atom based cut-off criterion.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

template<typename t_interaction_spec, typename t_perturbation_spec>
inline
interaction::Standard_Pairlist_Algorithm2<t_interaction_spec, t_perturbation_spec>::
Standard_Pairlist_Algorithm2()
  : interaction::Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>()
{
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm2<t_interaction_spec, t_perturbation_spec>::
prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  m_cutoff_short2 = sim.param().pairlist.cutoff_short * sim.param().pairlist.cutoff_short;
  m_cutoff_long2 = sim.param().pairlist.cutoff_long * sim.param().pairlist.cutoff_long;

  if (!t_interaction_spec::do_atomic_cutoff)
    calc_cg_cog(topo, conf, sim);
  
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm2<t_interaction_spec, t_perturbation_spec>::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim, 
       Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
       size_t begin, size_t end, size_t stride)
{
  DEBUG(7, "standard pairlist update");
  const double update_start = util::now();

  Periodicity_type periodicity(conf.current().box);

  // empty the pairlist
  for(size_t i=0; i<topo.num_atoms(); ++i)
    nbs.pairlist()[i].clear();

  if(t_perturbation_spec::do_perturbation){
    // and the perturbed pairlist
    for(size_t i=0; i<topo.num_atoms(); ++i)
      nbs.perturbed_pairlist()[i].clear();

  }

  // exclusion iterators
  std::set<int>::const_reverse_iterator e, e_to;

  size_t cg1 = 0, cg2;
  int cg1_atom, cg2_atom, cg2_end;

  const size_t num_cg = topo.num_chargegroups();
  const size_t num_solute_cg = topo.num_solute_chargegroups();
  
  math::Vec r;

  std::vector<int>::const_iterator
    cg1_begin = topo.chargegroups().begin(),
    cg1_end = topo.chargegroups().begin() + 1,
    cg_end = topo.chargegroups().end(),
    cg2_begin;

  // loop over all chargegroups ==> cg1
  for( ; cg1_end != cg_end; ++cg1_begin, ++cg1_end, ++cg1){
      
    //==================================================
    // add the intra chargegroup pairs
    if (cg1 < num_solute_cg){
      for(cg1_atom = *cg1_begin; cg1_atom < *cg1_end; ++cg1_atom){
	for(cg2_atom = cg1_atom + 1; cg2_atom < *cg1_end; ++cg2_atom){
	    
	  e = topo.all_exclusion(cg1_atom).rbegin();
	  e_to = topo.all_exclusion(cg1_atom).rend();
	    
	  for( ; e != e_to; ++e){
	    if (cg2_atom >  *e) goto not_excluded_intra;
	    if (cg2_atom == *e) goto excluded_intra;
	  }
	not_excluded_intra:
	  DEBUG(9, "shortrange pair " << cg1_atom << " - " << cg2_atom);
	  nbs.add_shortrange_pair(topo, conf, sim, cg1_atom, cg2_atom);
	excluded_intra:
	  {}
	}
      }
    }
      
    // intra chargegroup pairs done
    //==================================================
      
    //==================================================
    // loop over all solute chargegroup 2
    cg2_begin = cg1_end;
    cg2 = cg1 + 1;
      
    for( ; cg2 < num_solute_cg; ++cg2_begin, ++cg2){
	
      // get distance
      periodicity.nearest_image(m_cog(cg1), m_cog(cg2), r);
	
      const double dist2 = r(0)*r(0) + r(1)*r(1) + r(2)*r(2);
      
      if (cg1 == 38 && cg2 == 3642){
	std::cout << "THE CULPRIT dist2 = " << dist2 << "\n\n";
	std::cout << m_cog(cg1) << "\n"
		  << m_cog(cg2) << "\n"
		  << r << "\n"
		  << "\n\n";
      }

      DEBUG(9, "cg1=" << cg1 << " cg2=" << cg2 << " dist2 = " << dist2);

      if (dist2 > m_cutoff_long2) continue;
	
      if (dist2 > m_cutoff_short2){
	// longrange, no exclusions...

	cg2_end = *(cg2_begin + 1);
	for(cg1_atom = *cg1_begin; cg1_atom < *cg1_end; ++cg1_atom){
	  for(cg2_atom = *cg2_begin; cg2_atom < cg2_end; ++cg2_atom){

	    DEBUG(9, "longrange pair " << cg1_atom << " - " << cg2_atom);
	    nbs.add_longrange_pair(topo, conf, sim, cg1_atom, cg2_atom, periodicity);
	  }
	}
      }
      else{
	// shortrange, exclusions
	cg2_end = *(cg2_begin + 1);
	for(cg1_atom = *cg1_begin; cg1_atom < *cg1_end; ++cg1_atom){
	  for(cg2_atom = *cg2_begin; cg2_atom < cg2_end; ++cg2_atom){
	    
	    e = topo.all_exclusion(cg1_atom).rbegin();
	    e_to = topo.all_exclusion(cg1_atom).rend();
	    
	    for( ; e != e_to; ++e){
	      if (cg2_atom >  *e) goto not_excluded_inter;
	      if (cg2_atom == *e) goto excluded_inter;
	    }
	  not_excluded_inter:
	    DEBUG(9, "shortrange pair " << cg1_atom << " - " << cg2_atom);
	    nbs.add_shortrange_pair(topo, conf, sim, cg1_atom, cg2_atom);
	  excluded_inter:
	    {}
	  }
	}
      }
	
    } // loop over solute chargegroups ==> cg2
      
      // solute chargegroup - solute chargegroup interactions done
      //==================================================

      //==================================================
      // solute/solvent - solvent chargegroup interactions

    for( ; cg2 < num_cg; ++cg2_begin, ++cg2){
	
      // get distance
      periodicity.nearest_image(m_cog(cg1), m_cog(cg2), r);
	
      const double dist2 = r(0)*r(0) + r(1)*r(1) + r(2)*r(2);
	
      if (dist2 > m_cutoff_long2) continue;
	
      if (dist2 > m_cutoff_short2){
	// longrange, no exclusions...

	cg2_end = *(cg2_begin + 1);
	for(cg1_atom = *cg1_begin; cg1_atom < *cg1_end; ++cg1_atom){
	  for(cg2_atom = *cg2_begin; cg2_atom < cg2_end; ++cg2_atom){

	    DEBUG(9, "longrange pair " << cg1_atom << " - " << cg2_atom);
	    nbs.add_longrange_pair(topo, conf, sim, cg1_atom, cg2_atom, periodicity);
	  }
	}
      }
      else{
	// shortrange, no exclusions (solvent)
	cg2_end = *(cg2_begin + 1);
	for(cg1_atom = *cg1_begin; cg1_atom < *cg1_end; ++cg1_atom){
	  for(cg2_atom = *cg2_begin; cg2_atom < cg2_end; ++cg2_atom){

	    DEBUG(9, "shortrange pair " << cg1_atom << " - " << cg2_atom);
	    nbs.add_shortrange_pair(topo, conf, sim, cg1_atom, cg2_atom);
	  }
	}
      }
	
    } // loop over solute chargegroups ==> cg2
      
      // solute/solvent chargegroup - solvent chargegroup interactions done
      //==================================================

  } // loop over chargegroups ==> cg1

  m_timing += util::now() - update_start;
  
  DEBUG(7, "pairlist done");
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Standard_Pairlist_Algorithm2<t_interaction_spec, t_perturbation_spec>::
calc_cg_cog(topology::Topology & topo,
	    configuration::Configuration & conf,
	    simulation::Simulation & sim)
{

  // first put the chargegroups into the box
  math::Periodicity<t_interaction_spec::boundary_type> 
    periodicity(conf.current().box);
  
  periodicity.put_chargegroups_into_box(conf, topo);

  const size_t num_cg = topo.num_chargegroups();
  const size_t num_solute_cg = topo.num_solute_chargegroups();
  
  m_cog.resize(num_cg);
  
  size_t cg = 0, atom = 0;
  std::vector<int>::const_iterator
    cg_end = topo.chargegroups().begin() + 1;

  math::VArray const & pos = conf.current().pos;

  for( ; cg < num_solute_cg; ++cg, ++cg_end){

    double x=0, y=0, z=0;
    size_t c = 0;
    
    for( ; atom < unsigned(*cg_end); ++atom, ++c){
      
      x += pos(atom)(0);
      y += pos(atom)(1);
      z += pos(atom)(2);
      
    }
    
    m_cog(cg)(0) = x / c;
    m_cog(cg)(1) = y / c;
    m_cog(cg)(2) = z / c;
    
  }

  --cg_end;
  for( ; cg < num_cg; ++cg, ++cg_end){

    m_cog(cg)(0) = pos(*cg_end)(0);
    m_cog(cg)(1) = pos(*cg_end)(1);
    m_cog(cg)(2) = pos(*cg_end)(2);
    
  }


  // print them out... (seems to be reasonable, no solvent check yet)
  std::cout.precision(12);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  
  std::cout << "center of geometries\n"
	    << "====================\n";
  
  for(int i=0; i < m_cog.size(); ++i){
    std::cout << "cg " << i << ":"
	      << "\t" << m_cog(i)(0)
	      << "\t" << m_cog(i)(1)
	      << "\t" << m_cog(i)(2)
	      << "\n";
  }


}

