/**
 *  grid_pairlist_algorithm.tcc
 * create an atomic pairlist with a
 * chargegroup or an atom based cut-off criterion.
 * using a grid.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include <util/debug.h>

template<typename t_interaction_spec, typename t_perturbation_spec>
inline
interaction::Grid_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
Grid_Pairlist_Algorithm()
  : interaction::Standard_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>()
{
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Grid_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
prepare(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim)
{
  
  // prepare the range filter (cutoff)
  set_cutoff(sim.param().pairlist.cutoff_short,
	     sim.param().pairlist.cutoff_long);

  // prepare the range filter (center of geometries)    
  prepare_cog(topo, conf, sim);
  DEBUG(7, "range filter prepared (cog)");
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Grid_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
update(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
       size_t begin, size_t end, size_t stride)
{
  DEBUG(7, "pairlist update");
   
  const double update_start = util::now();

  // empty the pairlist
  for(size_t i=0; i<topo.num_atoms(); ++i)
    nbs.pairlist()[i].clear();

  if(t_perturbation_spec::do_perturbation){
    // and the perturbed pairlist
    for(size_t i=0; i<topo.num_atoms(); ++i)
      nbs.perturbed_pairlist()[i].clear();

  }

  if(t_perturbation_spec::do_perturbation){
    // and the perturbed pairlist
    nbs.perturbed_pairlist().clear();
    nbs.perturbed_pairlist().resize(topo.num_atoms());
  }
  
  DEBUG(7, "pairlist(s) resized");
  
  
  Periodicity_type periodicity(conf.current().box);


  DEBUG(7, "create a grid");
  Chargegroup_Grid_type 
    a_grid(periodicity, 0.5*m_cutoff_short, m_cutoff_long);

  DEBUG(7, "grid the cog's");
  grid_cog(topo, conf, sim, a_grid);
  
  int num_cells[3];
  int the_cell[3];
  
  a_grid.num_cells(num_cells);
  
  std::vector<size_t>::const_iterator cg_st, cg_to;

  // loop over all cells
  DEBUG(8, "num_cell: " << num_cells[0] << " " 
	<< num_cells[1] << " " << num_cells[2]);

  int x, y, z;

  // cells in x
  for(x=0; x < num_cells[0]; ++x){

    the_cell[0] = x;
    // cells in y
    for(y=0; y < num_cells[1]; ++y){
      the_cell[1] = y;
      // cells in z
      for(z=0; z < num_cells[2]; ++z){
	the_cell[2] = z;

	DEBUG(10, "cc " << std::setw(3) << the_cell[0] 
	      << " | " << std::setw(3) << the_cell[1] 
	      << " | " << std::setw(3) << the_cell[2]);

	DEBUG(10, "number of cg (of this grid point): " 
	      << a_grid.grid()[the_cell[0]][the_cell[1]][the_cell[2]].size());
	      
	// iterator over the chargegroups in the cell
	cg_st = a_grid.grid()[the_cell[0]][the_cell[1]][the_cell[2]].begin();
	cg_to = a_grid.grid()[the_cell[0]][the_cell[1]][the_cell[2]].end();

	// intra-cell interaction
	intra_cell(topo, conf, sim, nbs, cg_st, cg_to, periodicity);

	// only in central box
	inter_cell<false>(topo, conf, sim, nbs, cg_st, cg_to, a_grid, the_cell,
	   periodicity);

	// and also the periodic images
	inter_cell<true>(topo, conf, sim, nbs, cg_st, cg_to, a_grid, the_cell,
	   periodicity);
      }
    }
  }
  
  // and that's it...
  m_timing += util::now() - update_start;

  DEBUG(7, "pairlist done");
}

template<typename t_interaction_spec, typename t_perturbation_spec>
inline void
interaction::Grid_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
intra_cell(topology::Topology & topo,
	   configuration::Configuration & conf,
	   simulation::Simulation & sim,
	   Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
	   std::vector<size_t>::const_iterator &cg_st, 
	   std::vector<size_t>::const_iterator &cg_to,
	   Periodicity_type const & periodicity)
{

  DEBUG(12, "intra cell");

  std::vector<size_t>::const_iterator cg_it, cg_it2;
  
  // do the intra - cell
  for(cg_it = cg_st; cg_it != cg_to; ++cg_it){
    
    topology::Chargegroup_Iterator cg1 =  
      topo.chargegroup_begin();
    cg1 += *cg_it;
    
    if (unsigned(**cg1) < topo.solute().num_atoms()){
      // not solvent, intra cg
      DEBUG(12, "intra cg: " << *cg_it);
      
      do_cg_interaction_intra(topo, conf, sim, nbs, cg1,
			      periodicity, 13);
    }
    
    // second chargegroups in the cell
    for(cg_it2 = cg_it + 1;  cg_it2 != cg_to; ++cg_it2){

      DEBUG(12, "intra cell: " << *cg_it << " -- " << *cg_it2);
      
      topology::Chargegroup_Iterator cg2 =  
	topo.chargegroup_begin();
      cg2 += *cg_it2;
      
      // ASSUME SHORTRANGE
      if (unsigned(**cg2) < topo.solute().num_atoms()){
	// exclusions! (because cg2 is not solvent)
	do_cg_interaction_excl(topo, conf, sim, nbs, cg1, cg2, periodicity, 13);
      }
      else{
	// no exclusions... (at least cg2 is solvent)
	do_cg_interaction(topo, conf, sim, nbs, cg1, cg2, periodicity, 13);
      }
    }
  }
}


template<typename t_interaction_spec, typename t_perturbation_spec>
template<bool periodic>
inline void
interaction::Grid_Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec>::
inter_cell(topology::Topology & topo,
	   configuration::Configuration & conf,
	   simulation::Simulation & sim,
	   Nonbonded_Set<t_interaction_spec, t_perturbation_spec> & nbs,
	   std::vector<size_t>::const_iterator &cg_st, 
	   std::vector<size_t>::const_iterator &cg_to,
	   Chargegroup_Grid_type & grid,
	   int cell[3],
	   Periodicity_type const & periodicity)
{
 
  DEBUG(12, "inter cell");
  
  std::vector<size_t>::const_iterator cg_it;
  
  int start = 13;
  int end = 13;
  
  if (periodic){
    // only half...
    start = 14;
    end = 26;
  }
  
  for(int pc=start; pc <= end; ++pc){

    DEBUG(13, "periodic copy: " << pc << " periodic= " << periodic);
    DEBUG(13, "shift cell: " << periodicity.shift(pc).cell[0] << " / "
	  << periodicity.shift(pc).cell[1] << " / " 
	  << periodicity.shift(pc).cell[2]);
    
    // shift the central cell
    int our_cell[3];
    for(int d=0; d<3; ++d)
      our_cell[d] = cell[d] + periodicity.shift(pc).cell[d];

    // get an iterator over the mask
    Cell_Cell_Iterator<t_interaction_spec::boundary_type, periodic>
      it(grid, our_cell);
    if (it.eol()) {
      DEBUG(12, "no interactions in range..."); 
      continue;
    }
    
    do{
      
      for(cg_it = cg_st; cg_it != cg_to; ++cg_it){
	DEBUG(12, "inter cell " << *cg_it << "  --  " << *it);
	
	// get the chargegroup iterators
	topology::Chargegroup_Iterator cg1 =  
	  topo.chargegroup_begin();
	
	topology::Chargegroup_Iterator cg2 =  
	  topo.chargegroup_begin();

	cg1 += *cg_it;
	cg2 += *it;
	
	if (!t_interaction_spec::do_atomic_cutoff){
	  // filter out interactions based on chargegroup distances
	  if (range_chargegroup_pair(topo, conf, sim, nbs, *cg_it, *it, cg1, cg2, 
				     pc, periodicity))
	    continue;
	}
	else{
	  throw std::string("grid based pairlist not implemented"
			    "for atomic cutoff ???");
	}
	
	DEBUG(13, "\tshortrange!");
	
	// SHORTRANGE
	if (unsigned(**cg2) >= topo.solute().num_atoms() || 
	    unsigned(**cg1) >= topo.solute().num_atoms()){
	  // no exclusions... (at least cg2 is solvent)
	  do_cg_interaction(topo, conf, sim, nbs, cg1, cg2, periodicity, pc);
	}
	else{
	  // exclusions! (because no cg is solvent)
	  if (*cg_it <= *it)
	    do_cg_interaction_excl(topo, conf, sim, nbs, cg1, cg2,
				   periodicity, pc);
	  else
	    do_cg_interaction_inv_excl(topo, conf, sim, nbs, cg1, cg2,
				       periodicity, pc);
	}

      }
    } while (++it);
    DEBUG(12, "out of loop");
    
  } // end of loop over periodic copies
  
}
