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

#include "../../debug.h"

template<typename t_simulation, typename t_nonbonded_spec>
inline
interaction::Grid_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>::
Grid_Pairlist_Algorithm()
  : interaction::Standard_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>()
{
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Grid_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>::
update(t_simulation &sim, t_nonbonded_interaction & nonbonded_interaction)
{
  DEBUG(7, "pairlist update");
   
  // empty the pairlist
  nonbonded_interaction.pairlist().clear();
  nonbonded_interaction.pairlist().resize(sim.topology().num_atoms());

  if(t_nonbonded_spec::do_perturbation){
    // and the perturbed pairlist
    nonbonded_interaction.perturbed_pairlist().clear();
    nonbonded_interaction.perturbed_pairlist().resize(sim.topology().num_atoms());
  }
  
  DEBUG(7, "pairlist(s) resized");
  
  // prepare the range filter (cutoff)
  set_cutoff(sim.nonbonded().cutoff_short(), sim.nonbonded().cutoff_long());
  
  // prepare the range filter (center of geometries)    
  prepare_cog(sim);
  DEBUG(7, "range filter prepared (cog)");

  DEBUG(7, "create a grid");
  Chargegroup_Grid<t_simulation> 
    a_grid(sim.system().periodicity(), 0.5*m_cutoff_short, m_cutoff_long);

  DEBUG(7, "grid the cog's");
  grid_cog(sim, a_grid);
  
  int num_cells[3];
  int the_cell[3];
  
  a_grid.num_cells(num_cells);
  
  std::vector<size_t>::const_iterator cg_st, cg_it, cg_it2, cg_to;

  // loop over all cells
  DEBUG(8, "num_cell: " << num_cells[0] << " " << num_cells[1] << " " << num_cells[2]);
  
  for(the_cell[0]=0; the_cell[0] < num_cells[0]; ++the_cell[0]){
    for(the_cell[1]=0; the_cell[1] < num_cells[1]; ++the_cell[1]){
      DEBUG(8, "num_cell: " << num_cells[0] << " " << num_cells[1] << " " << num_cells[2]);
      for(the_cell[2]=0; the_cell[2] < num_cells[2]; ++the_cell[2]){
	
	DEBUG(8, "cc " << std::setw(3) << the_cell[0] 
	      << " | " << std::setw(3) << the_cell[1] 
	      << " | " << std::setw(3) << the_cell[2]);

	DEBUG(8, "number of cg (of this grid point): " 
	      << a_grid.grid()[the_cell[0]][the_cell[1]][the_cell[2]].size());
	      
	cg_st = a_grid.grid()[the_cell[0]][the_cell[1]][the_cell[2]].begin();
	cg_to = a_grid.grid()[the_cell[0]][the_cell[1]][the_cell[2]].end();

	intra_cell(sim, nonbonded_interaction, cg_st, cg_to);

	// only in central box
	inter_cell<t_nonbonded_interaction, false>
	  (sim, nonbonded_interaction, cg_st, cg_to, a_grid, the_cell);

	// and also the periodic images
	inter_cell<t_nonbonded_interaction, true>
	  (sim, nonbonded_interaction, cg_st, cg_to, a_grid, the_cell);
      }
    }
  }
  
  // and that's it...
  DEBUG(7, "pairlist done");
}

template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction>
inline void
interaction::Grid_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>::
intra_cell(t_simulation &sim, t_nonbonded_interaction & nonbonded_interaction,
	   std::vector<size_t>::const_iterator &cg_st, 
	   std::vector<size_t>::const_iterator &cg_to)
{

  DEBUG(8, "intra cell");

  std::vector<size_t>::const_iterator cg_it, cg_it2;
  
  // do the intra - cell
  for(cg_it = cg_st; cg_it != cg_to; ++cg_it){
    
    simulation::chargegroup_iterator cg1 =  
      sim.topology().chargegroup_begin();
    cg1 += *cg_it;
    
    if (unsigned(**cg1) < sim.topology().solute().num_atoms()){
      // not solvent, intra cg
      DEBUG(8, "intra cg: " << *cg_it);
      
      do_cg_interaction_intra(sim, nonbonded_interaction, cg1);	  
    }
    
    // second chargegroups in the cell
    for(cg_it2 = cg_it + 1;  cg_it2 != cg_to; ++cg_it2){

      DEBUG(8, "intra cell: " << *cg_it << " -- " << *cg_it2);
      
      simulation::chargegroup_iterator cg2 =  
	sim.topology().chargegroup_begin();
      cg2 += *cg_it2;
      
      // ASSUME SHORTRANGE
      if (unsigned(**cg2) < sim.topology().solute().num_atoms()){
	// exclusions! (because cg2 is not solvent)
	do_cg_interaction_excl(sim, nonbonded_interaction, cg1, cg2);
      }
      else{
	// no exclusions... (at least cg2 is solvent)
	do_cg_interaction(sim, nonbonded_interaction, cg1, cg2);
      }
    }
  }
}


template<typename t_simulation, typename t_nonbonded_spec>
template<typename t_nonbonded_interaction, bool periodic>
inline void
interaction::Grid_Pairlist_Algorithm<t_simulation, t_nonbonded_spec>::
inter_cell(t_simulation &sim, t_nonbonded_interaction & nonbonded_interaction,
	   std::vector<size_t>::const_iterator &cg_st, 
	   std::vector<size_t>::const_iterator &cg_to,
	   Chargegroup_Grid<t_simulation> &grid,
	   int cell[3])
{
 
  DEBUG(8, "inter cell");
  
  std::vector<size_t>::const_iterator cg_it;
  
  int start = 13;
  int end = 13;
  
  if (periodic){
    // only half...
    start = 14;
    end = 26;
  }
  
  for(int pc=start; pc <= end; ++pc){

    DEBUG(11, "periodic copy: " << pc << " periodic= " << periodic);
    DEBUG(11, "shift cell: " << sim.system().periodicity().shift(pc).cell[0] << " / "
	  << sim.system().periodicity().shift(pc).cell[1] << " / " 
	  << sim.system().periodicity().shift(pc).cell[2]);
    
    // shift the central cell
    int our_cell[3];
    for(int d=0; d<3; ++d)
      our_cell[d] = cell[d] + sim.system().periodicity().shift(pc).cell[d];

    // get an iterator over the mask
    Cell_Cell_Iterator<t_simulation, periodic> it(grid, our_cell);
    if (it.eol()) {
      DEBUG(8, "no interactions in range..."); 
      continue;
    }
    
    do{
      
      for(cg_it = cg_st; cg_it != cg_to; ++cg_it){
	DEBUG(8, "inter cell " << *cg_it << "  --  " << *it);
	
	// get the chargegroup iterators
	simulation::chargegroup_iterator cg1 =  
	  sim.topology().chargegroup_begin();
	
	simulation::chargegroup_iterator cg2 =  
	  sim.topology().chargegroup_begin();

	cg1 += *cg_it;
	cg2 += *it;
	
	if (!t_nonbonded_spec::do_atomic_cutoff){
	  // filter out interactions based on chargegroup distances
	  if (range_chargegroup_pair(sim, nonbonded_interaction,
				     *cg_it, *it, cg1, cg2, 
				     sim.system().periodicity().shift(pc)))
	    continue;
	}
	
	DEBUG(11, "\tshortrange!");
	
	// SHORTRANGE
	if (unsigned(**cg2) >= sim.topology().solute().num_atoms() || 
	    unsigned(**cg1) >= sim.topology().solute().num_atoms()){
	  // no exclusions... (at least cg2 is solvent)
	  do_cg_interaction(sim, nonbonded_interaction, cg1, cg2);
	}
	else{
	  // exclusions! (because no cg is solvent)
	  if (*cg_it <= *it)
	    do_cg_interaction_excl(sim, nonbonded_interaction, cg1, cg2);
	  else
	    do_cg_interaction_inv_excl(sim, nonbonded_interaction, cg1, cg2);
	}

      }
    } while (++it);
    DEBUG(8, "out of loop");
    
  } // end of loop over periodic copies
  
}
