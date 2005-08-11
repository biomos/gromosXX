/**
 * @file grid_pairlist_algorithm.cc
 * grid pairlist algorithm
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <interaction/interaction_types.h>
#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <interaction/nonbonded/interaction/nonbonded_parameter.h>

#include <interaction/nonbonded/interaction/nonbonded_term.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_term.h>

#include <interaction/nonbonded/interaction/nonbonded_innerloop.h>
#include <interaction/nonbonded/interaction/perturbed_nonbonded_innerloop.h>

#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>
#include <interaction/nonbonded/pairlist/grid_pairlist_algorithm.h>

#include <interaction/nonbonded/interaction_spec.h>
#include <interaction/nonbonded/innerloop_template.h>

#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::Grid_Pairlist_Algorithm::Grid_Pairlist_Algorithm()
  : interaction::Pairlist_Algorithm()
{
}

int interaction::Grid_Pairlist_Algorithm::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::ostream & os,
 bool quiet)
{

  if (!sim.param().pairlist.grid){
    io::messages.add("Grid_Pairlist_Algorithm",
		     "wrong input parameters",
		     io::message::error);
    return 1;
  }
  
  if (conf.boundary_type != math::rectangular){
    io::messages.add("Grid Pairlist Algorithm",
		     "only implemented for rectangular boundary conditions",
		     io::message::error);
    return 1;
  }
  
  if (sim.param().pcouple.virial != math::molecular_virial &&
      sim.param().pcouple.virial != math::no_virial){
    io::messages.add("Grid Pairlist Algorithm",
		     "only implemented for molecular virial",
		     io::message::error);
    return 1;
  }

  if (sim.param().multicell.multicell){
    io::messages.add("Grid Pairlist Algorithm",
		     "not compatible with MULTICELL simulations",
		     io::message::error);
  }

  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);

  grid_properties(topo, conf, sim);

  if (!quiet){
    std::cout << "GridPairlistAlgorithm\n"
	      << "\tcells             " 
	      << std::setw(10) << m_grid.Na
	      << std::setw(10) << m_grid.Nb
	      << std::setw(10) << m_grid.Nc << "\n"
	      << "\textenced cells    "
	      << std::setw(10) << m_grid.Na_ex
	      << std::setw(10) << m_grid.Nb_ex
	      << std::setw(10) << m_grid.Nc_ex << "\n"
	      << "\tcell size         "
	      << std::setw(10) << m_grid.a
	      << std::setw(10) << m_grid.b
	      << std::setw(10) << m_grid.c << "\n";

    const int Ncell = m_grid.Na * m_grid.Nb * m_grid.Nc;
    const int N = m_grid.Na * m_grid.Nb;

    const double P = topo.num_chargegroups();
    const double V = math::volume(conf.current().box, conf.boundary_type);
    // const double Vcell = m_grid.a * m_grid.b * m_grid.c;

    const double Vcut = 4.0 / 3.0 * math::Pi * m_cutoff_long_2 * m_cutoff_long;
    const double Pcut = P * Vcut / V;
    
    const double Pcell = P / Ncell;
    const double Player = P / m_grid.Nc;
    
    std::cout << "\tparticles / cell    " << std::setw(10) << Pcell << "\n"
	      << "\tparticles / layer   " << std::setw(10) << Player << "\n";

    // mask size:
    int Nmask = 0;
    for(int z=0; z<=m_grid.mask_z; ++z){
      for(unsigned int y=0; y<m_grid.mask[z].size(); y+=2){
	Nmask += m_grid.mask[z][y+1] - m_grid.mask[z][y] - 1;
      }
    }
    
    std::cout << "\tcells in mask       " << std::setw(10) << Nmask << "\n"
	      << "\tpairs               " << std::setw(10) << 0.5 * P * P << "\n"
	      << "\tpairs (grid)        " << std::setw(10) << P * Nmask * Pcell << "\n"
	      << "\tpairs (cutoff)      " << std::setw(10) << 0.5 * P * Pcut << "\n";

    // just for fun, already try this
    if (prepare_grid(topo, conf, sim)){
      std::cout << "\terror during grid preparation!\n";
      io::messages.add("Grid_Pairlist_Algorithm",
		       "error during grid preparation",
		       io::message::error);
      return 1;
    }
    // collapse_grid();

    int occupied = 0;
    for(int z=0; z<m_grid.Nc; ++z){
      for(int i=0; i<N; ++i){
	if (m_grid.count[z][i])
	  ++occupied;
      }
    }

    std::cout << "\toccupied            " << std::setw(10) << occupied << "\n"
	      << "\tparticle / occ cell " << std::setw(10) << double(P) / occupied << "\n";
    

    // print_mask();

    std::cout << "END\n";
  }

  return 0;
}

/**
 * calculate grid properties,
 * put chargegroups on grid,
 * including the center of geometries
 */
int interaction::Grid_Pairlist_Algorithm::prepare
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(7, "standard pairlist algorithm : prepare");
  
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
  
  // first put the chargegroups into the box
  // _prepare_cog(conf, topo);

  grid_properties(topo, conf, sim);

  if (prepare_grid(topo, conf, sim)){
    return 1;
  }

  collapse_grid();

  // print_grid();

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// pairlist update
////////////////////////////////////////////////////////////////////////////////
void interaction::Grid_Pairlist_Algorithm::update
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 unsigned int begin,
 unsigned int end,
 unsigned int stride
 )
{
  if (sim.param().pairlist.atomic_cutoff){
    // see standard_pairlist_algorithm_atomic.cc
    io::messages.add("Grid based pairlist with atomic cutoff not implemented",
        "GridPairlistAlgorithm",
        io::message::critical);
    assert(false);
  }
  else{
    SPLIT_INNERLOOP(_update, topo, conf, sim, storage,
		    pairlist, begin, end, stride);
  }    
}

template<typename t_interaction_spec>
void interaction::Grid_Pairlist_Algorithm::_update
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 unsigned int begin,
 unsigned int end,
 unsigned int stride
 )
{
  Nonbonded_Innerloop<t_interaction_spec> innerloop(*m_param);
  innerloop.init(sim);

  // empty the pairlist
  for(unsigned int i=0; i<topo.num_atoms(); ++i)
    pairlist[i].clear();
  
  const double update_start = util::now();

  std::vector<Grid::Particle> p_plane;
  std::vector<int> cell_start;

  const int c_ex = (m_grid.Nc_ex - m_grid.Nc) / 2;
  const int N = m_grid.Na * m_grid.Nb;
  const int num_solute_cg = topo.num_solute_chargegroups();

  for(int z=begin; z < m_grid.Nc + c_ex; z+=stride){

    prepare_plane(z, p_plane, cell_start);
    // print_plane(z, p_plane, cell_start);

    // std::cerr << "plane " << z << " printed..." << std::endl;
    
    if (z < m_grid.Nc){
      // do self interaction
      const int i_first = m_grid.cell_start[z][0];
      const int i_last  = m_grid.cell_start[z][N];

      int base = -1;
      int start = -1;

      for(int i=i_first; i < i_last; ++i){

	if (m_grid.p_cell[z][i].i < num_solute_cg){ // self interaction
	  DEBUG(8, "self " << m_grid.p_cell[z][i].i);
	  for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		a_to = topo.chargegroups()[m_grid.p_cell[z][i].i + 1];
	      a1 < a_to; ++a1){
	    for(int a2 = a1+1; a2 < a_to; ++a2){
	      if (excluded_solute_pair(topo, a1, a2))
		continue;
	      pairlist[a1].push_back(a2);
	    }
	  }
	}
	
	if (base != m_grid.cell_index[z][i]){
	  base = m_grid.cell_index[z][i];
	  start = i;
	}
	else{ // more than 1 cg in the cell
	  if (m_grid.p_cell[z][i].i < num_solute_cg){ // solute - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      if (m_grid.p_cell[z][j].i < num_solute_cg){ // solute - solute

		const int ii = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][i].i : m_grid.p_cell[z][j].i;
		const int jj = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][j].i : m_grid.p_cell[z][i].i;

		for(int a1 = topo.chargegroups()[ii],
		      a_to = topo.chargegroups()[ii + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[jj],
			a2_to = topo.chargegroups()[jj + 1];
		      a2 < a2_to; ++a2){
		    
		    if (excluded_solute_pair(topo, a1, a2))
		      continue;
		    pairlist[a1].push_back(a2);
		  }
		}
	      }
	      else{ // solute - solvent
		for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
			a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		      a2 < a2_to; ++a2){
		    pairlist[a1].push_back(a2);
		  }
		}

	      }
	    }
	  }
	  else{ // solvent - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		    a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		  a1 < a_to; ++a1){
		for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
		      a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		    a2 < a2_to; ++a2){
		  pairlist[a2].push_back(a1);
		}
	      }
	    }
	  }
	}
      } // loop over plane
    } // self interaction

    ////////////////////////////////////////////////////////////////////////////////
    // 
    // INTER CELL INTERACTIONS
    //
    ////////////////////////////////////////////////////////////////////////////////

    // loop over all mask levels (inside box)
    for(int mask_z=0; mask_z < int(m_grid.mask.size()); ++mask_z){
      if (z - mask_z < 0) break;
      if (z - mask_z >= m_grid.Nc) continue;
      
      const int i_level = z - mask_z;

      const int i_first = m_grid.cell_start[i_level][0];
      const int i_last  = m_grid.cell_start[i_level][N];

      // loop over all particles in this level
      for(int i=i_first; i < i_last; ++i){
	
	assert(m_grid.cell_index[i_level].size() > unsigned(i));
	const int base = m_grid.cell_index[i_level][i];
	
	if (m_grid.p_cell[i_level][i].i < num_solute_cg){ // solute - ?
	  
	  assert(m_grid.mask.size() > unsigned(mask_z));
	  for(unsigned int m=0; m<m_grid.mask[mask_z].size(); m+=2){
	    
	    assert(cell_start.size() > unsigned(base + m_grid.mask[mask_z][m]));
	  
	    const int first = cell_start[base + m_grid.mask[mask_z][m]];
	    const int last  = cell_start[base + m_grid.mask[mask_z][m+1]];
	    
	    for(int j=first; j<last; ++j){
	      
	      const double d2 = 
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) *
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) +
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) *
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) +
		(m_grid.p_cell[i_level][i].z - p_plane[j].z) *
		(m_grid.p_cell[i_level][i].z - p_plane[j].z);

	      if (d2 > m_cutoff_long_2) continue;

	      // no self interaction here...
	      DEBUG(8, "inter cell: " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      assert(m_grid.p_cell[i_level][i].i != p_plane[j].i);

	      if (p_plane[j].i < num_solute_cg){ 
		// solute - solute

		if (d2 > m_cutoff_short_2){

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						       m_shift_vector[p_plane[j].shift_index]);
		    }
		  }
		}
		else{
		  
		  const int ii = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    m_grid.p_cell[i_level][i].i : p_plane[j].i;
		  const int jj = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    p_plane[j].i : m_grid.p_cell[i_level][i].i;
		  
		  DEBUG(8, "rewritten: " << ii
			<< " - " << jj);
		  
		  for(int a1 = topo.chargegroups()[ii],
			a_to = topo.chargegroups()[ii + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[jj],
			  a2_to = topo.chargegroups()[jj + 1];
			a2 < a2_to; ++a2){
		      
		      if (excluded_solute_pair(topo, a1, a2))
			continue;
		      pairlist[a1].push_back(a2);
		    }
		  }
		}
	      }
	      else{ // solute - solvent
		if (d2 > m_cutoff_short_2){

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						       m_shift_vector[p_plane[j].shift_index]);
		    }
		  }
		}
		else{

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      pairlist[a1].push_back(a2);
		    }
		  }
		}
		
	      }
	      
	    } // j in mask row
	  }  // mask
	} // solute - ?
	else{ // solvent - ?
	  assert(m_grid.mask.size() > unsigned(mask_z));
	  for(unsigned int m=0; m<m_grid.mask[mask_z].size(); m+=2){
	    
	    assert(cell_start.size() > unsigned(base + m_grid.mask[mask_z][m]));
	  
	    const int first = cell_start[base + m_grid.mask[mask_z][m]];
	    const int last  = cell_start[base + m_grid.mask[mask_z][m+1]];
	    
	    for(int j=first; j<last; ++j){

	      const double d2 = 
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) *
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) +
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) *
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) +
		(m_grid.p_cell[i_level][i].z - p_plane[j].z) *
		(m_grid.p_cell[i_level][i].z - p_plane[j].z);

	      if (d2 > m_cutoff_long_2) continue;

	      DEBUG(8, "inter (solvent - ?): " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      DEBUG(8, "i_level=" << i_level << " d2=" << d2);
	      
	      if (d2 > m_cutoff_short_2){
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){

		    // maybe for fast solvent loop, check cg of j...
		    innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						     m_shift_vector[p_plane[j].shift_index]);
		  }
		}
	      }
	      else{
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){
		    pairlist[a2].push_back(a1);
		  }
		}
	      }
	      
	    }
	    
	  }

	} // solvent - ?

	// std::cout << "\n";
	
      } // particle in the current mask plane

    } // mask levels

  } // the extended planes

  m_timing += util::now() - update_start;
}

////////////////////////////////////////////////////////////////////////////////
// perturbed pairlist update
////////////////////////////////////////////////////////////////////////////////

void interaction::Grid_Pairlist_Algorithm::update_perturbed
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 interaction::Pairlist & perturbed_pairlist,
 unsigned int begin,
 unsigned int end, 
 unsigned int stride
 )
{
  SPLIT_PERT_INNERLOOP(_update_perturbed,
		       topo, conf, sim,
		       storage, pairlist, perturbed_pairlist,
		       begin, end, stride);
}


template<typename t_interaction_spec, typename t_perturbation_details>
void interaction::Grid_Pairlist_Algorithm::_update_perturbed
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::Storage & storage,
 interaction::Pairlist & pairlist,
 interaction::Pairlist & perturbed_pairlist,
 unsigned int begin,
 unsigned int end, 
 unsigned int stride
 )
{
  math::Periodicity<t_interaction_spec::boundary_type> periodicity(conf.current().box);

  Nonbonded_Innerloop<t_interaction_spec> innerloop(*m_param);
  innerloop.init(sim);

  Perturbed_Nonbonded_Innerloop<t_interaction_spec, t_perturbation_details>
    perturbed_innerloop(*m_param);
  perturbed_innerloop.init(sim);
  perturbed_innerloop.set_lambda(topo.lambda(), topo.lambda_exp());

  // empty the pairlist
  assert(pairlist.size() == topo.num_atoms());
  assert(perturbed_pairlist.size() == topo.num_atoms());
  for(unsigned int i=0; i<topo.num_atoms(); ++i){
    pairlist[i].clear();
    perturbed_pairlist[i].clear();
  }
  
  const double update_start = util::now();

  std::vector<Grid::Particle> p_plane;
  std::vector<int> cell_start;

  const int c_ex = (m_grid.Nc_ex - m_grid.Nc) / 2;
  const int N = m_grid.Na * m_grid.Nb;
  const int num_solute_cg = topo.num_solute_chargegroups();

  for(int z=begin; z < m_grid.Nc + c_ex; z+=stride){

    prepare_plane(z, p_plane, cell_start);
    // print_plane(z, p_plane, cell_start);
    
    if (z < m_grid.Nc){
      // do self interaction
      const int i_first = m_grid.cell_start[z][0];
      const int i_last  = m_grid.cell_start[z][N];

      int base = -1;
      int start = -1;

      for(int i=i_first; i < i_last; ++i){

	if (m_grid.p_cell[z][i].i < num_solute_cg){ // self interaction
	  DEBUG(8, "self " << m_grid.p_cell[z][i].i);
	  for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		a_to = topo.chargegroups()[m_grid.p_cell[z][i].i + 1];
	      a1 < a_to; ++a1){
	    for(int a2 = a1+1; a2 < a_to; ++a2){
	      if (excluded_solute_pair(topo, a1, a2))
		continue;
	      if (insert_pair<t_perturbation_details>
		  (topo, pairlist, perturbed_pairlist,
		   a1, a2, sim.param().perturbation.scaled_only))
		;
	      else if (insert_pair<t_perturbation_details>
		       (topo, pairlist, perturbed_pairlist,
			a2, a1, sim.param().perturbation.scaled_only))
		;
	      else
		pairlist[a1].push_back(a2);
	    }
	  }
	}
	
	if (base != m_grid.cell_index[z][i]){
	  base = m_grid.cell_index[z][i];
	  start = i;
	}
	else{ // more than 1 cg in the cell
	  if (m_grid.p_cell[z][i].i < num_solute_cg){ // solute - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      if (m_grid.p_cell[z][j].i < num_solute_cg){ // solute - solute

		const int ii = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][i].i : m_grid.p_cell[z][j].i;
		const int jj = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][j].i : m_grid.p_cell[z][i].i;

		for(int a1 = topo.chargegroups()[ii],
		      a_to = topo.chargegroups()[ii + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[jj],
			a2_to = topo.chargegroups()[jj + 1];
		      a2 < a2_to; ++a2){
		    
		    if (excluded_solute_pair(topo, a1, a2))
		      continue;

		    if (insert_pair<t_perturbation_details>
			(topo, pairlist, perturbed_pairlist,
			 a1, a2, sim.param().perturbation.scaled_only))
		      ;
		    else if (insert_pair<t_perturbation_details>
			     (topo, pairlist, perturbed_pairlist,
			      a2, a1, sim.param().perturbation.scaled_only))
		      ;
		    else
		      pairlist[a1].push_back(a2);
		  }
		}
	      }
	      else{ // solute - solvent
		for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
			a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		      a2 < a2_to; ++a2){

		    if (insert_pair<t_perturbation_details>
			(topo, pairlist, perturbed_pairlist,
			 a1, a2, sim.param().perturbation.scaled_only))
		      ;
		    else
		      pairlist[a1].push_back(a2);
		  }
		}

	      }
	    }
	  }
	  else{ // solvent - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		    a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		  a1 < a_to; ++a1){
		for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
		      a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		    a2 < a2_to; ++a2){

		  if (insert_pair<t_perturbation_details>
		      (topo, pairlist, perturbed_pairlist,
		       a2, a1, sim.param().perturbation.scaled_only))
		    ;
		  else
		    pairlist[a2].push_back(a1);
		}
	      }
	    }
	  }
	}
      } // loop over plane
    } // self interaction

    ////////////////////////////////////////////////////////////////////////////////
    // 
    // INTER CELL INTERACTIONS
    //
    ////////////////////////////////////////////////////////////////////////////////

    // loop over all mask levels (inside box)
    for(int mask_z=0; mask_z < int(m_grid.mask.size()); ++mask_z){
      if (z - mask_z < 0) break;
      if (z - mask_z >= m_grid.Nc) continue;
      
      const int i_level = z - mask_z;

      const int i_first = m_grid.cell_start[i_level][0];
      const int i_last  = m_grid.cell_start[i_level][N];

      // loop over all particles in this level
      for(int i=i_first; i < i_last; ++i){
	
	assert(m_grid.cell_index[i_level].size() > unsigned(i));
	const int base = m_grid.cell_index[i_level][i];
	
	if (m_grid.p_cell[i_level][i].i < num_solute_cg){ // solute - ?
	  
	  assert(m_grid.mask.size() > unsigned(mask_z));
	  for(unsigned int m=0; m<m_grid.mask[mask_z].size(); m+=2){
	    
	    assert(cell_start.size() > unsigned(base + m_grid.mask[mask_z][m]));
	  
	    const int first = cell_start[base + m_grid.mask[mask_z][m]];
	    const int last  = cell_start[base + m_grid.mask[mask_z][m+1]];
	    
	    for(int j=first; j<last; ++j){
	      
	      const double d2 = 
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) *
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) +
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) *
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) +
		(m_grid.p_cell[i_level][i].z - p_plane[j].z) *
		(m_grid.p_cell[i_level][i].z - p_plane[j].z);

	      if (d2 > m_cutoff_long_2) continue;

	      // no self interaction here...
	      DEBUG(8, "inter cell: " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      assert(m_grid.p_cell[i_level][i].i != p_plane[j].i);

	      if (p_plane[j].i < num_solute_cg){ 
		// solute - solute

		if (d2 > m_cutoff_short_2){

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      if (calculate_pair<t_interaction_spec, t_perturbation_details>
			  (topo, conf, storage, innerloop, perturbed_innerloop, a1, a2,
			   m_shift_vector[p_plane[j].shift_index],
			   sim.param().perturbation.scaled_only))
			;
		      // careful: shift index needs to change here !!!
		      else if (calculate_pair<t_interaction_spec, t_perturbation_details>
			       (topo, conf, storage, innerloop, perturbed_innerloop, a2, a1,
				m_reverse_shift_vector[p_plane[j].shift_index],
				sim.param().perturbation.scaled_only))
			;
		      else
			innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
							 m_shift_vector[p_plane[j].shift_index]);
		    }
		  }
		}
		else{
		  // exclusions => order
		  const int ii = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    m_grid.p_cell[i_level][i].i : p_plane[j].i;
		  const int jj = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    p_plane[j].i : m_grid.p_cell[i_level][i].i;
		  
		  DEBUG(8, "rewritten: " << ii
			<< " - " << jj);
		  
		  for(int a1 = topo.chargegroups()[ii],
			a_to = topo.chargegroups()[ii + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[jj],
			  a2_to = topo.chargegroups()[jj + 1];
			a2 < a2_to; ++a2){
		      
		      if (excluded_solute_pair(topo, a1, a2))
			continue;

		      if (insert_pair<t_perturbation_details>
			  (topo, pairlist, perturbed_pairlist,
			   a1, a2, sim.param().perturbation.scaled_only))
			;
		      else if (insert_pair<t_perturbation_details>
			       (topo, pairlist, perturbed_pairlist,
				a2, a1, sim.param().perturbation.scaled_only))
			;
		      else
			pairlist[a1].push_back(a2);
		    }
		  }
		}
	      }
	      else{ // solute - solvent
		if (d2 > m_cutoff_short_2){

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      if (calculate_pair<t_interaction_spec, t_perturbation_details>
			  (topo, conf, storage, innerloop, perturbed_innerloop, a1, a2,
			   m_shift_vector[p_plane[j].shift_index],
			   sim.param().perturbation.scaled_only))
			;
		      else
			innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
							 m_shift_vector[p_plane[j].shift_index]);
		    }
		  }
		}
		else{

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){

		      if (insert_pair<t_perturbation_details>
			  (topo, pairlist, perturbed_pairlist,
			   a1, a2, sim.param().perturbation.scaled_only))
			;
		      else
			pairlist[a1].push_back(a2);
		    }
		  }
		}
		
	      }
	      
	    } // j in mask row
	  }  // mask
	} // solute - ?
	else{ // solvent - ?
	  assert(m_grid.mask.size() > unsigned(mask_z));
	  for(unsigned int m=0; m<m_grid.mask[mask_z].size(); m+=2){
	    
	    assert(cell_start.size() > unsigned(base + m_grid.mask[mask_z][m]));
	  
	    const int first = cell_start[base + m_grid.mask[mask_z][m]];
	    const int last  = cell_start[base + m_grid.mask[mask_z][m+1]];
	    
	    for(int j=first; j<last; ++j){

	      const double d2 = 
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) *
		(m_grid.p_cell[i_level][i].x - p_plane[j].x) +
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) *
		(m_grid.p_cell[i_level][i].y - p_plane[j].y) +
		(m_grid.p_cell[i_level][i].z - p_plane[j].z) *
		(m_grid.p_cell[i_level][i].z - p_plane[j].z);

	      if (d2 > m_cutoff_long_2) continue;

	      DEBUG(8, "inter (solvent - ?): " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      DEBUG(8, "i_level=" << i_level << " d2=" << d2);
	      
	      if (d2 > m_cutoff_short_2){
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){

		    // maybe for fast solvent loop, check cg of j...

		    // careful : shift index needs to change!
		    if (calculate_pair<t_interaction_spec, t_perturbation_details>
			(topo, conf, storage, innerloop, perturbed_innerloop, a2, a1,
			 m_reverse_shift_vector[p_plane[j].shift_index],
			 sim.param().perturbation.scaled_only))
		      ;
		    else
		      innerloop.lj_crf_innerloop_shift(topo, conf, a1, a2, storage,
						       m_shift_vector[p_plane[j].shift_index]);
		  }
		}
	      }
	      else{
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){

		    if (insert_pair<t_perturbation_details>
			(topo, pairlist, perturbed_pairlist,
			 a2, a1, sim.param().perturbation.scaled_only))
		      ;
		    else
		      pairlist[a2].push_back(a1);
		  }
		}
	      }
	      
	    }
	    
	  }

	} // solvent - ?

      } // particle in the current mask plane

    } // mask levels

  } // the extended planes

  m_timing += util::now() - update_start;
}


template<typename t_perturbation_details>
inline bool interaction::Grid_Pairlist_Algorithm::insert_pair
(
 topology::Topology & topo,
 interaction::Pairlist & pairlist,
 interaction::Pairlist & perturbed_pairlist,
 int a1, int a2,
 bool scaled_only
 )
{
  if (topo.is_perturbed(a1)){
    if (t_perturbation_details::do_scaling && scaled_only){
      // ok, only perturbation if it is a scaled pair...
      std::pair<int, int> 
	energy_group_pair(topo.atom_energy_group(a1),
			  topo.atom_energy_group(a2));
      
      if (topo.energy_group_scaling().count(energy_group_pair))
	perturbed_pairlist[a1].push_back(a2);
      else
	pairlist[a1].push_back(a2);
    } // scaling
    else{
      perturbed_pairlist[a1].push_back(a2);
    }
    return true;
  }

  return false;
}


template<typename t_interaction_spec, typename t_perturbation_details>
inline bool interaction::Grid_Pairlist_Algorithm::calculate_pair
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 interaction::Storage & storage,
 Nonbonded_Innerloop<t_interaction_spec> & innerloop,
 Perturbed_Nonbonded_Innerloop
 <t_interaction_spec, t_perturbation_details> & perturbed_innerloop,
 int a1, int a2,
 math::Vec const & shift,
 bool scaled_only
 )
{

  if (topo.is_perturbed(a1)){
    if (t_perturbation_details::do_scaling && scaled_only){
      // ok, only perturbation if it is a scaled pair...
      std::pair<int, int> 
	energy_group_pair(topo.atom_energy_group(a1),
			  topo.atom_energy_group(a2));
	      
      if (topo.energy_group_scaling().count(energy_group_pair))
	perturbed_innerloop.
	  perturbed_lj_crf_innerloop_shift(topo, conf, a1, a2, storage, shift);
      else
	innerloop.
	  lj_crf_innerloop_shift(topo, conf, a1, a2, storage, shift);
    } // scaling
    else{
      perturbed_innerloop.
	perturbed_lj_crf_innerloop_shift(topo, conf, a1, a2, storage, shift);
    }

    return true;
  }
  return false;
}


bool interaction::Grid_Pairlist_Algorithm
::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  assert(i<j);
  
  std::set<int>::reverse_iterator
    e = topo.all_exclusion(i).rbegin(),
    e_to = topo.all_exclusion(i).rend();

  for( ; e != e_to; ++e){
    if (j > unsigned(*e)) break;
    if (j == unsigned(*e)){
      DEBUG(11, "\texcluded");
      return true;
    }
      
  }
  DEBUG(12, "\tnot excluded");
  return false;
}
