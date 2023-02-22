/**
 * @file extended_grid_pairlist_algorithm.cc
 * extended grid pairlist algorithm
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"

#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/extended_grid_pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"

#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::Extended_Grid_Pairlist_Algorithm::Extended_Grid_Pairlist_Algorithm()
  : interaction::Failing_Pairlist_Algorithm()
{
}

interaction::Extended_Grid_Pairlist_Algorithm::~Extended_Grid_Pairlist_Algorithm()
{
  if (fallback_algorithm != NULL)
    delete fallback_algorithm;
}

int interaction::Extended_Grid_Pairlist_Algorithm::init
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 std::ostream & os,
 bool quiet)
{
  if (!sim.param().pairlist.grid || !sim.param().pairlist.grid_size){
    io::messages.add("wrong input parameters",
                     "Grid_Pairlist_Algorithm",
		     io::message::error);
    return 1;
  }
  
  if (conf.boundary_type != math::rectangular){
    io::messages.add("only implemented for rectangular boundary conditions",
                     "Grid_Pairlist_Algorithm",
		     io::message::error);
    return 1;
  }
  
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);

  grid_properties(topo, conf, sim);

  if (!quiet){
    
    os.precision(4);
    os.setf(std::ios::fixed, std::ios::floatfield);
    
    os << "\textended grid pairlist algorithm\n"
       << "\t\tcells             "
       << std::setw(10) << m_grid.Na
       << std::setw(10) << m_grid.Nb
       << std::setw(10) << m_grid.Nc << "\n"
       << "\t\textenced cells    "
       << std::setw(10) << m_grid.Na_ex
       << std::setw(10) << m_grid.Nb_ex
       << std::setw(10) << m_grid.Nc_ex << "\n"
       << "\t\tcell size         "
       << std::setw(10) << m_grid.a
       << std::setw(10) << m_grid.b
       << std::setw(10) << m_grid.c << "\n";

    const int Ncell = m_grid.Na * m_grid.Nb * m_grid.Nc;
    const int N = m_grid.Na * m_grid.Nb;

    double P = 0.0;
    if(sim.param().pairlist.atomic_cutoff)
      P = topo.num_atoms();
    else
      P = topo.num_chargegroups();

    const double V = math::volume(conf.current().box, conf.boundary_type);
    // const double Vcell = m_grid.a * m_grid.b * m_grid.c;

    const double Vcut = 4.0 / 3.0 * math::Pi * m_cutoff_long_2 * m_cutoff_long;
    const double Pcut = P * Vcut / V;
    
    const double Pcell = P / Ncell;
    const double Player = P / m_grid.Nc;
    
    os << "\t\tparticles / cell    " << std::setw(10) << Pcell << "\n"
       << "\t\tparticles / layer   " << std::setw(10) << Player << "\n";

    // mask size:
    int Nmask = 0;
    for(int z=0; z<=m_grid.mask_z; ++z){
      for(unsigned int y=0; y<m_grid.mask[z].size(); y+=2){
	assert(int(m_grid.mask.size()) > z);
	assert(m_grid.mask[z].size() > y+1);
	Nmask += m_grid.mask[z][y+1] - m_grid.mask[z][y] - 1;
      }
    }
    
    os << "\t\tcells in mask       " << std::setw(10) << Nmask << "\n"
       << "\t\tpairs               " << std::setw(10) << 0.5 * P * P << "\n"
       << "\t\tpairs (grid)        " << std::setw(10) << P * Nmask * Pcell << "\n"
       << "\t\tpairs (cutoff)      " << std::setw(10) << 0.5 * P * Pcut << "\n";

    // just for fun, already try this
    if (prepare_grid(topo, conf, sim)){
      os << "\t\terror during grid preparation!\n";
      io::messages.add("error during grid preparation. Use standard pairlist instead.",
                       "Grid_Pairlist_Algorithm",
		       io::message::error);
      return 1;
    }
    // collapse_grid();

    int occupied = 0;
    for(int z=0; z<m_grid.Nc; ++z){
      for(int i=0; i<N; ++i){
	assert(int(m_grid.count.size()) > z &&
	       int(m_grid.count[z].size()) > i);
	if (m_grid.count[z][i])
	  ++occupied;
      }
    }

    if (!quiet){
      os << "\t\toccupied            " << std::setw(10) << occupied << "\n"
	 << "\t\tparticle / occ cell " << std::setw(10) << double(P) / occupied << "\n";
      
      // print_mask();

      // os << "END\n";
    }
  }

  fallback_algorithm = new Standard_Pairlist_Algorithm();
  fallback_algorithm->init(topo, conf, sim, os, true);

  return 0;
}

/**
 * calculate grid properties,
 * put chargegroups on grid,
 * including the center of geometries
 */
int interaction::Extended_Grid_Pairlist_Algorithm::prepare
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  DEBUG(7, "grid pairlist algorithm : prepare");
  timer().start_subtimer("pairlist prepare");
  
  set_cutoff(sim.param().pairlist.cutoff_short, 
	     sim.param().pairlist.cutoff_long);
		 
  // first put the chargegroups into the box
  // _prepare_cog(conf, topo);

  grid_properties(topo, conf, sim);

  failed = false;

  if (prepare_grid(topo, conf, sim)){
    timer().stop_subtimer("pairlist prepare");
    std::ostringstream msg;
    msg << "At step " << sim.steps() << ": Could not prepare grid. "
            "Falling back to standard algoritm for this step.";
    io::messages.add(msg.str(), "Extended_Grid_Pairlist", io::message::notice);
    failed = true;
    return fallback_algorithm->prepare(topo, conf, sim);
  }

  collapse_grid();

  // print_grid();
  timer().stop_subtimer("pairlist prepare");

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// pairlist update
////////////////////////////////////////////////////////////////////////////////
void interaction::Extended_Grid_Pairlist_Algorithm::update
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 unsigned int begin,
 unsigned int end,
 unsigned int stride
 )
{
  if (failed) {
    fallback_algorithm->update(topo, conf, sim, pairlist, begin, end, stride);
    return;
  }

  if (sim.param().pairlist.atomic_cutoff){
    _update_atomic(topo, conf, sim, pairlist, begin, end, stride);
  }
  else{
    _update(topo, conf, sim, pairlist, begin, end, stride);
  }    
}

void interaction::Extended_Grid_Pairlist_Algorithm::_update
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 unsigned int begin,
 unsigned int end,
 unsigned int stride
 )
{
  // empty the pairlist
  pairlist.clear();
  timer().start_subtimer("pairlist");
  
  DEBUG(10,"Extended_Grid_Pairlist_Algorithm::_update");
  std::vector<Grid::Particle> p_plane;
  std::vector<int> cell_start;

  const int c_ex = (m_grid.Nc_ex - m_grid.Nc) / 2;
  const int N = m_grid.Na * m_grid.Nb;
  const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_solute_atoms = topo.num_solute_atoms();
	const simulation::qmmm_enum qmmm = sim.param().qmmm.qmmm;

#ifdef HAVE_LIBCUDART  
  const bool no_cuda = sim.param().innerloop.method != simulation::sla_cuda;
#endif

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
		// Exclude all QM atoms from self-interaction
		if (!qmmm || !topo.is_qm(topo.chargegroup(m_grid.p_cell[z][i].i))) {
	  DEBUG(8, "self " << m_grid.p_cell[z][i].i);
	  for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		a_to = topo.chargegroups()[m_grid.p_cell[z][i].i + 1];
	      a1 < a_to; ++a1){
	    for(int a2 = a1+1; a2 < a_to; ++a2){
	      if (excluded_solute_pair(topo, a1, a2))
		continue;
	      pairlist.solute_short[a1].push_back(a2);
	    }
	  }
	}
		else {
			DEBUG(9, "Skipping CG " << m_grid.p_cell[z][i].i);
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

		if (Pairlist_Algorithm::qm_excluded(
			topo, qmmm, topo.chargegroup(ii), topo.chargegroup(jj))) {
			DEBUG(9, "Skipping CGs " << ii << " " << jj);
			continue;
		}
		for(int a1 = topo.chargegroups()[ii],
		      a_to = topo.chargegroups()[ii + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[jj],
			a2_to = topo.chargegroups()[jj + 1];
		      a2 < a2_to; ++a2){
		    if (excluded_solute_pair(topo, a1, a2))
		      continue;
		    pairlist.solute_short[a1].push_back(a2);
		  }
		}
	      }
	      else{ // solute - solvent
				if (Pairlist_Algorithm::qm_excluded(
						topo, qmmm
						, topo.chargegroup(m_grid.p_cell[z][i].i)
						, topo.chargegroup(m_grid.p_cell[z][j].i))) {
					DEBUG(9, "Skipping CGs " << m_grid.p_cell[z][i].i << "-" << m_grid.p_cell[z][j].i);
					continue;
				}
		for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
			a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		      a2 < a2_to; ++a2){
		    pairlist.solute_short[a1].push_back(a2);
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
		    a2 < a2_to; ++a2) {
                  if (a2 >= num_solute_atoms) { //solvent-solvent
#ifdef HAVE_LIBCUDART
                    // do only add the atoms to the pairlist if we are not doing
                    // CUDA
                    if (no_cuda)
                      pairlist.solvent_short[a2].push_back(a1);
#else
                    pairlist.solvent_short[a2].push_back(a1);
#endif
                  } else {
										if (Pairlist_Algorithm::qm_excluded(topo, qmmm, a1, a2)) {
											DEBUG(9, "Skipping pair: " << a1 << "-" << a2);
											continue;
										}
                    pairlist.solute_short[a2].push_back(a1);
			}
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

		if (Pairlist_Algorithm::qm_excluded(topo, qmmm
				, topo.chargegroup(m_grid.p_cell[i_level][i].i)
				, topo.chargegroup(p_plane[j].i))) {
			DEBUG(9, "Skipping pairs: (" << topo.chargegroup(m_grid.p_cell[i_level][i].i)
								  << "-" << topo.chargegroup(m_grid.p_cell[i_level][i].i + 1) - 1
								  << ")-(" << topo.chargegroup(p_plane[j].i)
								  << "-" << topo.chargegroup(p_plane[j].i + 1) - 1
								  << ")");
			continue;
		}

		if (d2 > m_cutoff_short_2){ // LONGRANGE

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      pairlist.solute_long[a1].push_back(a2);
		    }
		  }
		}
		else{ // SHORTRANGE
		  
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
		      pairlist.solute_short[a1].push_back(a2);
		    }
		  }
		}
	      }
	      else{ // solute - solvent
		if (Pairlist_Algorithm::qm_excluded(topo, qmmm
				, topo.chargegroup(m_grid.p_cell[i_level][i].i)
				, topo.chargegroup(p_plane[j].i))) {
			DEBUG(9, "Skipping pairs: (" << topo.chargegroup(m_grid.p_cell[i_level][i].i)
								  << "-" << topo.chargegroup(m_grid.p_cell[i_level][i].i + 1) - 1
								  << ")-(" << topo.chargegroup(p_plane[j].i)
								  << "-" << topo.chargegroup(p_plane[j].i + 1) - 1
								  << ")");
			continue;
		}
		if (d2 > m_cutoff_short_2){ // LONGRANGE

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      
		      pairlist.solute_long[a1].push_back(a2);
		    }
		  }
		}
		else{ // SHORTRANGE

		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
		      pairlist.solute_short[a1].push_back(a2);
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
	      
	      if (d2 > m_cutoff_short_2){ // LONGRANGE
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2) {
                    if (a2 >= num_solute_atoms) { //solvent-solvent
#ifdef HAVE_LIBCUDART
                      if (no_cuda)
                        pairlist.solvent_long[a2].push_back(a1);
#else
                      pairlist.solvent_long[a2].push_back(a1);
#endif
                    } else

											if (Pairlist_Algorithm::qm_excluded(topo, qmmm, a1, a2)) {
												DEBUG(9, "Skipping pair: " << a1 << "-" << a2);
												continue;
											}
                        pairlist.solute_long[a2].push_back(a1);
		  }
		}
	      }
	      else{ // SHORTRANGE
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2) {
                    if (a2 >= num_solute_atoms) {
#ifdef HAVE_LIBCUDART
                      if (no_cuda)
                        pairlist.solvent_short[a2].push_back(a1);
#else
                      pairlist.solvent_short[a2].push_back(a1);
#endif
                    } else
										if (Pairlist_Algorithm::qm_excluded(topo, qmmm, a1, a2)) {
												DEBUG(9, "Skipping pair: " << a1 << "-" << a2);
												continue;
											}
                      pairlist.solute_short[a2].push_back(a1);
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
  timer().stop_subtimer("pairlist");
}

void interaction::Extended_Grid_Pairlist_Algorithm::_update_atomic
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 unsigned int begin,
 unsigned int end,
 unsigned int stride
 )
{
  // empty the pairlist
  pairlist.clear();
  timer().start_subtimer("pairlist");
  
  DEBUG(10,"Extended_Grid_Pairlist_Algorithm::_update_atomic");
  std::vector<Grid::Particle> p_plane;
  std::vector<int> cell_start;

  const int c_ex = (m_grid.Nc_ex - m_grid.Nc) / 2;
  const int N = m_grid.Na * m_grid.Nb;
  //const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_solute_atoms = topo.num_solute_atoms();

#ifdef HAVE_LIBCUDART  
  const bool no_cuda = sim.param().innerloop.method != simulation::sla_cuda;
#endif

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

	if (base != m_grid.cell_index[z][i]){
	  base = m_grid.cell_index[z][i];
	  start = i;
	}
	else{ // more than 1 atom in the cell
	  if (m_grid.p_cell[z][i].i < num_solute_atoms){ // solute - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      if (m_grid.p_cell[z][j].i < num_solute_atoms){ // solute - solute

		const int ii = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][i].i : m_grid.p_cell[z][j].i;
		const int jj = (m_grid.p_cell[z][i].i < m_grid.p_cell[z][j].i) ? 
		  m_grid.p_cell[z][j].i : m_grid.p_cell[z][i].i;

		if(excluded_solute_pair(topo,ii,jj))
		  continue;
		pairlist.solute_short[ii].push_back(jj);
	      }
	      else{ // solute - solvent
                pairlist.solute_short[m_grid.p_cell[z][i].i].push_back(m_grid.p_cell[z][j].i);
	      }
	    }
	  }
	  else{ // solvent - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(8, "intra cell " << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
              if (m_grid.p_cell[z][j].i >= num_solute_atoms) { //solvent-solvent
		if(!excluded_solvent_pair(topo, m_grid.p_cell[z][j].i, m_grid.p_cell[z][i].i)){
		  
#ifdef HAVE_LIBCUDART
		  // do only add the atoms to the pairlist if we are not doing
		  // CUDA
		  if (no_cuda)
		    pairlist.solvent_short[m_grid.p_cell[z][j].i].push_back(m_grid.p_cell[z][i].i);
#else
		  pairlist.solvent_short[m_grid.p_cell[z][j].i].push_back(m_grid.p_cell[z][i].i);
#endif
		}
		
	      } else
		pairlist.solute_short[m_grid.p_cell[z][j].i].push_back(m_grid.p_cell[z][i].i);
	      
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
	
	if (m_grid.p_cell[i_level][i].i < num_solute_atoms){ // solute - ?
	  
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

	      if (p_plane[j].i < num_solute_atoms){ 
		// solute - solute

		if (d2 > m_cutoff_short_2){ // LONGRANGE

		  pairlist.solute_long[m_grid.p_cell[i_level][i].i].push_back(p_plane[j].i);
		}
		else{ // SHORTRANGE
		  
		  const int ii = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    m_grid.p_cell[i_level][i].i : p_plane[j].i;
		  const int jj = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    p_plane[j].i : m_grid.p_cell[i_level][i].i;
		  
		  DEBUG(8, "rewritten: " << ii
			<< " - " << jj);
		  
		  if (excluded_solute_pair(topo, ii, jj))
		    continue;
		  pairlist.solute_short[ii].push_back(jj);
		}
	      }
	      else{ // solute - solvent
		if (d2 > m_cutoff_short_2){ // LONGRANGE

		  pairlist.solute_long[m_grid.p_cell[i_level][i].i].push_back(p_plane[j].i);
		}
		else{ // SHORTRANGE

		  pairlist.solute_short[m_grid.p_cell[i_level][i].i].push_back(p_plane[j].i);
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
	      
	      if (d2 > m_cutoff_short_2){ // LONGRANGE
		if (p_plane[j].i >= num_solute_atoms) { //solvent-solvent
		  if(!excluded_solvent_pair(topo, p_plane[j].i, m_grid.p_cell[i_level][i].i)){
#ifdef HAVE_LIBCUDART
		  if (no_cuda)
		    pairlist.solvent_long[p_plane[j].i].push_back(m_grid.p_cell[i_level][i].i);
#else
		  pairlist.solvent_long[p_plane[j].i].push_back(m_grid.p_cell[i_level][i].i);
#endif
		  }
		} else
		  pairlist.solute_long[p_plane[j].i].push_back(m_grid.p_cell[i_level][i].i);
	      }
	      else{ // SHORTRANGE
		if (p_plane[j].i >= num_solute_atoms) {
		  if(!excluded_solvent_pair(topo, p_plane[j].i, m_grid.p_cell[i_level][i].i)){
#ifdef HAVE_LIBCUDART
		    if (no_cuda)
		      pairlist.solvent_short[p_plane[j].i].push_back(m_grid.p_cell[i_level][i].i);
#else
		    pairlist.solvent_short[p_plane[j].i].push_back(m_grid.p_cell[i_level][i].i);
#endif
		  }
		} else
		  pairlist.solute_short[p_plane[j].i].push_back(m_grid.p_cell[i_level][i].i);
	      }
	      
	    }
	    
	  }

	} // solvent - ?

	// std::cout << "\n";
	
      } // particle in the current mask plane

    } // mask levels

  } // the extended planes
  timer().stop_subtimer("pairlist");
}

////////////////////////////////////////////////////////////////////////////////
// perturbed pairlist update
////////////////////////////////////////////////////////////////////////////////

void interaction::Extended_Grid_Pairlist_Algorithm::update_perturbed
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 interaction::PairlistContainer & perturbed_pairlist,
 unsigned int begin,
 unsigned int end, 
 unsigned int stride
 ) {
  if (failed) {
    DEBUG(10,"Extended_Grid_Pairlist_Algorithm::update_perturbed. FALLBACK");
    fallback_algorithm->update_perturbed(topo, conf, sim, pairlist, perturbed_pairlist, begin, end, stride);

    return;
  }
  _update_perturbed(topo, conf, sim, pairlist, perturbed_pairlist, begin, end, stride);
}

void interaction::Extended_Grid_Pairlist_Algorithm::_update_perturbed
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 interaction::PairlistContainer & pairlist,
 interaction::PairlistContainer & perturbed_pairlist,
 unsigned int begin,
 unsigned int end, 
 unsigned int stride
 )
{
  // check whether we do scaling && scaling only
  bool scaled_only = (sim.param().perturbation.scaling && sim.param().perturbation.scaled_only);
 
// empty the pairlist
  assert(pairlist.size() == topo.num_atoms());
  assert(perturbed_pairlist.size() == topo.num_atoms());

  pairlist.clear();
  perturbed_pairlist.clear();
  timer().start_subtimer("perturbed pairlist");

  std::vector<Grid::Particle> p_plane;
  std::vector<int> cell_start;

  const int c_ex = (m_grid.Nc_ex - m_grid.Nc) / 2;
  const int N = m_grid.Na * m_grid.Nb;
  const int num_solute_cg = topo.num_solute_chargegroups();
  const int num_solute_atoms = topo.num_solute_atoms();
  
  DEBUG(10, "Extended_Grid_Pairlist_Algorithm::_update_perturbed");
  
#ifdef HAVE_LIBCUDART  
  const bool no_cuda = sim.param().innerloop.method != simulation::sla_cuda;
#endif
  
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
          DEBUG(10,"INTRACELL");
	if (m_grid.p_cell[z][i].i < num_solute_cg){ // self interaction
	  DEBUG(11, "self interaction" << m_grid.p_cell[z][i].i);
          DEBUG(10,"\tsolute-solute (intracell)");
	  for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		a_to = topo.chargegroups()[m_grid.p_cell[z][i].i + 1];
	      a1 < a_to; ++a1){
	    for(int a2 = a1+1; a2 < a_to; ++a2){
              DEBUG(10,"\t\ta1 " << a1 << " a2 " << a2);
	      if (excluded_solute_pair(topo, a1, a2))
		continue;
              
	      if (insert_pair
		  (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
		   a1, a2, scaled_only)){
                DEBUG(10,"\t\t\tadded a2 to a1 solute shortrange pairlist");  
              }
	      else if (insert_pair
		       (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
			a2, a1, scaled_only)){
                DEBUG(10,"\t\t\tadded a1 to a2 solute shortrange pairlist");  
              }
	      else{
		pairlist.solute_short[a1].push_back(a2);
                DEBUG(10,"\t\t\tadd a2 to a1 solute shortrange pairlist");
              }
	    }
	  }
	}
	
	if (base != m_grid.cell_index[z][i]){
          DEBUG(10,"\tbase!=cell index (intracell)");
	  base = m_grid.cell_index[z][i];
	  start = i;
	}
	else{ // more than 1 cg in the cell
          DEBUG(10,"more than 1 cg in the cell");
	  if (m_grid.p_cell[z][i].i < num_solute_cg){ // solute - ?
            DEBUG(10,"\tsolute - ? (intracell)");
	    for(int j=start; j<i; ++j){
	      DEBUG(11, "intra cell solute" << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      
	      if (m_grid.p_cell[z][j].i < num_solute_cg){ // solute - solute
                DEBUG(10,"\t\tsolute-solute (intracell)");
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
		    DEBUG(10,"\t\t\ta1 " << a1 << " a2 " << a2);
		    if (excluded_solute_pair(topo, a1, a2))
		      continue;
                    
		    if (insert_pair
			(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
			 a1, a2, scaled_only)){
                       DEBUG(10,"\t\t\tadded a2 to a1 solute shortrange pairlist");   
                    }
		    else if (insert_pair
			     (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
			      a2, a1, scaled_only)){
                        DEBUG(10,"\t\t\tadded a1 to a2 solute shortrange pairlist");
                    }
		    else{
		      pairlist.solute_short[a1].push_back(a2);
                      DEBUG(10,"\t\t\tadd a2 to a1 solute shortrange pairlist");
                    }
		  }
		}
	      }
	      else{ // solute - solvent
                DEBUG(10,"\t\tsolute-solvent (intracell)");
		for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
			a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		      a2 < a2_to; ++a2){
                    DEBUG(10,"\t\t\ta1 " << a1 << " a2 " << a2);
		    if (insert_pair
			(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
			 a1, a2, scaled_only)){
                      DEBUG(10,"\t\t\tadded a2 to a1 solute shortrange pairlist");   
                    }
		    else{
		      pairlist.solute_short[a1].push_back(a2);
                      DEBUG(10,"\t\t\tadd a2 to a1 solute shortrange pairlist");
                    }
		  }
		}

	      }
	    }
	  }
	  else{ // solvent - ?
	    for(int j=start; j<i; ++j){
	      DEBUG(11, "intra cell solvent" << m_grid.p_cell[z][i].i << " - " 
		    << m_grid.p_cell[z][j].i);
	      DEBUG(10,"\tsolvent-? (intracell)");
	      for(int a1 = topo.chargegroups()[m_grid.p_cell[z][i].i],
		    a_to = topo.chargegroups()[m_grid.p_cell[z][i].i+1];
		  a1 < a_to; ++a1){
		for(int a2 = topo.chargegroups()[m_grid.p_cell[z][j].i],
		      a2_to = topo.chargegroups()[m_grid.p_cell[z][j].i + 1];
		    a2 < a2_to; ++a2){
                  DEBUG(10,"\t\ta1 " << a1 << " a2 " << a2);
		  if (insert_pair //only a2 can be solute
		      (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
		       a2, a1, scaled_only)){
		    DEBUG(10,"\t\tadded a1 to a2 solute shortrange pairlist");
                  }
                  else if (a2 >= num_solute_atoms){ //solvent-solvent
#ifdef HAVE_LIBCUDART
                    // do only add the atoms to the pairlist if we are not doing
                    // CUDA
                    if (no_cuda)
                      pairlist.solvent_short[a2].push_back(a1);
#else
                    pairlist.solvent_short[a2].push_back(a1);
#endif
                    DEBUG(10,"\t\tadd a1 to a2 solvent shortrange pairlist")
                  }
                  else{ //solvent-solute
                    pairlist.solute_short[a2].push_back(a1);
                    DEBUG(10,"\t\tadd a1 to a2 solute shortrange pairlist");
                  }
                  
		}
	      }
	    }
	  } //solvent - ?
	}// more than 1 cg in the cell
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
        
	DEBUG(10,"INTERCELL");
	if (m_grid.p_cell[i_level][i].i < num_solute_cg){ // solute - ?
          DEBUG(10,"\tsolute-? (intercell)")
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
	      DEBUG(11, "inter cell: " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      assert(m_grid.p_cell[i_level][i].i != p_plane[j].i);

	      if (p_plane[j].i < num_solute_cg){ 
		// solute - solute
	      DEBUG(10,"\t\tsolute-solute (intercell)");
                if (d2 > m_cutoff_short_2){ // LONGRANGE
                  DEBUG(10,"\t\t\tlongrange (solute-solute intercell)");
                  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
                          a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
                  a1 < a_to; ++a1){
                    for(int a2 = topo.chargegroups()[p_plane[j].i],
                            a2_to = topo.chargegroups()[p_plane[j].i + 1];
                    a2 < a2_to; ++a2){
                      DEBUG(10,"\t\t\ta1 " << a1 << " a2 " << a2);
                      if (insert_pair
                              (topo, pairlist.solute_long, perturbed_pairlist.solute_long,
                              a1, a2, scaled_only)){
                          DEBUG(10,"\t\t\t\tperturbed: added a2 to a1 solute longrange pairlist");
                      }
                      // careful: shift index needs to change here !!!
                      else if (insert_pair
                              (topo, pairlist.solute_long, perturbed_pairlist.solute_long,
                              a2, a1, scaled_only)){
                          DEBUG(10,"\t\t\t\tperturbed: added a1 to a2 solute longrange pairlist");
                      }
                      else {
                        pairlist.solute_long[a1].push_back(a2);
                        DEBUG(10,"\t\t\t\tadd a2 to a1 solute longrange pairlist");
                      }
                    }
                  }
                }
		else{ // SHORTRANGE
                  DEBUG(10,"\t\t\tshortrange (solute-solute intercell)");
		  // exclusions => order
		  const int ii = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    m_grid.p_cell[i_level][i].i : p_plane[j].i;
		  const int jj = (m_grid.p_cell[i_level][i].i < p_plane[j].i) ? 
		    p_plane[j].i : m_grid.p_cell[i_level][i].i;
		  
		  DEBUG(8, "\t\t\trewritten: " << ii
			<< " - " << jj);
		  
		  for(int a1 = topo.chargegroups()[ii],
			a_to = topo.chargegroups()[ii + 1];
                        a1 < a_to; ++a1){
                    for(int a2 = topo.chargegroups()[jj],
                            a2_to = topo.chargegroups()[jj + 1];
                            a2 < a2_to; ++a2){
                      DEBUG(10,"\t\t\ta1 " << a1 << " a2 " << a2);                                   
                      if (excluded_solute_pair(topo, a1, a2))
                        continue;
                                            
                      if (insert_pair
                              (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
                              a1, a2, scaled_only)){
                         DEBUG(10,"\t\t\t\tperturbed: added a2 to a1 solute shortrange pairlist");
                      }
                      else if (insert_pair
                              (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
                              a2, a1, scaled_only)){
                         DEBUG(10,"\t\t\t\tperturbed: added a1 to a2 solute shortrange pairlist");
                      }
                      else{
                        pairlist.solute_short[a1].push_back(a2);
                        DEBUG(10,"\t\t\t\tadd a2 to a1 solute shortrange pairlist");
                      }
                    }
		  }
		}
	      }
              else{ // solute - solvent
                DEBUG(10,"\t\tsolute-solvent (intercell)");
                if (d2 > m_cutoff_short_2){ // LONGRANGE
                  DEBUG(10,"\t\t\tlongrange (solute-solvent intercell)");
                  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
                          a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
                          a1 < a_to; ++a1){
                    for(int a2 = topo.chargegroups()[p_plane[j].i],
                            a2_to = topo.chargegroups()[p_plane[j].i + 1];
                            a2 < a2_to; ++a2){
                      DEBUG(10,"\t\t\ta1 " << a1 << " a2 " << a2);
                      if (insert_pair
                              (topo, pairlist.solute_long, perturbed_pairlist.solute_long,
                              a1, a2, scaled_only)){
                          DEBUG(10,"\t\t\t\tperturbed: added a2 to a1 solute longrange pairlist");
                      }
                      else {
                        pairlist.solute_long[a1].push_back(a2);
                        DEBUG(10,"\t\t\t\tadd a2 to a1 solute longrange pairlist");
                      }
                    }
                  }
                }
                else{ // SHORTRANGE
                  DEBUG(10,"\t\t\tshortrange (solute-solvent intercell)");
		  for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
			a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		    for(int a2 = topo.chargegroups()[p_plane[j].i],
			  a2_to = topo.chargegroups()[p_plane[j].i + 1];
			a2 < a2_to; ++a2){
                      DEBUG(10,"\t\t\ta1 " << a1 << " a2 " << a2);
		      if (insert_pair
			  (topo, pairlist.solute_short, perturbed_pairlist.solute_short,
			   a1, a2, scaled_only)){
                          DEBUG(10,"\t\t\t\tperturbed: added a2 to a1 solute shortrange pairlist");
                      }
		      else{
			pairlist.solute_short[a1].push_back(a2);
                        DEBUG(10,"\t\t\t\tadd a2 to a1 solute shortrange pairlist");
                      }
		    }
		  }
		}
		
	      }
	      
	    } // j in mask row
	  }  // mask
	} // solute - ?
	else{ // i=solvent - j=?
          DEBUG(10,"\tsolvent - ? (intercell)");
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

	      DEBUG(11, "inter (solvent - ?): " << m_grid.p_cell[i_level][i].i << " - " << p_plane[j].i);
	      DEBUG(11, "i_level=" << i_level << " d2=" << d2);
	      
	      if (d2 > m_cutoff_short_2){ // LONGRANGE
                DEBUG(10,"\t\tlongrange (solvent - ? intercell)");
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		      a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
		       a2_to = topo.chargegroups()[p_plane[j].i + 1];
		       a2 < a2_to; ++a2){
		    // maybe for fast solvent loop, check cg of j...
                    DEBUG(10,"\t\ta1 " << a1 << " a2 " << a2);
		    // careful : shift index needs to change!
		    if (insert_pair
			(topo, pairlist.solute_long, perturbed_pairlist.solute_long,
			 a2, a1, scaled_only)){
		       DEBUG(10,"\t\t\tperturbed: a1 added to a2 solute longrange pairlist");
                    }
                    else if (a2 >= num_solute_atoms){ //solvent -solvent //BUG: a1 >= num_solute_atoms && a1 >= num_solute_atoms --martina
#ifdef HAVE_LIBCUDART
                    // do only add the atoms to the pairlist if we are not doing
                    // CUDA
                    if (no_cuda)
                      pairlist.solvent_long[a1].push_back(a2);
#else
                      pairlist.solvent_long[a1].push_back(a2);
#endif
                      DEBUG(10,"\t\t\tsolvent-solvent: add a2 to a1 solvent longrange pairlist");
                    }
                    else{ //if (a1 >= num_solute_atoms && a2 < num_solute_atoms){ //solvent-solute
                      pairlist.solute_long[a2].push_back(a1);
                      DEBUG(10,"\t\t\tsolvent-solute: add a1 to a2 solute longrange pairlist");
                    }
                    /*else{
                       io::messages.add("error during perturbed grid update. Use standard pairlist instead.",
                       "Grid_Pairlist_Algorithm",
		       io::message::error);
                       return;
                    }*/
		  }//for
		}
	      }
	      else{ // SHORTRANGE
                DEBUG(10,"\t\tshortrange (solvent - ? intercell)");
		for(int a1 = topo.chargegroups()[m_grid.p_cell[i_level][i].i],
		      a_to = topo.chargegroups()[m_grid.p_cell[i_level][i].i + 1];
		    a1 < a_to; ++a1){
		  for(int a2 = topo.chargegroups()[p_plane[j].i],
			a2_to = topo.chargegroups()[p_plane[j].i + 1];
		      a2 < a2_to; ++a2){
                    DEBUG(10,"\t\ta1 " << a1 << " a2 " << a2);
		    if (insert_pair
			(topo, pairlist.solute_short, perturbed_pairlist.solute_short,
			 a2, a1, scaled_only)){ //only a2 could be solute (=perturbed). therefore no need to check for a1
		       DEBUG(10,"\t\t\tperturbed: a1 added to a2 solute shortrange pairlist");
                    }
                    else if (a2 >= num_solute_atoms){ //solvent-solvent
#ifdef HAVE_LIBCUDART
                    // do only add the atoms to the pairlist if we are not doing
                    // CUDA
                      if (no_cuda)
                        pairlist.solvent_short[a2].push_back(a1);
#else
                      pairlist.solvent_short[a2].push_back(a1);
#endif
                      DEBUG(10,"\t\t\tsolvent-solvent: add a1 to a2 solvent shortrange pairlist");
                    }
                    else{ //if (a1 >= num_solute_atoms && a2 < num_solute_atoms){ //solvent-solute
                      pairlist.solute_short[a2].push_back(a1);
                      DEBUG(10,"\t\t\tsolvent-solute: add a1 to a2 solute shortrange pairlist");
                    }
		  }
		}
	      }
	      
	    }
	    
	  }

	} // solvent - ?

      } // particle in the current mask plane

    } // mask levels

  } // the extended planes
  timer().stop_subtimer("perturbed pairlist");
}

inline bool interaction::Extended_Grid_Pairlist_Algorithm::insert_pair
(
 topology::Topology & topo,
 interaction::Pairlist & pairlist,
 interaction::Pairlist & perturbed_pairlist,
 int a1, int a2,
 bool scaled_only
 )
{
  DEBUG(10,"\t\t\t\tinsert perturbed pair?");
  if (topo.is_perturbed(a1) || topo.is_eds_perturbed(a1)){
    if (scaled_only){
        DEBUG(10,"\t\t\t\t\tscaled")
      // ok, only perturbation if it is a scaled pair...
      std::pair<int, int> 
	energy_group_pair(topo.atom_energy_group(a1),
			  topo.atom_energy_group(a2));
      
      if (topo.energy_group_scaling().count(energy_group_pair)){
	perturbed_pairlist[a1].push_back(a2);
        DEBUG(10,"\t\t\t\t\t\tadd a2 to perturbed a1 pairlist");
      }
      else{
	pairlist[a1].push_back(a2);
        DEBUG(10,"\t\t\t\t\t\tadd a2 to a1 pairlist");
      }
    } // scaling
    else{
      perturbed_pairlist[a1].push_back(a2);
      DEBUG(10,"\t\t\t\t\tadd a2 to perturbed a1 pairlist");
    }
    return true;
  }
  DEBUG(10,"\t\t\t\tno insertion");
  return false;
}

bool interaction::Extended_Grid_Pairlist_Algorithm
::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  DEBUG(10,"\t\t\t\texcluded solute pair?");
  return topo.all_exclusion(i).is_excluded(j);
}

bool interaction::Extended_Grid_Pairlist_Algorithm
::excluded_solvent_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  DEBUG(10,"\t\t\t\texcluded solvent pair?");
  
  // determine what molecule we are in
  int ii = i - topo.num_solute_atoms();
  int jj = j - topo.num_solute_atoms();
  unsigned int solv_i=0;
  unsigned int solv_j=0;
  //DEBUG(10, "ii: " << ii << " jj: " << jj);
  //DEBUG(10, "num_solvents " << topo.num_solvents());
  
  while(solv_i < topo.num_solvents() && ii > topo.num_solvent_atoms(solv_i)){
    ii -= topo.num_solvent_atoms(solv_i);
    ++solv_i;
  }
  while(solv_j < topo.num_solvents() && jj > topo.num_solvent_atoms(solv_j)){
    jj -= topo.num_solvent_atoms(solv_j);
    ++solv_j;
  }
  
  if(solv_i!=solv_j) return false;
  //DEBUG(10, "solv_i: " << solv_i << " solv_j: " << solv_j);
  const int num_solv_at = topo.num_solvent_atoms(solv_i) / topo.num_solvent_molecules(solv_i);

  int m_i = ii / num_solv_at;
  int m_j = jj / num_solv_at;
  
  if(m_i == m_j) return true;
  else return false;
}

