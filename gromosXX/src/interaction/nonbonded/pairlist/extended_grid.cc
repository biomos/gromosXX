/**
 * @file extended_grid.cc
 * grid routines : prepare a grid for the Grid_Pairlist_Algorithm
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction_types.h"
#include "../../../math/periodicity.h"
#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/extended_grid_pairlist_algorithm.h"

#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

////////////////////////////////////////////////////////////////////////////////
// put chargegroups into box and on the grid
// the grid will be sparse...
////////////////////////////////////////////////////////////////////////////////

/**
 * put the chargegroups into the box and on the grid
 */
int interaction::Extended_Grid_Pairlist_Algorithm::prepare_grid
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  const int Nab = m_grid.Na * m_grid.Nb;
  
  // this is enough for small cells (when only about one particle fits in anyway)
  const int space = m_grid.Pcell * 3 + m_grid.shift_space;

  /*
  std::cout << "grid: Nab=" << Nab << " space=" << space
	    << " shift_space=" << m_grid.shift_space
	    << " Pcell=" << m_grid.Pcell << std::endl;
  */

  m_grid.p_cell.resize(m_grid.Nc);
  m_grid.cell_start.resize(m_grid.Nc);
  m_grid.count.resize(m_grid.Nc);
  m_grid.cell_index.resize(m_grid.Nc);
  
  for(int z=0; z<m_grid.Nc; ++z){
    m_grid.p_cell[z].resize(Nab * space);
    m_grid.count[z].assign(Nab, 0);
    m_grid.cell_start[z].resize(Nab+1);
    
    int b = 0;
    for(int i=0; i<Nab; ++i, b+=space)
      m_grid.cell_start[z][i] = b;
  }
  
  math::Periodicity<math::rectangular> periodicity(conf.current().box);

  math::VArray &pos = conf.current().pos;
  math::Vec v, v_box, trans;

  if(!sim.param().pairlist.atomic_cutoff){
      
    topology::Chargegroup_Iterator cg_it = topo.chargegroup_begin(),
      cg_to = topo.chargegroup_end();

    // solute chargegroups...
    unsigned int i = 0;
    for( ; i < topo.num_solute_chargegroups(); ++cg_it, ++i){
      // cog
      cg_it.cog(pos, v);

      // gather on first atom...
      v_box = v;
      periodicity.put_into_box(v_box);
      trans = v_box - v;

      // now grid the cg
       
      const int x = int((v_box(0) + 0.5 * abs(conf.current().box(0))) / m_grid.a);
      const int y = int((v_box(1) + 0.5 * abs(conf.current().box(1))) / m_grid.b);
      const int z = int((v_box(2) + 0.5 * abs(conf.current().box(2))) / m_grid.c);
  
      const int c = y * m_grid.Na + x;

      if(!(m_grid.cell_start.size() > unsigned(z) &&
           m_grid.cell_start[z].size() > unsigned(c)))
        return 1;
    
      if(!(m_grid.p_cell.size() > unsigned(z) &&
	   m_grid.p_cell[z].size() > unsigned(m_grid.cell_start[z][c] + m_grid.count[z][c])))
        return 1;

      if(!(m_grid.count.size() > unsigned(z) &&
	   m_grid.count[z].size() > unsigned(c)))
        return 1;
    
      m_grid.p_cell[z][m_grid.cell_start[z][c] + m_grid.count[z][c]] =
        Grid::Particle(i, v_box);

      ++m_grid.count[z][c];
    
      // atoms in a chargegroup
      topology::Atom_Iterator at_it = cg_it.begin(),
        at_to = cg_it.end();
      for( ; at_it != at_to; ++at_it){
        pos(*at_it) += trans;

      } // loop over atoms
    } // loop over solute cg's
 
    // solvent chargegroups
    for( ; cg_it != cg_to; ++cg_it, ++i){

      // cog is first atom
      v = pos(**cg_it);
      v_box = v;
      periodicity.put_into_box(v_box);
      trans = v_box - v;

      // now grid the cg
      
      const int x = int((v_box(0) + 0.5 * abs(conf.current().box(0))) / m_grid.a);
      const int y = int((v_box(1) + 0.5 * abs(conf.current().box(1))) / m_grid.b);
      const int z = int((v_box(2) + 0.5 * abs(conf.current().box(2))) / m_grid.c);
    
      const int c = y * m_grid.Na + x;

      if(!(m_grid.cell_start.size() > unsigned(z) &&
	   m_grid.cell_start[z].size() > unsigned(c)))
        return 1;
      if(!(m_grid.p_cell.size() > unsigned(z))) 
        return 1;

      if (m_grid.p_cell[z].size() <= unsigned(m_grid.cell_start[z][c] + m_grid.count[z][c])){

      /*
      std::cout << "ERROR: " << " z=" << z
		<< " cell_start=" << m_grid.cell_start[z][c]
		<< " count=" << m_grid.count[z][c]
		<< std::endl;
      */
      // don't add particle to grid. not enough space!
      }
      else{
        if(!(m_grid.p_cell[z].size() > unsigned(m_grid.cell_start[z][c] + m_grid.count[z][c])))
          return 1;
        if(!(m_grid.count.size() > unsigned(z) &&
	     m_grid.count[z].size() > unsigned(c)))
          return 1;

        m_grid.p_cell[z][m_grid.cell_start[z][c] + m_grid.count[z][c]] = 
	  Grid::Particle(i, v_box);
      }
      ++m_grid.count[z][c];
    
      // loop over the atoms
      topology::Atom_Iterator at_it = cg_it.begin(),
        at_to = cg_it.end();
      for( ; at_it != at_to; ++at_it){
        pos(*at_it) += trans;
      } // atoms
    } // solvent cg's
  } // charge-group based cutoff
  else {  // atomic cutoff
       
    // we can do all atoms in one go ...
    unsigned int i = 0;
    for( ; i < topo.num_atoms(); ++i){

      v=pos(i);
      v_box = v;
      periodicity.put_into_box(v_box);
      trans = v_box - v;

      // now grid the atom
       
      const int x = int((v_box(0) + 0.5 * abs(conf.current().box(0))) / m_grid.a);
      const int y = int((v_box(1) + 0.5 * abs(conf.current().box(1))) / m_grid.b);
      const int z = int((v_box(2) + 0.5 * abs(conf.current().box(2))) / m_grid.c);
  
      const int c = y * m_grid.Na + x;

      if(!(m_grid.cell_start.size() > unsigned(z) &&
           m_grid.cell_start[z].size() > unsigned(c)))
        return 1;
    
      if(!(m_grid.p_cell.size() > unsigned(z) &&
	   m_grid.p_cell[z].size() > unsigned(m_grid.cell_start[z][c] + m_grid.count[z][c])))
        return 1;

      if(!(m_grid.count.size() > unsigned(z) &&
	   m_grid.count[z].size() > unsigned(c)))
        return 1;
    
      m_grid.p_cell[z][m_grid.cell_start[z][c] + m_grid.count[z][c]] =
        Grid::Particle(i, v_box);

      ++m_grid.count[z][c];
    
      // move the atom back
      pos(i) += trans;

    } // loop over atoms
  } // atomic cutoff
 
  
  
  // check that there was enough space
  int max = 0;
  for(int z=0; z < m_grid.Nc; ++z){
    for(int c=0; c<Nab; ++c){
      if(!(m_grid.count.size() > unsigned(z) &&
	     m_grid.count[z].size() > unsigned(c)))
        return 1;
      
      if (m_grid.count[z][c] > space){
	if (m_grid.count[z][c] - space > max){
	  max = m_grid.count[z][c] - space;
	  std::cout << "grid[" << z << "][" << c << "] count = " << m_grid.count[z][c]
		    << " space = " << space << " new max = " << max << std::endl;
	}
      }
    }    
  }

  if (max){
    std::string ptype;
    if(sim.param().pairlist.atomic_cutoff)
      ptype ="atoms";
    else
      ptype = "chargegroups";
    
    std::cout << "not enough space to put " << ptype << " into cells:\n\t"
	      << "available space " << space << "\n\t"
	      << "additional space required " << max << "\n\n" << std::endl;
    
    io::messages.add("Not enough space to put "+ptype+" into cells!",
		     "Grid Pairlist Algorithm",
		     io::message::notice);

    m_grid.shift_space += max + 1;
    
    return prepare_grid(topo, conf, sim);
  }
  
  return 0;
}

/**
 * The name of this method is slightly misleading, since it
 * takes the grid generated by prepare_grid (a wasteful representation
 * on a non-extended grid), and generates a compact (memory-wise) 
 * representation on an extended grid. So yea, it collapses the grid,
 * but it also extends it. ;-) Also generates a particle-index
 * to (extended) cell-index map.
 */
void interaction::Extended_Grid_Pairlist_Algorithm::collapse_grid()
{
  const int N = m_grid.Na * m_grid.Nb;

  for(int z=0; z<m_grid.Nc; ++z){
    int p_cell_index = m_grid.count[z][0];

    // starts out as the cell index of the first cell in the extended grid.
    int cell_index_ex = (m_grid.Nb_ex - m_grid.Nb) / 2 * m_grid.Na_ex +
      (m_grid.Na_ex - m_grid.Na) / 2;

    m_grid.cell_index[z].clear();
    
    // need the index for all particles...

    // for any particle index pi we want the cell index ci,
    // _in an extended grid_.

    // this loops over the particles in the first cell of the plane
    for(int j=0; j<m_grid.count[z][0]; ++j)
      m_grid.cell_index[z].push_back(cell_index_ex);
    
    // loop over all the other cells (i.e. except the first)
    for(int i=1; i < N; ++i){

      // this points to the contents of the i-th cell in 
      // (non-collapsed part of) the particle vector:
      const int start = m_grid.cell_start[z][i];
      // the cell start index now becomes the number of particles we have
      // already "collapsed", since the subsequent particles will be inserted here
      m_grid.cell_start[z][i] = p_cell_index;

      // the current cell index in an "extended grid" formalism
      ++cell_index_ex;
      if (i % m_grid.Na == 0)
	cell_index_ex += m_grid.Na_ex - m_grid.Na;

      // loop over all particles in this cell.
      // the thing to notice about p_cell_index is that it only gets incremented 
      // for non-empty cells.
      for(int j=0; j < m_grid.count[z][i]; ++j, ++p_cell_index){

        // start is the beginning of cell i in the non-dense formalism,
        // p_cell_index is the number of particles in this plane that we have
        // already "densified".
        //std::cout << "c_g: " << z << " " << i << " " << p_cell_index << " " << start + j << std::endl;

	m_grid.p_cell[z][p_cell_index] = m_grid.p_cell[z][start + j];
        // this now contains the grid cell index for this particle
        // _on an extended grid_.
	m_grid.cell_index[z].push_back(cell_index_ex);

      } // loop over particles in cell

    } // loop over cells in plane

    // sentinel...
    m_grid.cell_start[z][N] = p_cell_index;

  } // loop over planes (z)
  
}

void interaction::Extended_Grid_Pairlist_Algorithm::grid_properties
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  const double s = sim.param().pairlist.grid_size;
  
  m_grid.Na = int(rint(abs(conf.current().box(0)) / s));
  m_grid.Nb = int(rint(abs(conf.current().box(1)) / s));
  m_grid.Nc = int(rint(abs(conf.current().box(2)) / s));
  
  m_grid.a =  abs(conf.current().box(0)) / m_grid.Na;
  m_grid.b =  abs(conf.current().box(1)) / m_grid.Nb;
  m_grid.c =  abs(conf.current().box(2))/ m_grid.Nc;  

  const int Ncell = m_grid.Na * m_grid.Nb * m_grid.Nc;
  
  double P = 0.0;
  if(sim.param().pairlist.atomic_cutoff)
    P = topo.num_atoms();
  else
    P = topo.num_chargegroups();
  // const double V = math::volume(conf.current().box, conf.boundary_type);
  // const double Vcell = m_grid.a * m_grid.b * m_grid.c;

  m_grid.Pcell = int(P / Ncell) + 1;

  m_grid.Na_ex = int(m_cutoff_long / m_grid.a);
  if (m_cutoff_long / m_grid.a > math::epsilon) ++m_grid.Na_ex;
  m_grid.Na_ex *= 2;
  m_grid.Na_ex += m_grid.Na;

  m_grid.Nb_ex = int(m_cutoff_long / m_grid.b);
  if (m_cutoff_long / m_grid.b > math::epsilon) ++m_grid.Nb_ex;
  m_grid.Nb_ex *= 2;
  m_grid.Nb_ex += m_grid.Nb;

  m_grid.Nc_ex = int(m_cutoff_long / m_grid.c);
  if (m_cutoff_long / m_grid.c > math::epsilon) ++m_grid.Nc_ex;
  m_grid.Nc_ex *= 2;
  m_grid.Nc_ex += m_grid.Nc;

  calculate_mask();

  // and the shift vectors
  m_shift_vector.clear();
  m_reverse_shift_vector.clear();

  for(int z=0; z<2; ++z){
    for(int y=-1; y<2; ++y){
      for(int x=-1; x<2; ++x){
	m_shift_vector.push_back(x * conf.current().box(0) + 
				 y * conf.current().box(1) + 
				 z * conf.current().box(2));

	m_reverse_shift_vector.push_back(-x * conf.current().box(0)
					 -y * conf.current().box(1)
					 -z * conf.current().box(2));
      }
    }
  }

}

/**
 * calculate the mask
 */
void interaction::Extended_Grid_Pairlist_Algorithm::calculate_mask()
{
  const double c2 = m_grid.c * m_grid.c;
  const double b2 = m_grid.b * m_grid.b;
  // const double a2 = m_grid.a * m_grid.a;

  m_grid.mask_z = int(m_cutoff_long / m_grid.c);
  if (m_cutoff_long / m_grid.c > math::epsilon) ++m_grid.mask_z;
  
  m_grid.mask.resize(m_grid.mask_z+1);
  

  double z_dist = 0.0, y_dist = 0.0;

  // special case of 0 plane
  {
    m_grid.mask[0].clear();
    int mask_y = int(sqrt(m_cutoff_long_2) / m_grid.b);
    if (sqrt(m_cutoff_long_2) / m_grid.b > math::epsilon) ++mask_y;
    
    for(int y=0; y <= mask_y; ++y){
      const int row = y * m_grid.Na_ex;

      if (y>1) y_dist = (y - 1) * (y - 1) * b2;
      else y_dist = 0.0;

      int mask_x = int(sqrt(m_cutoff_long_2 - y_dist) / m_grid.a);
      if (sqrt(m_cutoff_long_2 - y_dist) / m_grid.a > math::epsilon) ++mask_x;
      // begin 'till one past end
      m_grid.mask[0].push_back(row - mask_x);
      m_grid.mask[0].push_back(row + mask_x + 1);
    }
    // don't do self interaction over mask...
    m_grid.mask[0][0] = 1;
  }

  for(int z=1; z<=m_grid.mask_z; ++z){
    m_grid.mask[z].clear();

    if (z>1) z_dist = (z - 1) * (z - 1) * c2;
    else z_dist = 0.0;
    
    int mask_y = int(sqrt(m_cutoff_long_2 - z_dist) / m_grid.b);
    if (sqrt(m_cutoff_long_2 - z_dist) / m_grid.b > math::epsilon) ++mask_y;
    
    for(int y=mask_y; y>=0; --y){
      const int row = -y * m_grid.Na_ex;

      if (y>1) y_dist = (y - 1) * (y - 1) * b2;
      else y_dist = 0.0;

      int mask_x = int(sqrt(m_cutoff_long_2 - y_dist - z_dist) / m_grid.a);
      if (sqrt(m_cutoff_long_2 - y_dist - z_dist) / m_grid.a > math::epsilon) ++mask_x;

      // begin 'till one past end
      assert(m_grid.mask.size() > unsigned(z));
      m_grid.mask[z].push_back(row - mask_x);
      m_grid.mask[z].push_back(row + mask_x + 1);
    }
  }

  for(int z=1; z <= m_grid.mask_z; ++z){
    int row = 0;
    for(int y = m_grid.mask[z].size() - 4; y >= 0; y-=2){
      row += 2 * m_grid.Na_ex;
      m_grid.mask[z].push_back(row + m_grid.mask[z][y]);
      m_grid.mask[z].push_back(row + m_grid.mask[z][y+1]);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// prepare a plane
////////////////////////////////////////////////////////////////////////////////
void interaction::Extended_Grid_Pairlist_Algorithm::prepare_plane
(
 int z,
 std::vector<Grid::Particle> & p_plane, 
 std::vector<int> & cell_start
 )
{
  int z_shift = 0;
  if (z >= m_grid.Nc){
    z_shift = 9;
    z -= m_grid.Nc;
  }

  // reserve enough space for anything...
  p_plane.resize(m_grid.p_cell[z].size() * 4);
  cell_start.resize(m_grid.Na_ex * m_grid.Nb_ex + 1);

  const int a_ex = (m_grid.Na_ex - m_grid.Na) / 2;
  const int b_ex = (m_grid.Nb_ex - m_grid.Nb) / 2;
  
  // index into the (newly constructed) plane
  int j = 0;
  int cj = 0;
  int cs = 0;
  
  /*
    --------------------    1 = i
    |j |            |  |    2 = i_ex_to
    |  |            |  |    3 = i_ex
    --------------------    4 = i_to
    |  |            |  |
    |  |            |  |    first copy i_ex -> i_to to j
    |  |1 2      3 4|  |    then  copy i -> i_to
    |  |            |  |    then  copy i -> i_ex_to
    --------------------
    |  |            |  |    then add Na to i, i_ex_to, i_ex, i_to
    |  |            |  |
    --------------------
  */


  // the upper extended area
  int i = (m_grid.Nb - b_ex) * m_grid.Na;
  int i_to = i + m_grid.Na;
  int i_ex = i_to - a_ex;
  int i_ex_to = i + a_ex;

  for(int e=0; e<b_ex; ++e){

    const int pi = m_grid.cell_start[z][i];
    const int pi_to = m_grid.cell_start[z][i_to];
    const int pi_ex = m_grid.cell_start[z][i_ex];
    const int pi_ex_to = m_grid.cell_start[z][i_ex_to];

    // shift the particles
    for(int p=pi_ex; p < pi_to; ++p, ++j){
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift, m_shift_vector);
    }
    for(int p=pi; p < pi_to; ++p, ++j){
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 1, m_shift_vector);
    }
    for(int p=pi; p < pi_ex_to; ++p, ++j){
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 2, m_shift_vector);
    }

    // adapt the cell_start
    for(int c=i_ex; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i_ex] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i_ex];
    for(int c=i; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i];
    for(int c=i; c < i_ex_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_ex_to] - m_grid.cell_start[z][i];

    i += m_grid.Na;
    i_to += m_grid.Na;
    i_ex += m_grid.Na;
    i_ex_to += m_grid.Na;
  }

  // the center area
  
  i = 0;
  i_to = i + m_grid.Na;
  i_ex = i_to - a_ex;
  i_ex_to = i + a_ex;

  for(int e=0; e<m_grid.Nb; ++e){

    const int pi = m_grid.cell_start[z][i];
    const int pi_to = m_grid.cell_start[z][i_to];
    const int pi_ex = m_grid.cell_start[z][i_ex];
    const int pi_ex_to = m_grid.cell_start[z][i_ex_to];
    
    // shift particles
    for(int p=pi_ex; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 3, m_shift_vector);
    for(int p=pi; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 4, m_shift_vector);
    for(int p=pi; p < pi_ex_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 5, m_shift_vector);

    // adapt the cell_start
    for(int c=i_ex; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i_ex] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i_ex];
    for(int c=i; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i];
    for(int c=i; c < i_ex_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_ex_to] - m_grid.cell_start[z][i];
    // std::cerr << "cj after upper (" << e << ") = " << cj << std::endl;

    i += m_grid.Na;
    i_to += m_grid.Na;
    i_ex += m_grid.Na;
    i_ex_to += m_grid.Na;
  }

  // and the final part
  i = 0;
  i_to = i + m_grid.Na;
  i_ex = i_to - a_ex;
  i_ex_to = i + a_ex;

  for(int e=0; e<b_ex; ++e){

    const int pi = m_grid.cell_start[z][i];
    const int pi_to = m_grid.cell_start[z][i_to];
    const int pi_ex = m_grid.cell_start[z][i_ex];
    const int pi_ex_to = m_grid.cell_start[z][i_ex_to];
    
    // shift particles
    for(int p=pi_ex; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 6, m_shift_vector);
    for(int p=pi; p < pi_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 7, m_shift_vector);
    for(int p=pi; p < pi_ex_to; ++p, ++j)
      p_plane[j].shift(m_grid.p_cell[z][p], z_shift + 8, m_shift_vector);

    // adapt the cell_start
    for(int c=i_ex; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i_ex] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i_ex];
    for(int c=i; c < i_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_to] - m_grid.cell_start[z][i];
    for(int c=i; c < i_ex_to; ++c, ++cj)
      cell_start[cj] = m_grid.cell_start[z][c] - m_grid.cell_start[z][i] + cs;
    cs += m_grid.cell_start[z][i_ex_to] - m_grid.cell_start[z][i];

    i += m_grid.Na;
    i_to += m_grid.Na;
    i_ex += m_grid.Na;
    i_ex_to += m_grid.Na;
  }

  // sentinel
  cell_start[cj] = cs;

}

//////////////////////////////////////////////////////////////////////
// DEBUG FUNCTIONS
//////////////////////////////////////////////////////////////////////

/**
 * check the grid and print it...
 */
void interaction::Extended_Grid_Pairlist_Algorithm::print_grid()
{
  std::cout << "THE GRID\n========\n";

  const int N = m_grid.Na * m_grid.Nb;
  int num_P = 0;
  int errors = 0;

  const double ha = 0.5 * m_grid.Na * m_grid.a;
  const double hb = 0.5 * m_grid.Nb * m_grid.b;
  const double hc = 0.5 * m_grid.Nc * m_grid.c;

  for(int z=0; z<m_grid.Nc; ++z){
    
    // std::cout << "plane " << z << ":\n";
    int ci = 0;
    
    for(int i=0; i<N; ++i){
      if (m_grid.count[z][i] != m_grid.cell_start[z][i+1] - m_grid.cell_start[z][i]){
	std::cout << "grid error: count and cell_start don't match!" << std::endl;
	++errors;
      }
      
      if (m_grid.count[z][i] != 0){
	
	
	for(int n=0; n<m_grid.count[z][i]; ++n, ++ci){

	  // look up extended index
	  const int ind = m_grid.cell_index[z][ci];
	  const int indy = ind / m_grid.Na_ex - (m_grid.Nb_ex - m_grid.Nb) / 2;
	  const int indx = ind % m_grid.Na_ex - (m_grid.Na_ex - m_grid.Na) / 2;
	  if (indy * m_grid.Na + indx != i){
	    std::cout << "grid error: extended index and index don't match" << std::endl;
	    ++errors;
	  }

	  Grid::Particle const & p = m_grid.p_cell[z][m_grid.cell_start[z][i] + n];
	  
	  if (p.shift_index != 0){
	    std::cout << "grid error: shift index not zero in central box" << std::endl;
	    ++errors;
	  }
	  
	  if (p.z + hc < z * m_grid.c || p.z + hc > (z+1) * m_grid.c ||
	      p.x + ha < indx * m_grid.a || p.x + ha > (indx+1) * m_grid.a ||
	      p.y + hb < indy * m_grid.b || p.y + hb > (indy+1) * m_grid.b){
	    std::cout << "grid error: particle cog not inside cell" << std::endl;
	    // std::cout << "x=" << p.x << " y=" << p.y << " z=" << p.z << "\n"
	    // << "indx*a=" << indx*m_grid.a << " indy*b=" << indy*m_grid.b
	    // << " z*c=" << z*m_grid.c << "\n";
	    ++errors;
	  }
	  
	  std::cout << "cell " << i << " (" << ind << " = [" << indx << ", " << indy << "]) ["
		    << m_grid.count[z][i] << "] ";

	  std::cout << p.i << "\n";
	  
	  ++num_P;
	}
	
	std::cout << "\n";
      }
    }
  }
  
  std::cout << "particles on grid: " << num_P << "\n";
  
  if (errors == 0){
    std::cout << "no errors detected\n";
  }
  else
    std::cout << errors << " errors detected!\n";
}


void interaction::Extended_Grid_Pairlist_Algorithm::print_mask()
{
  std::cout << "\tmask\n";
  for(int z=0; z <= m_grid.mask_z; ++z){
    std::cout << "\n\tplane " << z << ":\n\t";

    for(unsigned int y=0; y < m_grid.mask[z].size(); y+=2){
      
      assert(m_grid.mask.size() > unsigned(z));
      assert(m_grid.mask[z].size() > y+1);

      std::cout << std::setw(5) << m_grid.mask[z][y] << " -> " 
		<< std::setw(5) << m_grid.mask[z][y+1] << "\t";

      if ((y + 2) % 10 == 0) std::cout << "\n\t";

    }
  }
  std::cout << "\n";
}


/**
 * check the grid and print it...
 */
void interaction::Extended_Grid_Pairlist_Algorithm::print_plane
(
 int z,
 std::vector<Grid::Particle> & p_plane, 
 std::vector<int> & cell_start
)
{
  std::cout << "PLANE " << std::setw(3) << z << "\n=========\n";

  std::cout.precision(3);
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  
  const int N = m_grid.Na_ex * m_grid.Nb_ex;
  int num_P = 0;

  // const double ha = 0.5 * m_grid.Na * m_grid.a;
  // const double hb = 0.5 * m_grid.Nb * m_grid.b;
  // const double hc = 0.5 * m_grid.Nc * m_grid.c;

  int ci = 0;
    
  for(int i=0; i<N; ++i){
    
    if (cell_start[i+1] - cell_start[i] != 0){
	
      std::cout << "cell " << i << " [" << cell_start[i+1] - cell_start[i] << "] = ";
      
      for(int n = cell_start[i]; n < cell_start[i+1]; ++n){
	
	Grid::Particle const & p = p_plane[n];
	
	std::cout << p.i << " (" << p.x << " | " << p.y << " | " << p.z << " : " << p.shift_index << ") ";
	++num_P;
      }
      
      std::cout << "\n";
	++ci;
    }
  }
  
  std::cout << "particles on plane: " << num_P << "\n";
}


