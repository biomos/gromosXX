/**
 * @file chargegroup_grid.tcc
 * create a grid for the
 * chargegroups (or anything...)
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../../debug.h"

template<typename t_simulation>
inline
interaction::Chargegroup_Grid<t_simulation>
::Chargegroup_Grid(math::Periodicity<
		   t_simulation::system_type::boundary_type> 
		   & periodicity,
		   double const size,
		   double const cutoff)
  : m_periodicity(periodicity)
{
  DEBUG(8, "initialize grid");
  DEBUG(8, "\tcells\t\tsize");
  for(size_t d=0; d<3; ++d){
    double const length = sqrt(dot(periodicity.box()(d), periodicity.box()(d)));
    m_num_cells[d] = int(rint(length / size));
    m_size(d) = length / m_num_cells[d];
    DEBUG(8, "\t" << m_num_cells[d] << "\t\t" << m_size(d));
  }  
  
  m_grid.resize(m_num_cells[0]);
  for(size_t d=0; d<m_num_cells[0]; ++d){
    m_grid[d].resize(m_num_cells[1]);
    for(size_t e=0; e<m_num_cells[1]; ++e)
      m_grid[d][e].resize(m_num_cells[2]);
  }

  // initialize the shift vectors
  {
    int index=0;
    for(int k=-1; k<2; ++k){
      for(int l=-1; l<2; ++l){
	for(int m=-1; m<2; ++m, ++index){

	  m_shift[index].cell[0] = k*m_num_cells[0];
	  m_shift[index].cell[1] = l*m_num_cells[1];
	  m_shift[index].cell[2] = m*m_num_cells[2];

	  m_shift[index].pos = 
	    k*m_periodicity.box()(0) +
	    l*m_periodicity.box()(1) +
	    m*m_periodicity.box()(2);

	}
      }
    }
  }

  // loop corrections
  // K direction
  double m_k = sqrt(dot(m_periodicity.box()(0), m_periodicity.box()(0)));
  // L x M
  math::Vec cr = math::cross(m_periodicity.box()(1), m_periodicity.box()(2));
  // unit (L x M)
  cr = cr / sqrt(math::dot(cr, cr));
  // projection of K on unit vector (L x M)
  int cells_k = int(cutoff / fabs(dot(m_periodicity.box()(0), cr)) * m_num_cells[0]) + 1;
  DEBUG(8, "cells inside cutoff (k) " << cells_k);

  // L direction
  double m_l = sqrt(dot(m_periodicity.box()(1), m_periodicity.box()(1)));
  // K x M
  cr = math::cross(m_periodicity.box()(2), m_periodicity.box()(0));
  // unit (K x M)
  cr = cr / sqrt(math::dot(cr, cr));
  // projection of L on unit vector (K x M)
  int cells_l = int(cutoff / fabs(dot(m_periodicity.box()(1), cr)) * m_num_cells[1]) + 1;
  DEBUG(8, "cells inside cutoff (l) " << cells_l);

  // M direction
  double m_m = sqrt(dot(m_periodicity.box()(2), m_periodicity.box()(2)));
  // K x L
  cr = math::cross(m_periodicity.box()(0), m_periodicity.box()(1));
  // unit (K x L)
  cr = cr / sqrt(math::dot(cr, cr));
  // projection of M on unit vector (K x L)
  int cells_m = int(cutoff / fabs(dot(m_periodicity.box()(2), cr)) * m_num_cells[2]) + 1;
  DEBUG(8, "cells inside cutoff (m) " << cells_m);

  // get the longest diagonal:
  math::Vec diag1 = m_periodicity.box()(0) / int(m_num_cells[0]) +
    m_periodicity.box()(1) / int(m_num_cells[1]) +
    m_periodicity.box()(2) / int(m_num_cells[2]);

  math::Vec diag2 = m_periodicity.box()(0) / int(m_num_cells[0]) +
    m_periodicity.box()(1) / int(m_num_cells[1]) -
    m_periodicity.box()(2) / int(m_num_cells[2]);

  math::Vec diag3 = m_periodicity.box()(0) / int(m_num_cells[0]) -
    m_periodicity.box()(1) / int(m_num_cells[1]) +
    m_periodicity.box()(2) / int(m_num_cells[2]);

  math::Vec diag4 = m_periodicity.box()(1) / int(m_num_cells[1]) -
    m_periodicity.box()(0) / int(m_num_cells[0]) +
    m_periodicity.box()(2) / int(m_num_cells[2]);

  double d1 = dot(diag1, diag1);
  double d2 = dot(diag2, diag2);
  double d3 = dot(diag3, diag3);
  double d4 = dot(diag4, diag4);

  DEBUG(10, "diag1: " << d1);
  DEBUG(10, "diag2: " << d2);
  DEBUG(10, "diag3: " << d3);
  DEBUG(10, "diag4: " << d4);
  
  d1 = (d1 > d2) ? d1 : d2;
  d3 = (d3 > d4) ? d3 : d4;
  d1 = (d1 > d3) ? d1 : d3;
  d1 = sqrt(d1);

  // generate a mask!
  
  int in = 0, out = 0;
  int first;
  
  // in direction of the lattice vectors!

  std::pair<int, int> px(-cells_k, cells_k);
  m_mask.k = px;
  m_mask.l.resize(2*cells_k+1);
  m_mask.m.resize(2*cells_k+1);
  
  for(int i = -cells_k; i <= cells_k; ++i){
    // m_mask.m[i+cells_k].resize(2*cells_l+1);

    // set an impossible first value
    m_mask.l[i+cells_k].first = cells_l+1;
    DEBUG(15, "reset m_mask.l[" << i + cells_k << "] = " << cells_l+1);

    for(int j = -cells_l; j <= cells_l; ++j){
     
      first = cells_m + 1;
      for(int k = -cells_m; k <= cells_m; ++k){

	math::Vec d =
	  i * m_periodicity.box()(0) / int(m_num_cells[0]) + 
	  j * m_periodicity.box()(1) / int(m_num_cells[1]) + 
	  k * m_periodicity.box()(2) / int(m_num_cells[2]);

	// worst case:
	// it is the longest diagonal closer
	double dist = sqrt(dot(d,d));
	dist -= d1;	

	if (dist > cutoff){
	  ++out;
	  DEBUG(11, std::setw(6) << i << std::setw(6) << j 
		<< std::setw(6) << k << std::setw(8) << "out!" 
		<< std::setw(20) << dist << ")");
	  if (first < k){
	    // we have a range!
	    std::pair<int, int> p(first, k-1);

	    assert(int(m_mask.l.size()) > i+cells_k);
	    
	    // m_mask.m[i+cells_k][j+cells_l] = p;
	    m_mask.m[i+cells_k].push_back(p);
	    m_mask.l[i+cells_k].second = j;

	    DEBUG(15, "mask.m[" << i+cells_k << "," << m_mask.m[i+cells_k].size() - 1
		  << "] = " << p.first << " .. " << p.second);
	    DEBUG(15, "mask.l[" << i+cells_k << "] = " <<
		  m_mask.l[i+cells_k].first << " .. " <<
		  m_mask.l[i+cells_k].second);

	    first = cells_m + 1;
	    break;
	  } // at the end of a range
	} // outside
	else{
	  // inside
	  ++in;
	  DEBUG(11, std::setw(6) << i << std::setw(6) << j 
		<< std::setw(6) << k << std::setw(8) << "in!" 
		<< std::setw(20) << dist << ")");
	  if (first > k){
	    // beginning of a range
	    first = k;
	    assert(int(m_mask.l.size()) > i+cells_k);
	    if (m_mask.l[i+cells_k].first > j){
	      // also the beginning of the l range
	      m_mask.l[i+cells_k].first = j;
	      DEBUG(15, "assign m_mask.l[" << i+cells_k << "] = " << j);
	    }
	    
	  } // beginning of a range

	} // inside

      } // loop over m

      if (first < cells_m + 1){
	// still open!
	std::pair<int, int> p(first, cells_m);

	// m_mask.m[i+cells_k][j+cells_l] = p;
	m_mask.m[i+cells_k].push_back(p);
	m_mask.l[i+cells_k].second = j;

	DEBUG(15, "closing open range");
	DEBUG(15, "mask.m[" << i+cells_k << "," << m_mask.m[i+cells_k].size()-1
	      << "] = " << p.first << " .. " << p.second);
	DEBUG(15, "mask.l[" << i+cells_k << "] = " <<
	      m_mask.l[i+cells_k].first << " .. " <<
	      m_mask.l[i+cells_k].second);

      }
      
    }
  }

  DEBUG(8, "ratio: " << in << " / " << out);
  print_mask();
  
}

template<typename t_simulation>
inline void
interaction::Chargegroup_Grid<t_simulation>
::print_mask()
{
  std::cout << "MASK\n";
  
  for(int k=m_mask.lower_k(); k<=m_mask.upper_k(); ++k){

    for(int l=m_mask.lower_l(k); l<=m_mask.upper_l(k); ++l){

      std::cout << std::setw(5) << k << std::setw(5) << l
		<< std::setw(5) << m_mask.lower_m(k, l)
		<< " .. "
		<< std::setw(5) << m_mask.upper_m(k,l)
		<< "\n";
    }
  }
  std::cout << "END\n";
  
}

template<typename t_simulation>
inline std::vector<std::vector<std::vector<std::vector<size_t> > > > &
interaction::Chargegroup_Grid<t_simulation>
::grid()
{
  return m_grid;
}

template<typename t_simulation>
inline typename interaction::Chargegroup_Grid<t_simulation>::shift_struct &
interaction::Chargegroup_Grid<t_simulation>
::shift(size_t const i)
{
  return m_shift[i];
}

template<typename t_simulation>
inline void
interaction::Chargegroup_Grid<t_simulation>
::num_cells(int cells[3])
{
  for(int d=0; d<3; ++d)
    cells[d] = m_num_cells[d];
}

template<typename t_simulation>
inline void
interaction::Chargegroup_Grid<t_simulation>
::add(math::Vec & v, size_t i)
{
  math::Vec n;
  m_periodicity.box_components(v, n);
  DEBUG(10, "point: " << v(0) << " | " << v(1) << " | " << v(2));
  DEBUG(11, "n:     " << n(0) << " | " << n(1) << " | " << n(2));
  DEBUG(10, "grid:  " << int((n[0]+0.5)*m_num_cells[0]) << " | "
	<< int((n[1]+0.5)*m_num_cells[1]) << " | "
	<< int((n[2]+0.5)*m_num_cells[2]));

  int nx = int((n[0]+0.5)*m_num_cells[0]);
  int ny = int((n[1]+0.5)*m_num_cells[1]);
  int nz = int((n[2]+0.5)*m_num_cells[2]);
  
  m_grid[nx][ny][nz].push_back(i);
  
}

template<typename t_simulation, bool periodic>
inline 
interaction::Cell_Cell_Iterator<t_simulation, periodic>
::Cell_Cell_Iterator(Chargegroup_Grid<t_simulation> & grid,
		     int cell[3])
  : m_grid(grid),
    m_eol(false),
    m_outside(false)
{
  m_cell[0] = cell[0];
  m_cell[1] = cell[1];
  m_cell[2] = cell[2];

  DEBUG(10, "cell: " << m_cell[0] << "\t" << m_cell[1] << "\t" << m_cell[2]);
  
  // initialize the iterator
  if (!periodic)
    k = 0;
  else 
    k = m_grid.m_mask.lower_k();

  if (k + m_cell[0] < 0) k = - m_cell[0];
  if (k > m_grid.m_mask.upper_k() || 
      k + m_cell[0] >= int(m_grid.m_num_cells[0])){
    m_eol = true;
    DEBUG(10, "returning from constructor (k) -- empty");
    return;
  }
  
  if (!periodic)
    l = 0;
  else
    l = m_grid.m_mask.lower_l(k);

  if (l + m_cell[1] < 0) l = - m_cell[1];
  if (l > m_grid.m_mask.upper_l(k) || 
      l + m_cell[1] >= int(m_grid.m_num_cells[1])){
    DEBUG(10, "need to ++ because of l");
    it = to - 1;
    // reinit l
    l = m_grid.m_mask.upper_l(k);
    // initialize m
    m = m_grid.m_mask.upper_m(k, l);
    
    DEBUG(10, "++ from constructor (l)");
    m_eol = ! (++(*this));
    DEBUG(10, "returning from l with m_eol=" << m_eol);
    return;
  }
    
  if (!periodic)
    m = 1;
  else 
    m = m_grid.m_mask.lower_m(k, l);

  if (m + m_cell[2] < 0) m = - m_cell[2];
  if (m > m_grid.m_mask.upper_m(k, l) || 
      m + m_cell[2] >= int(m_grid.m_num_cells[2])){
    it = to - 1;
    // reinit m
    m = m_grid.m_mask.upper_m(k, l);

    DEBUG(10, "++ from constructor (m)");
    m_eol = ! (++(*this));
    DEBUG(10, "returning from m with m_eol=" << m_eol);
    return;
  }
  
  // everything in range
  DEBUG(10, "trying to initialize iterator with " << k << " | " << l << " | " << m);

  it = m_grid.grid()[k+m_cell[0]][l+m_cell[1]][m+m_cell[2]].begin();
  to = m_grid.grid()[k+m_cell[0]][l+m_cell[1]][m+m_cell[2]].end();
  
  if (it == to){
    DEBUG(10, "empty beginning");
    --it;
    DEBUG(10, "++ from constructor (empty beginning)");
    m_eol = ! (++(*this));
    return;
  }
  

}

template<typename t_simulation, bool periodic>
inline bool
interaction::Cell_Cell_Iterator<t_simulation, periodic>
::operator++()
{
  DEBUG(10, "++ with k = " << k << " l = " << l << " m = " << m);

  ++it;
  while(it == to){
    // go to next bin
    ++m;
    DEBUG(10, "++m = " << m);

    while (m > m_grid.m_mask.upper_m(k, l) ||
	   m + m_cell[2] >= int(m_grid.m_num_cells[2]) ){

      ++l;
      DEBUG(10, "++l = " << l);

      while(l > m_grid.m_mask.upper_l(k) ||
	    l + m_cell[1] >= int(m_grid.m_num_cells[1]) ){
	++k;
	DEBUG(10, "++k = " << k);

	if (k > m_grid.m_mask.upper_k() || 
	    k + m_cell[0] >= int(m_grid.m_num_cells[0]) ){
	  // end of the list reached!
	  DEBUG(8, "end reached");
	  return false;
	}
	// k is ok
	// we had a ++k, so k must be >0 (if !periodic)
	l = m_grid.m_mask.lower_l(k);
	DEBUG(10, "new l = " << l);
	if (l + m_cell[1] < 0) {
	  l = - m_cell[1];
	  DEBUG(10, "corrected l = " << l);
	}
	
      }
      // l is ok
      // we had an l++ so it should be ok...
      m = m_grid.m_mask.lower_m(l, k);
      DEBUG(10, "new m = " << m);
      if (m + m_cell[2] < 0){
	m = - m_cell[2];
	DEBUG(10, "corrected m = " << m);
      }
      
    }
    // m is ok
    DEBUG(10, "trying for chargegroups in " << k + m_cell[0] 
	  << " | " << l + m_cell[1] << " | " << m + m_cell[2]);
    DEBUG(10, "with the k = " << k 
	  << " l = " << l << " m = " << m);

    it = m_grid.grid()[k+m_cell[0]][l+m_cell[1]][m+m_cell[2]].begin();
    to = m_grid.grid()[k+m_cell[0]][l+m_cell[1]][m+m_cell[2]].end();

  }
      
  DEBUG(8, "found cg " << *it << " in " << k+m_cell[0] << " | " << l+m_cell[1]
	<< " | " << m+m_cell[2]);

  return true;
  
}

template<typename t_simulation, bool periodic>
inline size_t
interaction::Cell_Cell_Iterator<t_simulation, periodic>
::operator*()
{
  return *it;
}

template<typename t_simulation, bool periodic>
inline bool
interaction::Cell_Cell_Iterator<t_simulation, periodic>
::eol()
{
  return m_eol;
}

