/**
 * @file chargegroup_grid.h
 * create a grid for the
 * chargegroups.
 */

#ifndef INCLUDED_CHARGEGROUP_GRID_H
#define INCLUDED_CHARGEGROUP_GRID_H

namespace interaction
{
  template<typename t_simulation, bool periodic=false>
  class Cell_Cell_Iterator;

  /**
   * @class Chargegroup_Grid
   */
  template<typename t_simulation>
  class Chargegroup_Grid
  {
  public:
    struct mask_struct
    {
      /**
       * mask: k
       */
      std::pair<int, int> k;
      /**
       * mask: l
       */
      std::vector<std::pair<int, int> > l;
      /**
       * mask: m
       */
      std::vector<std::vector<std::pair<int, int> > > m;

      // accessors...
      int lower_k()const 
      {
	return k.first;
      }
      int upper_k()const
      {
	return k.second;
      }
      int lower_l(int const k)const
      {
	return l[k-lower_k()].first;
      }
      int upper_l(int const k)const
      {
	return l[k-lower_k()].second;
      }
      int lower_m(int const k, int const l)const
      {
	return m[k-lower_k()][l-lower_l(k)].first;
      }
      int upper_m(int const k, int const l)const
      {
	return m[k-lower_k()][l-lower_l(k)].second;
      }
      std::pair<int, int> & pair_k()
      {
	return k;
      }
      std::pair<int, int> & pair_l(int const k)
      {
	return l[k-lower_k()];
      }
      std::pair<int, int> & pair_m(int const k, int const l)
      {
	return m[k-lower_k()][l-lower_l(k)];
      }
      
    };

    /**
     * Constructor.
     */
    Chargegroup_Grid(
		     t_simulation 
		     & sim,
		     double const size,
		     double const cutoff);

    /**
     * the (chargegroup) grid.
     */
    std::vector<std::vector<std::vector<std::vector<size_t> > > > & grid();

    /**
     * add a point to the grid
     */
    void add(math::Vec & v, size_t i);
 
    /**
     * print the mask
     */
    void print_mask();
   
    /**
     * the number of cells in each direction
     */
    void num_cells(int cells[3]);
    
  protected:
    /**
     * the periodicity class (a reference)
     */
    math::Periodicity<t_simulation::system_type::boundary_type>
    & m_periodicity;
    /**
     * the grid.
     */
    std::vector<std::vector<std::vector<std::vector<size_t> > > > m_grid;
    /**
     * the real sizes of the cells.
     */
    math::Vec m_size;
    /**
     * the number of cells in every direction.
     */
    size_t m_num_cells[3];
    /**
     * the MASK!
     */
    mask_struct m_mask;
    
    friend class Cell_Cell_Iterator<t_simulation, true>;
    friend class Cell_Cell_Iterator<t_simulation, false>;
    
  };
  
  /**
   * @class Cell_Cell_Iterator
   */
  template<typename t_simulation, bool central>
  class Cell_Cell_Iterator
  {
  public:
    Cell_Cell_Iterator(Chargegroup_Grid<t_simulation> & grid, 
		       int cell[3]);
    bool eol();
    
    bool operator++();
    
    size_t operator*();

  protected:
    Chargegroup_Grid<t_simulation> & m_grid;
    bool m_eol;
    int m_cell[3];
    int k;
    int l;
    int m;
    bool m_outside;

    std::vector<size_t>::const_iterator it;
    std::vector<size_t>::const_iterator to;

  };
  
}

#include "chargegroup_grid.tcc"

#endif
