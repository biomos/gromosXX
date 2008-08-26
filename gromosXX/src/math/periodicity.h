/**
 * @file periodicity.h
 * periodic boundary conditions (triclinic)
 */

#ifndef INCLUDED_PERIODICITY_H
#define INCLUDED_PERIODICITY_H

#include "boundary_implementation.h"

namespace math
{

  /**
   * @class Periodicity
   * the periodic boundary condition functions.
   */
  template<boundary_enum b>
  class Periodicity : public Boundary_Implementation<b>
  {
  public:
    /**
     * Constructor.
     * If b is any no specific code will be generated,
     * otherwise one can specify optimized code.
     */
    Periodicity(math::Box const & bb);
    /**
     * puts a vector into the box (centered at (0, 0, 0).
     */
    void put_into_box(Vec &v)const;
    /**
     * puts a vector into the box centered at (Kx/2, Ly/2, Mz/2).
     */
    void put_into_positive_box(Vec &v)const;
    /**
     * put chargegroups into the box.
     */
    void put_chargegroups_into_box(configuration::Configuration & conf,
				   topology::Topology const & topo )const;
    
    /**
     * put chargegroups into the box and save the lattice shifts.
     */
    void put_chargegroups_into_box_saving_shifts(
                                   configuration::Configuration & conf,
				   topology::Topology const & topo )const;


    void gather_chargegroups(configuration::Configuration & conf, 
			     topology::Topology const & topo)const;

    void gather_molecules_into_box(configuration::Configuration & conf, 
				   topology::Topology const & topo)const;
    
  private:
    
  };
  
} // math

// inline functions
#include "periodicity.cc"

#endif
  
  
    
