/**
 * @file perturbed_shake.h
 * the perturbed shake algorithm.
 */

#ifndef INCLUDED_PERTURBED_SHAKE_H
#define INCLUDED_PERTURBED_SHAKE_H

namespace algorithm
{
  /**
   * @class Perturbed_Shake
   * implements perturbed SHAKE.
   */
  template <typename t_simulation>
  class Perturbed_Shake
    : public Shake<t_simulation>
  {
  public:
    typedef t_simulation simulation_type;
    /**
     * Constructor.
     */
    Perturbed_Shake(double const tolerance = 0.000001, 
		    int const max_iteractions = 1000);
    /**
     * Destructor.
     */
    virtual ~Perturbed_Shake();

    /**
     * initialization
     */
    void init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
	      io::InInput &input);
    
    /**
     * shake solute.
     */
    int solute(typename simulation_type::topology_type &topo,
	       typename simulation_type::system_type &sys,
	       double dt);
    /**
     * add all bonds to the solute constraint vector and
     * remove them from the bond vector.
     */
    virtual void
    add_bond_length_constraints(typename t_simulation::topology_type &topo);
    /**
     * add bonds connecting an atom of type iac to the
     * constraint vector and remove from the bond vector.
     */
    virtual void
    add_bond_length_constraints(int iac, std::vector<int> const &atom_iac,
				typename t_simulation::topology_type &topo);
    
    /**
     * add bonds connecting an atom of mass mass to the
     * constraint vector and remove from the bond vector.
     */
    virtual void
    add_bond_length_constraints(double mass,
				math::SArray const &atom_mass,
				typename t_simulation::topology_type &topo);
  };
  
} // algorithm

// template methods
#include "perturbed_shake.tcc"

#endif

  
  
