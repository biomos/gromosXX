/**
 * @file shake.h
 * the shake algorithm.
 */

#ifndef INCLUDED_SHAKE_H
#define INCLUDED_SHAKE_H

namespace algorithm
{
  /**
   * @class Shake
   * implements the shake algorithm.
   */
  template<typename t_simulation>
  class Shake
  {
  public:
    typedef t_simulation simulation_type;
    /**
     * Constructor.
     */
    Shake(double const tolerance = 0.000001, int const max_iterations = 1000);
    /**
     * Destructor.
     */
    virtual ~Shake();
    
    /**
     * initialization.
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
     * shake solvent.
     */
    int solvent(typename simulation_type::topology_type &topo,
		typename simulation_type::system_type &sys,
		double dt);

    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * add bond type.
     */
    void add_bond_type(interaction::bond_type_struct s);
    /**
     * add bond type.
     */
    void add_bond_type(double K, double r0);

    /**
     * add all bonds to the solute constraint vector and
     * remove them from the bond vector.
     */
    virtual void
    add_bond_length_constraints(typename t_simulation::
				topology_type &topo);
    /**
     * add bonds connecting an atom of type iac to the
     * constraint vector and remove from the bond vector.
     */
    virtual void 
    add_bond_length_constraints(int iac, std::vector<int> 
				const &atom_iac,
				typename t_simulation::
				topology_type &topo);
    
    /**
     * add bonds connecting an atom of mass mass to the
     * constraint vector and remove from the bond vector.
     */
    virtual void 
    add_bond_length_constraints(double mass,
				math::SArray const &atom_mass,
				typename t_simulation::
				topology_type &topo);

  protected:

    template<typename t_distance_struct>
    bool _shake(typename simulation_type::topology_type const &topo,
		typename simulation_type::system_type &sys,
		int const first, 
		std::vector<bool> &skip_now,
		std::vector<bool> &skip_next,
		// std::vector<simulation::compound::distance_constraint_struct>
		std::vector<t_distance_struct>
		& constr, 
		bool do_constraint_force = false, size_t force_offset = 0);

    double m_tolerance;
    const int max_iterations;
    
    std::vector<interaction::bond_type_struct> m_bond_parameter;
    
  };
  
} //algorithm

// template methods
#include "shake.tcc"

#endif
  
