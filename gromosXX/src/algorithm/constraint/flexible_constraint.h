/**
 * @file flexible_constraint.h
 * the flexible shake algorithm
 */

#ifndef INCLUDED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_FLEXIBLE_CONSTRAINT_H

namespace algorithm
{
  /**
   * @class Flexible_Constraint
   * calculates the flexible constraint distance
   */
  template<typename t_simulation>
  class Flexible_Constraint :
    public Shake<t_simulation> 
  {
  public:
    /**
     * Constructor.
     */
    Flexible_Constraint(double const tolerance = 0.000001,
			int const max_iterations = 1000);

    /**
     * initialization.
     */
    void init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
	      io::InInput &input);

    /**
     * shake solute.
     * @override
     */
    int solute(typename simulation_type::topology_type &topo,
	       typename simulation_type::system_type &sys,
	       double dt);

    // Accessor
    std::vector<double> const & r0()const;
    std::vector<double> & r0();
    std::vector<double> const & vel()const;
    std::vector<double> & vel();
    std::vector<double> const & K()const;
    std::vector<double> & K();
    
    /**
     * add all bonds to the flexible solute constraint vector and
     * remove them from the bond vector.
     */
    void add_bond_length_constraints(simulation::Solute &solute);
    /**
     * add bonds connecting an atom of type iac to the
     * flexible constraint vector and remove from the bond vector.
     */
    void add_bond_length_constraints(int iac,
				     std::vector<int> const &atom_iac,
				     simulation::Solute &solute);
    
    /**
     * add bonds connecting an atom of mass mass to the
     * flexible constraint vector and remove from the bond vector.
     */
    void add_bond_length_constraints(double mass,
				     math::SArray const &atom_mass,
				     simulation::Solute &solute);

  private:
    void calc_distance(typename simulation_type::topology_type const &topo, 
		       typename simulation_type::system_type &sys, 
		       int const first,
		       std::vector<simulation::compound::
		       distance_constraint_struct>
		       & constr, double const dt);

    std::vector<double> force_on_constraint;
    double Ekin;
    double Epot;

    std::vector<double> m_r0;
    std::vector<double> m_vel;
    std::vector<double> m_K;
    
    int m_lfcon;

  };
  
  
} // algorithm

#include "flexible_constraint.tcc"

#endif
