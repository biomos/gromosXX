/**
 * @file perturbed_flexible_constraint.h
 * the flexible shake algorithm
 */

#ifndef INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H
#define INCLUDED_PERTURBED_FLEXIBLE_CONSTRAINT_H

namespace algorithm
{
  /**
   * @class Perturbed_Flexible_Constraint
   * calculates the flexible constraint distance of perturbed constraints...
   */
  template<typename t_simulation, interaction::virial_enum do_virial>
  class Perturbed_Flexible_Constraint :
    public Flexible_Constraint<t_simulation, do_virial> 
  {
  public:
    /**
     * Constructor.
     */
    Perturbed_Flexible_Constraint(double const tolerance = 0.000001,
			int const max_iterations = 1000);
    
    /**
     * initialization.
     */
    void init(t_simulation &sim, io::Argument &args, io::InTopology &topo,
	      io::InInput &input);

    /**
     * shake solute.slightly different than in class Flexible_Constraint
     * @override
     */
     int solute(typename simulation_type::topology_type &topo,
          typename simulation_type::system_type &sys,
          double dt);

    // Accessor
    std::vector<double> const & A_K()const;
    std::vector<double> & A_K();
    std::vector<double> const & B_K()const;
    std::vector<double> & B_K();
    
    /**
     * add all bonds to the flexible solute constraint vector and
     * remove them from the bond vector.
     */
    virtual void
    add_bond_length_constraints(typename simulation_type::topology_type &topo);
    /**
     * add bonds connecting an atom of type iac to the
     * flexible constraint vector and remove from the bond vector.
     */
    virtual void
    add_bond_length_constraints(int iac,
				std::vector<int> const &atom_iac,
				typename simulation_type::topology_type &topo);
    
    /**
     * add bonds connecting an atom of mass mass to the
     * flexible constraint vector and remove from the bond vector.
     */
    virtual void
    add_bond_length_constraints(double mass,
				math::SArray const &atom_mass,
				typename simulation_type::topology_type &topo);
    
  private:
    int m_lfcon;
    std::vector<double> m_A_K;
    std::vector<double> m_B_K;
    
  };
  
  
} // algorithm

#include "perturbed_flexible_constraint.tcc"

#endif
