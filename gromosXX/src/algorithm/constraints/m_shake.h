/**
 * @file m_shake.h
 * the m_shake algorithm.
 */

#ifndef INCLUDED_M_SHAKE_H
#define INCLUDED_M_SHAKE_H

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm
{
  /**
   * @class M_Shake
   * implements the m_shake algorithm.
   */
  class M_Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    M_Shake(double const tolerance = 0.000001,
	  int const max_iterations = 1000,
	  std::string const name = "M_Shake");

    /**
     * Destructor.
     */
    virtual ~M_Shake();
        
    /**
     * apply m_shake.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const & tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const & max_iterations()const {return m_max_iterations;}

    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const
    {
      return m_parameter;
    }
    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter()
    {
      return m_parameter;
    }
    /**
     * accessor to the constrained atoms
     */
    std::set<unsigned int> & constrained_atoms() {
      return m_constrained_atoms;
    }
     /**
     * accessor to the constrained atoms
     */
    const std::set<unsigned int> & constrained_atoms() const {
      return m_constrained_atoms;
    }

    /**
     * initialize startup positions and velocities
     * if required.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations
     */
    const int m_max_iterations;
    /**
     * bond parameter
     */
    std::vector<interaction::bond_type_struct> m_parameter;
    /**
     * the atoms that are involved in the contraints
     */
    std::set<unsigned int> m_constrained_atoms;
    /** 
     * rank and size for parallelization
     */
    int m_rank, m_size;

    template<math::boundary_enum B, math::virial_enum V>
    int m_shake_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     int first,
     double dt,
     std::vector<topology::two_body_term_struct> const & constr,
     math::GenericMatrix<double> const & factor,
     math::GenericVec<double> const & constr_length2,
     math::GenericVec<math::Vec> const & dist_old
     );
    
    template<math::boundary_enum B, math::virial_enum V>
    void solvent
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     double dt, int const max_iterations,
     int & error
     );

  };
  
} //algorithm

#endif
