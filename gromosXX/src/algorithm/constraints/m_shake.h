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
     * rank and size for parallelization
     */
    int m_rank, m_size;
    /**
     * the factor matrix
     */
    math::Matrix factor;
    /**
     * the constraint lengths squared
     */
    math::Vec constr_length2;
    /**
     * inverted masses
     */
    math::Vec mass_i;

    inline
    int m_shake_molecule
    (
     configuration::Configuration & conf,
     bool & convergence,
     int first,
     double dt2i,
     const std::vector<topology::two_body_term_struct> & constr,
     math::GenericVec<math::Vec> const & dist_old,
     bool do_virial
     );
    
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
