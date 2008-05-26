/**
 * @file shake.h
 * the shake algorithm.
 */

#ifndef INCLUDED_SHAKE_H
#define INCLUDED_SHAKE_H

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm
{
  /**
   * @class Shake
   * implements the shake algorithm.
   */
  class Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Shake(double const tolerance = 0.000001, 
	  int const max_iterations = 1000,
	  std::string const name = "Shake");

    /**
     * Destructor.
     */
    virtual ~Shake();
        
    /**
     * apply shake.
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
    double const tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const max_iterations()const {return m_max_iterations;}

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

    template<math::boundary_enum B, math::virial_enum V>
    int shake_iteration
    (
     topology::Topology const &topo,
     configuration::Configuration & conf,
     bool & convergence,
     int first,
     std::vector<bool> &skip_now,
     std::vector<bool> &skip_next,
     std::vector<topology::two_body_term_struct> const & constr,
     double dt,
     math::Periodicity<B> const & periodicity
     );
    
    
    template<math::boundary_enum B, math::virial_enum V>
    int dih_constr_iteration
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     bool & convergence,
     std::vector<bool> & skip_now,
     std::vector<bool> & skip_next,
     math::Periodicity<B> const & periodicity
     );


    template<math::boundary_enum B, math::virial_enum V>
    void solute
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation const & sim,
     int const max_iterations,
     int & error
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
