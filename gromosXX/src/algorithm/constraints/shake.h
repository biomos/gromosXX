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
  template<math::virial_enum do_virial>
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
		     bool quiet = false);

    /**
     * print out timing results.
     */
    virtual void print_timing(std::ostream & os);

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
     * time spent for solvent
     */
    double m_solvent_timing;

    template<math::boundary_enum b>
    int shake_iteration(topology::Topology const &topo,
			configuration::Configuration & conf,
			bool & convergence,
			int const first,
			std::vector<bool> &skip_now,
			std::vector<bool> &skip_next,
			std::vector<topology::two_body_term_struct> const & constr,
			double const dt,
			math::Periodicity<b> const & periodicity,
			bool do_constraint_force = false, size_t force_offset = 0);
    
    template<math::boundary_enum b>
    int solute(topology::Topology const & topo,
	       configuration::Configuration & conf,
	       double dt, int const max_iterations);
    
    template<math::boundary_enum b>
    int solvent(topology::Topology const & topo,
		configuration::Configuration & conf,
		double dt, int const max_iterations);
    

  };
  
} //algorithm

// template methods
#include "shake.tcc"

#endif
