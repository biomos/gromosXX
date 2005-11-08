/**
 * @file stochastic.h
 * stochastic dynamics
 */

#ifndef INCLUDED_STOCHASTIC_H
#define INCLUDED_STOCHASTIC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace algorithm
{
  /**
   * @class Stochastic_Dynamics_Pos
   * implements stochastic dynamics (langevin equations)
   * position calculation
   */
  class Stochastic_Dynamics_Pos : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Stochastic_Dynamics_Pos();

    /**
     * Destructor.
     */
    virtual ~Stochastic_Dynamics_Pos() { gsl_rng_free(m_rng); }
    
    /**
     * Stochastic Dynamics: calculate positions
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);
    
    /**
     * get random number generator
     */
    gsl_rng * rng() { return m_rng; }

  private:
    /**
     * calculate friction coefficients
     */
    int calc_friction_coeff(topology::Topology &topo, 
			    configuration::Configuration &conf,
			    simulation::Simulation &sim);
    
    /**
     * random number generator
     */
    gsl_rng * m_rng;

  };

  /**
   * @class Stochastic_Dynamics_Pos
   * implements stochastic dynamics (langevin equations)
   * add stochastic integrals
   */
  class Stochastic_Dynamics_Int : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Stochastic_Dynamics_Int(gsl_rng * rng) : Algorithm("Stochastic_Dynamics_Int"), m_rng(rng) {}

    /**
     * Destructor.
     */
    virtual ~Stochastic_Dynamics_Int(){}
    
    /**
     * Stochastic Dynamics: add stochastic integrals
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

  private:
    gsl_rng * m_rng;
  };

} // algorithm

#endif

