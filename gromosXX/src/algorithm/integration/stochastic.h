/**
 * @file stochastic.h
 * stochastic dynamics
 */

#ifndef INCLUDED_STOCHASTIC_H
#define INCLUDED_STOCHASTIC_H

#include <gromosXX/math/random.h>

namespace simulation {
  class Parameter;
}

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
    Stochastic_Dynamics_Pos(const simulation::Parameter& param);

    /**
     * Destructor.
     */
    virtual ~Stochastic_Dynamics_Pos() { delete m_rng; }
    
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
    math::RandomGenerator * rng() { return m_rng; }
    
    /**
     * @struct RandomVectors
     * storage for random vectors
     */
    struct RandomVectors {
      math::VArray vrand1;
      math::VArray vrand2;
      math::VArray vrand3;
      math::VArray vrand4;
      
      /**
       * resize the vectors
       */
      void resize(unsigned int size) {
        vrand1.resize(size);
        vrand2.resize(size);
        vrand3.resize(size);
        vrand4.resize(size);
      }
    };
    /**
     * accessor to the random vectos
     */
    algorithm::Stochastic_Dynamics_Pos::RandomVectors& random_vectors()
    { return m_vrand; }
    
    /**
     * Stochastic Dynamics: calculate positions
     */
    template<math::boundary_enum B>
    int _apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);
    /**
     * calculate friction coefficients
     */
    template<math::boundary_enum B>
    int calc_friction_coeff(topology::Topology &topo, 
			    configuration::Configuration &conf,
			    simulation::Simulation &sim);
    
    /**
     * random number generator
     */
    math::RandomGenerator * m_rng;
    /**
     * random vectors
     */
    RandomVectors m_vrand;
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
    Stochastic_Dynamics_Int(math::RandomGenerator * rng,
                            Stochastic_Dynamics_Pos::RandomVectors *vrand) : 
        Algorithm("Stochastic_Dynamics_Int"), m_rng(rng),
        m_vrand(vrand) {}

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

  protected:
    /**
     * pointer to the random number generator of the
     * Stochastic_Dynamics_Pos algorithm
     */
    math::RandomGenerator * m_rng;
    /**
     * pointer to the random vectors of the Stochastic_Dynamics_Pos algorithm
     */
    Stochastic_Dynamics_Pos::RandomVectors *m_vrand;
    
    /**
     * Stochastic Dynamics: add stochastic integrals
     */
    template<math::boundary_enum B>
    int _apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);
  };

} // algorithm

#endif

