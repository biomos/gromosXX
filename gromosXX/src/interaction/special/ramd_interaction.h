/**
 * @file ramd_interaction.h
 * Random Acceleration MD interaction
 */

#ifndef INCLUDED_RAMD_INTERACTION_H
#define INCLUDED_RAMD_INTERACTION_H

#include <gsl/gsl_rng.h>

namespace interaction
{
  /**
   * @class RAMD_Interaction
   * calculates the Random Acceleration MD interaction
   */
  class RAMD_Interaction : 
    public Interaction
  {
  public:
    /**
     * Constructor.
     */
    RAMD_Interaction() : Interaction("RAMD") {

      // just copied from Monte Carlo!!
      gsl_rng_env_setup();
      const gsl_rng_type * rng_type = gsl_rng_default;
      
      // get the random number generator
      m_rng = gsl_rng_alloc(rng_type);
    }
    
    /**
     * Destructor.
     */
    virtual ~RAMD_Interaction() 
    {
      gsl_rng_free(m_rng);
    }

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false);

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);
    
  private:
    /**
     * random number generator (GSL)
     */
    gsl_rng * m_rng;
    
  };
  
} // interaction

#endif
