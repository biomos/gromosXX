/**
 * @file monte_carlo.h
 * monte carlo move
 */

#ifndef INCLUDED_MONTE_CARLO_H
#define INCLUDED_MONTE_CARLO_H

#include <gsl/gsl_rng.h>

namespace interaction
{
  class Forcefield;
}

namespace algorithm
{
  /**
   * @class Monte_Carlo
   * implements a Monte Carlo move
   */
  class Monte_Carlo : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Monte_Carlo(interaction::Forcefield *ff) : Algorithm("MonteCarlo"), m_ff(*ff)
    {
      // enable control via environment variables
      gsl_rng_env_setup();
      const gsl_rng_type * rng_type = gsl_rng_default;
      
      // get the random number generator
      m_rng = gsl_rng_alloc(rng_type);
    }

    /**
     * Destructor.
     */
    virtual ~Monte_Carlo()
    {
      gsl_rng_free(m_rng);
    }
    
    /**
     * Monte Carlo step
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
		     bool quiet = false) 
    {
      if (!quiet)
	os << "\tMonte Carlo move\nEND\n";

      // set seed if not set by environment variable
      if (gsl_rng_default_seed == 0)
	gsl_rng_set (m_rng, sim.param().start.ig);

      for(unsigned int i = 1; i < sim.param().multibath.multibath.bath_index().size(); ++i){
        if(sim.param().multibath.multibath.bath(i).temperature !=
                sim.param().multibath.multibath.bath(0).temperature){
          io::messages.add("Chemical MONTECARLO only possible if all baths "
                           " have the same temperature.",
                           "Monte_Carlo", io::message::error);
          return -1;
        }
      }
      
      return 0;
    };

  private:
    /**
     * try a Monte-Carlo move
     */
    int move(topology::Topology & topo,
	     configuration::Configuration & conf,
	     simulation::Simulation & sim);
    /**
     * accept or reject new configuration
     */
    int accept(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim);
    
    /**
     * random number generator (GSL)
     */
    gsl_rng * m_rng;
    
    /**
     * Forcefield
     */
    interaction::Forcefield & m_ff;
    
    /**
     * lambda before MC move (for rejected moves)
     */
    double m_lambda_old;

  };
 
} // algorithm

#endif

