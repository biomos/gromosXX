/**
 * @file pairlist_algorithm.h
 * base class for the pairlist(s) generating
 * algorithms.
 */

#ifndef INCLUDED_PAIRLIST_ALGORITHM_H
#define INCLUDED_PAIRLIST_ALGORITHM_H

namespace interaction
{
  class Storage;
  class Pairlist;
  
  /**
   * @class Pairlist_Algorithm
   * creates a pairlist.
   */
  class Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Pairlist_Algorithm()
      : m_param(NULL),
	m_timing(0.0) {}

    /**
     * destructor.
     */
    virtual ~Pairlist_Algorithm() {}

    int init(Nonbonded_Parameter * param) { m_param = param; return 0; }
    
    /**
     * prepare the pairlist(s).
     */
    virtual void prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim) = 0;
    /**
     * update the pairlist
     */
    virtual void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation &sim,
			interaction::Storage &storage,
			interaction::Pairlist &pairlist,
			unsigned int begin, unsigned int end, 
			unsigned int stride) = 0;

    /**
     * update the pairlist, separating perturbed and non-perturbed interactions
     */
    virtual void update_perturbed(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  interaction::Storage & storage,
				  interaction::Pairlist & pairlist,
				  interaction::Pairlist & perturbed_pairlist,
				  unsigned int begin, unsigned int end, 
				  unsigned int stride) = 0;
    
    /**
     * print timing
     */
    virtual void print_timing(std::ostream & os) {}
    
    /**
     * timing accessor.
     */
    double timing() { return m_timing; }

  protected:
    /**
     * nonbonded parameters (needed to construct the Innerloop).
     */
    Nonbonded_Parameter * m_param;
    /**
     * timing information.
     */
    double m_timing;

  };
} // interaction

#endif
