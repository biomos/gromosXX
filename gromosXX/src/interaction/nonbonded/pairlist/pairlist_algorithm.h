/**
 * @file pairlist_algorithm.h
 * base class for the pairlist(s) generating
 * algorithms.
 */

#ifndef INCLUDED_PAIRLIST_ALGORITHM_H
#define INCLUDED_PAIRLIST_ALGORITHM_H

#include <algorithm.h>

namespace interaction
{
  class Storage;
  class Pairlist;
  class Nonbonded_Parameter;
  
  /**
   * @class Pairlist_Algorithm
   * creates a pairlist.
   */
  class Pairlist_Algorithm : public algorithm::Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Pairlist_Algorithm()
      : algorithm::Algorithm("PairlistAlgorithm"),
	m_param(NULL)
    {}

    /**
     * destructor.
     */
    virtual ~Pairlist_Algorithm() {}

    void set_parameter(Nonbonded_Parameter * param) { m_param = param; }
    
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
     * no output... (done from forcefield)
     */
    virtual void print_timing(std::ostream & os) {}
    
    /**
     * timing accessor.
     * so needs to be accessed...
     */
    double timing() { return m_timing; }

    /**
     * apply the algorithm
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim) 
    {
      std::cerr << "don't call apply on a pairlist algorithm -- use update" << std::endl;
      assert(false);
      return 1;
    }

  protected:
    /**
     * nonbonded parameters (needed to construct the Innerloop).
     */
    Nonbonded_Parameter * m_param;

  };
} // interaction

#endif
