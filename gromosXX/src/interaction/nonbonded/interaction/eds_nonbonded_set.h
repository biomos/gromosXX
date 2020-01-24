/**
 * @file eds_nonbonded_set.h
 * the (eds-perturbed) non bonded interactions for a set of atoms:
 * Lennard-Jones and Coulomb interactions
 * eds-perturbed and non-eds-perturbed.
 */

#ifndef INCLUDED_EDS_NONBONDED_SET_H
#define INCLUDED_EDS_NONBONDED_SET_H

namespace interaction
{
  /**
   * @class Eds_Nonbonded_Set
   * calculates the (eds-perturbed) nonbonded interactions.
   */
  class Eds_Nonbonded_Set : public Nonbonded_Set
  {
  public:    
    /**
     * Constructor.
     */
    Eds_Nonbonded_Set(Pairlist_Algorithm & pairlist_alg,
			    Nonbonded_Parameter & param,
			    int rank, int num_threads);

    /**
     * Destructor
     */
    virtual ~Eds_Nonbonded_Set() {}
    
    /**
     * initialize some things
     */
    virtual int init(topology::Topology const & topo,
		     configuration::Configuration const & conf,
		     simulation::Simulation const & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);
    
    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * store results in configuration
     */
    virtual int update_configuration(topology::Topology const & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation const & sim);
   
    PairlistContainer & perturbed_pairlist() { return m_perturbed_pairlist; }
    PairlistContainer const & perturbed_pairlist()const { return m_perturbed_pairlist; }

  protected:

    PairlistContainer m_perturbed_pairlist;

    Eds_Nonbonded_Outerloop m_eds_outerloop;
    
  };
  
} // interaction

#endif
