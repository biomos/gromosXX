/**
 * @file nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_NONBONDED_INTERACTION_H
#define INCLUDED_NONBONDED_INTERACTION_H

namespace interaction
{
  /**
   * @class Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  template<typename t_interaction_spec>
  class Nonbonded_Interaction : 
    public Interaction,
    public Nonbonded_Base,
    public Storage,
    public t_interaction_spec::nonbonded_innerloop_type
  {
  public:    

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Nonbonded_Interaction();
    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual int calculate_interactions(topology::Topology & topo,
				       configuration::Configuration & conf,
				       simulation::Simulation & sim);

    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     size_t const i, size_t const j);
    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     size_t const i, size_t const j,
			     int pc);

    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    size_t const i, size_t const j,
			    Periodicity_type const & periodicity);

    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    size_t const i, size_t const j,
			    Periodicity_type const & periodicity, int pc);
    
    // ACCESSORS
    /**
     * pairlist accessor.
     */
    Pairlist & pairlist() { return m_pairlist; }
    /**
     * const pairlist accessor.
     */
    Pairlist const & pairlist()const { return m_pairlist; }
    /**
     * perturbed pairlist accessor.
     */
    Pairlist & perturbed_pairlist() { return m_perturbed_pairlist; }
    /**
     * const perturbed pairlist accessor.
     */
    Pairlist const & perturbed_pairlist()const { return m_perturbed_pairlist; }

  protected:
    /**
     * size the arrays of storage.
     */
    void initialize(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim);
    
    /**
     * calculate the interactions.
     */
    void do_interactions(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation & sim,
			 std::vector<std::vector<size_t> > const & pairlist);
    /*    
			 Pairlist::iterator it, 
			 Pairlist::iterator to);
    */

    /**
     * calculate the 1,4-interactions.
     */
    void do_14_interactions(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    void do_RF_excluded_interactions(topology::Topology & topo,
				     configuration::Configuration & conf,
				     simulation::Simulation & sim);
 
    /**
     * the (shortrange) pairlist.
     */
    Pairlist m_pairlist;
    /**
     * the perturbed (shortrange) pairlist.
     */
    Pairlist m_perturbed_pairlist;
    
    /**
     * the pairlist update algorithm.
     */
    typename t_interaction_spec::pairlist_algorithm_type m_pairlist_algorithm;
    
  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
