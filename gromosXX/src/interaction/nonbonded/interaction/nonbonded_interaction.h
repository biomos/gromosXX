/**
 * @file nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_NONBONDED_INTERACTION_H
#define INCLUDED_NONBONDED_INTERACTION_H

namespace interaction
{

  template<typename t_interaction_spec, bool perturbed>
  class Pairlist_Algorithm;

  /**
   * @class Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  template<typename t_interaction_spec, bool perturbed>
  class Nonbonded_Interaction : 
    public Interaction,
    public Nonbonded_Parameter
  {
  public:    

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Nonbonded_Interaction(Pairlist_Algorithm<t_interaction_spec, perturbed> *pa);
    
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
     * Pairlist algorithm
     */
    Pairlist_Algorithm<t_interaction_spec, perturbed> & pairlist_algorithm()
    {
      return *m_pairlist_algorithm;
    }
    /**
     * printing function for timings
     */  
    virtual void print_timing(std::ostream & os);
      
    /**
     * size the arrays of storage.
     */
    void initialize(topology::Topology const & topo,
		    configuration::Configuration const & conf,
		    simulation::Simulation const & sim);

  protected:
    /**
     * the pairlist update algorithm.
     */
    Pairlist_Algorithm<t_interaction_spec, perturbed> *m_pairlist_algorithm;
    /**
     * a vector of nonbonded sets
     */
    std::vector<Nonbonded_Set<t_interaction_spec, perturbed> > m_nonbonded_set;
  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
