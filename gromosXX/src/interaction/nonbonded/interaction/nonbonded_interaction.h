/**
 * @file nonbonded_interaction.h
 * the non bonded interactions:
 * Lennard-Jones and Coulomb interactions.
 */

#ifndef INCLUDED_NONBONDED_INTERACTION_H
#define INCLUDED_NONBONDED_INTERACTION_H

namespace interaction
{

  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Pairlist_Algorithm;

  /**
   * @class Nonbonded_Interaction
   * calculates the nonbonded interactions.
   */
  template<typename t_interaction_spec, typename t_perturbation_spec>
  class Nonbonded_Interaction : 
    public Interaction,
    public Nonbonded_Parameter
  {
  public:    

    typedef math::Periodicity<t_interaction_spec::boundary_type> Periodicity_type;

    /**
     * Constructor.
     */
    Nonbonded_Interaction(Pairlist_Algorithm<t_interaction_spec, 
			  t_perturbation_spec> *pa);
    
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
    Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec> & pairlist_algorithm()
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

    /**
     * accessor to the longrange timing
     */
    double & longrange_timing() { return m_longrange_timing; }
    
  protected:
    /**
     * the pairlist update algorithm.
     */
    Pairlist_Algorithm<t_interaction_spec, t_perturbation_spec> *m_pairlist_algorithm;
    /**
     * a vector of nonbonded sets
     */
    std::vector<Nonbonded_Set<t_interaction_spec, t_perturbation_spec> > m_nonbonded_set;
    /**
     * longrange timing.
     */
    double m_longrange_timing;

    /**
     * number of omp threads.
     */
    int m_omp_num_threads;

  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
