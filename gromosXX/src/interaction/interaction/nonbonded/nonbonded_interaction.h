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
  template<typename t_simulation, typename t_nonbonded_spec>
  class Nonbonded_Interaction : 
    public Interaction<t_simulation>,
    public Nonbonded_Base,
    public Storage,
    public t_nonbonded_spec::nonbonded_innerloop_type
  {
  public:    
    /**
     * Constructor.
     * @param sim where to store forces and energies
     * (and virial contribution).
     */
    Nonbonded_Interaction(t_simulation &sim);
    
    /**
     * Destructor.
     */
    virtual ~Nonbonded_Interaction();

    /**
     * calculate the interactions.
     */
    virtual void calculate_interactions(t_simulation &sim);

    /**
     * add a shortrange interaction.
     */
    void add_shortrange_pair(t_simulation const &sim,
			     size_t const i, size_t const j);
    /**
     * add a longrange interaction.
     */
    void add_longrange_pair(t_simulation & sim,
			    size_t const i, size_t const j);

    // ACCESSORS
    /**
     * pairlist accessor.
     */
    Pairlist & pairlist();
    /**
     * const pairlist accessor.
     */
    Pairlist const & pairlist()const;
    /**
     * perturbed pairlist accessor.
     */
    Pairlist & perturbed_pairlist();
    /**
     * const perturbed pairlist accessor.
     */
    Pairlist const & perturbed_pairlist()const;

  protected:
    /**
     * size the arrays of storage.
     */
    void initialize(t_simulation &sim);
    
    /**
     * calculate the interactions.
     */
    void do_interactions(t_simulation &sim,
			 Pairlist::iterator it, 
			 Pairlist::iterator to);

    /**
     * calculate the 1,4-interactions.
     */
    void do_14_interactions(t_simulation &sim);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    void do_RF_excluded_interactions(t_simulation &sim);
 
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
    typename t_nonbonded_spec::pairlist_algorithm_type m_pairlist_algorithm;
    
  };
  
} // interaction

// template methods
#include "nonbonded_interaction.tcc"

#endif
