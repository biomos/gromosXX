/**
 * pairlist_algorithm.tcc
 * base class for the pairlist(s) generating
 * algorithms.
 */

#ifndef INCLUDED_PAIRLIST_ALGORITHM_H
#define INCLUDED_PAIRLIST_ALGORITHM_H

namespace interaction
{
  /**
   * @class Pairlist_Algorithm
   * creates a pairlist.
   */
  template<typename t_simulation, typename t_nonbonded_spec>
  class Pairlist_Algorithm
  {
  public:
    /**
     * Constructor.
     */
    Pairlist_Algorithm() {}
    
    /**
     * update the pairlist(s).
     */
    template<typename t_nonbonded_interaction>
    void update(t_simulation const &sim, t_nonbonded_interaction &nonbonded_interaction){};
    
  protected:
  };
} // interaction

#endif
