/**
 * @file chargegroup_range_pairlist_algorithm.h
 * a chargegroup based range pairlist algorithm
 */

#ifndef INCLUDED_CHARGEGROUP_RANGE_PAIRLIST_ALGORITHM_H
#define INCLUDED_CHARGEGROUP_RANGE_PAIRLIST_ALGORITHM_H

namespace interaction
{
  /**
   * @class Chargegroup_Range_Pairlist_Algorithm
   * the chargegroup based range pairlist algorithm
   */
  template<typename t_simulation,
    typename t_filter>
    class Chargegroup_Range_Pairlist_Algorithm
    : public Basic_Pairlist_Algorithm<t_simulation, t_filter>
    {
    public:
      /**
       * Constructor.
       */
      Chargegroup_Range_Pairlist_Algorithm(std::vector<std::vector<unsigned int> > 
					   &pairlist, Nonbonded_Base &base);
    
      /**
       * update the pairlist.
       */
      void update(t_simulation &sim);
    
    protected:

      void do_cg_interaction(simulation::chargegroup_iterator cg1,
			     simulation::chargegroup_iterator cg2);
    
      void do_cg_interaction_excl(t_simulation &sim,
				  simulation::chargegroup_iterator cg1,
				  simulation::chargegroup_iterator cg2);

      void do_cg_interaction_intra(t_simulation &sim,
				   simulation::chargegroup_iterator cg1);

    };
  
} // interaction

#include "chargegroup_range_pairlist_algorithm.tcc"

#endif
