/**
 * @file basic_pairlist_algorithm.h
 * the basic pairlist algorithm.
 */

#ifndef INCLUDED_BASIC_PAIRLIST_ALGORITHM_H
#define INCLUDED_BASIC_PAIRLIST_ALGORITHM_H

namespace interaction
{
  typedef std::vector<std::vector<size_t> > basic_pairlist_type;

  /**
   * @class Basic_Pairlist_Algorithm
   * the basic pairlist building algorithm.
   */
  template<typename t_simulation,
    typename t_filter = Basic_Filter<t_simulation, interaction::Storage> >
    class Basic_Pairlist_Algorithm
    {
      public:
      /**
       * Constructor.
       */
      Basic_Pairlist_Algorithm(basic_pairlist_type &pairlist,
			       basic_pairlist_type &perturbed_pl,
			       Nonbonded_Base &base);
      /**
       * update the pairlist.
       */
      void update(t_simulation &sim);
      
      /**
       * filter accessor.
       */
      t_filter & filter();
      
    protected:
      basic_pairlist_type & m_pairlist;
      basic_pairlist_type & m_perturbed_pairlist;
      t_filter m_filter;
      
    };
  
} // interaction

#include "basic_pairlist_algorithm.tcc"

#endif
