/**
 * @file basic_pairlist_algorithm.h
 * the basic pairlist algorithm.
 */

#ifndef INCLUDED_BASIC_PAIRLIST_ALGORITHM_H
#define INCLUDED_BASIC_PAIRLIST_ALGORITHM_H

namespace interaction
{
  /**
   * @class Basic_Pairlist_Algorithm
   * the basic pairlist building algorithm.
   */
  template<typename t_simulation,
    typename t_filter = Basic_Filter<t_simulation> >
    class Basic_Pairlist_Algorithm
    {
      public:
      /**
       * Constructor.
       */
      Basic_Pairlist_Algorithm(std::vector<std::vector<unsigned int> > &pairlist,
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
      std::vector<std::vector<unsigned int> > & m_pairlist;
      t_filter m_filter;
      
    };
  
} // interaction

#include "basic_pairlist_algorithm.tcc"

#endif
