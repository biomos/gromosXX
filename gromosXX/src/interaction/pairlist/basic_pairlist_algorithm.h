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
  class Basic_Pairlist_Algorithm :
    public t_filter
  {
  public:
    /**
     * Constructor.
     */
    Basic_Pairlist_Algorithm(std::vector<std::vector<unsigned int> > &pairlist);
    /**
     * update the pairlist.
     */
    void update(t_simulation &sim);

  protected:
    std::vector<std::vector<unsigned int> > & m_pairlist;
    
  };
  
} // interaction

#include "basic_pairlist_algorithm.tcc"

#endif
