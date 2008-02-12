/**
 * @file pairlist.h
 * the pairlist class.
 */

#ifndef INCLUDED_PAIRLIST_H
#define INCLUDED_PAIRLIST_H

namespace interaction
{
  /**
   * @class Pairlist
   * holds a Pairlist.
   * very easy implementation that just uses standard vectors.
   */
  class Pairlist :
    public std::vector< std::vector<unsigned int> >
  {
  public:
    /**
     * Constructor.
     */
    Pairlist() {}

  protected:
    
  };
  
  /**
   * sort and print the pairlist.
   */
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl);
  
} // interaction

#endif
