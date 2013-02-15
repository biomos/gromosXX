/**
 * @file pressuregroup_iterator.h
 * iterator over the pressure groups.
 */

#ifndef INCLUDED_PRESSUREGROUP_ITERATOR
#define INCLUDED_PRESSUREGROUP_ITERATOR

namespace topology
{
  /**
   * @class Pressuregroup_Iterator
   * iterator over the pressure groups (PRESSUREGROUPS block)
   */
  class Pressuregroup_Iterator 
    : public Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>
  {
  public:
    Pressuregroup_Iterator(std::vector<unsigned int>::const_iterator mol_it)
      : Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>(mol_it) {};
  };
  
} // topology

#endif

  
