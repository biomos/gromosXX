/**
 * @file temperaturegroup_iterator.h
 * iterator over the temperature groups.
 */

#ifndef INCLUDED_TEMPERATUREGROUP_ITERATOR
#define INCLUDED_TEMPERATUREGROUP_ITERATOR

namespace topology
{
  /**
   * @class Temperaturegroup_Iterator
   * iterator over the temperature groups (TEMPERATUREGROUPS block)
   */
  class Temperaturegroup_Iterator 
    : public Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>
  {
  public:
    Temperaturegroup_Iterator(std::vector<unsigned int>::const_iterator mol_it)
      : Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>(mol_it) {};
  };
  
} // topology

#endif

  
