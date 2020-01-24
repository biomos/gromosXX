/**
 * @file molecule_iterator.h
 * iterator over the molecules.
 */

#ifndef INCLUDED_MOLECULE_ITERATOR
#define INCLUDED_MOLECULE_ITERATOR

namespace topology
{
  /**
   * @class Molecule_Iterator
   * iterator over the molecules (SUBMOLECULES block)
   */
  class Molecule_Iterator 
    : public Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>
  {
  public:
    Molecule_Iterator(std::vector<unsigned int>::const_iterator mol_it)
      : Atomgroup_Iterator<std::vector<unsigned int>::const_iterator>(mol_it) {};
    

  private:
  };
  
} // topology

#endif

  
