/**
 * @file molecule_iterator.h
 * iterator over the molecules.
 */

#ifndef INCLUDED_MOLECULE_ITERATOR
#define INCLUDED_MOLECULE_ITERATOR

namespace simulation
{
  /**
   * @class Molecule_Iterator
   * iterator over the molecules (SUBMOLECULES block)
   */
  class Molecule_Iterator 
    : public Atom_Group_Iterator<std::vector<size_t>::const_iterator>
  {
  public:
    Molecule_Iterator(std::vector<size_t>::const_iterator mol_it)
      : Atom_Group_Iterator<std::vector<size_t>::const_iterator>(mol_it) {};
    

  private:
  };
  
} // simulation

#endif

  
