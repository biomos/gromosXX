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
    
    void com(math::VArray &pos, math::VArray &vel, math::SArray &mass,
	     math::Vec &com_pos, math::Vec &com_e_kin)
    {
      com_pos = 0.0;
      com_e_kin = 0.0;
      double m = 0.0;
      
      for(Atom_Iterator it=begin(), to=end(); it!=to; ++it){
	m += mass(*it);
	com_pos += mass(*it) * pos(*it);
	com_e_kin += mass(*it) * vel(*it);
      }
      com_pos /= m;
      com_e_kin = 0.5 * com_e_kin * com_e_kin / m;
    }

  private:
  };
  
} // simulation

#endif

  
