/**
 * @file chargegroup_iterator.h
 * iterator over the chargegroups.
 */

#ifndef INCLUDED_CHARGEGROUP_ITERATOR_H
#define INCLUDED_CHARGEGROUP_ITERATOR_H

namespace simulation
{
  /**
   * @class chargegroup_iterator
   * iterates over the chargegroups.
   */
  class chargegroup_iterator
    : public Atom_Group_Iterator<std::vector<int>::const_iterator>
  {
  public:
    chargegroup_iterator(std::vector<int>::const_iterator cg_it)
      : Atom_Group_Iterator<std::vector<int>::const_iterator>(cg_it)
    {
    }

    void cog(math::VArray &pos, math::Vec &v)
    {
      v = 0.0;
      for(Atom_Iterator it=begin(), to=end(); it!=to; ++it)
	v += pos(*it);
      v /= num_atoms();
    }
      
  private:

  };
  
}

#endif
  
