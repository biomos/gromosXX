/**
 * @file chargegroup_iterator.h
 * iterator over the chargegroups.
 */

#ifndef INCLUDED_CHARGEGROUP_ITERATOR_H
#define INCLUDED_CHARGEGROUP_ITERATOR_H

namespace topology
{
  /**
   * @class Chargegroup_Iterator
   * iterates over the chargegroups.
   */
  class Chargegroup_Iterator
    : public Atomgroup_Iterator<std::vector<int>::const_iterator>
  {
  public:
    Chargegroup_Iterator(std::vector<int>::const_iterator cg_it)
      : Atomgroup_Iterator<std::vector<int>::const_iterator>(cg_it)
    {
    }

    void cog(const math::VArray &pos, math::Vec &v)const
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
  
