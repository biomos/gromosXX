/**
 * @file chargegroup_iterator.h
 * iterator over the chargegroups.
 */

#ifndef INCLUDED_CHARGEGROUP_ITERATOR_H
#define INCLUDED_CHARGEGROUP_ITERATOR_H

#undef MODULE
#undef SUBMODULE
#define MODULE topology
#define SUBMODULE topology


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

    Chargegroup_Iterator()
    {
    }

    void cog(const math::VArray &pos, math::Vec &v)const
    {
      DEBUG(10, "cog: " << *begin() << " - " << *end() << " of " << pos.size());
      v = 0.0;
      for(Atom_Iterator it=begin(), to=end(); it!=to; ++it){
	assert(pos.size() > *it);
	v += pos(int(*it));
      }
      
      v /= num_atoms();
    }
      
  private:

  };
  
}

#endif
  
