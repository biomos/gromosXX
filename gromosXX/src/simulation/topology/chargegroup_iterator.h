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
  {
  public:
    class atom_iterator
    {
    public:
      atom_iterator(size_t atom)
	: m_atom(atom)
      {
      }
      bool operator==(atom_iterator &it)
      {
	return m_atom == it.m_atom;
      }
      bool operator!=(atom_iterator &it)
      {
	return m_atom != it.m_atom;
      }
      void operator++()
      {
	++m_atom;
      }
      size_t operator*()
      {
	return m_atom;
      }
    private:
      size_t m_atom;
    };

    chargegroup_iterator(std::vector<int>::const_iterator cg_it)
      : m_cg_it(cg_it)
    {
    }

    bool operator==(chargegroup_iterator &it)
    {
      return m_cg_it == it.m_cg_it;
    }
    bool operator!=(chargegroup_iterator &it)
    {
      return !(*this == it);
    }
    void operator++()
    {
      ++m_cg_it;
    }
    std::vector<int>::const_iterator operator*()
    {
      return m_cg_it;
    }
    atom_iterator begin()
    {
      return atom_iterator(*m_cg_it);
    }
    atom_iterator end()
    {
      return atom_iterator(*(m_cg_it+1));
    }
    int num_atoms()
    {
      return *(m_cg_it+1) - *m_cg_it;
    }
    void cog(math::VArray &pos, math::Vec &v)
    {
      v = 0.0;
      for(atom_iterator it=begin(), to=end(); it!=to; ++it)
	v += pos(*it);
      v /= num_atoms();
    }
      
  private:
    std::vector<int>::const_iterator m_cg_it;
  };
  
}

#endif
  
