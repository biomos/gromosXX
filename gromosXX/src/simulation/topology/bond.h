/**
 * @file bond.h
 * the bond topology class.
 */

#ifndef INCLUDED_BOND_H
#define INCLUDED_BOND_H

namespace simulation
{
  /**
   * @class bond
   * holds bond information.
   */
  class bond
  {
  public:
    /**
     * @struct bond_struct
     * bond information.
     */
    struct bond_struct
    {
      int i;
      int j;
      int type;
    };

    /**
     * @class iterator
     * iterator over the bonds.
     */
    class iterator
    {
    public:
      iterator(std::vector<bond_struct> &bi);
      bool eol();
      bool neol();
      void operator++();
      int i();
      int j();
      int type();
    protected:
      std::vector<bond_struct>::const_iterator m_bond_it;
      std::vector<bond_struct>::const_iterator m_bond_end;
    };
      
    iterator begin();
    void add(int i, int j, int type);

  private:      
    std::vector<bond_struct> m_bond_information;
      
  };
	  
  
} // simulation

#include "bond.tcc"

#endif
