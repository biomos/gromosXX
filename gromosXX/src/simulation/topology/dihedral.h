/**
 * @file dihedral.h
 * the dihedral topology class.
 */

#ifndef INCLUDED_DIHEDRAL_H
#define INCLUDED_DIHEDRAL_H

namespace simulation
{
  /**
   * @class Dihedral
   * holds dihedral information.
   */
  class Dihedral
  {
  public:
    /**
     * @struct dihedral_struct
     * dihedral information.
     */
    struct dihedral_struct
    {
      int i;
      int j;
      int k;
      int l;
      int type;
    };

    /**
     * @class iterator
     * iterator over the dihedrals.
     */
    class iterator
    {
    public:
      iterator(std::vector<dihedral_struct> &bi);
      bool eol();
      bool neol();
      void operator++();
      int i();
      int j();
      int k();
      int l();
      int type();
    protected:
      std::vector<dihedral_struct>::const_iterator 
	m_dihedral_it;
      std::vector<dihedral_struct>::const_iterator 
	m_dihedral_end;
    };
      
    iterator begin();
    void add(int i, int j, int k, int l, int type);

  private:      
    std::vector<dihedral_struct> m_dihedral_information;
      
  };
	  
  
} // simulation

#include "dihedral.tcc"

#endif
