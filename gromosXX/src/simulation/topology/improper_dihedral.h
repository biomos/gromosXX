/**
 * @file improper_dihedral.h
 * the improper dihedral topology class.
 */

#ifndef INCLUDED_IMPROPER_DIHEDRAL_H
#define INCLUDED_IMPROPER_DIHEDRAL_H

namespace simulation
{
  /**
   * @class Angle
   * holds angle information.
   */
  class Improper_dihedral
  {
  public:
    /**
     * @struct angle_struct
     * angle information.
     */
    struct improper_dihedral_struct
    {
      int i;
      int j;
      int k;
      int l;
      int type;
    };

    /**
     * @class iterator
     * iterator over the improper dihedrals.
     */
    class iterator
    {
    public:
      iterator(std::vector<improper_dihedral_struct> &bi);
      bool eol();
      bool neol();
      void operator++();
      int i();
      int j();
      int k();
      int l();
      int type();
    protected:
      std::vector<improper_dihedral_struct>::const_iterator 
	m_improper_dihedral_it;
      std::vector<improper_dihedral_struct>::const_iterator 
	m_improper_dihedral_end;
    };
      
    iterator begin();
    void add(int i, int j, int k, int l, int type);

  private:      
    std::vector<improper_dihedral_struct> m_improper_dihedral_information;
      
  };
	  
  
} // simulation

#include "improper_dihedral.tcc"

#endif
