/**
 * @file virtual_atom.h
 * 
 */

#ifndef INCLUDED_VIRTUAL_ATOM_H
#define INCLUDED_VIRTUAL_ATOM_H

namespace util
{
  /**
   * @class Virtual_Atom
   * 
   */

  class Virtual_Atom
  {
    public:
    Virtual_Atom();
    
    Virtual_Atom(int type, std::vector<int> atom, double dish = 0.1, 
		 double disc = 0.153,int orientation = 0);

    math::Vec pos(math::VArray const & position)const;

    void force(math::VArray const & position, math::Vec const f, math::VArray & force) const;
  
    int atom(int i)const { return m_atom[i]; }

  private:
    int m_type;
    std::vector<int> m_atom;
    double m_dish, m_disc;
    int m_orientation;
      
  };
  
  
}

#endif
