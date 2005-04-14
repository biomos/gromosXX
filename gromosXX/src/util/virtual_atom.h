/**
 * @file virtual_atom.h
 * 
 */

#ifndef INCLUDED_VIRTUAL_ATOM_H
#define INCLUDED_VIRTUAL_ATOM_H

namespace configuration
{
  class Configuration;
}

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

    math::Vec pos(configuration::Configuration & conf)const;
    void force(configuration::Configuration & conf, math::Vec const f)const;
    int atom(int i)const { return m_atom[i]; }

  private:

    template<math::boundary_enum B>
    void _pos(math::VArray const & position, math::Box const & box, 
	      math::Vec & p)const;
    template<math::boundary_enum B>
    void _force(math::VArray const & position, math::Box const & box, 
		math::Vec const & f, math::VArray & force)const;

    int m_type;
    std::vector<int> m_atom;
    double m_dish, m_disc;
    int m_orientation;
      
  };
  
  
}

#endif
