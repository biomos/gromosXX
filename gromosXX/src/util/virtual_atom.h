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
  enum virtual_type { va_explicit = 0,
		      va_CH1 = 1,
		      va_aromatic = 2,
		      va_CH2 = 3,
		      va_stereo_CH2 = 4,
		      va_stereo_CH3 = 5,
		      va_CH3 = 6,
		      va_ring = 7,
		      va_NH2 = 8,
		      va_3CH3 = 9,
		      va_cog = 10,
		      va_com = 11
  };
  
  /**
   * @class Virtual_Atom
   * get position of a virtual or pseudo atom based on 
   * positions of specified other atoms.
   * distribute forces on virtual atom over the atoms used
   * to calculate the virtual site.
   * 
   */

  class Virtual_Atom
  {
    public:
    /**
     * Constructor
     */
    Virtual_Atom();
    /**
     * Constructor
     */
    Virtual_Atom(virtual_type type, std::vector<int> atom, double dish = 0.1, 
		 double disc = 0.153,int orientation = 0);
    /**
     * calculate the position of the virtual atom
     */
    math::Vec pos(configuration::Configuration & conf)const;
    /**
     * distribute force f of virtual atom on the real atoms.
     */
    void force(configuration::Configuration & conf, math::Vec const f)const;
    /**
     * real atom accessor
     */
    int atom(int i)const { assert(i >= 0 && i < int(m_atom.size())); return m_atom[i]; } 
    /**
     * number of atoms that define virtual atom
     */
    int size()const { return m_atom.size(); }
    
  private:
    /**
     * calculate the position of the virtual site
     */
    template<math::boundary_enum B>
    void _pos(math::VArray const & position, math::Box const & box, 
	      math::Vec & p)const;
    /**
     * distribute the force of the virtual site on the real atoms
     */
    template<math::boundary_enum B>
    void _force(math::VArray const & position, math::Box const & box, 
		math::Vec const & f, math::VArray & force)const;
    /**
     * type of the virtual atom
     */
    virtual_type m_type;
    /**
     * atoms to specify the virtual site
     */
    std::vector<int> m_atom;
    /**
     * C-H bond length
     */
    double m_dish;
    /**
     * C-C bond length
     */
    double m_disc;
    /**
     * orientation
     */
    int m_orientation;
  };
}

#endif
