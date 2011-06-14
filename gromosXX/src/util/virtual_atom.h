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

namespace topology
{
  class Topology;
}

namespace util
{
/**
   * @enum virtual_type
   * holds the virtual atom type
   */
  enum virtual_type {
    /**
     * 0: explicit atom
     */
    va_explicit = 0,
    /**
     * 1: aliphatic CH1
     */
    va_CH1 = 1,
    /**
     * 2: CH1 (aromatic)
     */
    va_aromatic = 2,
    /**
     * 3: CH2 (non-stereospecific)
     */
    va_CH2 = 3,
    /**
     * 4: CH2 (stereospecific)
     */
    va_stereo_CH2 = 4,
    /**
     * 5: CH3
     */
    va_stereo_CH3 = 5,
    /**
     * 6: CH3 (non-stereospecific, Val, Leu)
     */
    va_CH3 = 6,
    /**
     * 7: 3CH3 (non-stereospecific)
     */
    va_3CH3 = 7,
    /**
     * -1: centre of geometry
     */
    va_cog = -1,
    /**
     * -2: centre of mass
     */
    va_com = -2
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
		 double disc = 0.153,  int orientation = 0);
    /**
     * calculate the position of the virtual atom
     */
    math::Vec pos(configuration::Configuration & conf, topology::Topology & topo)const;
    /**
     * distribute force f of virtual atom on the real atoms.
     */
    void force(configuration::Configuration & conf, topology::Topology & topo, math::Vec const f)const;
    /**
     * distribute force f of virtual atom on the real atoms (for eds).
     */
    void force(configuration::Configuration & conf, topology::Topology  & topo, math::Vec const f, math::VArray & force)const;
    /**
     * real atom accessor
     */
    int atom(int i)const { assert(i >= 0 && i < int(m_atom.size())); return m_atom[i]; } 
    /**
     * number of atoms that define virtual atom
     */
    int size()const { return m_atom.size(); }
    /**
     * accessor to the type
     */
    virtual_type type() const { return m_type; }
    
  private:
    /**
     * calculate the position of the virtual site
     */
    template<math::boundary_enum B>
    void _pos(math::VArray const & position, topology::Topology const & topo, math::Box const & box, 
	      math::Vec & p)const;
    /**
     * distribute the force of the virtual site on the real atoms
     */
    template<math::boundary_enum B>
    void _force(math::VArray const & position, topology::Topology const & topo,  math::Box const & box, 
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
