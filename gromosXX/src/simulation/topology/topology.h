/**
 * @file topology.h
 * the topology class
 */

#ifndef INCLUDED_TOPOLOGY_H
#define INCLUDED_TOPOLOGY_H

namespace simulation
{
  /**
   * @class topology
   * holds the topological information of
   * the simulated system
   * @sa simulation::simulation
   * @sa simulation::system
   */
  class topology
  {
  public:
    /**
     * Constructor
     */
    explicit topology();

    /**
     * masses accessor
     */
    math::SArray &mass();

    /**
     * masses const accessor
     */
    math::SArray const & mass()const;

    /**
     * @class bond
     * holds bond information.
     */
    class bond
    {
    private:
      /**
       * @struct bond
       * bond information.
       */
      struct bond_struct
      {
	int i;
	int j;
	int type;
      };
    public:
      /**
       * @class iterator
       * iterator over the bonds.
       */
      class iterator
      {
      public:
	iterator(std::vector<bond_struct> &bi);
	bool eol();
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
	
    /**
     * bond accessor.
     */
    bond & bonds();

    /**
     * set the number of solute atoms
     */
    void num_solute_atoms(size_t atoms);
    
    /**
     * get the number of solute atoms
     */
    size_t num_solute_atoms()const;
    
  private:
    /**
     * the number of solute atoms.
     */
    size_t m_solute_atoms;
    
    /**
     * the atom masses.
     */
    math::SArray m_mass;
    
    /**
     * the bonds.
     */
    bond m_bonds;
    
  }; // topology
  
} // simulation

// inline method definitions
#include "topology.tcc"

#endif
