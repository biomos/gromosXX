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
     * charge accessor
     */
    math::SArray &charge();
    
    /**
     * charge const accessor
     */
    math::SArray const & charge()const;
    
    /**
     * bond accessor.
     */
    bond & bonds();

    /**
     * soluteatom accessor.
     */
    soluteatom & soluteatoms();

    /**
     * set the capacity of solute atoms
     */
    void solute_atoms_capacity(size_t atoms);

    /**
     * set the number of solute atoms
     */
    void num_solute_atoms(size_t atoms);
    
    /**
     * get the number of solute atoms
     */
    size_t num_solute_atoms()const;

    void add_solute_atom(std::string name, int residue_nr, int iac,
			 double mass, double charge, bool chargegroup,
			 std::set<int> exclusions,
			 std::set<int> one_four_pairs);
    
  private:
    /**
     * the number of solute atoms.
     */
    size_t m_num_solute_atoms;
    
    /**
     * the soluteatoms.
     */
    soluteatom m_soluteatoms;

    /**
     * the atom masses.
     */
    math::SArray m_mass;
    
    /**
     * the atom charges.
     */
     math::SArray m_charge;

    /**
     * the atom exclusions.
     */
    std::vector< std::set<int> > m_exclusion;
    
    /**
     * the atom 1-4 interactions.
     */
    std::vector< std::set<int> > m_one_four_pair;
    
    /**
     * atom exclusions and 1-4 interactions.
     */
    std::vector< std::set<int> > m_all_exclusion;
    
    /**
     * the chargegroups.
     */
    std::vector<int> m_chargegroup;
    
    /**
     * the bonds.
     */
    bond m_bonds;
    
  }; // topology
  
} // simulation

// inline method definitions
#include "topology.tcc"

#endif
