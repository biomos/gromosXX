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
     * integer atom code accessor.
     */
    int iac(int i);
    
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
     * solvent accessor.
     * support for multiple solvents.
     */
    solvent & solvents(size_t i);

    /**
     * add a solvent.
     */
    void add_solvent(solvent solv);

    /**
     * add solvent to the simulation.
     * @param solv the solvent (multiple solvents).
     * @param num_molecules the number of solvent molecules to add.
     */
    void solvate(int solv, int num_molecules);
    
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
    
    /**
     * residue names.
     */
    std::vector<std::string> & residue_name();

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
     * the solvents (multiple solvent).
     */
    std::vector<solvent> m_solvents;
    
    /**
     * the integer atom code.
     */
    std::vector<int> m_iac;

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
    
    /**
     * residue names.
     */
    std::vector<std::string> m_residue_name;
    
  }; // topology
  
} // simulation

// inline method definitions
#include "topology.tcc"

#endif
