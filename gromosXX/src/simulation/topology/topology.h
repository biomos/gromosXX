/**
 * @file topology.h
 * the topology class
 */

#ifndef INCLUDED_TOPOLOGY_H
#define INCLUDED_TOPOLOGY_H

namespace simulation
{
  /**
   * @class Topology
   * holds the topological information of
   * the simulated system
   * @sa simulation::simulation
   * @sa simulation::system
   */
  class Topology
  {
  public:
    /**
     * Constructor
     */
    explicit Topology();

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
     * soluteatom accessor.
     */
    Solute & solute();

    /**
     * const solute accessor.
     */
    Solute const & solute()const;
    
    /**
     * number of solvents.
     */
    size_t num_solvents()const;
    
    /**
     * solvent accessor.
     * support for multiple solvents.
     */
    Solvent & solvent(size_t i);

    /**
     * add a solvent.
     */
    void add_solvent(Solvent solv);

    /**
     * add solvent to the simulation.
     * @param solv the solvent (multiple solvents).
     * @param num_molecules the number of solvent molecules to add.
     */
    void solvate(size_t solv, size_t num_molecules);
    
    /**
     * set the capacity of solute atoms
     */
    void resize(size_t atoms);

    /**
     * get the number of solute atoms
     */
    size_t num_solute_atoms()const;

    /**
     * get the total number of solvent atoms.
     */
    size_t num_solvent_atoms()const;

    /**
     * get the number of solvent molecules.
     */
    size_t num_solvent_molecules(size_t i)const;
    
    /**
     * get the number of solvent atoms.
     */
    size_t num_solvent_atoms(size_t i)const;

    /**
     * add a solute atom.
     */
    void add_solute_atom(std::string name, int residue_nr, int iac,
			 double mass, double charge, bool chargegroup,
			 std::set<int> exclusions,
			 std::set<int> one_four_pairs);
    
    /**
     * residue names.
     */
    std::vector<std::string> & residue_name();

    /**
     * all exclusions for atom i.
     */
    std::set<int> & all_exclusion(size_t const i);
    /**
     * 1,4 pairs of atom i.
     */
    std::set<int> & one_four_pair(size_t const i);
    /**
     * the number of chargegroups present.
     */
    size_t num_chargegroups();
    /**
     * the number of solute chargegroups.
     */
    size_t num_solute_chargegroups();
    /**
     * iterator over the chargegrops
     */
    chargegroup_iterator chargegroup_begin();
    /**
     * end of the chargegroup iterator.
     */
    chargegroup_iterator chargegroup_end();

  private:
    /**
     * the soluteatoms.
     */
    Solute m_solute;

    /**
     * the number of solvent molecules.
     */
    std::vector<size_t> m_num_solvent_molecules;
    
    /**
     * the number of solvent atoms.
     * vector for multiple solvents.
     */
    std::vector<size_t> m_num_solvent_atoms;
    
    /**
     * the solvents (multiple solvent).
     */
    std::vector<Solvent> m_solvent;
    
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
     * the number of solute chargegroups.
     */
    size_t m_num_solute_chargegroups;
    
    /**
     * residue names (solute and solvent).
     */
    std::vector<std::string> m_residue_name;
    
  }; // topology
  
} // simulation

// inline method definitions
#include "topology.tcc"

#endif
