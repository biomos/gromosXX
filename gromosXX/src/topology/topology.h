/**
 * @file topology.h
 * the topology class
 */

#ifndef INCLUDED_TOPOLOGY_H
#define INCLUDED_TOPOLOGY_H

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

namespace simulation
{
  class Multibath;
}

namespace topology
{
  /**
   * @class Topology
   * holds the topological information of
   * the simulated system
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
    int iac(size_t const i)const {assert(i < m_iac.size()); return m_iac[i];}
    
    /**
     * masses accessor
     */
    math::SArray &mass() {return m_mass;}

    /**
     * masses const accessor
     */
    math::SArray const & mass()const {return m_mass;}

    /**
     * charge accessor
     */
    math::SArray &charge() {return m_charge;}
    
    /**
     * charge const accessor
     */
    math::SArray const & charge()const{return m_charge;}
    
    /**
     * solute accessor.
     */
    Solute & solute() {return m_solute;}

    /**
     * perturbed solute accessor.
     */
    Perturbed_Solute & perturbed_solute() {return m_perturbed_solute;}

    /**
     * const solute accessor.
     */
    Solute const & solute()const{return m_solute;}

    /**
     * const perturbed solute accessor.
     */
    Perturbed_Solute const & perturbed_solute()const{return m_perturbed_solute;}
    
    /**
     * number of solvents.
     */
    size_t num_solvents()const {return m_num_solvent_molecules.size();}
    
    /**
     * solvent accessor.
     * support for multiple solvents.
     */
    Solvent & solvent(size_t i) {assert(i < m_solvent.size()); return m_solvent[i];}
    /**
     * const solvent accessor.
     * support for multiple solvents.
     */
    Solvent const & solvent(size_t i)const {assert(i < m_solvent.size()); return m_solvent[i];}

    /**
     * add a solvent.
     */
    void add_solvent(Solvent solv) {m_solvent.push_back(solv);}

    /**
     * add solvent to the simulation.
     * @param solv the solvent (multiple solvents).
     * @param num_molecules the number of solvent molecules to add.
     */
    void solvate(size_t solv, size_t num_molecules);
    
    /**
     * set the capacity of solute atoms
     */
    void resize(size_t const atoms);

    /**
     * get the total number of atoms.
     */
    size_t num_atoms()const {return num_solute_atoms() + num_solvent_atoms();}

    /**
     * get the number of solute atoms
     */
    size_t num_solute_atoms()const{return solute().num_atoms();}

    /**
     * get the total number of solvent atoms.
     */
    size_t num_solvent_atoms()const;

    /**
     * get the number of solvent molecules.
     */
    size_t num_solvent_molecules(size_t i)const {
      assert(i < m_num_solvent_molecules.size());
      return m_num_solvent_molecules[i];
    }
    
    /**
     * get the number of solvent atoms.
     */
    size_t num_solvent_atoms(size_t i)const{
      assert(i<m_num_solvent_atoms.size());
      return m_num_solvent_atoms[i];
    }

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
    std::vector<std::string> & residue_names() {return m_residue_name;}
    /**
     * const residue names.
     */
    std::vector<std::string> const & residue_names()const {return m_residue_name;}

    /**
     * all exclusions for atom i. Exclusions and 1,4 interactions.
     */
    std::set<int> & all_exclusion(size_t const i){
      assert(i < m_all_exclusion.size());
      return m_all_exclusion[i];
    }

    /**
     * const all exclusions for atom i. Exclusions and 1,4 interactions.
     */
    std::set<int> const & all_exclusion(size_t const i)const{
      assert(i < m_all_exclusion.size());
      return m_all_exclusion[i];
    }

    /**
     * exclusions for atom i.
     */
    std::set<int> & exclusion(size_t const i){
      assert(i < m_exclusion.size());
      return m_exclusion[i];
    }
    /**
     * exclusions
     */
    std::vector<std::set<int> > & exclusion() {return m_exclusion;}
    
    /**
     * 1,4 pairs of atom i.
     */
    std::set<int> & one_four_pair(size_t const i){
      assert(i < m_one_four_pair.size());
      return m_one_four_pair[i];
    }
    /**
     * 1,4 pairs 
     */
    std::vector<std::set<int> > & one_four_pair(){return m_one_four_pair;}
    
    /**
     * the number of chargegroups present.
     */
    size_t num_chargegroups()const {return m_chargegroup.size()-1;}

    /**
     * the number of solute chargegroups.
     */
    size_t num_solute_chargegroups()const {return m_num_solute_chargegroups;}
    /**
     * iterator over the chargegrops
     */
    Chargegroup_Iterator chargegroup_begin()const {
      return Chargegroup_Iterator(m_chargegroup.begin());}

    /**
     * end of the chargegroup iterator.
     */
    Chargegroup_Iterator chargegroup_end()const{
      return Chargegroup_Iterator(m_chargegroup.end()-1);}

    /**
     * the molecules.
     */
    std::vector<size_t> & molecules(){ return m_molecule;}
    /**
     * iterator over the molecules.
     */
    Molecule_Iterator molecule_begin(){
      return Molecule_Iterator(m_molecule.begin());};
    /**
     * end of molecule iterator.
     */
    Molecule_Iterator molecule_end(){
      return Molecule_Iterator(m_molecule.end());}

    /**
     * const energy group accessor.
     */
    std::vector<size_t> const & energy_groups()const {return m_energy_group;}

    /**
     * energy group accessor.
     */
    std::vector<size_t> & energy_groups() {return m_energy_group;}

    /**
     * const energy group of atoms accessor.
     */
    std::vector<size_t> const & atom_energy_group()const{
      return m_atom_energy_group;
    }

    /**
     * energy group of atoms accessor.
     */
    std::vector<size_t> & atom_energy_group(){
      return m_atom_energy_group;
    }
  
    /**
     * energy group of atom accessor
     */
    const size_t atom_energy_group(size_t i)const{
      assert(i < m_atom_energy_group.size());
      return m_atom_energy_group[i];
    };

    /**
     * calculate constraint degrees of freedom
     */
    void calculate_constraint_dof(simulation::Multibath & mb)const;
    
    /**
     * check state
     */
    int check_state()const;

    /**
     * return lambda
     */
    const double lambda()const {return m_lambda;}
    /**
     * set lambda
     */

    void lambda(double const l){m_lambda = l;}
    
    /**
     * return nlam (exponent to lambda)
     */
    const int lambda_exp()const {return m_lambda_exp;}

    /**
     * set lambda exponent.
     */
    void lambda_exp(int const l) { m_lambda_exp = l; }

    /**
     * is the atom perturbed?
     */
    bool is_perturbed(size_t const i)const {
      assert(i < m_is_perturbed.size()); 
      return m_is_perturbed[i];
    }
    /**
     * is the atom perturbed?
     */
    std::vector<bool> & is_perturbed() { return m_is_perturbed;}
    /**
     * recalculate lambda dependent properties.
     */
    void update_for_lambda();

    /**
     * interaction matrix for scaled interactions with energy groups
     */
    std::map<std::pair<int, int>, std::pair<double, double> > const &
    energy_group_scaling()const { return m_energy_group_scaling;}

    /**
     * interaction matrix for scaled interactions with energy groups
     */
    std::map<std::pair<int, int>, std::pair<double, double> > &
    energy_group_scaling() { return m_energy_group_scaling;}
    
  private:
    /**
     * the solute.
     */
    Solute m_solute;

    /**
     * the perturbed solute.
     */
    Perturbed_Solute m_perturbed_solute;

    /**
     * is the atom perturbed?
     */
    std::vector<bool> m_is_perturbed;
    
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
     * the molecules.
     */
    std::vector<size_t> m_molecule;

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

    /**
     * energy groups.
     */
    std::vector<size_t> m_energy_group;
    
    /**
     * energy group of atom
     */
    std::vector<size_t> m_atom_energy_group;
    
    /**
     * lambda.
     */
    double m_lambda;
    
    /**
     * lambda exponent
     */
    int m_lambda_exp;

    /**
     * interaction matrix for scaled interactions
     */
    std::map<std::pair<int, int>, std::pair<double, double> > 
    m_energy_group_scaling;
    
  }; // topology
  
} // topology

#endif
