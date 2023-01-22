/**
 * @file topology.h
 * the topology class
 */

#ifndef INCLUDED_TOPOLOGY_H
#define INCLUDED_TOPOLOGY_H

#include "core/core.h"
#include "solute.h"
#include "solvent.h"
#include "perturbed_atom.h"
#include "perturbed_solute.h"
#include "sd.h"
#include "exclusions.h"
#include "../interaction/interaction_types.h"
#include "../util/virtual_atom.h"

namespace simulation
{
  class Multibath;
  class Simulation;
  class Parameter;
}

namespace util {
  class LE_Coordinate;
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
     * copy constructor
     * can multiply topologies
     */
    explicit Topology(Topology const & topo, int mul_solute = 1, int mul_solvent = -1);

    /**
     * Destructor
     */
    ~Topology();

    /**
     * integer atom code accessor.
     */
    int iac(unsigned int const i)const {assert(i < m_iac.size()); return m_iac[i];}

    /**
     * integer atom code accessor
     */
    std::vector<int> &iac() {return m_iac;}

    /**
     * integer atom code const accessor
     */
    std::vector<int> const & iac()const {return m_iac;}

    /**
     * masses accessor
     */
    math::SArray &mass() {return m_mass;}

    /**
     * masses const accessor
     */
    math::SArray const & mass()const {return m_mass;}

    /**
     * mass of atom i
     */
    double mass(int i)const { return m_mass(i); }
    /**
     * masses accessor
     */
    math::SArray &inverse_mass() {return m_inverse_mass;}

    /**
     * masses const accessor
     */
    math::SArray const & inverse_mass()const {return m_inverse_mass;}

    /**
     * mass of atom i
     */
    double inverse_mass(int i)const { return m_inverse_mass(i); }

    /**
     * charge accessor
     */
    math::SArray &charge() {return m_charge;}

    /**
     * charge const accessor
     */
    math::SArray const & charge()const{return m_charge;}

    /**
     * charge accessor
     */
    double charge(int i)const { return m_charge(i); }

    /**
     * calculates the square of the sum of the charges
     *
     * @f[S^2 = (\sum_{i=1}^{N_q}\,q_i)@f]
     */
    double squared_sum_charges() const;

    /**
     * calculates the sum of the squared charges
     *
     * @f[\tilde{S}^2 = \sum_{i=1}^{N_q}\,q_i^2@f]
     */
    double sum_squared_charges() const;

    /**
     * polarisation charge accessor
     */
    math::SArray &coscharge() {return m_coscharge;}

    /**
     * polarisation charge const accessor
     */
    double coscharge(int i)const { return m_coscharge(i); }

    /**
     * polarisability accessor
     */
    math::SArray &polarisability() {return m_polarisability;}

    /**
     * polarisability const accessor
     */
    double polarisability(int i)const {return m_polarisability(i);}

    /**
     * polarisability damping E0 accessor
     */
    math::SArray &damping_level() {return m_damping_level;}

    /**
     * polarisability damping E0 const accessor
     */
    double damping_level(int i)const {return m_damping_level(i);}

    /**
     * polarisability damping p accessor
     */
    math::SArray &damping_power() {return m_damping_power;}

    /**
     * polarisability damping p const accessor
     */
    double damping_power(int i)const {return m_damping_power(i);}
   /**
     * polarisability gamma accessor
     */
    math::SArray &gamma() {return m_gamma;}

    /**
     * polarisability gamma const accessor
     */
    double gamma(int i)const {return m_gamma(i);}

    /**
     * polarisability gamma atm j accessor
     */
    std::vector<int> &gamma_j() {return m_gamma_j;}

    /**
     * polarisability gamma atom j const accessor
     */
    int gamma_j(int i)const {return m_gamma_j[i];}

    /**
     * polarisability gamma atm j accessor
     */
    std::vector<int> &gamma_k() {return m_gamma_k;}

    /**
     * polarisability gamma atom j const accessor
     */
    int gamma_k(int i)const {return m_gamma_k[i];}

    /**
     * coarse grain factor accessor
     */
    std::vector<int> &cg_factor() {return m_cg_factor;}

    /**
     * coarse grain factor const accessor
     */
    int cg_factor(int i)const {return m_cg_factor[i];}

    /**
     * total coarse grain factor
     */
    double & tot_cg_factor() {return m_tot_cg_factor;}

    /**
     * stochastic dynamics const accessor
     */
    stochastic_struct const & stochastic()const { return m_stochastic; }

    /**
     * stochastic dynamics accessor
     */
    stochastic_struct & stochastic() { return m_stochastic; }

    /**
     * solute accessor.
     */
    Solute & solute() {return m_solute;}

    /**
     * perturbed solute accessor.
     */
    Perturbed_Solute & perturbed_solute() {return m_perturbed_solute;}

    /**
     * eds-perturbed solute accessor.
     */
    EDS_Perturbed_Solute & eds_perturbed_solute() {return m_eds_perturbed_solute;}

    /**
     * virtual solute accessor.
     */
    util::Virtual_Atoms_Group & virtual_atoms_group() {return m_virtual_atoms_group;}

    /**
     * const solute accessor.
     */
    Solute const & solute()const{return m_solute;}

    /**
     * const perturbed solute accessor.
     */
    Perturbed_Solute const & perturbed_solute()const{return m_perturbed_solute;}

    /**
     * const eds-perturbed solute accessor.
     */
    EDS_Perturbed_Solute const & eds_perturbed_solute()const{return m_eds_perturbed_solute;}

    /**
     * virtual solute accessor.
     */
    util::Virtual_Atoms_Group const & virtual_atoms_group()const{return m_virtual_atoms_group;}

    /**
     * number of atom types.
     */
    int num_atomtype()const {return m_num_atomtype;}

    /**
     * set number of atom types.
     */
    void set_num_atomtype(int num) {m_num_atomtype = num;}

    /**
     * number of solvents.
     */
    unsigned int num_solvents()const {return unsigned(m_num_solvent_molecules.size());}

    /**
     * number of bond types.
     */
    int num_bondtype()const {return m_num_bondtype;}
    
    /**
     * set number of bond types.
     */
    void set_num_bondtype(int num) {m_num_bondtype = num;}
    
    /**
     * number of angle types.
     */
    int num_angletype()const {return m_num_angletype;}
    
    /**
     * set number of angle types.
     */
    void set_num_angletype(int num) {m_num_angletype = num;}
    
    /**
     * number of improper dihedral types.
     */
    int num_impropertype()const {return m_num_impropertype;}
    
    /**
     * set number of improper dihedral types.
     */
    void set_num_impropertype(int num) {m_num_impropertype = num;}
    
    /**
     * solvent accessor.
     * support for multiple solvents.
     */
    Solvent & solvent(unsigned int i) {assert(i < m_solvent.size()); return m_solvent[i];}
    /**
     * const solvent accessor.
     * support for multiple solvents.
     */
    Solvent const & solvent(unsigned int i)const
    {assert(i < m_solvent.size()); return m_solvent[i];}

    /**
     * add a solvent.
     */
    void add_solvent(Solvent solv) {m_solvent.push_back(solv);}

    /**
     * add solvent to the simulation.
     * @param solv the solvent (multiple solvents).
     * @param num_molecules the number of solvent molecules to add.
     */
    void solvate(unsigned int solv, unsigned int num_molecules);

    /**
     * change the number of solvents in the simulation.
     * @param solv the solvent (has to be 0).
     * @param num_molecules the number of solvent molecules to add.
     */
    void resolvate(unsigned int solv, unsigned int num_molecules);

    /**
     * set the capacity of solute atoms
     */
    void resize(unsigned int const atoms);

    /**
     * get the total number of atoms.
     */
    unsigned int num_atoms()const {return num_solute_atoms() + num_solvent_atoms();}

    /**
     * get the number of solute atoms
     */
    unsigned int num_solute_atoms()const{return solute().num_atoms();}

    /**
     * get the total number of solvent atoms.
     */
    unsigned int num_solvent_atoms()const;

    /**
     * get the number of solvent molecules.
     */
    unsigned int num_solvent_molecules(unsigned int i)const {
      assert(i < m_num_solvent_molecules.size());
      return m_num_solvent_molecules[i];
    }

    /**
     * get the number of solvent atoms.
     */
    unsigned int num_solvent_atoms(unsigned int i)const{
      assert(i<m_num_solvent_atoms.size());
      return m_num_solvent_atoms[i];
    }

    /**
     * add a solute atom.
     */
    void add_solute_atom(std::string name, int residue_nr, int iac,
			 double mass, double charge, bool chargegroup,
			 topology::excl_cont_t::value_type exclusions,
			 topology::excl_cont_t::value_type one_four_pairs);

    /**
     * residue names.
     */
    std::vector<std::string> & residue_names() {return m_residue_name;}
    /**
     * const residue names.
     */
    std::vector<std::string> const & residue_names()const {return m_residue_name;}

    /**
     * atom names.
     */
    std::map<std::string, int> & atom_names() {return m_atom_name;}
    /**
     * const atom names.
     */
    std::map<std::string, int> const & atom_names()const {return m_atom_name;}

    /**
     * bond types
     */
    std::vector<interaction::bond_type_struct> & bond_types_harm() {return m_bond_types_harm;}
    std::vector<interaction::bond_type_struct> & bond_types_quart() {return m_bond_types_quart;}
    std::vector<interaction::bond_type_struct> const & bond_types_harm() const {return m_bond_types_harm;}
    std::vector<interaction::bond_type_struct> const & bond_types_quart() const {return m_bond_types_quart;}

    /**
     * bond angle types
     */
    std::vector<interaction::angle_type_struct> & angle_types_harm() {return m_angle_types_harm;}
    std::vector<interaction::angle_type_struct> & angle_types_cosharm() {return m_angle_types_cosharm;}
    std::vector<interaction::angle_type_struct> const & angle_types_harm() const {return m_angle_types_harm;}
    std::vector<interaction::angle_type_struct> const & angle_types_cosharm() const {return m_angle_types_cosharm;}

    /**
     * dihedral angle types
     */
    std::vector<interaction::dihedral_type_struct> & dihedral_types() {return m_dihedral_types;}
    std::vector<interaction::dihedral_type_struct> const & dihedral_types() const {return m_dihedral_types;}

    /**
     * improper dihedral angle types
     */
    std::vector<interaction::improper_dihedral_type_struct> & impdihedral_types() {return m_impdihedral_types;}
    std::vector<interaction::improper_dihedral_type_struct> const & impdihedral_types() const {return m_impdihedral_types;}

    /**
     * virtual atom types
     */
    std::vector<interaction::virtual_atom_type_struct> & virtual_atom_types() {return m_virtual_atom_types;}
    std::vector<interaction::virtual_atom_type_struct> const & virtual_atom_types() const {return m_virtual_atom_types;}
	

    /**
     * all exclusions for atom i. Exclusions, 1,4 interactions and Lennard-Jones exceptions
     */
    excl_cont_t::value_type & all_exclusion(unsigned int const i){
      assert(i < m_all_exclusion.size());
      return m_all_exclusion[i];
    }

    /**
     * const all exclusions for atom i. Exclusions, 1,4 interactions and Lennard-Jones exceptions
     */
    excl_cont_t::value_type const & all_exclusion(unsigned int const i)const{
      assert(i < m_all_exclusion.size());
      return m_all_exclusion[i];
    }

    /**
     * QM all exclusions for atom i. QM Exclusions, QM 1,4 interactions and QM Lennard-Jones exceptions
     */
    excl_cont_t & qm_all_exclusion(){
      return m_qm_all_exclusion;
    }

    /**
     * QM all exclusions for atom i. QM Exclusions, QM 1,4 interactions and QM Lennard-Jones exceptions
     */
    excl_cont_t::value_type & qm_all_exclusion(unsigned int const i){
      assert(i < m_qm_all_exclusion.size());
      return m_qm_all_exclusion[i];
    }

    /**
     * const QM all exclusions for atom i. QM Exclusions, QM 1,4 interactions and QM Lennard-Jones exceptions
     */
    excl_cont_t::value_type const & qm_all_exclusion(unsigned int const i)const{
      assert(i < m_qm_all_exclusion.size());
      return m_qm_all_exclusion[i];
    }

    /**
     * update all exclusions - updates all_exclusion container after modification of exclusions, 1,4-pairs and LJ exceptions
     */
    void update_all_exclusion();

    /**
     * update chargegroup exclusion from atomic exclusion
     */
    void update_chargegroup_exclusion();

    /**
     * chargegroup exclusions for chargegroup i.
     */
    std::set<int> & chargegroup_exclusion(unsigned int i){
      assert(i < m_chargegroup_exclusion.size());
      return m_chargegroup_exclusion[i];
    }

    /**
     * chargegroup exclusions for chargegroup i.
     */
    std::set<int> const & chargegroup_exclusion(unsigned int i)const{
      assert(i < m_chargegroup_exclusion.size());
      return m_chargegroup_exclusion[i];
    }

    /**
     * exclusions for atom i.
     */
    excl_cont_t::value_type & exclusion(unsigned int const i){
      assert(i < m_exclusion.size());
      return m_exclusion[i];
    }
    /**
     * const exclusions for atom i.
     */
    excl_cont_t::value_type const & exclusion(unsigned int const i)const{
      assert(i < m_exclusion.size());
      return m_exclusion[i];
    }

    /**
     * exclusions
     */
    excl_cont_t & exclusion() {return m_exclusion;}

    /**
     * 1,4 pairs of atom i.
     */
    excl_cont_t::value_type & one_four_pair(unsigned int const i){
      assert(i < m_one_four_pair.size());
      return m_one_four_pair[i];
    }
    /**
     * 1,4 pairs
     */
    excl_cont_t & one_four_pair(){return m_one_four_pair;}

    /**
     * QM exclusions for atom i.
     */
    excl_cont_t::value_type & qm_exclusion(unsigned int const i){
      assert(i < m_qm_exclusion.size());
      return m_qm_exclusion[i];
    }
    /**
     * const QM exclusions for atom i.
     */
    excl_cont_t::value_type const & qm_exclusion(unsigned int const i)const{
      assert(i < m_qm_exclusion.size());
      return m_qm_exclusion[i];
    }

    /**
     * QM exclusions
     */
    excl_cont_t & qm_exclusion() {return m_qm_exclusion;}

    /**
     * QM 1,4 pairs of atom i.
     */
    excl_cont_t::value_type & qm_one_four_pair(unsigned int const i){
      assert(i < m_qm_one_four_pair.size());
      return m_qm_one_four_pair[i];
    }
    /**
     * QM 1,4 pairs
     */
    excl_cont_t & qm_one_four_pair(){return m_qm_one_four_pair;}

    /**
     * chargegroup accessor.
     */
    std::vector<int> const & chargegroups()const { return m_chargegroup; }

    int chargegroup(int i)const
    {
      assert(int(m_chargegroup.size()) > i);
      return m_chargegroup[i];
    }

    /**
     * the number of chargegroups present.
     */
    unsigned int num_chargegroups()const {return unsigned(m_chargegroup.size())-1;}

    /**
     * the number of solute chargegroups.
     */
    unsigned int num_solute_chargegroups()const {return m_num_solute_chargegroups;}

    /**
     * the number of solute (sub)molecules
     */
    unsigned int num_solute_molecules()const {return m_num_solute_molecules;}

    /**
     * the number of solute (sub)molecules
     */
    unsigned int & num_solute_molecules() {return m_num_solute_molecules;}

    /**
     * the number of solute temperature groups (const)
     */
    unsigned int num_solute_temperature_groups()const {return m_num_solute_temperature_groups;}

    /**
     * the number of solute temperature groups
     */
    unsigned int & num_solute_temperature_groups() {return m_num_solute_temperature_groups;}

    /**
     * the number of solute pressure groups (const)
     */
    unsigned int num_solute_pressure_groups()const {return m_num_solute_pressure_groups;}

    /**
     * the number of solute temperature groups
     */
    unsigned int & num_solute_pressure_groups() {return m_num_solute_pressure_groups;}

    /**
     * iterator over the chargegrops
     */
    Chargegroup_Iterator chargegroup_begin()const {
      return Chargegroup_Iterator(m_chargegroup.begin());}

    /**
     * iterator on a specified chargegroup
     */
    Chargegroup_Iterator chargegroup_it(const unsigned int i)const
    {
      return Chargegroup_Iterator(m_chargegroup.begin()+i);
    }

    /**
     * end of the chargegroup iterator.
     */
    Chargegroup_Iterator chargegroup_end()const{
      return Chargegroup_Iterator(m_chargegroup.end()-1);}

    /**
     * the molecules.
     */
    std::vector<unsigned int> & molecules(){ return m_molecule;}
    /**
     * the molecules (const)
     */
    std::vector<unsigned int> const & molecules()const{ return m_molecule;}
    /**
     * iterator over the molecules.
     */
    Molecule_Iterator molecule_begin()const{
      return Molecule_Iterator(m_molecule.begin());};
    /**
     * end of molecule iterator.
     */
    Molecule_Iterator molecule_end()const{
      return Molecule_Iterator(m_molecule.end()-1);}

    /**
     * the temperature groups
     */
    std::vector<unsigned int> & temperature_groups(){ return m_temperature_group;}
    /**
     * the temperature groups (const)
     */
    std::vector<unsigned int> const & temperature_groups()const{ return m_temperature_group;}
    /**
     * iterator over the temperature groups
     */
    Temperaturegroup_Iterator temperature_group_begin()const{
      return Temperaturegroup_Iterator(m_temperature_group.begin());};
    /**
     * end of temperature groups iterator
     */
    Temperaturegroup_Iterator temperature_group_end()const{
      return Temperaturegroup_Iterator(m_temperature_group.end()-1);}

    /**
     * the pressure groups
     */
    std::vector<unsigned int> & pressure_groups(){ return m_pressure_group;}
    /**
     * the pressure groups (const)
     */
    std::vector<unsigned int> const & pressure_groups()const{ return m_pressure_group;}
    /**
     * iterator over the pressure groups
     */
    Pressuregroup_Iterator pressure_group_begin()const{
      return Pressuregroup_Iterator(m_pressure_group.begin());};
    /**
     * end of pressure groups iterator
     */
    Pressuregroup_Iterator pressure_group_end()const{
      return Pressuregroup_Iterator(m_pressure_group.end()-1);}

    /**
     * const energy group accessor.
     */
    std::vector<unsigned int> const & energy_groups()const {return m_energy_group;}

    /**
     * energy group accessor.
     */
    std::vector<unsigned int> & energy_groups() {return m_energy_group;}

    /**
     * const energy group of atoms accessor.
     */
    std::vector<unsigned int> const & atom_energy_group()const{
      return m_atom_energy_group;
    }

    /**
     * energy group of atoms accessor.
     */
    std::vector<unsigned int> & atom_energy_group(){
      return m_atom_energy_group;
    }

    /**
     * energy group of atom accessor
     */
    unsigned int atom_energy_group(unsigned int i)const{
      assert(i < m_atom_energy_group.size());
      return m_atom_energy_group[i];
    };

    /**
     * calculate constraint degrees of freedom
     */
    void calculate_constraint_dof(simulation::Multibath & mb,
				  bool rottrans_constraints,
                                  bool position_constraints,
                                  bool dih_constraints,
                                  bool ang_constraints)const;

    /**
     * check state
     */
    int check_state()const;

    /**
     * return lambda
     */
    double lambda()const {return m_lambda;}

    /**
     * return lambda (Bcast lambda to slaves)
     */
    double & lambda() {return m_lambda;}

    /**
     * set lambda
     */
    void lambda(double const l){ m_old_lambda = m_lambda; m_lambda = l;}

    /**
     * old lambda
     */
    double old_lambda()const { return m_old_lambda; }

    /**
     * return nlam (exponent to lambda)
     */
    int lambda_exp()const {return m_lambda_exp;}

    /**
     * set lambda exponent.
     */
    void lambda_exp(int const l) { m_lambda_exp = l; }

    /**
     * is the atom perturbed?
     */
    bool is_perturbed(unsigned int const i)const {
      if (i >= num_solute_atoms()) return false;

      assert(i < m_is_perturbed.size());
      return m_is_perturbed[i];
    }
    /**
     * is the atom perturbed?
     */
    std::vector<bool> & is_perturbed() { return m_is_perturbed;}

    /**
     * is the atom eds-perturbed?
     */
    bool is_eds_perturbed(unsigned int const i)const {
      if (i >= num_solute_atoms()) return false;

      assert(i < m_is_eds_perturbed.size());
      return m_is_eds_perturbed[i];
    }
    /**
     * is the atom eds-perturbed?
     */
    std::vector<bool> & is_eds_perturbed() { return m_is_eds_perturbed;}

    /**
     * ORIOL_GAMD
     * list of acceleration groups of the gaussian accelerated md atoms
     */
    std::vector<int> & gamd_accel_group(){ return m_gamd_accel_group;}

    /** 
     * ORIOL_GAMD
     * Map of interactions and where should they be saved based on the groups that define them
     */
    std::map< std::vector<unsigned int>, unsigned int> & get_gamd_interaction_pairs(){ return m_gamd_interaction_pairs;}
    /** in which acceleration group is the atom (index)
     * returns 0 if it does not exist
     */
     unsigned int gamd_interaction_group(std::vector<unsigned int> key)const{
       std::map< std::vector<unsigned int>, unsigned int>::const_iterator searched = m_gamd_interaction_pairs.find(key);
       if (searched != m_gamd_interaction_pairs.end()) {
          return m_gamd_interaction_pairs.at(key);
       } else {
         return 0;
       }
     };

    /**
     * in which acceleration group is the atom (index)
     * returns 0 if not accelerated
     **/
    int gamd_accel_group(unsigned int const i)const {
      //assert(i < m_gamd_accel_group.size());
      if (i < m_gamd_accel_group.size()){
        return m_gamd_accel_group[i];
      }
      else {
        return 0;
      }
    }
    /**
     * is the atom polarisable?
     */
    bool is_polarisable(unsigned int const i)const {

      assert(i < m_is_polarisable.size());
      return m_is_polarisable[i];
    }
    /**
     * is the atom polarisable?
     */
    std::vector<bool> & is_polarisable() { return m_is_polarisable;}

    /**
     * is the atom coarse-grained?
     */
    bool is_coarse_grained(unsigned int const i)const {

      assert(i < m_is_coarse_grained.size());
      return m_is_coarse_grained[i];
    }
    /**
     * is the atom coarse-grained?
     */
    std::vector<bool> & is_coarse_grained() { return m_is_coarse_grained;}

    /**
     * is the atom QM? - accessor
     */
    unsigned is_qm(const unsigned i)const {
      assert(i < m_is_qm.size());
      return m_is_qm[i];
    }
    
    /**
     * is the atom QM? - mutator
     */
    unsigned& is_qm(const unsigned i) {
      assert(i < m_is_qm.size());
      return m_is_qm[i];
    }

    /**
     * is the atom in the QM buffer? - accessor
     */
    int is_qm_buffer(const unsigned i)const {
      assert(i < m_is_qm_buffer.size());
      return m_is_qm_buffer[i];
    }
    
    /**
     * is the atom in the QM buffer? - mutator
     */
    int& is_qm_buffer(const unsigned i) {
      assert(i < m_is_qm_buffer.size());
      return m_is_qm_buffer[i];
    }

    /**
     * is the atom in the adaptive QM buffer? - accessor
     */
    bool is_adaptive_qm_buffer(const unsigned i)const {
      assert(i < m_is_qm_buffer.size());
      return m_is_qm_buffer[i] > 0;
    }

    /**
     * QM delta charge accessor
     */
    math::SArray &qm_delta_charge() {return m_qm_delta_charge;}

    /**
     * QM delta charge const accessor
     */
    math::SArray const & qm_delta_charge()const {return m_qm_delta_charge;}

    /**
     * QM delta charge mutator
     */
    double &qm_delta_charge(unsigned i) { return m_qm_delta_charge[i]; }

    /**
     * QM delta charge accessor
     */
    double qm_delta_charge(unsigned i) const { return m_qm_delta_charge[i]; }

    /**
     * QM atomic number accessor
     */
    const unsigned& qm_atomic_number(unsigned i) const {
      assert(i < m_qm_atomic_number.size());
      return m_qm_atomic_number[i]; 
    }

    /**
     * QM atomic number mutator
     */
    unsigned& qm_atomic_number(unsigned i) {
      assert(i < m_qm_atomic_number.size());
      return m_qm_atomic_number[i]; 
    }
    
    /**
     * QM atomic number vector mutator
     */
    const std::vector<unsigned>& qm_atomic_number() const {
      return m_qm_atomic_number; 
    }

    /**
     * QM atomic number mutator
     */
    std::vector<unsigned> & qm_atomic_number() { return m_qm_atomic_number; }

    /**
     * QM to MM link accessor
     */
    const std::set< std::pair<unsigned,unsigned> > & qmmm_link() const { return m_qmmm_link; }

    /**
     * QM to MM link mutator
     */
    std::set< std::pair<unsigned,unsigned> > & qmmm_link() { return m_qmmm_link; }

    /**
     * if atoms are linked
     */
    bool are_linked(unsigned qmi, unsigned mmi) {
      return m_qmmm_link.count(std::make_pair( qmi, mmi ));
    }

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

    /**
     * matrix for changed lambda dependence for interactions
     * between energy groups  energy groups
     */
    //std::map<std::pair<int, int>, std::pair<int, double> > const &
    //energy_group_lambdadep()const { return m_energy_group_lambdadep;}

    /**
     * matrix for changed lambda dependence for interactions
     * between energy groups  energy groups
     */
    //std::map<std::pair<int, int>, std::pair<int, double> > &
    //energy_group_lambdadep() { return m_energy_group_lambdadep;}

    /**
     * lambda prime's for different lambda dependencies.
     * (so that they do not have to be recalculated for every interaction.
     */
    //std::vector<double> & lambda_prime() { return m_lambda_prime; }

    /**
     * d lambda prime / d lambda derivative.
     * (to store the d Energy / d lambda derivative in one energy class as
     * sum d E / d lambda prime * d lambda prime / d lambda
     * for all lambda prime).
     */
    //std::vector<double> & lambda_prime_derivative()
    //{ return m_lambda_prime_derivative; }

    /**
     * alpha parameters for the perturbed energy derivatives.
     */
    //std::vector<double> const & perturbed_energy_derivative_alpha()const
    //{ return m_perturbed_energy_derivative_alpha; }

    /**
     * alpha parameters for the perturbed energy derivatives.
     */
    //std::vector<double> & perturbed_energy_derivative_alpha()
    //{ return m_perturbed_energy_derivative_alpha; }

    /**
     * individual lambda values
     */
    std::vector< std::vector< double > >  & individual_lambda(int i)
    {
      return m_individual_lambda[i];
    }
    /**
     * individual lambda values as a const
     */
    std::vector< std::vector< double > >  const & individual_lambda(int i) const
    {
      return m_individual_lambda[i];
    }

    /**
     * individual lambda derivative
     */
    std::vector< std::vector< double > >  & individual_lambda_derivative(int i)
    {
      return m_individual_lambda_derivative[i];
    }
    /**
     * individual lambda derivative as a const
     */
    std::vector< std::vector< double > >  const & individual_lambda_derivative(int i)const
    {
      return m_individual_lambda_derivative[i];
    }

    /**
     * position restraints accessor.
     */
    std::vector<position_restraint_struct> & position_restraints()
    {
      return m_position_restraint;
    }

    /**
     * const position restraints accessor.
     */
    std::vector<position_restraint_struct> const & position_restraints()const
    {
      return m_position_restraint;
    }

    /**
     * xray restraints accessor.
     */
    std::vector<xray_restraint_struct> & xray_restraints()
    {
      return m_xray_restraint;
    }

    /**
     * const xray restraints accessor.
     */
    std::vector<xray_restraint_struct> const & xray_restraints()const
    {
      return m_xray_restraint;
    }

    /**
     * xray R-free HKLs accessor.
     */
    std::vector<xray_restraint_struct> & xray_rfree()
    {
      return m_xray_rfree;
    }

    /**
     * const xray R-free accessor.
     */
    std::vector<xray_restraint_struct> const & xray_rfree()const
    {
      return m_xray_rfree;
    }

    /**
     * xray element accessor.
     */
    std::vector<std::string> & xray_elements()
    {
      return m_xray_element;
    }
    /**
     * xray constant element accessor.
     */
    std::vector<std::string> const & xray_elements()const
    {
      return m_xray_element;
    }
    /**
     * xray solvent element accessor.
     */
    std::vector<std::string> & xray_solvelements()
    {
      return m_xray_solvelement;
    }
    /**
     * xray constant solvent element accessor.
     */
    std::vector<std::string> const & xray_solvelements()const
    {
      return m_xray_solvelement;
    }

    /**
     * xray B factor accessor.
     */
    std::vector<double> & xray_b_factors()
    {
      return m_xray_b_factor;
    }
    /**
     * xray constant B factor accessor.
     */
    std::vector<double> const & xray_b_factors()const
    {
      return m_xray_b_factor;
    }
    /**
     * xray solvent B factor accessor.
     */
    std::vector<double> & xray_solv_b_factors()
    {
      return m_xray_solv_b_factor;
    }
    /**
     * xray constant B factor accessor.
     */
    std::vector<double> const & xray_solv_b_factors()const
    {
      return m_xray_solv_b_factor;
    }

    /**
     * xray occupancy accessor.
     */
    std::vector<double> & xray_occupancies()
    {
      return m_xray_occupancy;
    }
    /**
     * xray constant B factor accessor.
     */
    std::vector<double> const & xray_occupancies()const
    {
      return m_xray_occupancy;
    }
    /**
     * xray solvent occupancy accessor.
     */
    std::vector<double> & xray_solv_occupancies()
    {
      return m_xray_solv_occupancy;
    }
    /**
     * xray constant occupancy accessor.
     */
    std::vector<double> const & xray_solv_occupancies()const
    {
      return m_xray_solv_occupancy;
    }
    /**
     * definition of the asymetric unit. Pointer to the first atom in every ASU
     */
    std::vector<unsigned int> & xray_asu()
    {
      return m_xray_asu;
    }
    /**
     * definition of the asymetric unit. Pointer to the first atom in every ASU. Const version
     */
    std::vector<unsigned int> const & xray_asu() const
    {
      return m_xray_asu;
    }
    /**
     * definition of the asymetric unit. Pointer to the first atom in every ASU
     */
    std::vector<unsigned int> & sym_asu()
    {
      return m_sym_asu;
    }
    /**
     * definition of the asymetric unit. Pointer to the first atom in every ASU. Const version
     */
    std::vector<unsigned int> const & sym_asu() const
    {
      return m_sym_asu;
    }
    /**
     * xray umbrella weight accessor.
     */
    std::vector<xray_umbrella_weight_struct> & xray_umbrella_weights()
    {
      return m_xray_umbrella_weights;
    }
    /**
     * xray umbrella weight accessor.
     */
    std::vector<xray_umbrella_weight_struct> const & xray_umbrella_weights()const
    {
      return m_xray_umbrella_weights;
    }

    /**
     * atom numbers of the atoms to be symmetry restraint.
     */
    std::vector<unsigned int> & xray_sym_restraints() {
      return m_xray_sym_restraints;
    }

    /**
     * atom numbers of the atoms to be symmetry restraint.
     */
    std::vector<unsigned int> const & xray_sym_restraints() const {
      return m_xray_sym_restraints;
    }

    /**
     * atom numbers of the atoms to be symmetry restraint.
     */
    std::vector<unsigned int> & sym_restraints() {
      return m_sym_restraints;
    }

    /**
     * atom numbers of the atoms to be symmetry restraint.
     */
    std::vector<unsigned int> const & sym_restraints() const {
      return m_sym_restraints;
    }

    /**
     * sasa parameter accessor.
     */
    std::vector<sasa_parameter_struct> & sasa_parameter()
    {
      return m_sasa_parameter;
    }
    /**
     * const sasa parameter accessor.
     */
    std::vector<sasa_parameter_struct> const & sasa_parameter()const
    {
      return m_sasa_parameter;
    }
    /**
     * get the sasa parameters for atom i.
     */
    sasa_parameter_struct const & sasa_parameter(unsigned int i)const
    {
      assert(i < m_sasa_parameter.size());
      return m_sasa_parameter[i];
    }
    /**
     * the total volume of all atoms
     */
    double & sasa_volume_tot()
    {
        return m_sasa_volume_tot;
    }
    /**
     * the total volume of all atoms
     */
    double const & sasa_volume_tot() const
    {
        return m_sasa_volume_tot;
    }
    /**
     * const distance restraints accessor.
     */
    std::vector<distance_restraint_struct> const & distance_restraints()const
    {
      return m_distance_restraint;
    }
    /**
     *  distance restraints accessor.
     */
    std::vector<distance_restraint_struct>  & distance_restraints()
    {
      return m_distance_restraint;
    }
    /**
     * const perturbed distance restraints accessor.
     */
    std::vector<perturbed_distance_restraint_struct> const & perturbed_distance_restraints()const
    {
      return m_perturbed_distance_restraint;
    }
    /**
     * perturbed distance restraints accessor.
     */
    std::vector<perturbed_distance_restraint_struct>  & perturbed_distance_restraints()
    {
      return m_perturbed_distance_restraint;
    }
    /**
     * const distancefield restraints accessor.
     */
    disfield_restraint_struct const & disfield_restraints()const
    {
      return m_disfield_restraint;
    }
    /**
     *  distancefield restraints accessor.
     */
    disfield_restraint_struct  & disfield_restraints()
    {
      return m_disfield_restraint;
    }
    /**
     * const perturbed distancefield restraints accessor.
     */
    perturbed_disfield_restraint_struct const & perturbed_disfield_restraints()const
    {
      return m_perturbed_disfield_restraint;
    }
    /**
     * perturbed distancefield restraints accessor.
     */
    perturbed_disfield_restraint_struct  & perturbed_disfield_restraints()
    {
      return m_perturbed_disfield_restraint;
    }
    /**
     * const eds distance restraints accessor.
     */
    std::vector<eds_distance_restraint_struct> const & eds_distance_restraints()const
    {
      return m_eds_distance_restraint;
    }
    /**
     *  distance restraints accessor.
     */
    std::vector<eds_distance_restraint_struct>  & eds_distance_restraints()
    {
      return m_eds_distance_restraint;
    }

    /**
     * const angle restraints accessor.
     */
    std::vector<angle_restraint_struct> const & angle_restraints()const
    {
      return m_angle_restraint;
    }
    /**
     * angle restraints accessor.
     */
    std::vector<angle_restraint_struct> & angle_restraints()
    {
      return m_angle_restraint;
    }
    /**
     * const perturbed angle restraints accessor.
     */
    std::vector<perturbed_angle_restraint_struct> const & perturbed_angle_restraints()const
    {
      return m_perturbed_angle_restraint;
    }
    /**
     * perturbed angle restraints accessor.
     */
    std::vector<perturbed_angle_restraint_struct> & perturbed_angle_restraints()
    {
      return m_perturbed_angle_restraint;
    }

    /**
     * const dihedral restraints accessor.
     */
    std::vector<dihedral_restraint_struct> const & dihedral_restraints()const
    {
      return m_dihedral_restraint;
    }
    /**
     * dihedral restraints accessor.
     */
    std::vector<dihedral_restraint_struct> & dihedral_restraints()
    {
      return m_dihedral_restraint;
    }
    /**
     * const perturbed dihedral restraints accessor.
     */
    std::vector<perturbed_dihedral_restraint_struct> const & perturbed_dihedral_restraints()const
    {
      return m_perturbed_dihedral_restraint;
    }
    /**
     * perturbed dihedral restraints accessor.
     */
    std::vector<perturbed_dihedral_restraint_struct> & perturbed_dihedral_restraints()
    {
      return m_perturbed_dihedral_restraint;
    }

    /**
     * j value restraints accessor.
     */
    std::vector<jvalue_restraint_struct> & jvalue_restraints()
    {
      return m_jvalue_restraint;
    }

    /**
     * const jvalue restraints accessor.
     */
    std::vector<jvalue_restraint_struct> const & jvalue_restraints()const
    {
      return m_jvalue_restraint;
    }

    /**
     * const order parameter restraints accessor.
     */
    std::vector<order_parameter_restraint_struct> const & order_parameter_restraints()const
    {
      return m_order_parameter_restraint;
    }
    /**
     * order parameter restraints accessor.
     */
    std::vector<order_parameter_restraint_struct> & order_parameter_restraints()
    {
      return m_order_parameter_restraint;
    }

    /**
     * RDC restraints accessor.
     */
    std::vector<std::vector<rdc_restraint_struct> > & rdc_restraints()
    {
      return m_rdc_restraint;
    }
    /**
     * const RDC restraints accessor.
     */
    std::vector<std::vector<rdc_restraint_struct> > const & rdc_restraints() const
    {
      return m_rdc_restraint;
    }


    /**
     * multi cell topology accessor
     */
    Topology & multicell_topo()
    {
      assert(m_multicell_topo != NULL);
      return *m_multicell_topo;
    }

    /**
     * multi cell topology const accessor
     */
    Topology const & multicell_topo()const
    {
      return *m_multicell_topo;
    }

    /**
     * virtual grains accessor
     */
    std::vector<virtual_grain_struct> & virtual_grains()
    {
      return m_virtual_grain;
    }

    /**
     * const virtual grains accessor
     */
    std::vector<virtual_grain_struct> const & virtual_grains()const
    {
      return m_virtual_grain;
    }

    /**
     * const accessor to local elevation coordinates
     */
    std::vector<util::LE_Coordinate*> const & le_coordinates()const
    {
      return m_le_coordinates;
    }
    /**
     * accessor to local elevation coordinates
     */
    std::vector<util::LE_Coordinate*> & le_coordinates()
    {
      return m_le_coordinates;
    }

    /**
     * initialise the topology.
     * - adjust submolecules if empty
     */
    void init(simulation::Simulation const & sim,
	      std::ostream & os = std::cout,
	      bool quiet = false);

    /**
     * direct neighbours of atom i.
     */
    std::set<unsigned int> & sasa_first_neighbour(unsigned int const i){
      assert(i < m_sasa_first_neighbour.size());
      return m_sasa_first_neighbour[i];
    }
    /**
     * direct neighbours list.
     */
    std::vector<std::set<unsigned int> > & sasa_first_neighbour(){return m_sasa_first_neighbour;}

    /**
     * second neighbours of atom i.
     */
    std::set<unsigned int> & sasa_second_neighbour(unsigned int const i){
      assert(i < m_sasa_second_neighbour.size());
      return m_sasa_second_neighbour[i];
    }
    /**
     * second neighbours list.
     */
    std::vector<std::set<unsigned int> > & sasa_second_neighbour(){return m_sasa_first_neighbour;}

    /**
     * third neighbours of atom i.
     */
    std::set<unsigned int> & sasa_third_neighbour(unsigned int const i){
      assert(i < m_sasa_third_neighbour.size());
      return m_sasa_third_neighbour[i];
    }
    /**
     * third neighbours list.
     */
    std::vector<std::set<unsigned int> > & sasa_third_neighbour(){return m_sasa_third_neighbour;}

    /**
     * higher neighbours of atom i.
     */
    std::set<unsigned int> & sasa_higher_neighbour(unsigned int const i){
      assert(i < m_sasa_higher_neighbour.size());
      return m_sasa_higher_neighbour[i];
    }
    /**
     * higher neighbours list.
     */
    std::vector<std::set<unsigned int> > & sasa_higher_neighbour(){return m_sasa_higher_neighbour;}

    /**
     * accessor to LJ exceptions
     */
    std::vector<lj_exception_struct> & lj_exceptions() { return m_lj_exceptions;}
    /**
     * const accessor to LJ exceptions
     */
    const std::vector<lj_exception_struct> & lj_exceptions() const { return m_lj_exceptions;}

    /**
     * accessor to QM LJ exceptions
     */
    std::vector<lj_exception_struct> & qm_lj_exceptions() { return m_qm_lj_exceptions;}
    /**
     * const accessor to QM LJ exceptions
     */
    const std::vector<lj_exception_struct> & qm_lj_exceptions() const { return m_qm_lj_exceptions;}

  private:
    /**
     * the number of atom types
     */
    int m_num_atomtype;
    /**
     * the number of bond types
     */
    int m_num_bondtype;
    /**
     * the number of angle types
     */
    int m_num_angletype;
    /**
     * the number of improper dihedral types
     */
    int m_num_impropertype;
    
    /**
     * the solute.
     */
    Solute m_solute;

    /**
     * the perturbed solute.
     */
    Perturbed_Solute m_perturbed_solute;

    /**
     * the eds-perturbed solute
     */
    EDS_Perturbed_Solute m_eds_perturbed_solute;

    /**
     * the virtual atoms
     */
    util::Virtual_Atoms_Group m_virtual_atoms_group;

    /**
     * is the atom perturbed?
     */
    std::vector<bool> m_is_perturbed;

    /**
     * is the atom eds-perturbed?
     */
    std::vector<bool> m_is_eds_perturbed;

    /**
     * ORIOL_GAMD
     **/
    std::vector<int> m_gamd_accel_group;
    std::map< std::vector<unsigned int>, unsigned int> m_gamd_interaction_pairs;

    /**
     * is the atom polarisable?
     */
    std::vector<bool> m_is_polarisable;

    /**
     * is the atom polarisable?
     */
    std::vector<bool> m_is_coarse_grained;

    /**
     * the number of solvent molecules.
     */
    std::vector<unsigned int> m_num_solvent_molecules;

    /**
     * the number of solvent atoms.
     * vector for multiple solvents.
     */
    std::vector<unsigned int> m_num_solvent_atoms;

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
     * the inverse of the atomic masses.
     */
    math::SArray m_inverse_mass;

    /**
     * the atom charges.
     */
    math::SArray m_charge;

    /**
     * stochastic dynamics variables
     */
    stochastic_struct m_stochastic;

    /**
     * the atom exclusions.
     */
    excl_cont_t m_exclusion;

    /**
     * the atom 1-4 interactions.
     */
    excl_cont_t m_one_four_pair;

    /**
     * atom exclusions and 1-4 interactions.
     */
    excl_cont_t m_all_exclusion;

    /**
     * chargegroup exclusions
     */
    std::vector<std::set<int> > m_chargegroup_exclusion;

    /**
     * the molecules.
     */
    std::vector<unsigned int> m_molecule;

    /**
     * the temperature groups
     */
    std::vector<unsigned int> m_temperature_group;

    /**
     * the pressure groups
     */
    std::vector<unsigned int> m_pressure_group;

    /**
     * the chargegroups.
     */
    std::vector<int> m_chargegroup;

    /**
     * the number of solute chargegroups.
     */
    unsigned int m_num_solute_chargegroups;

    /**
     * the number of solute molecules.
     */
    unsigned int m_num_solute_molecules;

    /**
     * the number of solute temperature groups
     */
    unsigned int m_num_solute_temperature_groups;

    /**
     * the number of solute pressure groups
     */
    unsigned int m_num_solute_pressure_groups;

    /**
     * residue names (solute and solvent).
     */
    std::vector<std::string> m_residue_name;

    /**
     * store all possible atom names (with their iac)
     */
    std::map<std::string, int> m_atom_name;

    /**
     * store all available bond types with harmonic/quartic force constant
     */
    std::vector<interaction::bond_type_struct> m_bond_types_harm;
    std::vector<interaction::bond_type_struct> m_bond_types_quart;

    /**
     * store all available angle types with harmonic/cosine harmonic force constant
     */
    std::vector<interaction::angle_type_struct> m_angle_types_harm;
    std::vector<interaction::angle_type_struct> m_angle_types_cosharm;

    /**
     * store all available dihedral types
     */
    std::vector<interaction::dihedral_type_struct> m_dihedral_types;

    /**
     * store all available improper dihedral types
     */
    std::vector<interaction::improper_dihedral_type_struct> m_impdihedral_types;

    /**
     * store all available virtual atom types
     */
    std::vector<interaction::virtual_atom_type_struct> m_virtual_atom_types;

    /**
     * energy groups.
     */
    std::vector<unsigned int> m_energy_group;

    /**
     * energy group of atom
     */
    std::vector<unsigned int> m_atom_energy_group;

    /**
     * lambda.
     */
    double m_lambda;

    /**
     * old lambda
     */
    double m_old_lambda;

    /**
     * lambda exponent
     */
    int m_lambda_exp;

    /**
     * interaction matrix for scaled interactions
     */
    std::map<std::pair<int, int>, std::pair<double, double> >
    m_energy_group_scaling;

    /**
     * interaction matrix for interactions with
     * changed lambda dependence
     */
    //std::map<std::pair<int, int>, std::pair<int, double> >
    //m_energy_group_lambdadep;

    /**
     * lambda primes for current lambda
     */
    //std::vector<double> m_lambda_prime;

    /**
     * lambda prime / lambda derivatives
     */
    //std::vector<double> m_lambda_prime_derivative;

    /**
     * alpha parameter values for the perturbed energy derivatives
     */
    //std::vector<double> m_perturbed_energy_derivative_alpha;
    /**
     * parameter values for the individual lambdas
     */
    struct individual_lambdas_struct{
      /**
       * a's
       */
      std::vector< std::map< std::pair<int, int>, double> > a;
      /**
       * b's
       */
      std::vector< std::map< std::pair<int, int>, double> > b;
       /**
       * c's
       */
      std::vector< std::map< std::pair<int, int>, double> > c;
       /**
       * d's
       */
      std::vector< std::map< std::pair<int, int>, double> > d;
      /**
       * e's
       */
      std::vector< std::map< std::pair<int, int>, double> > e;
    } /** individual_lambdas_struct */ m_individual_lambda_parameters;

    /**
     * lambda-values for individual lambdas
     */
    std::vector< std::vector< std::vector< double > > > m_individual_lambda;

    /**
     * lambda-derivative for individual lambdas
     */
    std::vector< std::vector< std::vector< double > > > m_individual_lambda_derivative;

    /**
     * position restraints / constraints
     */
    std::vector<position_restraint_struct> m_position_restraint;

    /**
     * distance restraints
     */
    std::vector<distance_restraint_struct> m_distance_restraint;
    /**
     * perturbed distance restraints
     */
    std::vector<perturbed_distance_restraint_struct> m_perturbed_distance_restraint;
    /**
     * distancefield restraints
     */
    disfield_restraint_struct m_disfield_restraint;
    /**
     * perturbed distancefield restraints
     */
    perturbed_disfield_restraint_struct m_perturbed_disfield_restraint;
    /**
     * eds distance restraints
     */
    std::vector<eds_distance_restraint_struct> m_eds_distance_restraint;
    /**
     * angle restraints
     */
    std::vector<angle_restraint_struct> m_angle_restraint;
    /**
     * perturbed angle restraints
     */
    std::vector<perturbed_angle_restraint_struct> m_perturbed_angle_restraint;
    /**
     * dihedral restraints
     */
    std::vector<dihedral_restraint_struct> m_dihedral_restraint;
    /**
     * perturbed dihedral restraints
     */
    std::vector<perturbed_dihedral_restraint_struct> m_perturbed_dihedral_restraint;
    /**
     * jvalue restraints / constraints
     */
    std::vector<jvalue_restraint_struct> m_jvalue_restraint;
    /**
     * rdc restraints
     */
    std::vector<std::vector<rdc_restraint_struct> > m_rdc_restraint;
    /**
     * xray restraints
     */
    std::vector<xray_restraint_struct> m_xray_restraint;
    /**
     * xray R-free hkls
     */
    std::vector<xray_restraint_struct> m_xray_rfree;
    /**
     * xray element name
     */
    std::vector<std::string> m_xray_element;
    /**
     * xray element name for solvent
     */
    std::vector<std::string> m_xray_solvelement;
    /**
     * xray B factor
     */
    std::vector<double> m_xray_b_factor;
    /**
     * xray solvent B factor
     */
    std::vector<double> m_xray_solv_b_factor;
    /**
     * xray occupancies
     */
    std::vector<double> m_xray_occupancy;
    /**
     * xray solvent occupancies
     */
    std::vector<double> m_xray_solv_occupancy;
    /**
     * xray symmetry periodic copy definition. i.e. pointer to the first atom in
     * every ASU in the unit cell
     */
    std::vector<unsigned int> m_xray_asu;
    /**
     * symmetry periodic copy definition. i.e. pointer to the first atom in
     * every ASU in the unit cell
     */
    std::vector<unsigned int> m_sym_asu;
    /**
     * xray umbrella weights
     */
    std::vector<xray_umbrella_weight_struct> m_xray_umbrella_weights;
    /**
     * the atoms for symmetry restrains
     */
    std::vector<unsigned int> m_xray_sym_restraints;
    /**
     * the atoms for symmetry restrains
     */
    std::vector<unsigned int> m_sym_restraints;
     /**
     * the last atom of roto-translational constraints
     */
    unsigned int m_rottrans_last_atom;
    /**
     * order parameter restraints
     */
    std::vector<order_parameter_restraint_struct> m_order_parameter_restraint;
    /**
     * expanded topology for multiple unit cell simulations
     */
    Topology * m_multicell_topo;

    /**
     * virtual grains
     */
    std::vector<virtual_grain_struct> m_virtual_grain;

    /**
     * the atom polarisabilities
     */
     math::SArray m_polarisability;

    /**
     * the polarisation cos-charge.
     */
     math::SArray m_coscharge;

    /**
     * the polarisability damping electric field offsef @f$ E_0 @f$
     */
     math::SArray m_damping_level;

    /**
     * the polarisabiliy damping power @f$ p @f$
     */
     math::SArray m_damping_power;

    /**
     * the polarisability damping electric field offsef @f$ E_0 @f$
     */
     math::SArray m_gamma;

    /**
     * the polarisabiliy offsite atom j
     */
     std::vector<int> m_gamma_j;
     /**
     * the polarisabiliy offsite atom k
     */
     std::vector<int> m_gamma_k;

     /**
     * the coarse grain factor
     */
     std::vector<int> m_cg_factor;
     /**
     * the total coarse grain factor
     */
     double m_tot_cg_factor;

    /**
     * local elevation coordinates
     */
     std::vector<util::LE_Coordinate*> m_le_coordinates;

    /**
     * sasa parameters
     */
    std::vector<sasa_parameter_struct> m_sasa_parameter;

    /**
     * total volume of all atoms
     */
    double m_sasa_volume_tot;

    /**
     * sasa neighbour lists
     */
    std::vector< std::set<unsigned int> > m_sasa_first_neighbour;
    std::vector< std::set<unsigned int> > m_sasa_second_neighbour;
    std::vector< std::set<unsigned int> > m_sasa_third_neighbour;
    std::vector< std::set<unsigned int> > m_sasa_higher_neighbour;

    /**
     * LJ exceptions
     */
    std::vector<lj_exception_struct> m_lj_exceptions;

    /**
     * Is the atom QM
     */
    std::vector<unsigned> m_is_qm;

    /**
     * Is the QM buffer (1: yes, 0: no, -1: temporarily disabled [adaptive buffer with cutoff])
     */
    std::vector<int> m_is_qm_buffer;

    /**
     * Delta-charges to be added to interactions between QM/buffer and MM atoms
     */
    math::SArray m_qm_delta_charge;

    /**
     * QM atomic number
     */
    std::vector<unsigned> m_qm_atomic_number;

    /**
     * links of QM to MM
     */
    std::set< std::pair<unsigned,unsigned> > m_qmmm_link;

    /**
     * the QMMM LJ exclusions.
     */
    excl_cont_t m_qm_exclusion;

    /**
     * the QMMM 1-4 interactions.
     */
    excl_cont_t m_qm_one_four_pair;
     
    /**
     * QMMM exclusions and 1-4 interactions.
     */
    excl_cont_t m_qm_all_exclusion;
     
     /**
     * QMMM LJ exceptions
     */
    std::vector<lj_exception_struct> m_qm_lj_exceptions;

  }; // topology

} // topology

#endif
