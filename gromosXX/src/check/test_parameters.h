/**
 * @file test_parameters.h
 * define parameters for testing
 */

#ifndef INCLUDED_TEST_PARAMETERS_H
#define INCLUDED_TEST_PARAMETERS_H

#include <cstdint>
#include "../stdheader.h"
#include "check.h"

namespace testing {

/**
 * @brief A struct that stores data regarding error that may occur during testing (also keeps the namespace testing clean)
 * 
 */
struct Error {

  /**
   * An enum that stores error codes that may occur during testing
   * 
   */
  enum error_codes {

    /**
     * code for success
     * 
     */
    error_success,

    /**
     * code for failure
     * 
     */
    error_failure,

    /**
     * error during parsing of arguments
     * 
     */
    error_parse_arguments,

    /**
     * error during reading of input files
     * 
     */
    error_read_input,

    /**
     * error during initialization of simulation
     * 
     */
    error_md_init,

    /**
     * error during run of simulation
     * 
     */
    error_md_run,

    /**
     * error caused by too low number of steps provided for simulation
     * 
     */
    error_md_steps

  };

};

/**
 * A struct that stores parameters for a test simulation
 * 
 */
struct Parameter {

  /**
   * @brief Construct a new Parameter object
   * 
   * @param name Name of the test
   * @param binary_name Name of the binary
   */
  Parameter(const std::string& name, const std::string& binary_name, const std::string& root) 
    : name_(name), binary_name_(binary_name), root_(root) {}

  /**
   * @brief Adds a feature (input, conf, topo, ...) and associated input file to the struct
   * 
   * @param feature Feature (input, conf, topo, ...)
   * @param file Name of the file (relative to root_)
   */
  void add_input(const std::string& feature, const std::string& file) {
    input_.push_back(std::make_pair(feature, TOP_SOURCE_DIR "/" + root_ + "/" + file));
  }

  /**
   * name of the test simulation
   * 
   */
  const std::string name_;

  /**
   * name of the binary
   * 
   */
  std::string binary_name_;

  /**
   * relative path to the root of the test file
   * 
   */
  const std::string root_;

  /**
   * vector of pairs of input features / file names
   * 
   */
  std::vector<std::pair<std::string, std::string>> input_;

};

/**
 * @brief A struct that stores keys for all expected test results (also keeps the namespace testing clean)
 * 
 */
struct Key {

  /**
   * @brief A collection of all keys for expected test results
   * 
   */
  enum keys {

    // energies that may be tested

    qm_energy_init,

    // "current" energies

    energy_total_current,

    energy_kinetic_total_current,

    energy_potential_total_current,

    energy_covalent_total_current,

    energy_bonds_total_current,

    energy_angles_total_current,

    energy_impropers_total_current,

    energy_dihedrals_total_current,

    energy_crossdihedrals_total_current,

    energy_nonbonded_total_current,

    energy_lennard_jones_total_current,

    energy_coulomb_reaction_field_total_current,

    lattice_total_current,

    lattice_sum_pair_total_current,

    lattice_sum_real_space_total_current,

    lattice_sum_k_reciprocal_space_total_current,

    lattice_sum_A_term_total_current,

    lattice_sum_self_total_current,

    lattice_sum_surface_total_current,

    polarisation_self_total_current,

    special_total_current,

    sasa_total_current,

    sasa_volume_total_current,

    constraints_total_current,

    distance_restraints_total_current,

    distancefield_restraints_total_current,

    dihedral_restraints_total_current,

    position_restraints_total_current,

    jvalue_restraints_total_current,

    xray_restraints_total_current,

    local_elevation_total_current,

    order_parameter_restraints_total_current,

    symmetry_restraints_total_current,

    eds_vmix_current,

    eds_vr_current,

    eds_emax_current,

    eds_emin_current,

    eds_glob_min_current,

    eds_glob_min_fluc_current,
    
    entropy_current,

    qmmm_total_current,

    bs_leus_energy_current,

    rdc_value_total_current,

    angle_restraints_total_current,

    // "old" energies

    energy_total_old,

    energy_kinetic_total_old,

    energy_potential_total_old,

    energy_covalent_total_old,

    energy_bonds_total_old,

    energy_angles_total_old,

    energy_impropers_total_old,

    energy_dihedrals_total_old,

    energy_crossdihedrals_total_old,

    energy_nonbonded_total_old,

    energy_lennard_jones_total_old,

    energy_coulomb_reaction_field_total_old,

    lattice_total_old,

    lattice_sum_pair_total_old,

    lattice_sum_real_space_total_old,

    lattice_sum_k_reciprocal_space_total_old,

    lattice_sum_A_term_total_old,

    lattice_sum_self_total_old,

    lattice_sum_surface_total_old,

    polarisation_self_total_old,

    special_total_old,

    sasa_total_old,

    sasa_volume_total_old,

    constraints_total_old,

    distance_restraints_total_old,

    distancefield_restraints_total_old,

    dihedral_restraints_total_old,

    position_restraints_total_old,

    jvalue_restraints_total_old,

    xray_restraints_total_old,

    local_elevation_total_old,

    order_parameter_restraints_total_old,

    symmetry_restraints_total_old,

    eds_vmix_old,

    eds_vr_old,

    eds_emax_old,

    eds_emin_old,

    eds_glob_min_old,

    eds_glob_min_fluc_old,
    
    entropy_old,

    qmmm_total_old,

    bs_leus_energy_old,

    rdc_value_total_old,

    angle_restraints_total_old,
    
    // parameters for QM/MM simulations

    qmmm,

    qm_ch,

    qm_software,

    qm_cutoff,

    qm_lj,

    qm_constraint,

    qm_mm_scale,

    // unit conversion factors

    unit_factor_length,

    unit_factor_energy,

    unit_factor_force,

    unit_factor_charge,

    // QM software binary name and input files

    binary,

    input_file,

    input_coordinate_file,

    input_pointcharges_file,

    output_file,

    output_gradient_file,

    output_mm_gradient_file,

    qm_worker_name,

    // QM zone parameters

    qm_zone_size_init,

    mm_zone_size_init,

    qm_zone_charge,

    qm_zone_spin_mult,

    qm_atom_charge_init,

    qm_atom_force_init,

    mm_atom_force_init,

    mm_atom_cos_force_init,

    // element symbols

    symbol_h,

    symbol_c,

    symbol_o,

    symbol_p,

    symbol_pd,

    // element numbers

    element_1,

    element_6,

    element_8,

    element_15,

    element_20,

    element_26,

    element_29,

    element_46

  };

};

/**
 * @brief A struct that stores expected test results
 * 
 */
struct Results {

  /**
   * @brief Results of type std::string
   * 
   */
  std::unordered_map<Key::keys, std::string> strs_;

  /**
   * @brief Results of type double
   * 
   */
  std::unordered_map<Key::keys, double> doubles_;

  /**
   * @brief Results of type int
   * 
   */
  std::unordered_map<Key::keys, int> ints_;

  /**
   * @brief A code for nullptr
   * 
   */
  const nullptr_t UNINITIALIZED = nullptr;

};

/**
 * @brief Stores parameters to test the mechanical embedding scheme with the Orca worker
 * 
 */
extern Parameter orca_parameter_mechanical;

/**
 * @brief Stores expected results for tests of the mechanical embedding scheme with the Orca worker
 * 
 */
extern Results orca_results_mechanical;

/**
 * @brief Stores parameters to test the electrostatic embedding scheme with the Orca worker
 * 
 */
extern Parameter orca_parameter_electrostatic;

/**
 * @brief Stores expected results for tests of the electrostatic embedding scheme with the Orca worker
 * 
 */
extern Results orca_results_electrostatic;

} // namespace testing

#endif /* INCLUDED_TEST_PARAMETERS_H */