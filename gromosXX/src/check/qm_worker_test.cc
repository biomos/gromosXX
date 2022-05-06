/**
 * @file qm_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"

#include "qm_worker_test.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

#include "../math/gmath.h"
#include "../math/transformation.h"
#include "../math/volume.h"

#include "../topology/topology.h"
#include "../util/template_split.h"

namespace testing {

QM_Worker_Test::QM_Worker_Test(const Parameter& parameter) : test_sim_(parameter) {}

void QM_Worker_Test::check_qmmm_parameter() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key::keys, int>& ints_res = results_.ints_;
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  EXPECT_EQ(param.qmmm.qmmm, ints_res[Key::qmmm]);
  EXPECT_EQ(param.qmmm.qm_ch, ints_res[Key::qm_ch]);
  EXPECT_EQ(param.qmmm.software, ints_res[Key::qm_software]);
  EXPECT_EQ(param.qmmm.cutoff, doubles_res[Key::qm_cutoff]);
  EXPECT_EQ(param.qmmm.qm_lj, ints_res[Key::qm_lj]);
  EXPECT_EQ(param.qmmm.qm_constraint, ints_res[Key::qm_constraint]);
  EXPECT_EQ(param.qmmm.mm_scale, doubles_res[Key::qm_mm_scale]);
}

void QM_Worker_Test::check_qm_zone_param() {
  // references to shorten the code
  const simulation::Parameter& param = test_sim_.sim().param();
  std::unordered_map<Key::keys, int>& ints_res = results_.ints_;
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // test the qm zone
  EXPECT_EQ(qm_zone_ptr->qm.size(), ints_res[Key::qm_zone_size_init]);
  EXPECT_EQ(qm_zone_ptr->mm.size(), ints_res[Key::mm_zone_size_init]);
  EXPECT_EQ(qm_zone_ptr->charge(), ints_res[Key::qm_zone_charge]);
  EXPECT_EQ(qm_zone_ptr->spin_mult(), ints_res[Key::qm_zone_spin_mult]);
  EXPECT_EQ(qm_zone_ptr->QM_energy(), doubles_res[Key::qm_energy_init]);
}

void QM_Worker_Test::check_qm_interaction_ptr() {
  // check if the qm interaction, qm_worker, and qm_zone have been initialized
  EXPECT_TRUE(qmmm_interaction_ptr != results_.UNINITIALIZED) << "qmmm_interaction not initialized";
  EXPECT_TRUE(qm_worker_ptr != results_.UNINITIALIZED) << "qm_worker not initialized";
  EXPECT_TRUE(qm_zone_ptr != results_.UNINITIALIZED) << "qm_zone not initialized";
}

void QM_Worker_Test::check_simulation_results() {
  check_simulation_results_energies();
  check_simulation_results_forces();
  check_simulation_results_velocities();
  check_simulation_results_positions();
  check_simulation_results_temperature_baths();
  check_simulation_results_bonded_terms();
  check_simulation_results_nonbonded_terms();
  check_simulation_results_special_terms();
  check_simulation_results_mass();
  check_simulation_results_temperature();
  check_simulation_results_volume();
  check_simulation_results_pressure();
}

void QM_Worker_Test::check_simulation_results_energies() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.total, doubles_res[Key::energy_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.kinetic_total, doubles_res[Key::energy_kinetic_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.potential_total, doubles_res[Key::energy_potential_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bonded_total, doubles_res[Key::energy_covalent_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_total, doubles_res[Key::energy_bonds_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_total, doubles_res[Key::energy_angles_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_total, doubles_res[Key::energy_impropers_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_total, doubles_res[Key::energy_dihedrals_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_total, doubles_res[Key::energy_crossdihedrals_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.nonbonded_total, doubles_res[Key::energy_nonbonded_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_total, doubles_res[Key::energy_lennard_jones_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_total, doubles_res[Key::energy_coulomb_reaction_field_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_pair_total, doubles_res[Key::lattice_sum_pair_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_realspace_total, doubles_res[Key::lattice_sum_real_space_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_kspace_total, doubles_res[Key::lattice_sum_k_reciprocal_space_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_a_term_total, doubles_res[Key::lattice_sum_A_term_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_self_total, doubles_res[Key::lattice_sum_self_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_surface_total, doubles_res[Key::lattice_sum_surface_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.self_total, doubles_res[Key::polarisation_self_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.special_total, doubles_res[Key::special_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_total, doubles_res[Key::sasa_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_volume_total, doubles_res[Key::sasa_volume_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.constraints_total, doubles_res[Key::constraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.distanceres_total, doubles_res[Key::distance_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.disfieldres_total, doubles_res[Key::distancefield_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihrest_total, doubles_res[Key::dihedral_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.posrest_total, doubles_res[Key::position_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.jvalue_total, doubles_res[Key::jvalue_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.xray_total, doubles_res[Key::xray_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.leus_total, doubles_res[Key::local_elevation_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.oparam_total, doubles_res[Key::order_parameter_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.symrest_total, doubles_res[Key::symmetry_restraints_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.eds_vmix, doubles_res[Key::eds_vmix_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.eds_vr, doubles_res[Key::eds_vr_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.eds_emax, doubles_res[Key::eds_emax_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.eds_emin, doubles_res[Key::eds_emin_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.eds_globmin, doubles_res[Key::eds_glob_min_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.eds_globminfluc, doubles_res[Key::eds_glob_min_fluc_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.entropy_term, doubles_res[Key::entropy_current] , epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.qm_total, doubles_res[Key::qmmm_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bsleus_total, doubles_res[Key::bs_leus_energy_current] , epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.rdc_total, doubles_res[Key::rdc_value_total_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angrest_total, doubles_res[Key::angle_restraints_total_current], epsilon_);
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.total, doubles_res[Key::energy_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.kinetic_total, doubles_res[Key::energy_kinetic_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.potential_total, doubles_res[Key::energy_potential_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bonded_total, doubles_res[Key::energy_covalent_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_total, doubles_res[Key::energy_bonds_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_total, doubles_res[Key::energy_angles_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_total, doubles_res[Key::energy_impropers_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_total, doubles_res[Key::energy_dihedrals_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_total, doubles_res[Key::energy_crossdihedrals_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.nonbonded_total, doubles_res[Key::energy_nonbonded_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_total, doubles_res[Key::energy_lennard_jones_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_total, doubles_res[Key::energy_coulomb_reaction_field_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_pair_total, doubles_res[Key::lattice_sum_pair_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_realspace_total, doubles_res[Key::lattice_sum_real_space_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_kspace_total, doubles_res[Key::lattice_sum_k_reciprocal_space_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_a_term_total, doubles_res[Key::lattice_sum_A_term_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_self_total, doubles_res[Key::lattice_sum_self_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_surface_total, doubles_res[Key::lattice_sum_surface_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.self_total, doubles_res[Key::polarisation_self_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.special_total, doubles_res[Key::special_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_total, doubles_res[Key::sasa_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_volume_total, doubles_res[Key::sasa_volume_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.constraints_total, doubles_res[Key::constraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.distanceres_total, doubles_res[Key::distance_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.disfieldres_total, doubles_res[Key::distancefield_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihrest_total, doubles_res[Key::dihedral_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.posrest_total, doubles_res[Key::position_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.jvalue_total, doubles_res[Key::jvalue_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.xray_total, doubles_res[Key::xray_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.leus_total, doubles_res[Key::local_elevation_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.oparam_total, doubles_res[Key::order_parameter_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.symrest_total, doubles_res[Key::symmetry_restraints_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.eds_vmix, doubles_res[Key::eds_vmix_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.eds_vr, doubles_res[Key::eds_vr_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.eds_emax, doubles_res[Key::eds_emax_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.eds_emin, doubles_res[Key::eds_emin_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.eds_globmin, doubles_res[Key::eds_glob_min_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.eds_globminfluc, doubles_res[Key::eds_glob_min_fluc_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.entropy_term, doubles_res[Key::entropy_old] , epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.qm_total, doubles_res[Key::qmmm_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bsleus_total, doubles_res[Key::bs_leus_energy_old] , epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.rdc_total, doubles_res[Key::rdc_value_total_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angrest_total, doubles_res[Key::angle_restraints_total_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_forces() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" forces
  EXPECT_NEAR(test_sim_.conf().current().force[0][0], doubles_res[Key::force_pos_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[0][1], doubles_res[Key::force_pos_0_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[0][2], doubles_res[Key::force_pos_0_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[1][0], doubles_res[Key::force_pos_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[1][1], doubles_res[Key::force_pos_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[1][2], doubles_res[Key::force_pos_1_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[200][0], doubles_res[Key::force_pos_200_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[200][1], doubles_res[Key::force_pos_200_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[200][2], doubles_res[Key::force_pos_200_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[210][0], doubles_res[Key::force_pos_210_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[210][1], doubles_res[Key::force_pos_210_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().force[210][2], doubles_res[Key::force_pos_210_2_current], epsilon_);
  // "old" forces
  EXPECT_NEAR(test_sim_.conf().old().force[0][0], doubles_res[Key::force_pos_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[0][1], doubles_res[Key::force_pos_0_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[0][2], doubles_res[Key::force_pos_0_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[1][0], doubles_res[Key::force_pos_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[1][1], doubles_res[Key::force_pos_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[1][2], doubles_res[Key::force_pos_1_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[200][0], doubles_res[Key::force_pos_200_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[200][1], doubles_res[Key::force_pos_200_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[200][2], doubles_res[Key::force_pos_200_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[210][0], doubles_res[Key::force_pos_210_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[210][1], doubles_res[Key::force_pos_210_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().force[210][2], doubles_res[Key::force_pos_210_2_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_velocities() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" velocities
  EXPECT_NEAR(test_sim_.conf().current().vel[0][0], doubles_res[Key::velocities_pos_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[0][1], doubles_res[Key::velocities_pos_0_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[0][2], doubles_res[Key::velocities_pos_0_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[1][0], doubles_res[Key::velocities_pos_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[1][1], doubles_res[Key::velocities_pos_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[1][2], doubles_res[Key::velocities_pos_1_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[200][0], doubles_res[Key::velocities_pos_200_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[200][1], doubles_res[Key::velocities_pos_200_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[200][2], doubles_res[Key::velocities_pos_200_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[210][0], doubles_res[Key::velocities_pos_210_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[210][1], doubles_res[Key::velocities_pos_210_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().vel[210][2], doubles_res[Key::velocities_pos_210_2_current], epsilon_);
  // "old" velocities
  EXPECT_NEAR(test_sim_.conf().old().vel[0][0], doubles_res[Key::velocities_pos_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[0][1], doubles_res[Key::velocities_pos_0_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[0][2], doubles_res[Key::velocities_pos_0_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[1][0], doubles_res[Key::velocities_pos_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[1][1], doubles_res[Key::velocities_pos_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[1][2], doubles_res[Key::velocities_pos_1_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[200][0], doubles_res[Key::velocities_pos_200_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[200][1], doubles_res[Key::velocities_pos_200_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[200][2], doubles_res[Key::velocities_pos_200_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[210][0], doubles_res[Key::velocities_pos_210_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[210][1], doubles_res[Key::velocities_pos_210_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().vel[210][2], doubles_res[Key::velocities_pos_210_2_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_positions() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  math::VArray positions_gathered_current;
  math::VArray positions_gathered_old;
  configuration::Configuration& conf = test_sim_.conf();
  topology::Topology& topo = test_sim_.topo();
  // helper function to gather positions - should be replaced with a new function in e.g. out_configuration.cc
  SPLIT_BOUNDARY(preprocess_positions, conf, topo, conf.current().pos, positions_gathered_current);
  SPLIT_BOUNDARY(preprocess_positions, conf, topo, conf.old().pos, positions_gathered_old);
  // "current" positions
  EXPECT_NEAR(positions_gathered_current[0][0], doubles_res[Key::positions_pos_0_0_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[0][1], doubles_res[Key::positions_pos_0_1_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[0][2], doubles_res[Key::positions_pos_0_2_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[1][0], doubles_res[Key::positions_pos_1_0_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[1][1], doubles_res[Key::positions_pos_1_1_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[1][2], doubles_res[Key::positions_pos_1_2_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[200][0], doubles_res[Key::positions_pos_200_0_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[200][1], doubles_res[Key::positions_pos_200_1_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[200][2], doubles_res[Key::positions_pos_200_2_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[210][0], doubles_res[Key::positions_pos_210_0_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[210][1], doubles_res[Key::positions_pos_210_1_current], epsilon_);
  EXPECT_NEAR(positions_gathered_current[210][2], doubles_res[Key::positions_pos_210_2_current], epsilon_);
  // "old" positions
  EXPECT_NEAR(positions_gathered_old[0][0], doubles_res[Key::positions_pos_0_0_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[0][1], doubles_res[Key::positions_pos_0_1_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[0][2], doubles_res[Key::positions_pos_0_2_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[1][0], doubles_res[Key::positions_pos_1_0_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[1][1], doubles_res[Key::positions_pos_1_1_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[1][2], doubles_res[Key::positions_pos_1_2_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[200][0], doubles_res[Key::positions_pos_200_0_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[200][1], doubles_res[Key::positions_pos_200_1_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[200][2], doubles_res[Key::positions_pos_200_2_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[210][0], doubles_res[Key::positions_pos_210_0_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[210][1], doubles_res[Key::positions_pos_210_1_old], epsilon_);
  EXPECT_NEAR(positions_gathered_old[210][2], doubles_res[Key::positions_pos_210_2_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_temperature_baths() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.kinetic_energy[0], doubles_res[Key::kinetic_total_bath_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.kinetic_energy[1], doubles_res[Key::kinetic_total_bath_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.com_kinetic_energy[0], doubles_res[Key::centre_of_mass_bath_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.com_kinetic_energy[1], doubles_res[Key::centre_of_mass_bath_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ir_kinetic_energy[0], doubles_res[Key::internal_rotational_bath_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ir_kinetic_energy[1], doubles_res[Key::internal_rotational_bath_1_current], epsilon_);
  // "old" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.kinetic_energy[0], doubles_res[Key::kinetic_total_bath_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.kinetic_energy[1], doubles_res[Key::kinetic_total_bath_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.com_kinetic_energy[0], doubles_res[Key::centre_of_mass_bath_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.com_kinetic_energy[1], doubles_res[Key::centre_of_mass_bath_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ir_kinetic_energy[0], doubles_res[Key::internal_rotational_bath_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ir_kinetic_energy[1], doubles_res[Key::internal_rotational_bath_1_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_bonded_terms() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_energy[0], doubles_res[Key::bond_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_energy[1], doubles_res[Key::bond_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_energy[0], doubles_res[Key::angle_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_energy[1], doubles_res[Key::angle_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_energy[0], doubles_res[Key::improper_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_energy[1], doubles_res[Key::improper_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_energy[0], doubles_res[Key::dihedral_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_energy[1], doubles_res[Key::dihedral_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_energy[0], doubles_res[Key::crossdihedral_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_energy[1], doubles_res[Key::crossdihedral_energy_group_1_current], epsilon_);
  // "old" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_energy[0], doubles_res[Key::bond_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_energy[1], doubles_res[Key::bond_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_energy[0], doubles_res[Key::angle_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_energy[1], doubles_res[Key::angle_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_energy[0], doubles_res[Key::improper_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_energy[1], doubles_res[Key::improper_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_energy[0], doubles_res[Key::dihedral_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_energy[1], doubles_res[Key::dihedral_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_energy[0], doubles_res[Key::crossdihedral_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_energy[1], doubles_res[Key::crossdihedral_energy_group_1_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_nonbonded_terms() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[0][0], doubles_res[Key::lennard_jones_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[1][0], doubles_res[Key::lennard_jones_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[1][1], doubles_res[Key::lennard_jones_group_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[0][0], doubles_res[Key::coulomb_reaction_field_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[1][0], doubles_res[Key::coulomb_reaction_field_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[1][1], doubles_res[Key::coulomb_reaction_field_group_1_1_current], epsilon_);
  // lattice sum energies cannot be tested like this as they are not calculated right now. values in .tre file correspond to totals
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[0][0], doubles_res[Key::lattice_sum_reciprocal_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[1][0], doubles_res[Key::lattice_sum_reciprocal_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[1][1], doubles_res[Key::lattice_sum_reciprocal_group_1_1_current], epsilon_);  
  // "old" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[0][0], doubles_res[Key::lennard_jones_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[1][0], doubles_res[Key::lennard_jones_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[1][1], doubles_res[Key::lennard_jones_group_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[0][0], doubles_res[Key::coulomb_reaction_field_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[1][0], doubles_res[Key::coulomb_reaction_field_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[1][1], doubles_res[Key::coulomb_reaction_field_group_1_1_old], epsilon_);
  // lattice sum energies cannot be tested like this as they are not calculated right now. values in .tre file correspond to totals
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[0][0], doubles_res[Key::lattice_sum_reciprocal_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[1][0], doubles_res[Key::lattice_sum_reciprocal_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[1][1], doubles_res[Key::lattice_sum_reciprocal_group_1_1_old], epsilon_); 
}

void QM_Worker_Test::check_simulation_results_special_terms() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // some values are not stored with Gromos, instead "0" is printed
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.constraints_energy[0], doubles_res[Key::contraints_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.posrest_energy[0], doubles_res[Key::pos_restraints_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.distanceres_energy[0], doubles_res[Key::dist_restraints_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.disfieldres_energy[0], doubles_res[Key::disfield_res_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angrest_energy[0], doubles_res[Key::angle_restraint_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihrest_energy[0], doubles_res[Key::dihedral_restraints_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_energy[0], doubles_res[Key::sasa_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_volume_energy[0], doubles_res[Key::sasa_vol_group_0_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.rdc_energy[0], doubles_res[Key::rdc_group_0_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_0_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.constraints_energy[1], doubles_res[Key::contraints_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.posrest_energy[1], doubles_res[Key::pos_restraints_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.distanceres_energy[1], doubles_res[Key::dist_restraints_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.disfieldres_energy[1], doubles_res[Key::disfield_res_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angrest_energy[1], doubles_res[Key::angle_restraint_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihrest_energy[1], doubles_res[Key::dihedral_restraints_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_energy[1], doubles_res[Key::sasa_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_volume_energy[1], doubles_res[Key::sasa_vol_group_1_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.rdc_energy[1], doubles_res[Key::rdc_group_1_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_1_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_1_current], epsilon_);
  // "old" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.constraints_energy[0], doubles_res[Key::contraints_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.posrest_energy[0], doubles_res[Key::pos_restraints_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.distanceres_energy[0], doubles_res[Key::dist_restraints_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.disfieldres_energy[0], doubles_res[Key::disfield_res_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angrest_energy[0], doubles_res[Key::angle_restraint_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihrest_energy[0], doubles_res[Key::dihedral_restraints_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_energy[0], doubles_res[Key::sasa_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_volume_energy[0], doubles_res[Key::sasa_vol_group_0_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.rdc_energy[0], doubles_res[Key::rdc_group_0_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_0_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.constraints_energy[1], doubles_res[Key::contraints_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.posrest_energy[1], doubles_res[Key::pos_restraints_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.distanceres_energy[1], doubles_res[Key::dist_restraints_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.disfieldres_energy[1], doubles_res[Key::disfield_res_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angrest_energy[1], doubles_res[Key::angle_restraint_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihrest_energy[1], doubles_res[Key::dihedral_restraints_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_energy[1], doubles_res[Key::sasa_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_volume_energy[1], doubles_res[Key::sasa_vol_group_1_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.rdc_energy[1], doubles_res[Key::rdc_group_1_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_1_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_1_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_mass() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // current
  EXPECT_NEAR(math::sum(test_sim_.topo().mass()), doubles_res[Key::mass_current], epsilon_);
  // old
  EXPECT_NEAR(math::sum(test_sim_.topo().mass()), doubles_res[Key::mass_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_temperature() {
  std::unordered_map<Key::keys, int>& ints_res = results_.ints_;
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  configuration::Energy& energies_current = test_sim_.conf().current().energies;
  configuration::Energy& energies_old = test_sim_.conf().old().energies;
  simulation::Multibath& baths = test_sim_.sim().multibath();
  // lambda for quick conversion
  auto energy_to_temperature = [](double kinetic_energy, double dof) { return 2 * kinetic_energy / math::k_Boltzmann / dof; };
  // current
  EXPECT_EQ(baths.size(), ints_res[Key::num_temperature_coupling_baths_current]);
  EXPECT_NEAR(energy_to_temperature(energies_current.kinetic_energy[0], baths[0].dof), doubles_res[Key::temperature_total_bath_0_current], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_current.kinetic_energy[1], baths[1].dof), doubles_res[Key::temperature_total_bath_1_current], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_current.com_kinetic_energy[0], baths[0].com_dof), doubles_res[Key::temperature_com_bath_0_current], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_current.com_kinetic_energy[1], baths[1].com_dof), doubles_res[Key::temperature_com_bath_1_current], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_current.ir_kinetic_energy[0], baths[0].ir_dof), doubles_res[Key::temperature_ir_bath_0_current], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_current.ir_kinetic_energy[1], baths[1].ir_dof), doubles_res[Key::temperature_ir_bath_1_current], epsilon_);
  EXPECT_NEAR(baths[0].scale, doubles_res[Key::temperature_scaling_factor_bath_0_current], epsilon_);
  EXPECT_NEAR(baths[1].scale, doubles_res[Key::temperature_scaling_factor_bath_1_current], epsilon_);
  // old
  EXPECT_EQ(baths.size(), ints_res[Key::num_temperature_coupling_baths_current]);
  EXPECT_NEAR(energy_to_temperature(energies_old.kinetic_energy[0], baths[0].dof), doubles_res[Key::temperature_total_bath_0_old], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_old.kinetic_energy[1], baths[1].dof), doubles_res[Key::temperature_total_bath_1_old], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_old.com_kinetic_energy[0], baths[0].com_dof), doubles_res[Key::temperature_com_bath_0_old], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_old.com_kinetic_energy[1], baths[1].com_dof), doubles_res[Key::temperature_com_bath_1_old], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_old.ir_kinetic_energy[0], baths[0].ir_dof), doubles_res[Key::temperature_ir_bath_0_old], epsilon_);
  EXPECT_NEAR(energy_to_temperature(energies_old.ir_kinetic_energy[1], baths[1].ir_dof), doubles_res[Key::temperature_ir_bath_1_old], epsilon_);
  EXPECT_NEAR(baths[0].scale, doubles_res[Key::temperature_scaling_factor_bath_0_old], epsilon_);
  EXPECT_NEAR(baths[1].scale, doubles_res[Key::temperature_scaling_factor_bath_1_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_volume() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // current
  EXPECT_NEAR(math::volume(test_sim_.conf().current().box, test_sim_.conf().boundary_type), doubles_res[Key::volume_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(0)(0), doubles_res[Key::box_k_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(0)(1), doubles_res[Key::box_k_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(0)(2), doubles_res[Key::box_k_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(1)(0), doubles_res[Key::box_l_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(1)(1), doubles_res[Key::box_l_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(1)(2), doubles_res[Key::box_l_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(2)(0), doubles_res[Key::box_m_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(2)(1), doubles_res[Key::box_m_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().box(2)(2), doubles_res[Key::box_m_2_current], epsilon_);
  // old
  EXPECT_NEAR(math::volume(test_sim_.conf().old().box, test_sim_.conf().boundary_type), doubles_res[Key::volume_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(0)(0), doubles_res[Key::box_k_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(0)(1), doubles_res[Key::box_k_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(0)(2), doubles_res[Key::box_k_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(1)(0), doubles_res[Key::box_l_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(1)(1), doubles_res[Key::box_l_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(1)(2), doubles_res[Key::box_l_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(2)(0), doubles_res[Key::box_m_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(2)(1), doubles_res[Key::box_m_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().box(2)(2), doubles_res[Key::box_m_2_old], epsilon_);
}

void QM_Worker_Test::check_simulation_results_pressure() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // current
  double pressure_current = (test_sim_.conf().current().pressure_tensor(0, 0)
  + test_sim_.conf().current().pressure_tensor(1, 1)
  + test_sim_.conf().current().pressure_tensor(2, 2)) / 3.0;
  // the virial is stored internally as just the outer product of positions and forces
  // so without the -0.5 prefactor
  test_sim_.conf().current().virial_tensor *= -0.5;
  double virial_current = (test_sim_.conf().current().virial_tensor(0, 0)
  + test_sim_.conf().current().virial_tensor(1, 1)
  + test_sim_.conf().current().virial_tensor(2, 2)) / 3.0;
  double molecular_kinetic_energy_current = (test_sim_.conf().current().kinetic_energy_tensor(0, 0)
  + test_sim_.conf().current().kinetic_energy_tensor(1, 1)
  + test_sim_.conf().current().kinetic_energy_tensor(2, 2)) / 3.0;
  EXPECT_NEAR(pressure_current, doubles_res[Key::pressure_current], epsilon_);
  EXPECT_NEAR(virial_current, doubles_res[Key::virial_current], epsilon_);
  EXPECT_NEAR(molecular_kinetic_energy_current, doubles_res[Key::molecular_kinetic_energy_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(0, 0), doubles_res[Key::pressure_tensor_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(0, 1), doubles_res[Key::pressure_tensor_0_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(0, 2), doubles_res[Key::pressure_tensor_0_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(1, 0), doubles_res[Key::pressure_tensor_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(1, 1), doubles_res[Key::pressure_tensor_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(1, 2), doubles_res[Key::pressure_tensor_1_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(2, 0), doubles_res[Key::pressure_tensor_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(2, 1), doubles_res[Key::pressure_tensor_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().pressure_tensor(2, 2), doubles_res[Key::pressure_tensor_2_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(0, 0), doubles_res[Key::virial_tensor_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(0, 1), doubles_res[Key::virial_tensor_0_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(0, 2), doubles_res[Key::virial_tensor_0_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(1, 0), doubles_res[Key::virial_tensor_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(1, 1), doubles_res[Key::virial_tensor_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(1, 2), doubles_res[Key::virial_tensor_1_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(2, 0), doubles_res[Key::virial_tensor_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(2, 1), doubles_res[Key::virial_tensor_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().virial_tensor(2, 2), doubles_res[Key::virial_tensor_2_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(0, 0), doubles_res[Key::molecular_kinetic_energy_tensor_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(0, 1), doubles_res[Key::molecular_kinetic_energy_tensor_0_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(0, 2), doubles_res[Key::molecular_kinetic_energy_tensor_0_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(1, 0), doubles_res[Key::molecular_kinetic_energy_tensor_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(1, 1), doubles_res[Key::molecular_kinetic_energy_tensor_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(1, 2), doubles_res[Key::molecular_kinetic_energy_tensor_1_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(2, 0), doubles_res[Key::molecular_kinetic_energy_tensor_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(2, 1), doubles_res[Key::molecular_kinetic_energy_tensor_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().kinetic_energy_tensor(2, 2), doubles_res[Key::molecular_kinetic_energy_tensor_2_2_current], epsilon_);
  // old
  double pressure_old = (test_sim_.conf().old().pressure_tensor(0, 0)
  + test_sim_.conf().old().pressure_tensor(1, 1)
  + test_sim_.conf().old().pressure_tensor(2, 2)) / 3.0;
  // the virial is stored internally as just the outer product of positions and forces
  // so without the -0.5 prefactor
  test_sim_.conf().old().virial_tensor *= -0.5;
  double virial_old = (test_sim_.conf().old().virial_tensor(0, 0)
  + test_sim_.conf().old().virial_tensor(1, 1)
  + test_sim_.conf().old().virial_tensor(2, 2)) / 3.0;
  double molecular_kinetic_energy_old = (test_sim_.conf().old().kinetic_energy_tensor(0, 0)
  + test_sim_.conf().old().kinetic_energy_tensor(1, 1)
  + test_sim_.conf().old().kinetic_energy_tensor(2, 2)) / 3.0;
  EXPECT_NEAR(pressure_old, doubles_res[Key::pressure_old], epsilon_);
  EXPECT_NEAR(virial_old, doubles_res[Key::virial_old], epsilon_);
  EXPECT_NEAR(molecular_kinetic_energy_old, doubles_res[Key::molecular_kinetic_energy_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(0, 0), doubles_res[Key::pressure_tensor_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(0, 1), doubles_res[Key::pressure_tensor_0_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(0, 2), doubles_res[Key::pressure_tensor_0_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(1, 0), doubles_res[Key::pressure_tensor_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(1, 1), doubles_res[Key::pressure_tensor_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(1, 2), doubles_res[Key::pressure_tensor_1_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(2, 0), doubles_res[Key::pressure_tensor_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(2, 1), doubles_res[Key::pressure_tensor_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().pressure_tensor(2, 2), doubles_res[Key::pressure_tensor_2_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(0, 0), doubles_res[Key::virial_tensor_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(0, 1), doubles_res[Key::virial_tensor_0_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(0, 2), doubles_res[Key::virial_tensor_0_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(1, 0), doubles_res[Key::virial_tensor_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(1, 1), doubles_res[Key::virial_tensor_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(1, 2), doubles_res[Key::virial_tensor_1_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(2, 0), doubles_res[Key::virial_tensor_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(2, 1), doubles_res[Key::virial_tensor_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().virial_tensor(2, 2), doubles_res[Key::virial_tensor_2_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(0, 0), doubles_res[Key::molecular_kinetic_energy_tensor_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(0, 1), doubles_res[Key::molecular_kinetic_energy_tensor_0_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(0, 2), doubles_res[Key::molecular_kinetic_energy_tensor_0_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(1, 0), doubles_res[Key::molecular_kinetic_energy_tensor_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(1, 1), doubles_res[Key::molecular_kinetic_energy_tensor_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(1, 2), doubles_res[Key::molecular_kinetic_energy_tensor_1_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(2, 0), doubles_res[Key::molecular_kinetic_energy_tensor_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(2, 1), doubles_res[Key::molecular_kinetic_energy_tensor_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().kinetic_energy_tensor(2, 2), doubles_res[Key::molecular_kinetic_energy_tensor_2_2_old], epsilon_);
}

}