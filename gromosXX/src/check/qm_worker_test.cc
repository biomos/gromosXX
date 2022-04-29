/**
 * @file qm_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"

#include "qm_worker_test.h"

#include "../interaction/qmmm/qmmm_interaction.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

QM_Worker_Test::QM_Worker_Test(const Parameter& parameter, const Results& results) : test_sim_(parameter), results_(results) {}

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

}