/**
 * @file peptide_tutorial_eq_test.cc
 */

#include "peptide_tutorial_eq_test.h"

namespace testing {

Peptide_Tutorial_Eq_Test::Peptide_Tutorial_Eq_Test()
  : Simulation_Test(Parameter("Peptide Tutorial - Equilibration Test", "tutorial_test", "src/check/data/tutorial/eq")) {
  init_parameters();
  init_results();
}

void Peptide_Tutorial_Eq_Test::SetUp() {
  // initialize test simulation
  int err = test_sim_.init_simulation();
  ASSERT_EQ(err, 0) << "Initialization of the simulation unsuccessful. Error code: " << err;
}

void Peptide_Tutorial_Eq_Test::TearDown() {}


void Peptide_Tutorial_Eq_Test::init_parameters() {
  test_sim_.parameter().add_input("topo", "peptide_2Cl_54a7.top");
  test_sim_.parameter().add_input("input", "eq_peptide_5.imd");
  test_sim_.parameter().add_input("conf", "eq_peptide_4.cnf");
  test_sim_.parameter().add_input("refpos", "peptide_2Cl_h2o.rpr");
  test_sim_.parameter().add_input("posresspec", "peptide_2Cl_h2o.por");
  test_sim_.parameter().add_input("trc", "eq_peptide_5_test.trc");
  test_sim_.parameter().add_input("tre", "eq_peptide_5_test.tre");
  test_sim_.parameter().add_input("trv", "eq_peptide_5_test.trv");
  test_sim_.parameter().add_input("trf", "eq_peptide_5_test.trf");
  test_sim_.parameter().add_input("fin", "eq_peptide_5_test_final.cnf");
}

void Peptide_Tutorial_Eq_Test::init_results() {
  init_results_energies();
  init_results_forces();
  init_results_velocities();
  init_results_positions();

  init_results_baths();
  init_results_bonded_terms();
  init_results_nonbonded_terms();
  init_results_special_terms();

  init_results_mass();
  init_results_temperature();
  init_results_volume();
  init_results_pressure();
}

void Peptide_Tutorial_Eq_Test::init_results_parameters() {

}

void Peptide_Tutorial_Eq_Test::init_results_energies() {
  // "current" energies (step 1)
  results_.doubles_[Key::energy_total_current]                           =  -3.768597878e+04;
  results_.doubles_[Key::energy_kinetic_total_current]                   =   5.580929992e+03;
  results_.doubles_[Key::energy_potential_total_current]                 =  -4.326690877e+04;
  results_.doubles_[Key::energy_covalent_total_current]                  =   2.124197300e+02;
  results_.doubles_[Key::energy_bonds_total_current]                     =   0.000000000e+00;
  results_.doubles_[Key::energy_angles_total_current]                    =   1.144685558e+02;
  results_.doubles_[Key::energy_impropers_total_current]                 =   3.273666900e+01;
  results_.doubles_[Key::energy_dihedrals_total_current]                 =   6.521450510e+01;
  results_.doubles_[Key::energy_crossdihedrals_total_current]            =   0.000000000e+00;
  results_.doubles_[Key::energy_nonbonded_total_current]                 =  -4.347932850e+04;
  results_.doubles_[Key::energy_lennard_jones_total_current]             =   7.774977418e+03;
  results_.doubles_[Key::energy_coulomb_reaction_field_total_current]    =  -5.125430592e+04;
  results_.doubles_[Key::lattice_total_current]                          =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_pair_total_current]                 =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_space_total_current]           =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_k_reciprocal_space_total_current]   =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_A_term_total_current]               =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_self_total_current]                 =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_surface_total_current]              =   0.000000000e+00;
  results_.doubles_[Key::polarisation_self_total_current]                =   0.000000000e+00;
  results_.doubles_[Key::special_total_current]                          =   0.000000000e+00;
  results_.doubles_[Key::sasa_total_current]                             =   0.000000000e+00;
  results_.doubles_[Key::sasa_volume_total_current]                      =   0.000000000e+00;
  results_.doubles_[Key::constraints_total_current]                      =   0.000000000e+00;
  results_.doubles_[Key::distance_restraints_total_current]              =   0.000000000e+00;
  results_.doubles_[Key::distancefield_restraints_total_current]         =   0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_total_current]              =   0.000000000e+00;
  results_.doubles_[Key::position_restraints_total_current]              =   0.000000000e+00;
  results_.doubles_[Key::jvalue_restraints_total_current]                =   0.000000000e+00;
  results_.doubles_[Key::xray_restraints_total_current]                  =   0.000000000e+00;
  results_.doubles_[Key::local_elevation_total_current]                  =   0.000000000e+00;
  results_.doubles_[Key::order_parameter_restraints_total_current]       =   0.000000000e+00;
  results_.doubles_[Key::symmetry_restraints_total_current]              =   0.000000000e+00;
  results_.doubles_[Key::eds_vmix_current]                               =   0.000000000e+00;
  results_.doubles_[Key::eds_vr_current]                                 =   0.000000000e+00;
  results_.doubles_[Key::eds_emax_current]                               =   0.000000000e+00;
  results_.doubles_[Key::eds_emin_current]                               =   0.000000000e+00;
  results_.doubles_[Key::eds_glob_min_current]                           =   0.000000000e+00;
  results_.doubles_[Key::eds_glob_min_fluc_current]                      =   0.000000000e+00;
  results_.doubles_[Key::entropy_current]                                =   0.000000000e+00;
  results_.doubles_[Key::qmmm_total_current]                             =   0.000000000e+00;
  results_.doubles_[Key::bs_leus_energy_current]                         =   0.000000000e+00;
  results_.doubles_[Key::rdc_value_total_current]                        =   0.000000000e+00;
  results_.doubles_[Key::angle_restraints_total_current]                 =   0.000000000e+00;
  // "old" energies (step 2)
  results_.doubles_[Key::energy_total_old]                               =  -3.765699411e+04;
  results_.doubles_[Key::energy_kinetic_total_old]                       =   5.605194523e+03;
  results_.doubles_[Key::energy_potential_total_old]                     =  -4.326218863e+04;
  results_.doubles_[Key::energy_covalent_total_old]                      =   2.175196817e+02;
  results_.doubles_[Key::energy_bonds_total_old]                         =   0.000000000e+00;
  results_.doubles_[Key::energy_angles_total_old]                        =   1.182049469e+02;
  results_.doubles_[Key::energy_impropers_total_old]                     =   3.448947132e+01;
  results_.doubles_[Key::energy_dihedrals_total_old]                     =   6.482526348e+01;
  results_.doubles_[Key::energy_crossdihedrals_total_old]                =   0.000000000e+00;
  results_.doubles_[Key::energy_nonbonded_total_old]                     =  -4.347970831e+04;
  results_.doubles_[Key::energy_lennard_jones_total_old]                 =   7.765395115e+03;
  results_.doubles_[Key::energy_coulomb_reaction_field_total_old]        =  -5.124510343e+04;
  results_.doubles_[Key::lattice_total_old]                              =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_pair_total_old]                     =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_space_total_old]               =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_k_reciprocal_space_total_old]       =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_A_term_total_old]                   =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_self_total_old]                     =   0.000000000e+00;
  results_.doubles_[Key::lattice_sum_surface_total_old]                  =   0.000000000e+00;
  results_.doubles_[Key::polarisation_self_total_old]                    =   0.000000000e+00;
  results_.doubles_[Key::special_total_old]                              =   0.000000000e+00;
  results_.doubles_[Key::sasa_total_old]                                 =   0.000000000e+00;
  results_.doubles_[Key::sasa_volume_total_old]                          =   0.000000000e+00;
  results_.doubles_[Key::constraints_total_old]                          =   0.000000000e+00;
  results_.doubles_[Key::distance_restraints_total_old]                  =   0.000000000e+00;
  results_.doubles_[Key::distancefield_restraints_total_old]             =   0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_total_old]                  =   0.000000000e+00;
  results_.doubles_[Key::position_restraints_total_old]                  =   0.000000000e+00;
  results_.doubles_[Key::jvalue_restraints_total_old]                    =   0.000000000e+00;
  results_.doubles_[Key::xray_restraints_total_old]                      =   0.000000000e+00;
  results_.doubles_[Key::local_elevation_total_old]                      =   0.000000000e+00;
  results_.doubles_[Key::order_parameter_restraints_total_old]           =   0.000000000e+00;
  results_.doubles_[Key::symmetry_restraints_total_old]                  =   0.000000000e+00;
  results_.doubles_[Key::eds_vmix_old]                                   =   0.000000000e+00;
  results_.doubles_[Key::eds_vr_old]                                     =   0.000000000e+00;
  results_.doubles_[Key::eds_emax_old]                                   =   0.000000000e+00;
  results_.doubles_[Key::eds_emin_old]                                   =   0.000000000e+00;
  results_.doubles_[Key::eds_glob_min_old]                               =   0.000000000e+00;
  results_.doubles_[Key::eds_glob_min_fluc_old]                          =   0.000000000e+00;
  results_.doubles_[Key::entropy_old]                                    =   0.000000000e+00;
  results_.doubles_[Key::qmmm_total_old]                                 =   0.000000000e+00;
  results_.doubles_[Key::bs_leus_energy_old]                             =   0.000000000e+00;
  results_.doubles_[Key::rdc_value_total_old]                            =   0.000000000e+00;
  results_.doubles_[Key::angle_restraints_total_old]                     =   0.000000000e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_forces() {
  // "current" forces (step 1)
  results_.doubles_[Key::force_pos_0_0_current]   =     26.649770601;
  results_.doubles_[Key::force_pos_0_1_current]   =    126.750300610;
  results_.doubles_[Key::force_pos_0_2_current]   =    136.112932052;
  results_.doubles_[Key::force_pos_1_0_current]   =    -69.262104105;
  results_.doubles_[Key::force_pos_1_1_current]   =    221.725038245;
  results_.doubles_[Key::force_pos_1_2_current]   =    131.187727391;
  results_.doubles_[Key::force_pos_200_0_current] =   -170.863895259;
  results_.doubles_[Key::force_pos_200_1_current] =    751.316517420;
  results_.doubles_[Key::force_pos_200_2_current] =    504.976378288;
  results_.doubles_[Key::force_pos_210_0_current] =    410.981734710;
  results_.doubles_[Key::force_pos_210_1_current] =   -200.499712637;
  results_.doubles_[Key::force_pos_210_2_current] =   -200.694592731;
  // "old" forces (step 2)
  results_.doubles_[Key::force_pos_0_0_old]       =    145.260022387;
  results_.doubles_[Key::force_pos_0_1_old]       =   -230.630429692;
  results_.doubles_[Key::force_pos_0_2_old]       =    154.433461281;
  results_.doubles_[Key::force_pos_1_0_old]       =   -329.462223752;
  results_.doubles_[Key::force_pos_1_1_old]       =     68.652711047;
  results_.doubles_[Key::force_pos_1_2_old]       =   -165.619297732;
  results_.doubles_[Key::force_pos_200_0_old]     =   -132.130520048;
  results_.doubles_[Key::force_pos_200_1_old]     =    778.254833186;
  results_.doubles_[Key::force_pos_200_2_old]     =    538.634191517;
  results_.doubles_[Key::force_pos_210_0_old]     =    419.567599987;
  results_.doubles_[Key::force_pos_210_1_old]     =   -199.156399741;
  results_.doubles_[Key::force_pos_210_2_old]     =   -176.754424569;
}

void Peptide_Tutorial_Eq_Test::init_results_velocities() {
  // "current" forces (step 3) -> step not written out in test by default (STEP = 2)
  // consequence of the leap-frog algorithm
  results_.doubles_[Key::velocities_pos_0_0_current]   =   0;
  results_.doubles_[Key::velocities_pos_0_1_current]   =   0;
  results_.doubles_[Key::velocities_pos_0_2_current]   =   0;
  results_.doubles_[Key::velocities_pos_1_0_current]   =   0;
  results_.doubles_[Key::velocities_pos_1_1_current]   =   0;
  results_.doubles_[Key::velocities_pos_1_2_current]   =   0;
  results_.doubles_[Key::velocities_pos_200_0_current] =   0;
  results_.doubles_[Key::velocities_pos_200_1_current] =   0;
  results_.doubles_[Key::velocities_pos_200_2_current] =   0;
  results_.doubles_[Key::velocities_pos_210_0_current] =   0;
  results_.doubles_[Key::velocities_pos_210_1_current] =   0;
  results_.doubles_[Key::velocities_pos_210_2_current] =   0;
  // "old" velocities (step 2)
  results_.doubles_[Key::velocities_pos_0_0_old]       =   0;
  results_.doubles_[Key::velocities_pos_0_1_old]       =   0;
  results_.doubles_[Key::velocities_pos_0_2_old]       =   0;
  results_.doubles_[Key::velocities_pos_1_0_old]       =   0;
  results_.doubles_[Key::velocities_pos_1_1_old]       =   0;
  results_.doubles_[Key::velocities_pos_1_2_old]       =   0;
  results_.doubles_[Key::velocities_pos_200_0_old]     =   0;
  results_.doubles_[Key::velocities_pos_200_1_old]     =   0;
  results_.doubles_[Key::velocities_pos_200_2_old]     =   0;
  results_.doubles_[Key::velocities_pos_210_0_old]     =   0;
  results_.doubles_[Key::velocities_pos_210_1_old]     =   0;
  results_.doubles_[Key::velocities_pos_210_2_old]     =   0;
}

void Peptide_Tutorial_Eq_Test::init_results_positions() {
  // "current" positions (step 1)
  results_.doubles_[Key::positions_pos_0_0_current]   =    0.374727460;
  results_.doubles_[Key::positions_pos_0_1_current]   =    3.016802118;
  results_.doubles_[Key::positions_pos_0_2_current]   =    2.996245342;
  results_.doubles_[Key::positions_pos_1_0_current]   =    0.217979094;
  results_.doubles_[Key::positions_pos_1_1_current]   =    3.062091094;
  results_.doubles_[Key::positions_pos_1_2_current]   =    2.980292802;
  results_.doubles_[Key::positions_pos_200_0_current] =    1.204936057;
  results_.doubles_[Key::positions_pos_200_1_current] =    0.730794425;
  results_.doubles_[Key::positions_pos_200_2_current] =    0.319094146;
  results_.doubles_[Key::positions_pos_210_0_current] =    3.150721755;
  results_.doubles_[Key::positions_pos_210_1_current] =    1.749879698;
  results_.doubles_[Key::positions_pos_210_2_current] =    2.560643665;
  //// "old" positions (step 2)
  results_.doubles_[Key::positions_pos_0_0_old]       =    0.373467647;
  results_.doubles_[Key::positions_pos_0_1_old]       =    3.021811870;
  results_.doubles_[Key::positions_pos_0_2_old]       =    2.995563795;
  results_.doubles_[Key::positions_pos_1_0_old]       =    0.220417764;
  results_.doubles_[Key::positions_pos_1_1_old]       =    3.063919446;
  results_.doubles_[Key::positions_pos_1_2_old]       =    2.980627788;
  results_.doubles_[Key::positions_pos_200_0_old]     =    1.204395605;
  results_.doubles_[Key::positions_pos_200_1_old]     =    0.730659967;
  results_.doubles_[Key::positions_pos_200_2_old]     =    0.319599441;
  results_.doubles_[Key::positions_pos_210_0_old]     =    3.150233351;
  results_.doubles_[Key::positions_pos_210_1_old]     =    1.751560907;
  results_.doubles_[Key::positions_pos_210_2_old]     =    2.562971299;
}

void Peptide_Tutorial_Eq_Test::init_results_baths() {
  // "current" (step 1)
  results_.doubles_[Key::kinetic_total_bath_0_current]       = 1.567102052e+02;
  results_.doubles_[Key::kinetic_total_bath_1_current]       = 5.424219787e+03;
  results_.doubles_[Key::centre_of_mass_bath_0_current]      = 9.277884352e+00;
  results_.doubles_[Key::centre_of_mass_bath_1_current]      = 2.609116032e+03;
  results_.doubles_[Key::internal_rotational_bath_0_current] = 1.474323208e+02;
  results_.doubles_[Key::internal_rotational_bath_1_current] = 2.815103755e+03;
  // "old" (step 2)
  results_.doubles_[Key::kinetic_total_bath_0_old]           = 1.527815394e+02;
  results_.doubles_[Key::kinetic_total_bath_1_old]           = 5.452412983e+03;
  results_.doubles_[Key::centre_of_mass_bath_0_old]          = 9.468183733e+00;
  results_.doubles_[Key::centre_of_mass_bath_1_old]          = 2.627749251e+03;
  results_.doubles_[Key::internal_rotational_bath_0_old]     = 1.433133557e+02;
  results_.doubles_[Key::internal_rotational_bath_1_old]     = 2.824663732e+03;
}

void Peptide_Tutorial_Eq_Test::init_results_bonded_terms() {
  // "current" (step 1)
  results_.doubles_[Key::bond_energy_group_0_current]               = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_1_current]               = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_2_current]               = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_3_current]               = 0.000000000e+00;

  results_.doubles_[Key::angle_energy_group_0_current]              = 1.144685558e+02;
  results_.doubles_[Key::angle_energy_group_1_current]              = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_2_current]              = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_3_current]              = 0.000000000e+00;

  results_.doubles_[Key::improper_energy_group_0_current]           = 3.273666900e+01;
  results_.doubles_[Key::improper_energy_group_1_current]           = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_2_current]           = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_3_current]           = 0.000000000e+00;

  results_.doubles_[Key::dihedral_energy_group_0_current]           = 6.521450510e+01;
  results_.doubles_[Key::dihedral_energy_group_1_current]           = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_2_current]           = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_3_current]           = 0.000000000e+00;

  results_.doubles_[Key::crossdihedral_energy_group_0_current]      = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_1_current]      = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_2_current]      = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_3_current]      = 0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::bond_energy_group_0_old]                   = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_1_old]                   = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_2_old]                   = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_3_old]                   = 0.000000000e+00;

  results_.doubles_[Key::angle_energy_group_0_old]                  = 1.182049469e+02;
  results_.doubles_[Key::angle_energy_group_1_old]                  = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_2_old]                  = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_3_old]                  = 0.000000000e+00;

  results_.doubles_[Key::improper_energy_group_0_old]               = 3.448947132e+01;
  results_.doubles_[Key::improper_energy_group_1_old]               = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_2_old]               = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_3_old]               = 0.000000000e+00;

  results_.doubles_[Key::dihedral_energy_group_0_old]               = 6.482526348e+01;
  results_.doubles_[Key::dihedral_energy_group_1_old]               = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_2_old]               = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_3_old]               = 0.000000000e+00;

  results_.doubles_[Key::crossdihedral_energy_group_0_old]          = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_1_old]          = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_2_old]          = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_3_old]          = 0.000000000e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_nonbonded_terms() {
  // "current" (step 1)
  results_.doubles_[Key::lennard_jones_group_0_0_current]                =  -6.778133283e+01;
  results_.doubles_[Key::lennard_jones_group_1_0_current]                =   3.185093680e+00;
  results_.doubles_[Key::lennard_jones_group_2_0_current]                =   5.663016353e+00;
  results_.doubles_[Key::lennard_jones_group_3_0_current]                =  -1.173779820e+02;
  results_.doubles_[Key::lennard_jones_group_1_1_current]                =   0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_2_1_current]                =  -6.072268390e-03;
  results_.doubles_[Key::lennard_jones_group_3_1_current]                =   3.388879535e+01;
  results_.doubles_[Key::lennard_jones_group_2_2_current]                =   0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_3_2_current]                =   3.190641029e+01;
  results_.doubles_[Key::lennard_jones_group_3_3_current]                =   7.885499489e+03;

  results_.doubles_[Key::coulomb_reaction_field_group_0_0_current]       =  -7.019423964e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_1_0_current]       =  -2.734854419e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_2_0_current]       =  -3.153413192e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_3_0_current]       =  -1.254261627e+03;
  results_.doubles_[Key::coulomb_reaction_field_group_1_1_current]       =  -7.382455923e+01;
  results_.doubles_[Key::coulomb_reaction_field_group_2_1_current]       =   6.731830466e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_3_1_current]       =  -4.464844477e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_2_2_current]       =  -7.382455923e+01;
  results_.doubles_[Key::coulomb_reaction_field_group_3_2_current]       =  -3.329436188e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_3_3_current]       =  -4.778892978e+04;

  results_.doubles_[Key::lattice_sum_real_group_0_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_2_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_1_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_2_1_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_1_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_2_2_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_2_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_3_current]             =  0.000000000e+00;

  results_.doubles_[Key::lattice_sum_reciprocal_group_0_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_2_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_1_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_2_1_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_1_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_2_2_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_2_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_3_current]       =  0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::lennard_jones_group_0_0_old]                =  -6.768001068e+01;
  results_.doubles_[Key::lennard_jones_group_1_0_old]                =   3.192585454e+00;
  results_.doubles_[Key::lennard_jones_group_2_0_old]                =   5.404158079e+00;
  results_.doubles_[Key::lennard_jones_group_3_0_old]                =  -1.174791683e+02;
  results_.doubles_[Key::lennard_jones_group_1_1_old]                =   0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_2_1_old]                =  -6.072268390e-03;
  results_.doubles_[Key::lennard_jones_group_3_1_old]                =   3.399752090e+01;
  results_.doubles_[Key::lennard_jones_group_2_2_old]                =   0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_3_2_old]                =   3.141642517e+01;
  results_.doubles_[Key::lennard_jones_group_3_3_old]                =   7.876549677e+03;

  results_.doubles_[Key::coulomb_reaction_field_group_0_0_old]       =  -7.037190949e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_1_0_old]       =  -2.721257313e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_2_0_old]       =  -3.178743149e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_3_0_old]       =  -1.245086074e+03;
  results_.doubles_[Key::coulomb_reaction_field_group_1_1_old]       =  -7.382455923e+01;
  results_.doubles_[Key::coulomb_reaction_field_group_2_1_old]       =   6.731830466e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_3_1_old]       =  -4.523458248e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_2_2_old]       =  -7.382455923e+01;
  results_.doubles_[Key::coulomb_reaction_field_group_3_2_old]       =  -3.269145127e+02;
  results_.doubles_[Key::coulomb_reaction_field_group_3_3_old]       =  -4.778612059e+04;

  results_.doubles_[Key::lattice_sum_real_group_0_0_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_0_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_2_0_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_0_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_1_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_2_1_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_1_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_2_2_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_2_old]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_3_3_old]             =  0.000000000e+00;

  results_.doubles_[Key::lattice_sum_reciprocal_group_0_0_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_0_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_2_0_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_0_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_1_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_2_1_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_1_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_2_2_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_2_old]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_3_3_old]       =  0.000000000e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_special_terms() {
  // "current" (step 1)
  results_.doubles_[Key::contraints_group_0_current]          = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_0_current]      = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_0_current]     = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_0_current]        = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_0_current] = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_0_current]                = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_0_current]            = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_0_current]              = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_0_current]                 = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_0_current]     = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_0_current]       = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_0_current]     = 0.000000000e+00;
  results_.doubles_[Key::contraints_group_1_current]          = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_1_current]      = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_1_current]     = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_1_current]        = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_1_current] = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_1_current]                = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_1_current]            = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_1_current]              = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_1_current]                 = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_1_current]     = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_1_current]       = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_1_current]     = 0.000000000e+00;
  results_.doubles_[Key::contraints_group_2_current]          = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_2_current]      = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_2_current]     = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_2_current]        = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_2_current] = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_2_current]                = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_2_current]            = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_2_current]              = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_2_current]                 = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_2_current]     = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_2_current]       = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_2_current]     = 0.000000000e+00;
  results_.doubles_[Key::contraints_group_3_current]          = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_3_current]      = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_3_current]     = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_3_current]        = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_3_current] = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_3_current]                = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_3_current]            = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_3_current]              = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_3_current]                 = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_3_current]     = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_3_current]       = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_3_current]     = 0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::contraints_group_0_old]              = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_0_old]          = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_0_old]         = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_0_old]            = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_0_old]                    = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_0_old]                = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_0_old]                  = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_0_old]                     = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_0_old]         = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_0_old]           = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_0_old]         = 0.000000000e+00;
  results_.doubles_[Key::contraints_group_1_old]              = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_1_old]          = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_1_old]         = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_1_old]            = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_1_old]                    = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_1_old]                = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_1_old]                  = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_1_old]                     = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_1_old]         = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_1_old]           = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_1_old]         = 0.000000000e+00;
  results_.doubles_[Key::contraints_group_2_old]              = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_2_old]          = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_2_old]         = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_2_old]            = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_2_old]                    = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_2_old]                = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_2_old]                  = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_2_old]                     = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_2_old]         = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_2_old]           = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_2_old]         = 0.000000000e+00;
  results_.doubles_[Key::contraints_group_3_old]              = 0.000000000e+00;
  results_.doubles_[Key::pos_restraints_group_3_old]          = 0.000000000e+00;
  results_.doubles_[Key::dist_restraints_group_3_old]         = 0.000000000e+00;
  results_.doubles_[Key::disfield_res_group_3_old]            = 0.000000000e+00;
  results_.doubles_[Key::dihedral_restraints_group_3_old]     = 0.000000000e+00;
  results_.doubles_[Key::sasa_group_3_old]                    = 0.000000000e+00;
  results_.doubles_[Key::sasa_vol_group_3_old]                = 0.000000000e+00;
  results_.doubles_[Key::jvalue_group_3_old]                  = 0.000000000e+00;
  results_.doubles_[Key::rdc_group_3_old]                     = 0.000000000e+00;
  results_.doubles_[Key::local_elevation_group_3_old]         = 0.000000000e+00;
  results_.doubles_[Key::path_integral_group_3_old]           = 0.000000000e+00;
  results_.doubles_[Key::angle_restraint_group_3_old]         = 0.000000000e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_mass() {
  // "current" (step 1)
  results_.doubles_[Key::mass_current] = 1.712372440e+04;
  // "old" (step 2)
  results_.doubles_[Key::mass_old]     = 1.712372440e+04;
}

void Peptide_Tutorial_Eq_Test::init_results_temperature() {
  // "current" (step 1)
  results_.ints_[Key::num_temperature_coupling_baths_current]       = 2;
  results_.doubles_[Key::temperature_total_bath_0_current]          = 2.548396555e+02;
  results_.doubles_[Key::temperature_total_bath_1_current]          = 2.396247238e+02;
  results_.doubles_[Key::temperature_com_bath_0_current]            = 2.479733473e+02;
  results_.doubles_[Key::temperature_com_bath_1_current]            = 2.304012597e+02;
  results_.doubles_[Key::temperature_ir_bath_0_current]             = 2.551387720e+02;
  results_.doubles_[Key::temperature_ir_bath_1_current]             = 2.485912636e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_current] = 1.001508485e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_current] = 1.002539348e+00;
  // "old" (step 2)
  results_.ints_[Key::num_temperature_coupling_baths_old]           = 2;
  results_.doubles_[Key::temperature_total_bath_0_old]              = 2.484509215e+02;
  results_.doubles_[Key::temperature_total_bath_1_old]              = 2.408702092e+02;
  results_.doubles_[Key::temperature_com_bath_0_old]                = 2.530595472e+02;
  results_.doubles_[Key::temperature_com_bath_1_old]                = 2.320466894e+02;
  results_.doubles_[Key::temperature_ir_bath_0_old]                 = 2.480107033e+02;
  results_.doubles_[Key::temperature_ir_bath_1_old]                 = 2.494354693e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_old]     = 1.002044754e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_old]     = 1.002493556e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_volume() {
  // "current" (step 1)
  results_.doubles_[Key::volume_current]  = 2.893999506e+01;
  results_.doubles_[Key::box_k_0_current] = 3.070196349e+00;
  results_.doubles_[Key::box_k_1_current] = 0.000000000e+00;
  results_.doubles_[Key::box_k_2_current] = 0.000000000e+00;
  results_.doubles_[Key::box_l_0_current] = 0.000000000e+00;
  results_.doubles_[Key::box_l_1_current] = 3.070196349e+00;
  results_.doubles_[Key::box_l_2_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_0_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_1_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_2_current] = 3.070196349e+00;
  // "old" (step 2)
  results_.doubles_[Key::volume_old]      = 2.893999506e+01;
  results_.doubles_[Key::box_k_0_old]     = 3.070196349e+00;
  results_.doubles_[Key::box_k_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_k_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_l_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_l_1_old]     = 3.070196349e+00;
  results_.doubles_[Key::box_l_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_2_old]     = 3.070196349e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_pressure() {
  // "current" (step 1)
  results_.doubles_[Key::pressure_current]                            =  0.000000000e+00;
  results_.doubles_[Key::virial_current]                              =  7.105427358e-15;
  results_.doubles_[Key::molecular_kinetic_energy_current]            =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_0_0_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_0_1_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_0_2_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_1_0_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_1_1_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_1_2_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_2_0_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_2_1_current]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_2_2_current]                 =  0.000000000e+00;
  results_.doubles_[Key::virial_tensor_0_0_current]                   = -1.064009562e+02;
  results_.doubles_[Key::virial_tensor_0_1_current]                   =  6.617571422e+01;
  results_.doubles_[Key::virial_tensor_0_2_current]                   =  2.514787226e+01;
  results_.doubles_[Key::virial_tensor_1_0_current]                   =  6.617571422e+01;
  results_.doubles_[Key::virial_tensor_1_1_current]                   =  4.446228639e+01;
  results_.doubles_[Key::virial_tensor_1_2_current]                   = -3.401802863e+01;
  results_.doubles_[Key::virial_tensor_2_0_current]                   =  2.514787226e+01;
  results_.doubles_[Key::virial_tensor_2_1_current]                   = -3.401802863e+01;
  results_.doubles_[Key::virial_tensor_2_2_current]                   =  6.193866985e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_current] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_current] =  0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::pressure_old]                            =  0.000000000e+00;
  results_.doubles_[Key::virial_old]                              = -1.421085472e-14;
  results_.doubles_[Key::molecular_kinetic_energy_old]            =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_0_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_0_1_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_0_2_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_1_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_1_1_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_1_2_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_2_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_2_1_old]                 =  0.000000000e+00;
  results_.doubles_[Key::pressure_tensor_2_2_old]                 =  0.000000000e+00;
  results_.doubles_[Key::virial_tensor_0_0_old]                   = -1.286827355e+02;
  results_.doubles_[Key::virial_tensor_0_1_old]                   =  9.400657447e+01;
  results_.doubles_[Key::virial_tensor_0_2_old]                   =  1.579647032e+00;
  results_.doubles_[Key::virial_tensor_1_0_old]                   =  9.400657447e+01;
  results_.doubles_[Key::virial_tensor_1_1_old]                   =  6.373544397e+01;
  results_.doubles_[Key::virial_tensor_1_2_old]                   = -3.618618363e+01;
  results_.doubles_[Key::virial_tensor_2_0_old]                   =  1.579647032e+00;
  results_.doubles_[Key::virial_tensor_2_1_old]                   = -3.618618363e+01;
  results_.doubles_[Key::virial_tensor_2_2_old]                   =  6.494729151e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_old] =  0.000000000e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_old] =  0.000000000e+00;
}

void Peptide_Tutorial_Eq_Test::check_parameter_init() {

}

void Peptide_Tutorial_Eq_Test::check_simulation_results_bonded_terms() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_energy[0], doubles_res[Key::bond_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_energy[1], doubles_res[Key::bond_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_energy[2], doubles_res[Key::bond_energy_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.bond_energy[3], doubles_res[Key::bond_energy_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_energy[0], doubles_res[Key::angle_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_energy[1], doubles_res[Key::angle_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_energy[2], doubles_res[Key::angle_energy_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angle_energy[3], doubles_res[Key::angle_energy_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_energy[0], doubles_res[Key::improper_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_energy[1], doubles_res[Key::improper_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_energy[2], doubles_res[Key::improper_energy_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.improper_energy[3], doubles_res[Key::improper_energy_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_energy[0], doubles_res[Key::dihedral_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_energy[1], doubles_res[Key::dihedral_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_energy[2], doubles_res[Key::dihedral_energy_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihedral_energy[3], doubles_res[Key::dihedral_energy_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_energy[0], doubles_res[Key::crossdihedral_energy_group_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_energy[1], doubles_res[Key::crossdihedral_energy_group_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_energy[2], doubles_res[Key::crossdihedral_energy_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crossdihedral_energy[3], doubles_res[Key::crossdihedral_energy_group_3_current], epsilon_);
  // "old" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_energy[0], doubles_res[Key::bond_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_energy[1], doubles_res[Key::bond_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_energy[2], doubles_res[Key::bond_energy_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.bond_energy[3], doubles_res[Key::bond_energy_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_energy[0], doubles_res[Key::angle_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_energy[1], doubles_res[Key::angle_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_energy[2], doubles_res[Key::angle_energy_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angle_energy[3], doubles_res[Key::angle_energy_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_energy[0], doubles_res[Key::improper_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_energy[1], doubles_res[Key::improper_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_energy[2], doubles_res[Key::improper_energy_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.improper_energy[3], doubles_res[Key::improper_energy_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_energy[0], doubles_res[Key::dihedral_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_energy[1], doubles_res[Key::dihedral_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_energy[2], doubles_res[Key::dihedral_energy_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihedral_energy[3], doubles_res[Key::dihedral_energy_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_energy[0], doubles_res[Key::crossdihedral_energy_group_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_energy[1], doubles_res[Key::crossdihedral_energy_group_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_energy[2], doubles_res[Key::crossdihedral_energy_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crossdihedral_energy[3], doubles_res[Key::crossdihedral_energy_group_3_old], epsilon_);
}

void Peptide_Tutorial_Eq_Test::check_simulation_results_nonbonded_terms() {
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // "current" energies
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[0][0], doubles_res[Key::lennard_jones_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[1][0], doubles_res[Key::lennard_jones_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[2][0], doubles_res[Key::lennard_jones_group_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[3][0], doubles_res[Key::lennard_jones_group_3_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[1][1], doubles_res[Key::lennard_jones_group_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[2][1], doubles_res[Key::lennard_jones_group_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[3][1], doubles_res[Key::lennard_jones_group_3_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[2][2], doubles_res[Key::lennard_jones_group_2_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[3][2], doubles_res[Key::lennard_jones_group_3_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.lj_energy[3][3], doubles_res[Key::lennard_jones_group_3_3_current], epsilon_);

  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[0][0], doubles_res[Key::coulomb_reaction_field_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[1][0], doubles_res[Key::coulomb_reaction_field_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[2][0], doubles_res[Key::coulomb_reaction_field_group_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[3][0], doubles_res[Key::coulomb_reaction_field_group_3_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[1][1], doubles_res[Key::coulomb_reaction_field_group_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[2][1], doubles_res[Key::coulomb_reaction_field_group_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[3][1], doubles_res[Key::coulomb_reaction_field_group_3_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[2][2], doubles_res[Key::coulomb_reaction_field_group_2_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[3][2], doubles_res[Key::coulomb_reaction_field_group_3_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.crf_energy[3][3], doubles_res[Key::coulomb_reaction_field_group_3_3_current], epsilon_);

  // lattice sum energies cannot be tested like this as they are not calculated right now. values in .tre file correspond to totals
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_2_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_3_current], epsilon_);

  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[0][0], doubles_res[Key::lattice_sum_reciprocal_group_0_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[1][0], doubles_res[Key::lattice_sum_reciprocal_group_1_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[2][0], doubles_res[Key::lattice_sum_reciprocal_group_2_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[3][0], doubles_res[Key::lattice_sum_reciprocal_group_3_0_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[1][1], doubles_res[Key::lattice_sum_reciprocal_group_1_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[2][1], doubles_res[Key::lattice_sum_reciprocal_group_2_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[3][1], doubles_res[Key::lattice_sum_reciprocal_group_3_1_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[2][2], doubles_res[Key::lattice_sum_reciprocal_group_2_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[3][2], doubles_res[Key::lattice_sum_reciprocal_group_3_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.ls_k_energy[3][3], doubles_res[Key::lattice_sum_reciprocal_group_3_3_current], epsilon_); 
  // "old" energies
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[0][0], doubles_res[Key::lennard_jones_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[1][0], doubles_res[Key::lennard_jones_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[2][0], doubles_res[Key::lennard_jones_group_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[3][0], doubles_res[Key::lennard_jones_group_3_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[1][1], doubles_res[Key::lennard_jones_group_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[2][1], doubles_res[Key::lennard_jones_group_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[3][1], doubles_res[Key::lennard_jones_group_3_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[2][2], doubles_res[Key::lennard_jones_group_2_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[3][2], doubles_res[Key::lennard_jones_group_3_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.lj_energy[3][3], doubles_res[Key::lennard_jones_group_3_3_old], epsilon_);

  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[0][0], doubles_res[Key::coulomb_reaction_field_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[1][0], doubles_res[Key::coulomb_reaction_field_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[2][0], doubles_res[Key::coulomb_reaction_field_group_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[3][0], doubles_res[Key::coulomb_reaction_field_group_3_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[1][1], doubles_res[Key::coulomb_reaction_field_group_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[2][1], doubles_res[Key::coulomb_reaction_field_group_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[3][1], doubles_res[Key::coulomb_reaction_field_group_3_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[2][2], doubles_res[Key::coulomb_reaction_field_group_2_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[3][2], doubles_res[Key::coulomb_reaction_field_group_3_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.crf_energy[3][3], doubles_res[Key::coulomb_reaction_field_group_3_3_old], epsilon_);

  // lattice sum energies cannot be tested like this as they are not calculated right now. values in .tre file correspond to totals
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_2_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_total, doubles_res[Key::lattice_sum_real_group_3_3_old], epsilon_);

  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[0][0], doubles_res[Key::lattice_sum_reciprocal_group_0_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[1][0], doubles_res[Key::lattice_sum_reciprocal_group_1_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[2][0], doubles_res[Key::lattice_sum_reciprocal_group_2_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[3][0], doubles_res[Key::lattice_sum_reciprocal_group_3_0_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[1][1], doubles_res[Key::lattice_sum_reciprocal_group_1_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[2][1], doubles_res[Key::lattice_sum_reciprocal_group_2_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[3][1], doubles_res[Key::lattice_sum_reciprocal_group_3_1_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[2][2], doubles_res[Key::lattice_sum_reciprocal_group_2_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[3][2], doubles_res[Key::lattice_sum_reciprocal_group_3_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.ls_k_energy[3][3], doubles_res[Key::lattice_sum_reciprocal_group_3_3_old], epsilon_);
}

void Peptide_Tutorial_Eq_Test::check_simulation_results_special_terms() {
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

  EXPECT_NEAR(test_sim_.conf().current().energies.constraints_energy[2], doubles_res[Key::contraints_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.posrest_energy[2], doubles_res[Key::pos_restraints_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.distanceres_energy[2], doubles_res[Key::dist_restraints_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.disfieldres_energy[2], doubles_res[Key::disfield_res_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angrest_energy[2], doubles_res[Key::angle_restraint_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihrest_energy[2], doubles_res[Key::dihedral_restraints_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_energy[2], doubles_res[Key::sasa_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_volume_energy[2], doubles_res[Key::sasa_vol_group_2_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_2_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.rdc_energy[2], doubles_res[Key::rdc_group_2_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_2_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_2_current], epsilon_);
  
  EXPECT_NEAR(test_sim_.conf().current().energies.constraints_energy[3], doubles_res[Key::contraints_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.posrest_energy[3], doubles_res[Key::pos_restraints_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.distanceres_energy[3], doubles_res[Key::dist_restraints_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.disfieldres_energy[3], doubles_res[Key::disfield_res_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.angrest_energy[3], doubles_res[Key::angle_restraint_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.dihrest_energy[3], doubles_res[Key::dihedral_restraints_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_energy[3], doubles_res[Key::sasa_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.sasa_volume_energy[3], doubles_res[Key::sasa_vol_group_3_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_3_current], epsilon_);
  EXPECT_NEAR(test_sim_.conf().current().energies.rdc_energy[3], doubles_res[Key::rdc_group_3_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_3_current], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_3_current], epsilon_);

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

  EXPECT_NEAR(test_sim_.conf().old().energies.constraints_energy[2], doubles_res[Key::contraints_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.posrest_energy[2], doubles_res[Key::pos_restraints_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.distanceres_energy[2], doubles_res[Key::dist_restraints_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.disfieldres_energy[2], doubles_res[Key::disfield_res_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angrest_energy[2], doubles_res[Key::angle_restraint_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihrest_energy[2], doubles_res[Key::dihedral_restraints_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_energy[2], doubles_res[Key::sasa_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_volume_energy[2], doubles_res[Key::sasa_vol_group_2_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_2_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.rdc_energy[2], doubles_res[Key::rdc_group_2_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_2_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_2_old], epsilon_);
  
  EXPECT_NEAR(test_sim_.conf().old().energies.constraints_energy[3], doubles_res[Key::contraints_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.posrest_energy[3], doubles_res[Key::pos_restraints_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.distanceres_energy[3], doubles_res[Key::dist_restraints_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.disfieldres_energy[3], doubles_res[Key::disfield_res_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.angrest_energy[3], doubles_res[Key::angle_restraint_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.dihrest_energy[3], doubles_res[Key::dihedral_restraints_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_energy[3], doubles_res[Key::sasa_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.sasa_volume_energy[3], doubles_res[Key::sasa_vol_group_3_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::jvalue_group_3_old], epsilon_);
  EXPECT_NEAR(test_sim_.conf().old().energies.rdc_energy[3], doubles_res[Key::rdc_group_3_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::local_elevation_group_3_old], epsilon_);
  EXPECT_NEAR(0, doubles_res[Key::path_integral_group_3_old], epsilon_);
}

TEST_F(Peptide_Tutorial_Eq_Test, check_init) {
  // check if input vile has been loaded correctly
  check_parameter_init(); 
}

TEST_F(Peptide_Tutorial_Eq_Test, check_simulation) {
  // run the simulation specified by the parameters object
  test_sim_.run_simulation();
  check_simulation_results();
}


} // namespace testing