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
  results_.doubles_[Key::kinetic_total_bath_0_current]       = 3.330665706e+02;
  results_.doubles_[Key::kinetic_total_bath_1_current]       = 7.393007732e+03;
  results_.doubles_[Key::centre_of_mass_bath_0_current]      = 4.160893165e+00;
  results_.doubles_[Key::centre_of_mass_bath_1_current]      = 1.184185934e+03;
  results_.doubles_[Key::internal_rotational_bath_0_current] = 3.289056774e+02;
  results_.doubles_[Key::internal_rotational_bath_1_current] = 6.208821799e+03;
  // "old" (step 2)
  results_.doubles_[Key::kinetic_total_bath_0_old]           = 4.305606395e+02;
  results_.doubles_[Key::kinetic_total_bath_1_old]           = 7.366368543e+03;
  results_.doubles_[Key::centre_of_mass_bath_0_old]          = 4.134134679e+00;
  results_.doubles_[Key::centre_of_mass_bath_1_old]          = 1.185360778e+03;
  results_.doubles_[Key::internal_rotational_bath_0_old]     = 4.264265048e+02;
  results_.doubles_[Key::internal_rotational_bath_1_old]     = 6.181007764e+03;
}

void Peptide_Tutorial_Eq_Test::init_results_bonded_terms() {
  // "current" (step 1)
  results_.doubles_[Key::bond_energy_group_0_current]               = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_1_current]               = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_0_current]              = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_1_current]              = 2.781467266e+03;
  results_.doubles_[Key::improper_energy_group_0_current]           = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_1_current]           = 1.247690747e+03;
  results_.doubles_[Key::dihedral_energy_group_0_current]           = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_1_current]           = 6.857307867e+02;
  results_.doubles_[Key::crossdihedral_energy_group_0_current]      = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_1_current]      = 0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::bond_energy_group_0_old]                   = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_1_old]                   = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_0_old]                  = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_1_old]                  = 2.812642458e+03;
  results_.doubles_[Key::improper_energy_group_0_old]               = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_1_old]               = 1.242791317e+03;
  results_.doubles_[Key::dihedral_energy_group_0_old]               = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_1_old]               = 6.879624537e+02;
  results_.doubles_[Key::crossdihedral_energy_group_0_old]          = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_1_old]          = 0.000000000e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_nonbonded_terms() {
  // "current" (step 1)
  results_.doubles_[Key::lennard_jones_group_0_0_current]                =  0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_1_0_current]                = -1.410201839e+02;
  results_.doubles_[Key::lennard_jones_group_1_1_current]                = -4.868703494e+03;
  results_.doubles_[Key::coulomb_reaction_field_group_0_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_1_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_1_1_current]       = -2.559496780e+04;
  results_.doubles_[Key::lattice_sum_real_group_0_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_1_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_0_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_1_current]       =  0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::lennard_jones_group_0_0_old]                    =  0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_1_0_old]                    = -1.411110293e+02;
  results_.doubles_[Key::lennard_jones_group_1_1_old]                    = -4.869270312e+03;
  results_.doubles_[Key::coulomb_reaction_field_group_0_0_old]           =  0.000000000e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_1_0_old]           =  0.000000000e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_1_1_old]           = -2.559545299e+04;
  results_.doubles_[Key::lattice_sum_real_group_0_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_1_old]                 =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_0_0_old]           =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_0_old]           =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_1_old]           =  0.000000000e+00;
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
}

void Peptide_Tutorial_Eq_Test::init_results_mass() {
  // "current" (step 1)
  results_.doubles_[Key::mass_current] = 1.991781000e+04;
  // "old" (step 2)
  results_.doubles_[Key::mass_old]     = 1.991781000e+04;
}

void Peptide_Tutorial_Eq_Test::init_results_temperature() {
  // "current" (step 1)
  results_.ints_[Key::num_temperature_coupling_baths_current]       = 2;
  results_.doubles_[Key::temperature_total_bath_0_current]          = 2.968677248e+02;
  results_.doubles_[Key::temperature_total_bath_1_current]          = 2.801836889e+02;
  results_.doubles_[Key::temperature_com_bath_0_current]            = 3.336222763e+02;
  results_.doubles_[Key::temperature_com_bath_1_current]            = 3.798019531e+02;
  results_.doubles_[Key::temperature_ir_bath_0_current]             = 2.963187086e+02;
  results_.doubles_[Key::temperature_ir_bath_1_current]             = 2.666979806e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_current] = 1.000019771e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_current] = 1.000014899e+00;
  // "old" (step 2)
  results_.ints_[Key::num_temperature_coupling_baths_old]           = 2;
  results_.doubles_[Key::temperature_total_bath_0_old]              = 3.837658875e+02;
  results_.doubles_[Key::temperature_total_bath_1_old]              = 2.791741040e+02;
  results_.doubles_[Key::temperature_com_bath_0_old]                = 3.314704712e+02;
  results_.doubles_[Key::temperature_com_bath_1_old]                = 3.801787590e+02;
  results_.doubles_[Key::temperature_ir_bath_0_old]                 = 3.841776103e+02;
  results_.doubles_[Key::temperature_ir_bath_1_old]                 = 2.655032371e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_old]     = 1.000023732e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_old]     = 1.000022414e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_volume() {
  // "current" (step 1)
  results_.doubles_[Key::volume_current]  = 5.091972755e+01;
  results_.doubles_[Key::box_k_0_current] = 3.706483096e+00;
  results_.doubles_[Key::box_k_1_current] = 0.000000000e+00;
  results_.doubles_[Key::box_k_2_current] = 0.000000000e+00;
  results_.doubles_[Key::box_l_0_current] = 0.000000000e+00;
  results_.doubles_[Key::box_l_1_current] = 3.706483096e+00;
  results_.doubles_[Key::box_l_2_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_0_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_1_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_2_current] = 3.706483096e+00;
  // "old" (step 2)
  results_.doubles_[Key::volume_old]      = 5.091976752e+01;
  results_.doubles_[Key::box_k_0_old]     = 3.706484065e+00;
  results_.doubles_[Key::box_k_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_k_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_l_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_l_1_old]     = 3.706484065e+00;
  results_.doubles_[Key::box_l_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_2_old]     = 3.706484065e+00;
}

void Peptide_Tutorial_Eq_Test::init_results_pressure() {
  // "current" (step 1)
  results_.doubles_[Key::pressure_current]                            =  1.774560087e+00;
  results_.doubles_[Key::virial_current]                              =  3.507492242e+02;
  results_.doubles_[Key::molecular_kinetic_energy_current]            =  3.959292823e+02;
  results_.doubles_[Key::pressure_tensor_0_0_current]                 = -3.087931850e+00;
  results_.doubles_[Key::pressure_tensor_0_1_current]                 =  8.257587370e+00;
  results_.doubles_[Key::pressure_tensor_0_2_current]                 = -1.101848872e+01;
  results_.doubles_[Key::pressure_tensor_1_0_current]                 =  1.307025805e+01;
  results_.doubles_[Key::pressure_tensor_1_1_current]                 = -2.212711333e+00;
  results_.doubles_[Key::pressure_tensor_1_2_current]                 = -1.464952771e+01;
  results_.doubles_[Key::pressure_tensor_2_0_current]                 = -5.037615305e+00;
  results_.doubles_[Key::pressure_tensor_2_1_current]                 = -6.298825539e+00;
  results_.doubles_[Key::pressure_tensor_2_2_current]                 =  1.062432344e+01;
  results_.doubles_[Key::virial_tensor_0_0_current]                   =  5.050716577e+02;
  results_.doubles_[Key::virial_tensor_0_1_current]                   = -1.937614474e+02;
  results_.doubles_[Key::virial_tensor_0_2_current]                   =  2.745860522e+02;
  results_.doubles_[Key::virial_tensor_1_0_current]                   = -3.162913873e+02;
  results_.doubles_[Key::virial_tensor_1_1_current]                   =  4.637087300e+02;
  results_.doubles_[Key::virial_tensor_1_2_current]                   =  3.558157852e+02;
  results_.doubles_[Key::virial_tensor_2_0_current]                   =  1.223138298e+02;
  results_.doubles_[Key::virial_tensor_2_1_current]                   =  1.432080454e+02;
  results_.doubles_[Key::virial_tensor_2_2_current]                   =  8.346728497e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_current] =  4.264533335e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_current] =  1.647560216e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_current] = -5.943169643e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_current] =  1.647560216e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_current] =  4.073734009e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_current] = -1.715919472e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_current] = -5.943169643e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_current] = -1.715919472e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_current] =  3.539611126e+02;
  // "old" (step 2)
  results_.doubles_[Key::pressure_old]                            =  1.581469292e+00;
  results_.doubles_[Key::virial_old]                              =  3.560431584e+02;
  results_.doubles_[Key::molecular_kinetic_energy_old]            =  3.963071828e+02;
  results_.doubles_[Key::pressure_tensor_0_0_old]                 = -3.120251261e+00;
  results_.doubles_[Key::pressure_tensor_0_1_old]                 =  8.399222892e+00;
  results_.doubles_[Key::pressure_tensor_0_2_old]                 = -1.109642505e+01;
  results_.doubles_[Key::pressure_tensor_1_0_old]                 =  1.295917267e+01;
  results_.doubles_[Key::pressure_tensor_1_1_old]                 = -2.700620433e+00;
  results_.doubles_[Key::pressure_tensor_1_2_old]                 = -1.490663574e+01;
  results_.doubles_[Key::pressure_tensor_2_0_old]                 = -5.347501517e+00;
  results_.doubles_[Key::pressure_tensor_2_1_old]                 = -6.223445946e+00;
  results_.doubles_[Key::pressure_tensor_2_2_old]                 =  1.056527957e+01;
  results_.doubles_[Key::virial_tensor_0_0_old]                   =  5.061933042e+02;
  results_.doubles_[Key::virial_tensor_0_1_old]                   = -1.971145922e+02;
  results_.doubles_[Key::virial_tensor_0_2_old]                   =  2.762657115e+02;
  results_.doubles_[Key::virial_tensor_1_0_old]                   = -3.132103835e+02;
  results_.doubles_[Key::virial_tensor_1_1_old]                   =  4.764302609e+02;
  results_.doubles_[Key::virial_tensor_1_2_old]                   =  3.623069342e+02;
  results_.doubles_[Key::virial_tensor_2_0_old]                   =  1.298987867e+02;
  results_.doubles_[Key::virial_tensor_2_1_old]                   =  1.412339313e+02;
  results_.doubles_[Key::virial_tensor_2_2_old]                   =  8.550591018e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_old] =  4.267520698e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_old] =  1.672864624e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_old] = -6.247980288e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_old] =  1.672864624e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_old] =  4.076727786e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_old] = -1.721427904e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_old] = -6.247980288e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_old] = -1.721427904e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_old] =  3.544966999e+02;
}

void Peptide_Tutorial_Eq_Test::check_parameter_init() {

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