/**
 * @file orca_worker_test_electrostatic.cc
 */

#include <gtest/gtest.h>

#include "orca_worker_test_electrostatic.h"
#include "../interaction/qmmm/orca_worker.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

Orca_Worker_Test_Electrostatic::Orca_Worker_Test_Electrostatic() : Orca_Worker_Test("Orca QM/MM Tests Electrostatic Embedding") {
  // initialize test files and expected results
  init_parameters();
  init_results();
}

void Orca_Worker_Test_Electrostatic::init_parameters() {
  test_sim_.parameter().add_input("topo", "md.top");
  test_sim_.parameter().add_input("input", "md_el.imd");
  test_sim_.parameter().add_input("conf", "md.cnf");
  test_sim_.parameter().add_input("qmmm", "md_el.qmmm");
  test_sim_.parameter().add_input("trc", "md_el.trc");
  test_sim_.parameter().add_input("tre", "md_el.tre");
  test_sim_.parameter().add_input("trv", "md_el.trv");
  test_sim_.parameter().add_input("trf", "md_el.trf");
  test_sim_.parameter().add_input("fin", "md_el_final.cnf");
}

void Orca_Worker_Test_Electrostatic::init_results() {
  // call helper functions
  // first test (initialization)
  init_results_parameters();
  init_results_binary_name();
  init_results_units();
  init_results_elements();
  init_results_files();
  init_results_qm_zone_init();
  // second test (energy comparison)
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

void Orca_Worker_Test_Electrostatic::init_results_parameters() {
  // parameters from the input file
  results_.ints_[Key::qmmm] = simulation::qmmm_electrostatic;
  results_.ints_[Key::qm_software] = simulation::qm_orca;
  results_.ints_[Key::qm_ch] = simulation::qm_ch_constant;
  results_.doubles_[Key::qm_cutoff] = 1.2;
  results_.ints_[Key::qm_lj] = simulation::qm_lj_off;
  results_.ints_[Key::qm_constraint] = simulation::qm_constr_off;
  results_.doubles_[Key::qm_mm_scale] = -1.0;
}

void Orca_Worker_Test_Electrostatic::init_results_files() {
  // input file names
  results_.strs_[Key::input_file] = "xphos_el.inp";
  results_.strs_[Key::input_coordinate_file] = "xphos_el.xyz";
  results_.strs_[Key::input_pointcharges_file] = "xphos_el.pc";
  results_.strs_[Key::output_file] = "xphos_el.out";
  results_.strs_[Key::output_gradient_file] = "xphos_el.engrad";
  results_.strs_[Key::output_mm_gradient_file] = "xphos_el.pcgrad";
}

void Orca_Worker_Test_Electrostatic::init_results_qm_zone_init() {
  // QM zone
  results_.ints_[Key::qm_zone_size_init] = 90;
  results_.ints_[Key::mm_zone_size_init] = 1194;
  results_.ints_[Key::qm_zone_charge] = 0;
  results_.ints_[Key::qm_zone_spin_mult] = 1;
  results_.doubles_[Key::qm_energy_init] = 0.0;
  results_.doubles_[Key::qm_atom_charge_init] = 0.0;
  results_.doubles_[Key::qm_atom_force_init] = 0.0;
  results_.doubles_[Key::mm_atom_force_init] = 0.0;
  results_.doubles_[Key::mm_atom_cos_force_init] = 0.0;
}

void Orca_Worker_Test_Electrostatic::init_results_energies() {
  // "current" energies (step 1)
  results_.doubles_[Key::energy_total_current]                           =  -9.426338718e+05;
  results_.doubles_[Key::energy_kinetic_total_current]                   =   7.726074303e+03;
  results_.doubles_[Key::energy_potential_total_current]                 =  -9.503599461e+05;
  results_.doubles_[Key::energy_covalent_total_current]                  =   4.714888800e+03;
  results_.doubles_[Key::energy_bonds_total_current]                     =   0.000000000e+00;
  results_.doubles_[Key::energy_angles_total_current]                    =   2.781467266e+03;
  results_.doubles_[Key::energy_impropers_total_current]                 =   1.247690747e+03;
  results_.doubles_[Key::energy_dihedrals_total_current]                 =   6.857307867e+02;
  results_.doubles_[Key::energy_crossdihedrals_total_current]            =   0.000000000e+00;
  results_.doubles_[Key::energy_nonbonded_total_current]                 =  -3.060469148e+04;
  results_.doubles_[Key::energy_lennard_jones_total_current]             =  -5.009723678e+03;
  results_.doubles_[Key::energy_coulomb_reaction_field_total_current]    =  -2.559496780e+04;
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
  results_.doubles_[Key::qmmm_total_current]                             =  -9.244701434e+05;
  results_.doubles_[Key::bs_leus_energy_current]                         =   0.000000000e+00;
  results_.doubles_[Key::rdc_value_total_current]                        =   0.000000000e+00;
  results_.doubles_[Key::angle_restraints_total_current]                 =   0.000000000e+00;
  // "old" energies (step 2)
  results_.doubles_[Key::energy_total_old]                               =  -9.426383998e+05;
  results_.doubles_[Key::energy_kinetic_total_old]                       =   7.796929182e+03;
  results_.doubles_[Key::energy_potential_total_old]                     =  -9.504353290e+05;
  results_.doubles_[Key::energy_covalent_total_old]                      =   4.743396228e+03;
  results_.doubles_[Key::energy_bonds_total_old]                         =   0.000000000e+00;
  results_.doubles_[Key::energy_angles_total_old]                        =   2.812642458e+03;
  results_.doubles_[Key::energy_impropers_total_old]                     =   1.242791317e+03;
  results_.doubles_[Key::energy_dihedrals_total_old]                     =   6.879624537e+02;
  results_.doubles_[Key::energy_crossdihedrals_total_old]                =   0.000000000e+00;
  results_.doubles_[Key::energy_nonbonded_total_old]                     =  -3.060583433e+04;
  results_.doubles_[Key::energy_lennard_jones_total_old]                 =  -5.010381342e+03;
  results_.doubles_[Key::energy_coulomb_reaction_field_total_old]        =  -2.559545299e+04;
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
  results_.doubles_[Key::qmmm_total_old]                                 =  -9.245728909e+05;
  results_.doubles_[Key::bs_leus_energy_old]                             =   0.000000000e+00;
  results_.doubles_[Key::rdc_value_total_old]                            =   0.000000000e+00;
  results_.doubles_[Key::angle_restraints_total_old]                     =   0.000000000e+00;
}

void Orca_Worker_Test_Electrostatic::init_results_forces() {
  // "current" forces (step 1)
  results_.doubles_[Key::force_pos_0_0_current]   =    436.158347386;
  results_.doubles_[Key::force_pos_0_1_current]   =   -372.599172224;
  results_.doubles_[Key::force_pos_0_2_current]   =   -357.771249851;
  results_.doubles_[Key::force_pos_1_0_current]   =   -853.160469825;
  results_.doubles_[Key::force_pos_1_1_current]   =   -642.435169900;
  results_.doubles_[Key::force_pos_1_2_current]   =   3652.335764283;
  results_.doubles_[Key::force_pos_200_0_current] =   -169.922327428;
  results_.doubles_[Key::force_pos_200_1_current] =   -120.092514533;
  results_.doubles_[Key::force_pos_200_2_current] =    161.500127483;
  results_.doubles_[Key::force_pos_210_0_current] =    239.213725210;
  results_.doubles_[Key::force_pos_210_1_current] =   -185.061101920;
  results_.doubles_[Key::force_pos_210_2_current] =    -51.118786278;
  // "old" forces (step 2)
  results_.doubles_[Key::force_pos_0_0_old]       =    322.794611469;
  results_.doubles_[Key::force_pos_0_1_old]       =    -61.321755176;
  results_.doubles_[Key::force_pos_0_2_old]       =   -287.529491621;
  results_.doubles_[Key::force_pos_1_0_old]       =   -894.924699449;
  results_.doubles_[Key::force_pos_1_1_old]       =   -631.250853775;
  results_.doubles_[Key::force_pos_1_2_old]       =   3173.369707558;
  results_.doubles_[Key::force_pos_200_0_old]     =   -175.874840911;
  results_.doubles_[Key::force_pos_200_1_old]     =   -130.311468863;
  results_.doubles_[Key::force_pos_200_2_old]     =    170.904652248;
  results_.doubles_[Key::force_pos_210_0_old]     =    245.611394915;
  results_.doubles_[Key::force_pos_210_1_old]     =   -197.900427118;
  results_.doubles_[Key::force_pos_210_2_old]     =    -48.298160058;
}

void Orca_Worker_Test_Electrostatic::init_results_velocities() {
  // "current" forces (step 3) -> step not written out in test by default (STEP = 2)
  // consequence of the leap-frog algorithm
  results_.doubles_[Key::velocities_pos_0_0_current]   =    0.130304736;
  results_.doubles_[Key::velocities_pos_0_1_current]   =   -1.093819539;
  results_.doubles_[Key::velocities_pos_0_2_current]   =   -0.333630857;
  results_.doubles_[Key::velocities_pos_1_0_current]   =    0.074216900;
  results_.doubles_[Key::velocities_pos_1_1_current]   =   -0.326119514;
  results_.doubles_[Key::velocities_pos_1_2_current]   =    0.777797065;
  results_.doubles_[Key::velocities_pos_200_0_current] =    0.726016950;
  results_.doubles_[Key::velocities_pos_200_1_current] =   -0.059555086;
  results_.doubles_[Key::velocities_pos_200_2_current] =    0.148548304;
  results_.doubles_[Key::velocities_pos_210_0_current] =    0.128325067;
  results_.doubles_[Key::velocities_pos_210_1_current] =    0.399790335;
  results_.doubles_[Key::velocities_pos_210_2_current] =   -0.750546025;
  // "old" velocities (step 2)
  results_.doubles_[Key::velocities_pos_0_0_old]       =   -0.029814731;
  results_.doubles_[Key::velocities_pos_0_1_old]       =   -1.063376043;
  results_.doubles_[Key::velocities_pos_0_2_old]       =   -0.190999184;
  results_.doubles_[Key::velocities_pos_1_0_old]       =    0.111469518;
  results_.doubles_[Key::velocities_pos_1_1_old]       =   -0.299833744;
  results_.doubles_[Key::velocities_pos_1_2_old]       =    0.645675963;
  results_.doubles_[Key::velocities_pos_200_0_old]     =    0.813240182;
  results_.doubles_[Key::velocities_pos_200_1_old]     =    0.005084874;
  results_.doubles_[Key::velocities_pos_200_2_old]     =    0.063770841;
  results_.doubles_[Key::velocities_pos_210_0_old]     =    0.006491142;
  results_.doubles_[Key::velocities_pos_210_1_old]     =    0.497946269;
  results_.doubles_[Key::velocities_pos_210_2_old]     =   -0.726571782;
}

void Orca_Worker_Test_Electrostatic::init_results_positions() {
  // "current" positions (step 1)
  results_.doubles_[Key::positions_pos_0_0_current]   =    3.159254531;
  results_.doubles_[Key::positions_pos_0_1_current]   =    0.681110048;
  results_.doubles_[Key::positions_pos_0_2_current]   =    1.245389912;
  results_.doubles_[Key::positions_pos_1_0_current]   =    3.105987491;
  results_.doubles_[Key::positions_pos_1_1_current]   =    0.772524833;
  results_.doubles_[Key::positions_pos_1_2_current]   =    1.272514968;
  results_.doubles_[Key::positions_pos_200_0_current] =    2.703872896;
  results_.doubles_[Key::positions_pos_200_1_current] =    0.769687228;
  results_.doubles_[Key::positions_pos_200_2_current] =    0.986185259;
  results_.doubles_[Key::positions_pos_210_0_current] =    2.804205664;
  results_.doubles_[Key::positions_pos_210_1_current] =    3.684085782;
  results_.doubles_[Key::positions_pos_210_2_current] =    2.657992100;
  //// "old" positions (step 2)
  results_.doubles_[Key::positions_pos_0_0_old]       =    3.159240450;
  results_.doubles_[Key::positions_pos_0_1_old]       =    0.680578538;
  results_.doubles_[Key::positions_pos_0_2_old]       =    1.245294739;
  results_.doubles_[Key::positions_pos_1_0_old]       =    3.106044039;
  results_.doubles_[Key::positions_pos_1_1_old]       =    0.772375118;
  results_.doubles_[Key::positions_pos_1_2_old]       =    1.272838139;
  results_.doubles_[Key::positions_pos_200_0_old]     =    2.704280224;
  results_.doubles_[Key::positions_pos_200_1_old]     =    0.769689971;
  results_.doubles_[Key::positions_pos_200_2_old]     =    0.986217403;
  results_.doubles_[Key::positions_pos_210_0_old]     =    2.804209643;
  results_.doubles_[Key::positions_pos_210_1_old]     =    3.684335719;
  results_.doubles_[Key::positions_pos_210_2_old]     =    2.657629509;
}

void Orca_Worker_Test_Electrostatic::init_results_baths() {
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

void Orca_Worker_Test_Electrostatic::init_results_bonded_terms() {
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

void Orca_Worker_Test_Electrostatic::init_results_nonbonded_terms() {
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

void Orca_Worker_Test_Electrostatic::init_results_special_terms() {
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

void Orca_Worker_Test_Electrostatic::init_results_mass() {
  // "current" (step 1)
  results_.doubles_[Key::mass_current] = 1.991781000e+04;
  // "old" (step 2)
  results_.doubles_[Key::mass_old]     = 1.991781000e+04;
}

void Orca_Worker_Test_Electrostatic::init_results_temperature() {
  // "current" (step 1)
  results_.ints_[Key::num_temperature_coupling_baths_current]       = 2;
  results_.doubles_[Key::temperature_total_bath_0_current]          = 2.968675272e+02;
  results_.doubles_[Key::temperature_total_bath_1_current]          = 2.801836889e+02;
  results_.doubles_[Key::temperature_com_bath_0_current]            = 3.336290580e+02;
  results_.doubles_[Key::temperature_com_bath_1_current]            = 3.798019531e+02;
  results_.doubles_[Key::temperature_ir_bath_0_current]             = 2.963184326e+02;
  results_.doubles_[Key::temperature_ir_bath_1_current]             = 2.666979806e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_current] = 1.000019771e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_current] = 1.000014899e+00;
  // "old" (step 2)
  results_.ints_[Key::num_temperature_coupling_baths_old]           = 2;
  results_.doubles_[Key::temperature_total_bath_0_old]              = 3.837655402e+02;
  results_.doubles_[Key::temperature_total_bath_1_old]              = 2.791741044e+02;
  results_.doubles_[Key::temperature_com_bath_0_old]                = 3.314835071e+02;
  results_.doubles_[Key::temperature_com_bath_1_old]                = 3.801787590e+02;
  results_.doubles_[Key::temperature_ir_bath_0_old]                 = 3.841771128e+02;
  results_.doubles_[Key::temperature_ir_bath_1_old]                 = 2.655032375e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_old]     = 1.000023732e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_old]     = 1.000022414e+00;
}

void Orca_Worker_Test_Electrostatic::init_results_volume() {
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

void Orca_Worker_Test_Electrostatic::init_results_pressure() {
  // "current" (step 1)
  results_.doubles_[Key::pressure_current]                            =  1.802914765e+00;
  results_.doubles_[Key::virial_current]                              =  3.500273180e+02;
  results_.doubles_[Key::molecular_kinetic_energy_current]            =  3.959292823e+02;
  results_.doubles_[Key::pressure_tensor_0_0_current]                 = -3.093682555e+00;
  results_.doubles_[Key::pressure_tensor_0_1_current]                 =  8.268917846e+00;
  results_.doubles_[Key::pressure_tensor_0_2_current]                 = -1.103667867e+01;
  results_.doubles_[Key::pressure_tensor_1_0_current]                 =  1.309011686e+01;
  results_.doubles_[Key::pressure_tensor_1_1_current]                 = -2.251838747e+00;
  results_.doubles_[Key::pressure_tensor_1_2_current]                 = -1.458671234e+01;
  results_.doubles_[Key::pressure_tensor_2_0_current]                 = -4.996534585e+00;
  results_.doubles_[Key::pressure_tensor_2_1_current]                 = -6.379765980e+00;
  results_.doubles_[Key::pressure_tensor_2_2_current]                 =  1.075426560e+01;
  results_.doubles_[Key::virial_tensor_0_0_current]                   =  5.052180699e+02;
  results_.doubles_[Key::virial_tensor_0_1_current]                   = -1.940499198e+02;
  results_.doubles_[Key::virial_tensor_0_2_current]                   =  2.750491660e+02;
  results_.doubles_[Key::virial_tensor_1_0_current]                   = -3.167969899e+02;
  results_.doubles_[Key::virial_tensor_1_1_current]                   =  4.647049086e+02;
  results_.doubles_[Key::virial_tensor_1_2_current]                   =  3.542165144e+02;
  results_.doubles_[Key::virial_tensor_2_0_current]                   =  1.212679202e+02;
  results_.doubles_[Key::virial_tensor_2_1_current]                   =  1.452687781e+02;
  results_.doubles_[Key::virial_tensor_2_2_current]                   =  8.015897545e+01;
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
  results_.doubles_[Key::pressure_old]                            =  1.613723372e+00;
  results_.doubles_[Key::virial_old]                              =  3.552220291e+02;
  results_.doubles_[Key::molecular_kinetic_energy_old]            =  3.963072392e+02;
  results_.doubles_[Key::pressure_tensor_0_0_old]                 = -3.129626137e+00;
  results_.doubles_[Key::pressure_tensor_0_1_old]                 =  8.388690086e+00;
  results_.doubles_[Key::pressure_tensor_0_2_old]                 = -1.110619490e+01;
  results_.doubles_[Key::pressure_tensor_1_0_old]                 =  1.299149964e+01;
  results_.doubles_[Key::pressure_tensor_1_1_old]                 = -2.664260532e+00;
  results_.doubles_[Key::pressure_tensor_1_2_old]                 = -1.487290677e+01;
  results_.doubles_[Key::pressure_tensor_2_0_old]                 = -5.280628076e+00;
  results_.doubles_[Key::pressure_tensor_2_1_old]                 = -6.148213173e+00;
  results_.doubles_[Key::pressure_tensor_2_2_old]                 =  1.063505679e+01;
  results_.doubles_[Key::virial_tensor_0_0_old]                   =  5.064319489e+02;
  results_.doubles_[Key::virial_tensor_0_1_old]                   = -1.968464027e+02;
  results_.doubles_[Key::virial_tensor_0_2_old]                   =  2.765144185e+02;
  results_.doubles_[Key::virial_tensor_1_0_old]                   = -3.140334005e+02;
  results_.doubles_[Key::virial_tensor_1_1_old]                   =  4.755045847e+02;
  results_.doubles_[Key::virial_tensor_1_2_old]                   =  3.614481177e+02;
  results_.doubles_[Key::virial_tensor_2_0_old]                   =  1.281961623e+02;
  results_.doubles_[Key::virial_tensor_2_1_old]                   =  1.393184302e+02;
  results_.doubles_[Key::virial_tensor_2_2_old]                   =  8.372955379e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_old] =  4.267520302e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_old] =  1.672867459e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_old] = -6.248016444e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_old] =  1.672867459e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_old] =  4.076728204e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_old] = -1.721436458e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_old] = -6.248016444e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_old] = -1.721436458e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_old] =  3.544968668e+02;
}

void Orca_Worker_Test_Electrostatic::check_qm_atoms_init() {
  // references to shorten the code
  std::unordered_map<Key::keys, int>& ints_res = results_.ints_;
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // test the atoms of the QM zone
  for (const auto& atom : qm_zone_ptr->qm) {
    EXPECT_EQ(atom.qm_charge, doubles_res[Key::qm_atom_charge_init]);
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force(i), doubles_res[Key::qm_atom_force_init]);
    }
    switch (atom.index) { // test some selected atoms
      case 0:  EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break;
      case 48: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_15]); break;
      case 49: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_46]); break;
      case 89: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break;
    }
  }
}

void Orca_Worker_Test_Electrostatic::check_mm_atoms_init() {
  // references to shorten the code
  std::unordered_map<Key::keys, int>& ints_res = results_.ints_;
  std::unordered_map<Key::keys, double>& doubles_res = results_.doubles_;
  // test the atoms of the MM zone
  for (const auto& atom : qm_zone_ptr->mm) {
    for (unsigned int i = 0; i < 3; ++i) {
      EXPECT_EQ(atom.force[i], doubles_res[Key::mm_atom_force_init]);
      EXPECT_EQ(atom.cos_force[i], doubles_res[Key::mm_atom_cos_force_init]);
    }
    switch (atom.index) { // test some selected atoms
      case 105:  EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break; 
      case 106:  EXPECT_EQ(atom.atomic_number, ints_res[Key::element_6]); break; 
      case 3232: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_8]); break;
      case 3233: EXPECT_EQ(atom.atomic_number, ints_res[Key::element_1]); break;
    }
  }
}

/**
 * Test if initialization is successful
 * 
 */
TEST_F(Orca_Worker_Test_Electrostatic, check_init) {
  // check if input vile has been loaded correctly
  check_parameter_init(); 
}

/**
 * Test if simulation runs successfully and matches expected energy terms
 * 
 */
TEST_F(Orca_Worker_Test_Electrostatic, check_simulation) {
  // run the simulation specified by the parameters object
  test_sim_.run_simulation();
  check_simulation_results();
}

} // namespace testing
