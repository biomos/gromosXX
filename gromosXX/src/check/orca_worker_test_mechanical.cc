/**
 * @file orca_worker_test_mechanical.cc
 */

#include <gtest/gtest.h>

#include "orca_worker_test_mechanical.h"
#include "../interaction/qmmm/orca_worker.h"
#include "../interaction/qmmm/qm_atom.h"
#include "../interaction/qmmm/mm_atom.h"

namespace testing {

Orca_Worker_Test_Mechanical::Orca_Worker_Test_Mechanical() : Orca_Worker_Test("Orca QM/MM Tests Mechanical Embedding") {
  // initialize test files and expected results
  init_parameters();
  init_results();
}

void Orca_Worker_Test_Mechanical::init_parameters() {
  test_sim_.parameter().add_input("topo", "md.top");
  test_sim_.parameter().add_input("input", "md_mc.imd");
  test_sim_.parameter().add_input("conf", "md.cnf");
  test_sim_.parameter().add_input("qmmm", "md_mc.qmmm");
  test_sim_.parameter().add_input("trc", "md_mc.trc");
  test_sim_.parameter().add_input("tre", "md_mc.tre");
  test_sim_.parameter().add_input("trv", "md_mc.trv");
  test_sim_.parameter().add_input("trf", "md_mc.trf");
  test_sim_.parameter().add_input("fin", "md_mc_final.cnf");
}

void Orca_Worker_Test_Mechanical::init_results() {
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

void Orca_Worker_Test_Mechanical::init_results_parameters() {
  // parameters from the input file
  results_.ints_[Key::qmmm] = simulation::qmmm_mechanical;
  results_.ints_[Key::qm_software] = simulation::qm_orca;
  results_.ints_[Key::qm_ch] = simulation::qm_ch_dynamic;
  results_.doubles_[Key::qm_cutoff] = 1.2;
  results_.ints_[Key::qm_lj] = simulation::qm_lj_off;
  results_.ints_[Key::qm_constraint] = simulation::qm_constr_off;
  results_.doubles_[Key::qm_mm_scale] = -1.0;
}

void Orca_Worker_Test_Mechanical::init_results_files() {
  // input file names
  results_.strs_[Key::input_file] = "xphos_mc.inp";
  results_.strs_[Key::input_coordinate_file] = "xphos_mc.xyz";
  results_.strs_[Key::input_pointcharges_file] = "xphos_mc.pc";
  results_.strs_[Key::output_file] = "xphos_mc.out";
  results_.strs_[Key::output_gradient_file] = "xphos_mc.engrad";
  results_.strs_[Key::output_mm_gradient_file] = "xphos_mc.pcgrad";
}

void Orca_Worker_Test_Mechanical::init_results_qm_zone_init() {
  // QM zone
  results_.ints_[Key::qm_zone_size_init] = 90;
  results_.ints_[Key::mm_zone_size_init] = 0;
  results_.ints_[Key::qm_zone_charge] = 0;
  results_.ints_[Key::qm_zone_spin_mult] = 1;
  results_.doubles_[Key::qm_energy_init] = 0.0;
  results_.doubles_[Key::qm_atom_charge_init] = 0.0;
  results_.doubles_[Key::qm_atom_force_init] = 0.0;
  results_.doubles_[Key::mm_atom_force_init] = 0.0;
  results_.doubles_[Key::mm_atom_cos_force_init] = 0.0;
}

void Orca_Worker_Test_Mechanical::init_results_energies() {
  // "current" energies (step 1)
  results_.doubles_[Key::energy_total_current]                           =  -9.427274939e+05;
  results_.doubles_[Key::energy_kinetic_total_current]                   =   7.725377958e+03;
  results_.doubles_[Key::energy_potential_total_current]                 =  -9.504528718e+05;
  results_.doubles_[Key::energy_covalent_total_current]                  =   4.714885375e+03;
  results_.doubles_[Key::energy_bonds_total_current]                     =   0.000000000e+00;
  results_.doubles_[Key::energy_angles_total_current]                    =   2.781463488e+03;
  results_.doubles_[Key::energy_impropers_total_current]                 =   1.247690511e+03;
  results_.doubles_[Key::energy_dihedrals_total_current]                 =   6.857313771e+02;
  results_.doubles_[Key::energy_crossdihedrals_total_current]            =   0.000000000e+00;
  results_.doubles_[Key::energy_nonbonded_total_current]                 =  -3.068805987e+04;
  results_.doubles_[Key::energy_lennard_jones_total_current]             =  -4.996215450e+03;
  results_.doubles_[Key::energy_coulomb_reaction_field_total_current]    =  -2.569184442e+04;
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
  results_.doubles_[Key::qmmm_total_current]                             =  -9.244796973e+05;
  results_.doubles_[Key::bs_leus_energy_current]                         =   0.000000000e+00;
  results_.doubles_[Key::rdc_value_total_current]                        =   0.000000000e+00;
  results_.doubles_[Key::angle_restraints_total_current]                 =   0.000000000e+00;
  // "old" energies (step 2)
  results_.doubles_[Key::energy_total_old]                               =  -9.427313974e+05;
  results_.doubles_[Key::energy_kinetic_total_old]                       =   7.795226609e+03;
  results_.doubles_[Key::energy_potential_total_old]                     =  -9.505266241e+05;
  results_.doubles_[Key::energy_covalent_total_old]                      =   4.743388758e+03;
  results_.doubles_[Key::energy_bonds_total_old]                         =   0.000000000e+00;
  results_.doubles_[Key::energy_angles_total_old]                        =   2.812634368e+03;
  results_.doubles_[Key::energy_impropers_total_old]                     =   1.242790246e+03;
  results_.doubles_[Key::energy_dihedrals_total_old]                     =   6.879641434e+02;
  results_.doubles_[Key::energy_crossdihedrals_total_old]                =   0.000000000e+00;
  results_.doubles_[Key::energy_nonbonded_total_old]                     =  -3.068880330e+04;
  results_.doubles_[Key::energy_lennard_jones_total_old]                 =  -4.996873197e+03;
  results_.doubles_[Key::energy_coulomb_reaction_field_total_old]        =  -2.569193010e+04;
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
  results_.doubles_[Key::qmmm_total_old]                                 =  -9.245812095e+05;
  results_.doubles_[Key::bs_leus_energy_old]                             =   0.000000000e+00;
  results_.doubles_[Key::rdc_value_total_old]                            =   0.000000000e+00;
  results_.doubles_[Key::angle_restraints_total_old]                     =   0.000000000e+00;
}

void Orca_Worker_Test_Mechanical::init_results_forces() {
  // "current" forces (step 1)
  results_.doubles_[Key::force_pos_0_0_current]   =    439.218233773;
  results_.doubles_[Key::force_pos_0_1_current]   =   -383.480941817;
  results_.doubles_[Key::force_pos_0_2_current]   =   -357.032960733;
  results_.doubles_[Key::force_pos_1_0_current]   =   -849.777182866;
  results_.doubles_[Key::force_pos_1_1_current]   =   -608.168893961;
  results_.doubles_[Key::force_pos_1_2_current]   =   3649.045083024;
  results_.doubles_[Key::force_pos_200_0_current] =   -170.894151393;
  results_.doubles_[Key::force_pos_200_1_current] =   -120.728440341;
  results_.doubles_[Key::force_pos_200_2_current] =    160.678090070;
  results_.doubles_[Key::force_pos_210_0_current] =    239.442654526;
  results_.doubles_[Key::force_pos_210_1_current] =   -185.211523222;
  results_.doubles_[Key::force_pos_210_2_current] =    -50.938980060;
  // "old" forces (step 2)
  results_.doubles_[Key::force_pos_0_0_old]       =    323.162208784;
  results_.doubles_[Key::force_pos_0_1_old]       =    -67.550663607;
  results_.doubles_[Key::force_pos_0_2_old]       =   -285.588989842;
  results_.doubles_[Key::force_pos_1_0_old]       =   -895.359865079;
  results_.doubles_[Key::force_pos_1_1_old]       =   -609.145393347;
  results_.doubles_[Key::force_pos_1_2_old]       =   3171.334370292;
  results_.doubles_[Key::force_pos_200_0_old]     =   -176.758981804;
  results_.doubles_[Key::force_pos_200_1_old]     =   -130.943962314;
  results_.doubles_[Key::force_pos_200_2_old]     =    170.042681989;
  results_.doubles_[Key::force_pos_210_0_old]     =    245.856723587;
  results_.doubles_[Key::force_pos_210_1_old]     =   -198.100479582;
  results_.doubles_[Key::force_pos_210_2_old]     =    -48.097047479;
}

void Orca_Worker_Test_Mechanical::init_results_velocities() {
  // "current" forces (step 3) -> step not written out in test by default (STEP = 2)
  // consequence of the leap-frog algorithm
  results_.doubles_[Key::velocities_pos_0_0_current]   =    0.134308379;
  results_.doubles_[Key::velocities_pos_0_1_current]   =   -1.109033388;
  results_.doubles_[Key::velocities_pos_0_2_current]   =   -0.332272542;
  results_.doubles_[Key::velocities_pos_1_0_current]   =    0.074573336;
  results_.doubles_[Key::velocities_pos_1_1_current]   =   -0.322044862;
  results_.doubles_[Key::velocities_pos_1_2_current]   =    0.777399550;
  results_.doubles_[Key::velocities_pos_200_0_current] =    0.724539015;
  results_.doubles_[Key::velocities_pos_200_1_current] =   -0.060163276;
  results_.doubles_[Key::velocities_pos_200_2_current] =    0.147616775;
  results_.doubles_[Key::velocities_pos_210_0_current] =    0.128485876;
  results_.doubles_[Key::velocities_pos_210_1_current] =    0.399433922;
  results_.doubles_[Key::velocities_pos_210_2_current] =   -0.750472237;
  // "old" velocities (step 2)
  results_.doubles_[Key::velocities_pos_0_0_old]       =   -0.025993535;
  results_.doubles_[Key::velocities_pos_0_1_old]       =   -1.075499701;
  results_.doubles_[Key::velocities_pos_0_2_old]       =   -0.190603424;
  results_.doubles_[Key::velocities_pos_1_0_old]       =    0.111844055;
  results_.doubles_[Key::velocities_pos_1_1_old]       =   -0.296679378;
  results_.doubles_[Key::velocities_pos_1_2_old]       =    0.645363120;
  results_.doubles_[Key::velocities_pos_200_0_old]     =    0.812200842;
  results_.doubles_[Key::velocities_pos_200_1_old]     =    0.004790434;
  results_.doubles_[Key::velocities_pos_200_2_old]     =    0.063266898;
  results_.doubles_[Key::velocities_pos_210_0_old]     =    0.006530257;
  results_.doubles_[Key::velocities_pos_210_1_old]     =    0.497689096;
  results_.doubles_[Key::velocities_pos_210_2_old]     =   -0.726597754;
}

void Orca_Worker_Test_Mechanical::init_results_positions() {
  // "current" positions (step 1)
  results_.doubles_[Key::positions_pos_0_0_current]   =    3.159255642;
  results_.doubles_[Key::positions_pos_0_1_current]   =    0.681106676;
  results_.doubles_[Key::positions_pos_0_2_current]   =    1.245389911;
  results_.doubles_[Key::positions_pos_1_0_current]   =    3.105987567;
  results_.doubles_[Key::positions_pos_1_1_current]   =    0.772525687;
  results_.doubles_[Key::positions_pos_1_2_current]   =    1.272514864;
  results_.doubles_[Key::positions_pos_200_0_current] =    2.703872582;
  results_.doubles_[Key::positions_pos_200_1_current] =    0.769687228;
  results_.doubles_[Key::positions_pos_200_2_current] =    0.986185198;
  results_.doubles_[Key::positions_pos_210_0_current] =    2.804205682;
  results_.doubles_[Key::positions_pos_210_1_current] =    3.684085704;
  results_.doubles_[Key::positions_pos_210_2_current] =    2.657992106;
    //// "old" positions (step 2)
  results_.doubles_[Key::positions_pos_0_0_old]       =    3.159243410;
  results_.doubles_[Key::positions_pos_0_1_old]       =    0.680569091;
  results_.doubles_[Key::positions_pos_0_2_old]       =    1.245294911;
  results_.doubles_[Key::positions_pos_1_0_old]       =    3.106044242;
  results_.doubles_[Key::positions_pos_1_1_old]       =    0.772377534;
  results_.doubles_[Key::positions_pos_1_2_old]       =    1.272837853;
  results_.doubles_[Key::positions_pos_200_0_old]     =    2.704279338;
  results_.doubles_[Key::positions_pos_200_1_old]     =    0.769689810;
  results_.doubles_[Key::positions_pos_200_2_old]     =    0.986217071;
  results_.doubles_[Key::positions_pos_210_0_old]     =    2.804209627;
  results_.doubles_[Key::positions_pos_210_1_old]     =    3.684335441;
  results_.doubles_[Key::positions_pos_210_2_old]     =    2.657629451;
}

void Orca_Worker_Test_Mechanical::init_results_baths() {
  // "current" (step 1)
  results_.doubles_[Key::kinetic_total_bath_0_current]       = 3.323082462e+02;
  results_.doubles_[Key::kinetic_total_bath_1_current]       = 7.393069712e+03;
  results_.doubles_[Key::centre_of_mass_bath_0_current]      = 4.168582379e+00;
  results_.doubles_[Key::centre_of_mass_bath_1_current]      = 1.184143946e+03;
  results_.doubles_[Key::internal_rotational_bath_0_current] = 3.281396639e+02;
  results_.doubles_[Key::internal_rotational_bath_1_current] = 6.208925767e+03;
  // "old" (step 2)
  results_.doubles_[Key::kinetic_total_bath_0_old]           = 4.287424173e+02;
  results_.doubles_[Key::kinetic_total_bath_1_old]           = 7.366484192e+03;
  results_.doubles_[Key::centre_of_mass_bath_0_old]          = 4.146932216e+00;
  results_.doubles_[Key::centre_of_mass_bath_1_old]          = 1.185289872e+03;
  results_.doubles_[Key::internal_rotational_bath_0_old]     = 4.245954851e+02;
  results_.doubles_[Key::internal_rotational_bath_1_old]     = 6.181194320e+03;
}

void Orca_Worker_Test_Mechanical::init_results_bonded_terms() {
  // "current" (step 1)
  results_.doubles_[Key::bond_energy_group_0_current]               = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_1_current]               = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_0_current]              = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_1_current]              = 2.781463488e+03;
  results_.doubles_[Key::improper_energy_group_0_current]           = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_1_current]           = 1.247690511e+03;
  results_.doubles_[Key::dihedral_energy_group_0_current]           = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_1_current]           = 6.857313771e+02;
  results_.doubles_[Key::crossdihedral_energy_group_0_current]      = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_1_current]      = 0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::bond_energy_group_0_old]                   = 0.000000000e+00;
  results_.doubles_[Key::bond_energy_group_1_old]                   = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_0_old]                  = 0.000000000e+00;
  results_.doubles_[Key::angle_energy_group_1_old]                  = 2.812634368e+03;
  results_.doubles_[Key::improper_energy_group_0_old]               = 0.000000000e+00;
  results_.doubles_[Key::improper_energy_group_1_old]               = 1.242790246e+03;
  results_.doubles_[Key::dihedral_energy_group_0_old]               = 0.000000000e+00;
  results_.doubles_[Key::dihedral_energy_group_1_old]               = 6.879641434e+02;
  results_.doubles_[Key::crossdihedral_energy_group_0_old]          = 0.000000000e+00;
  results_.doubles_[Key::crossdihedral_energy_group_1_old]          = 0.000000000e+00;
}

void Orca_Worker_Test_Mechanical::init_results_nonbonded_terms() {
  // "current" (step 1)
  results_.doubles_[Key::lennard_jones_group_0_0_current]                =  0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_1_0_current]                = -1.275111904e+02;
  results_.doubles_[Key::lennard_jones_group_1_1_current]                = -4.868704259e+03;
  results_.doubles_[Key::coulomb_reaction_field_group_0_0_current]       = -9.161646950e+01;
  results_.doubles_[Key::coulomb_reaction_field_group_1_0_current]       = -5.259865857e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_1_1_current]       = -2.559496809e+04;
  results_.doubles_[Key::lattice_sum_real_group_0_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_0_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_1_current]             =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_0_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_0_current]       =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_1_current]       =  0.000000000e+00;
  // "old" (step 2)
  results_.doubles_[Key::lennard_jones_group_0_0_old]                    =  0.000000000e+00;
  results_.doubles_[Key::lennard_jones_group_1_0_old]                    = -1.276017905e+02;
  results_.doubles_[Key::lennard_jones_group_1_1_old]                    = -4.869271407e+03;
  results_.doubles_[Key::coulomb_reaction_field_group_0_0_old]           = -9.072299234e+01;
  results_.doubles_[Key::coulomb_reaction_field_group_1_0_old]           = -5.750324629e+00;
  results_.doubles_[Key::coulomb_reaction_field_group_1_1_old]           = -2.559545679e+04;
  results_.doubles_[Key::lattice_sum_real_group_0_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_0_old]                 =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_real_group_1_1_old]                 =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_0_0_old]           =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_0_old]           =  0.000000000e+00;
  results_.doubles_[Key::lattice_sum_reciprocal_group_1_1_old]           =  0.000000000e+00;
}

void Orca_Worker_Test_Mechanical::init_results_special_terms() {
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

void Orca_Worker_Test_Mechanical::init_results_mass() {
  // "current" (step 1)
  results_.doubles_[Key::mass_current] = 1.991781000e+04;
  // "old" (step 2)
  results_.doubles_[Key::mass_old]     = 1.991781000e+04;
}

void Orca_Worker_Test_Mechanical::init_results_temperature() {
  // "current" (step 1)
  results_.ints_[Key::num_temperature_coupling_baths_current]       = 2;
  results_.doubles_[Key::temperature_total_bath_0_current]          = 2.961916207e+02;
  results_.doubles_[Key::temperature_total_bath_1_current]          = 2.801860379e+02;
  results_.doubles_[Key::temperature_com_bath_0_current]            = 3.342455953e+02;
  results_.doubles_[Key::temperature_com_bath_1_current]            = 3.797884863e+02;
  results_.doubles_[Key::temperature_ir_bath_0_current]             = 2.956283140e+02;
  results_.doubles_[Key::temperature_ir_bath_1_current]             = 2.667024465e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_current] = 1.000019788e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_current] = 1.000014899e+00;
  // "old" (step 2)
  results_.ints_[Key::num_temperature_coupling_baths_old]           = 2;
  results_.doubles_[Key::temperature_total_bath_0_old]              = 3.821449298e+02;
  results_.doubles_[Key::temperature_total_bath_1_old]              = 2.791784873e+02;
  results_.doubles_[Key::temperature_com_bath_0_old]                = 3.325096402e+02;
  results_.doubles_[Key::temperature_com_bath_1_old]                = 3.801560172e+02;
  results_.doubles_[Key::temperature_ir_bath_0_old]                 = 3.825275064e+02;
  results_.doubles_[Key::temperature_ir_bath_1_old]                 = 2.655112510e+02;
  results_.doubles_[Key::temperature_scaling_factor_bath_0_old]     = 1.000023817e+00;
  results_.doubles_[Key::temperature_scaling_factor_bath_1_old]     = 1.000022414e+00;
}

void Orca_Worker_Test_Mechanical::init_results_volume() {
// "current" (step 1)
  results_.doubles_[Key::volume_current]  = 5.091972555e+01;
  results_.doubles_[Key::box_k_0_current] = 3.706483047e+00;
  results_.doubles_[Key::box_k_1_current] = 0.000000000e+00;
  results_.doubles_[Key::box_k_2_current] = 0.000000000e+00;
  results_.doubles_[Key::box_l_0_current] = 0.000000000e+00;
  results_.doubles_[Key::box_l_1_current] = 3.706483047e+00;
  results_.doubles_[Key::box_l_2_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_0_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_1_current] = 0.000000000e+00;
  results_.doubles_[Key::box_m_2_current] = 3.706483047e+00;
  // "old" (step 2)
  results_.doubles_[Key::volume_old]      = 5.091976255e+01;
  results_.doubles_[Key::box_k_0_old]     = 3.706483945e+00;
  results_.doubles_[Key::box_k_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_k_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_l_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_l_1_old]     = 3.706483945e+00;
  results_.doubles_[Key::box_l_2_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_0_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_1_old]     = 0.000000000e+00;
  results_.doubles_[Key::box_m_2_old]     = 3.706483945e+00;
}

void Orca_Worker_Test_Mechanical::init_results_pressure() {
// "current" (step 1)
  results_.doubles_[Key::pressure_current]                            =  1.647660615e+00;
  results_.doubles_[Key::virial_current]                              =  3.539725488e+02;
  results_.doubles_[Key::molecular_kinetic_energy_current]            =  3.959217620e+02;
  results_.doubles_[Key::pressure_tensor_0_0_current]                 = -3.766366594e+00;
  results_.doubles_[Key::pressure_tensor_0_1_current]                 =  8.580545820e+00;
  results_.doubles_[Key::pressure_tensor_0_2_current]                 = -1.088919914e+01;
  results_.doubles_[Key::pressure_tensor_1_0_current]                 =  1.365516138e+01;
  results_.doubles_[Key::pressure_tensor_1_1_current]                 = -2.602074418e+00;
  results_.doubles_[Key::pressure_tensor_1_2_current]                 = -1.486827679e+01;
  results_.doubles_[Key::pressure_tensor_2_0_current]                 = -8.642592868e+00;
  results_.doubles_[Key::pressure_tensor_2_1_current]                 = -2.951667672e+00;
  results_.doubles_[Key::pressure_tensor_2_2_current]                 =  1.131142286e+01;
  results_.doubles_[Key::virial_tensor_0_0_current]                   =  5.223153755e+02;
  results_.doubles_[Key::virial_tensor_0_1_current]                   = -2.019800080e+02;
  results_.doubles_[Key::virial_tensor_0_2_current]                   =  2.712750175e+02;
  results_.doubles_[Key::virial_tensor_1_0_current]                   = -3.311790237e+02;
  results_.doubles_[Key::virial_tensor_1_1_current]                   =  4.736239941e+02;
  results_.doubles_[Key::virial_tensor_1_2_current]                   =  3.613946584e+02;
  results_.doubles_[Key::virial_tensor_2_0_current]                   =  2.140767301e+02;
  results_.doubles_[Key::virial_tensor_2_1_current]                   =  5.799942561e+01;
  results_.doubles_[Key::virial_tensor_2_2_current]                   =  6.597827676e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_current] =  4.264241988e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_current] =  1.647951115e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_current] = -5.962498400e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_current] =  1.647951115e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_current] =  4.073755365e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_current] = -1.714962828e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_current] = -5.962498400e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_current] = -1.714962828e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_current] =  3.539655505e+02;
  // "old" (step 2)
  results_.doubles_[Key::pressure_old]                                =  1.406391830e+00;
  results_.doubles_[Key::virial_old]                                  =  3.604852268e+02;
  results_.doubles_[Key::molecular_kinetic_energy_old]                =  3.962917958e+02;
  results_.doubles_[Key::pressure_tensor_0_0_old]                     = -3.797616893e+00;
  results_.doubles_[Key::pressure_tensor_0_1_old]                     =  8.707589398e+00;
  results_.doubles_[Key::pressure_tensor_0_2_old]                     = -1.100016851e+01;
  results_.doubles_[Key::pressure_tensor_1_0_old]                     =  1.352268475e+01;
  results_.doubles_[Key::pressure_tensor_1_1_old]                     = -3.175406759e+00;
  results_.doubles_[Key::pressure_tensor_1_2_old]                     = -1.505435219e+01;
  results_.doubles_[Key::pressure_tensor_2_0_old]                     = -8.842460616e+00;
  results_.doubles_[Key::pressure_tensor_2_1_old]                     = -3.209503002e+00;
  results_.doubles_[Key::pressure_tensor_2_2_old]                     =  1.119219914e+01;
  results_.doubles_[Key::virial_tensor_0_0_old]                       =  5.233793923e+02;
  results_.doubles_[Key::virial_tensor_0_1_old]                       = -2.049575108e+02;
  results_.doubles_[Key::virial_tensor_0_2_old]                       =  2.737765356e+02;
  results_.doubles_[Key::virial_tensor_1_0_old]                       = -3.275492667e+02;
  results_.doubles_[Key::virial_tensor_1_1_old]                       =  4.885229471e+02;
  results_.doubles_[Key::virial_tensor_1_2_old]                       =  3.660865035e+02;
  results_.doubles_[Key::virial_tensor_2_0_old]                       =  2.188415487e+02;
  results_.doubles_[Key::virial_tensor_2_1_old]                       =  6.451804943e+01;
  results_.doubles_[Key::virial_tensor_2_2_old]                       =  6.955334083e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_0_old]     =  4.266925171e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_1_old]     =  1.673668148e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_0_2_old]     = -6.286448795e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_0_old]     =  1.673668148e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_1_old]     =  4.076774680e+02;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_1_2_old]     = -1.719551595e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_0_old]     = -6.286448795e+00;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_1_old]     = -1.719551595e+01;
  results_.doubles_[Key::molecular_kinetic_energy_tensor_2_2_old]     =  3.545054022e+02;
}

void Orca_Worker_Test_Mechanical::check_qm_atoms_init() {
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

void Orca_Worker_Test_Mechanical::check_mm_atoms_init() {
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
TEST_F(Orca_Worker_Test_Mechanical, check_init) {
  // check if input vile has been loaded correctly
  check_parameter_init(); 
}

/**
 * Test if simulation runs successfully and matches expected energy terms
 * 
 */
TEST_F(Orca_Worker_Test_Mechanical, check_simulation) {
  // run the simulation specified by the parameters object
  test_sim_.run_simulation();
  check_simulation_results();
}

} // namespace testing