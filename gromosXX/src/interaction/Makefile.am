AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/interaction -I$(top_srcdir)/src/interaction/nonbonded -I$(top_srcdir)/src/interaction/nonbonded/pairlist -I$(top_srcdir)/src/interaction/nonbonded/interaction -I$(top_srcdir)/src/interaction/bonded -I$(top_srcdir)/src/interaction/forcefield -I$(top_srcdir)/src/interaction/special -I$(top_srcdir)/src/interaction/special/xray -I$(top_srcdir)/src/interaction/special/qmmm -I$(top_srcdir)/src/algorithm 


xxincludedir = $(includedir)/md++/interaction
xxinclude_HEADERS =	interaction.h\
	interaction_types.h\
	molecular_virial_interaction.h

xxintffdir = $(includedir)/md++/interaction/forcefield
xxintff_HEADERS = forcefield/forcefield.h\
	forcefield/create_forcefield.h

xxintbondeddir = $(includedir)/md++/interaction/bonded
xxintbonded_HEADERS = bonded/create_bonded.h\
	bonded/quartic_bond_interaction.h\
	bonded/harmonic_bond_interaction.h\
	bonded/cg_bond_interaction.h\
	bonded/angle_interaction.h\
	bonded/harm_angle_interaction.h\
	bonded/improper_dihedral_interaction.h\
	bonded/dihedral_interaction.h\
	bonded/dihedral_new_interaction.h\
	bonded/crossdihedral_interaction.h\
	bonded/perturbed_quartic_bond_interaction.h\
	bonded/perturbed_harmonic_bond_interaction.h\
	bonded/perturbed_cg_bond_interaction.h\
	bonded/perturbed_angle_interaction.h\
	bonded/perturbed_improper_dihedral_interaction.h\
	bonded/perturbed_dihedral_new_interaction.h\
	bonded/perturbed_dihedral_interaction.h

xxintnbdir = $(includedir)/md++/interaction/nonbonded
xxintnb_HEADERS = nonbonded/create_nonbonded.h\
	nonbonded/interaction_spec.h\
	nonbonded/innerloop_template.h

xxintnbintdir = $(includedir)/md++/interaction/nonbonded/interaction
xxintnbint_HEADERS = nonbonded/interaction/nonbonded_parameter.h\
	nonbonded/interaction/nonbonded_term.h\
	nonbonded/interaction/nonbonded_term.cc\
	nonbonded/interaction/nonbonded_innerloop.h\
	nonbonded/interaction/nonbonded_innerloop.cc\
	nonbonded/interaction/solvent_innerloop.cc\
	nonbonded/interaction/nonbonded_outerloop.h\
	nonbonded/interaction/perturbed_nonbonded_term.h\
	nonbonded/interaction/perturbed_nonbonded_term.cc\
	nonbonded/interaction/perturbed_nonbonded_innerloop.h\
	nonbonded/interaction/perturbed_nonbonded_innerloop.cc\
	nonbonded/interaction/perturbed_nonbonded_outerloop.h\
	nonbonded/interaction/perturbed_nonbonded_pair.h\
	nonbonded/interaction/eds_nonbonded_term.h\
	nonbonded/interaction/eds_nonbonded_term.cc\
	nonbonded/interaction/eds_nonbonded_innerloop.h\
	nonbonded/interaction/eds_nonbonded_innerloop.cc\
	nonbonded/interaction/eds_nonbonded_outerloop.h\
	nonbonded/interaction/nonbonded_set_interface.h\
	nonbonded/interaction/nonbonded_set.h\
	nonbonded/interaction/perturbed_nonbonded_set.h\
	nonbonded/interaction/eds_nonbonded_set.h\
	nonbonded/interaction/nonbonded_interaction.h\
	nonbonded/interaction/omp_nonbonded_interaction.h\
	nonbonded/interaction/mpi_nonbonded_master.h\
	nonbonded/interaction/mpi_nonbonded_slave.h\
	nonbonded/interaction/storage.h\
	nonbonded/interaction/latticesum.h\
	nonbonded/interaction/cuda_nonbonded_set.h\
	nonbonded/interaction/spc_table.h

xxintnbpairlistdir = $(includedir)/md++/interaction/nonbonded/pairlist
xxintnbpairlist_HEADERS = nonbonded/pairlist/pairlist.h\
	nonbonded/pairlist/pairlist_algorithm.h\
	nonbonded/pairlist/standard_pairlist_algorithm.h\
	nonbonded/pairlist/extended_grid_pairlist_algorithm.h\
	nonbonded/pairlist/grid_cell_pairlist.h

xxintnbspecialdir = $(includedir)/md++/interaction/special
xxintnbspecial_HEADERS = special/create_special.h\
	special/position_restraint_interaction.h\
	special/distance_restraint_interaction.h\
	special/distance_field_interaction.h\
	special/dihedral_restraint_interaction.h\
	special/perturbed_distance_restraint_interaction.h\
	special/perturbed_distance_field_interaction.h\
	special/eds_distance_restraint_interaction.h\
	special/bs_interaction.h\
	special/perturbed_dihedral_restraint_interaction.h\
	special/jvalue_restraint_interaction.h\
	special/pscale.h\
	special/xray_restraint_interaction.h\
	special/local_elevation_interaction.h\
	special/electric_field_interaction.h\
	special/adde_reweighting.h\
	special/nemd.h\
	special/order_parameter_restraint_interaction.h\
	special/rdc_restraint_interaction.h\
	special/qmmm_interaction.h\
	special/symmetry_restraint_interaction.h

xxintspecialxraydir = $(includedir)/md++/interaction/special/xray
xxintspecialxray_HEADERS = special/xray/dens.h\
	special/xray/sf.h\
	special/xray/bfactor.h

xxintspecialqmmmdir = $(includedir)/md++/interaction/special/qmmm
xxintspecialqmmm_HEADERS = special/qmmm/mm_atom.h\
	special/qmmm/qm_storage.h\
	special/qmmm/qm_worker.h\
	special/qmmm/mndo_worker.h\
	special/qmmm/turbomole_worker.h\
	special/qmmm/gathering.h

noinst_LTLIBRARIES = libgrointeraction.la

libgrointeraction_la_SOURCES =interaction.cc\
	molecular_virial_interaction.cc\
	forcefield/forcefield.cc\
	forcefield/create_forcefield.cc\
	nonbonded/create_nonbonded.cc\
	nonbonded/pairlist/pairlist.cc\
	nonbonded/pairlist/standard_pairlist_algorithm.cc\
	nonbonded/pairlist/standard_pairlist_algorithm_atomic.cc\
	nonbonded/pairlist/extended_grid.cc\
	nonbonded/pairlist/extended_grid_pairlist_algorithm.cc\
	nonbonded/pairlist/grid_cell_pairlist.cc\
	nonbonded/interaction/nonbonded_interaction.cc\
	nonbonded/interaction/omp_nonbonded_interaction.cc\
	nonbonded/interaction/mpi_nonbonded_master.cc\
	nonbonded/interaction/mpi_nonbonded_slave.cc\
	nonbonded/interaction/nonbonded_set.cc\
	nonbonded/interaction/nonbonded_outerloop.cc\
	nonbonded/interaction/perturbed_nonbonded_set.cc\
	nonbonded/interaction/perturbed_nonbonded_outerloop.cc\
	nonbonded/interaction/perturbed_nonbonded_pair.cc\
	nonbonded/interaction/latticesum.cc\
	nonbonded/interaction/eds_nonbonded_set.cc\
	nonbonded/interaction/eds_nonbonded_outerloop.cc\
	nonbonded/interaction/cuda_nonbonded_set.cc\
	bonded/create_bonded.cc\
	bonded/quartic_bond_interaction.cc\
	bonded/perturbed_quartic_bond_interaction.cc\
	bonded/harmonic_bond_interaction.cc\
	bonded/perturbed_harmonic_bond_interaction.cc\
	bonded/cg_bond_interaction.cc\
	bonded/perturbed_cg_bond_interaction.cc\
	bonded/angle_interaction.cc\
	bonded/harm_angle_interaction.cc\
	bonded/perturbed_angle_interaction.cc\
	bonded/improper_dihedral_interaction.cc\
	bonded/perturbed_improper_dihedral_interaction.cc\
	bonded/dihedral_interaction.cc\
	bonded/dihedral_new_interaction.cc\
	bonded/crossdihedral_interaction.cc\
	bonded/perturbed_dihedral_interaction.cc\
	bonded/perturbed_dihedral_new_interaction.cc\
	special/create_special.cc\
	special/position_restraint_interaction.cc\
	special/distance_restraint_interaction.cc\
	special/distance_field_interaction.cc\
	special/dihedral_restraint_interaction.cc\
	special/perturbed_distance_restraint_interaction.cc\
	special/perturbed_distance_field_interaction.cc\
	special/eds_distance_restraint_interaction.cc\
	special/bs_interaction.cc\
	special/perturbed_dihedral_restraint_interaction.cc\
	special/jvalue_restraint_interaction.cc\
	special/pscale.cc\
	special/xray_restraint_interaction.cc\
	special/local_elevation_interaction.cc\
	special/electric_field_interaction.cc\
	special/adde_reweighting.cc\
	special/nemd.cc\
	special/xray/dens.cc\
	special/xray/sf.cc\
	special/xray/bfactor.cc\
	special/order_parameter_restraint_interaction.cc\
	special/rdc_restraint_interaction.cc\
	special/qmmm/mm_atom.cc\
	special/qmmm/qm_worker.cc\
	special/qmmm/mndo_worker.cc\
	special/qmmm/turbomole_worker.cc\
	special/qmmm/gathering.cc\
	special/qmmm_interaction.cc\
	special/symmetry_restraint_interaction.cc

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
