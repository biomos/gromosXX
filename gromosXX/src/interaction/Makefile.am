INCLUDES = -I$(top_srcdir)/src

include $(top_srcdir)/src/interaction/nonbonded/split/Makefile.inc

xxincludedir = $(includedir)/gromosXX

xxinclude_HEADERS =		interaction.h \
				interaction_types.h \
				forcefield/forcefield.h \
				forcefield/create_forcefield.h \
				bonded/create_bonded.h \
				bonded/quartic_bond_interaction.h \
				bonded/quartic_bond_interaction.tcc \
				bonded/harmonic_bond_interaction.h \
				bonded/harmonic_bond_interaction.tcc \
				bonded/angle_interaction.h \
				bonded/angle_interaction.tcc \
				bonded/improper_dihedral_interaction.h \
				bonded/improper_dihedral_interaction.tcc \
				bonded/dihedral_interaction.h \
				bonded/dihedral_interaction.tcc \
				bonded/perturbed_quartic_bond_interaction.h \
				bonded/perturbed_quartic_bond_interaction.tcc \
				bonded/perturbed_harmonic_bond_interaction.h \
				bonded/perturbed_harmonic_bond_interaction.tcc \
				bonded/perturbed_angle_interaction.h \
				bonded/perturbed_angle_interaction.tcc \
				bonded/perturbed_improper_dihedral_interaction.h \
				bonded/perturbed_improper_dihedral_interaction.tcc \
				bonded/perturbed_dihedral_interaction.h \
				bonded/perturbed_dihedral_interaction.tcc \
				nonbonded/create_nonbonded.h \
				nonbonded/interaction_spec.h \
				nonbonded/interaction/nonbonded_parameter.h \
				nonbonded/interaction/nonbonded_term.h \
				nonbonded/interaction/nonbonded_term.tcc \
				nonbonded/interaction/nonbonded_innerloop.h \
				nonbonded/interaction/nonbonded_innerloop.tcc \
				nonbonded/interaction/nonbonded_outerloop.h \
				nonbonded/interaction/nonbonded_outerloop.tcc \
				nonbonded/interaction/perturbed_nonbonded_term.h \
				nonbonded/interaction/perturbed_nonbonded_term.tcc \
				nonbonded/interaction/perturbed_nonbonded_innerloop.h \
				nonbonded/interaction/perturbed_nonbonded_innerloop.tcc \
				nonbonded/interaction/perturbed_nonbonded_outerloop.h \
				nonbonded/interaction/perturbed_nonbonded_outerloop.tcc \
				nonbonded/interaction/perturbed_nonbonded_pair.h \
				nonbonded/interaction/perturbed_nonbonded_pair.tcc \
				nonbonded/interaction/nonbonded_set.h \
				nonbonded/interaction/nonbonded_set.tcc \
				nonbonded/interaction/nonbonded_interaction.h \
				nonbonded/interaction/nonbonded_interaction.tcc \
				nonbonded/interaction/storage.h \
				nonbonded/filter/chargegroup_grid.h \
				nonbonded/filter/chargegroup_grid.tcc \
				nonbonded/filter/exclusion_filter.h \
				nonbonded/filter/exclusion_filter.tcc \
				nonbonded/filter/filter.h \
				nonbonded/filter/range_filter.h \
				nonbonded/filter/range_filter.tcc \
				nonbonded/pairlist/grid_pairlist_algorithm.h \
				nonbonded/pairlist/grid_pairlist_algorithm.tcc \
				nonbonded/pairlist/pairlist.h \
				nonbonded/pairlist/pairlist.tcc \
				nonbonded/pairlist/pairlist_algorithm.h \
				nonbonded/pairlist/standard_pairlist_algorithm.h \
				nonbonded/pairlist/standard_pairlist_algorithm.tcc \
				special/create_special.h \
				special/position_restraint_interaction.h \
				special/position_restraint_interaction.tcc

noinst_LTLIBRARIES =		libgrointeraction.la

libgrointeraction_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libgrointeraction_la_SOURCES =	interaction.cc \
				forcefield/forcefield.cc \
				forcefield/create_forcefield.cc \
				bonded/create_bonded.cc \
				nonbonded/create_nonbonded.cc \
				special/create_special.cc \
				$(split_files)


bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

