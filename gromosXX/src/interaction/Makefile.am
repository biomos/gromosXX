INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		interaction.h \
                                interaction/nonbonded/storage.h \
				forcefield/forcefield.h \
				forcefield/forcefield.tcc \
				interaction/interaction.h \
				interaction/nonbonded/nonbonded_interaction.h \
				interaction/nonbonded/nonbonded_interaction.tcc \
				interaction/nonbonded/nonbonded_innerloop.h \
				interaction/nonbonded/nonbonded_innerloop.tcc \
				interaction/nonbonded/perturbed_nonbonded_interaction.h \
				interaction/nonbonded/perturbed_nonbonded_interaction.tcc \
				interaction/nonbonded/perturbed_nonbonded_innerloop.h \
				interaction/nonbonded/perturbed_nonbonded_innerloop.tcc \
				interaction/harmonic_bond_interaction.h \
				interaction/harmonic_bond_interaction.tcc \
				interaction/quartic_bond_interaction.h \
				interaction/quartic_bond_interaction.tcc \
                                interaction/angle_interaction.h \
                                interaction/angle_interaction.tcc \
                                interaction/improper_dihedral_interaction.h \
                                interaction/improper_dihedral_interaction.tcc \
				interaction/dihedral_interaction.h \
				interaction/dihedral_interaction.tcc \
				interaction/perturbed_quartic_bond_interaction.h \
				interaction/perturbed_quartic_bond_interaction.tcc \
				interaction/perturbed_harmonic_bond_interaction.h \
				interaction/perturbed_harmonic_bond_interaction.tcc \
				interaction/perturbed_angle_interaction.h \
				interaction/perturbed_angle_interaction.tcc \
				interaction/perturbed_improper_dihedral_interaction.h \
				interaction/perturbed_improper_dihedral_interaction.tcc \
				interaction/perturbed_dihedral_interaction.h \
				interaction/perturbed_dihedral_interaction.tcc \
				interaction/nonbonded/nonbonded_base.h \
				interaction/nonbonded/nonbonded_base.tcc \
				pairlist/pairlist.h \
				pairlist/pairlist.tcc \
				pairlist/pairlist_algorithm.h \
				pairlist/standard_pairlist_algorithm.h \
				pairlist/standard_pairlist_algorithm.tcc \
				pairlist/grid_pairlist_algorithm.h \
				pairlist/grid_pairlist_algorithm.tcc \
				pairlist/filter/filter.h \
				pairlist/filter/exclusion_filter.h \
				pairlist/filter/exclusion_filter.tcc \
				pairlist/filter/chargegroup_grid.h \
				pairlist/filter/chargegroup_grid.tcc \
				pairlist/filter/range_filter.h \
				pairlist/filter/range_filter.tcc \
				pairlist/filter/perturbation_filter.h \
				pairlist/filter/perturbation_filter.tcc


lib_LTLIBRARIES	= 		libinteraction.la

libinteraction_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libinteraction_la_SOURCES =	interaction.cc

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
