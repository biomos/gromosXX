INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		interaction.h \
                                interaction/storage.h \
				forcefield/forcefield.h \
				forcefield/forcefield.tcc \
				interaction/interaction.h \
				interaction/nonbonded_interaction.h \
				interaction/nonbonded_interaction.tcc \
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
				interaction/perturbed_angle_interaction.h \
				interaction/perturbed_angle_interaction.tcc \
				interaction/nonbonded_base.h \
				interaction/nonbonded_base.tcc \
				interaction/nonbonded_inner_loop.h \
				interaction/nonbonded_inner_loop.tcc \
				interaction/nonbonded_inner_loop_virial.h \
				interaction/nonbonded_inner_loop_virial.tcc \
				pairlist/basic_pairlist.h \
				pairlist/basic_pairlist.tcc \
				pairlist/basic_pairlist_algorithm.h \
				pairlist/basic_pairlist_algorithm.tcc \
				pairlist/basic_filter.h \
				pairlist/basic_filter.tcc \
				pairlist/perturbation_filter.h \
				pairlist/perturbation_filter.tcc \
				pairlist/twinrange_filter.h \
				pairlist/twinrange_filter.tcc \
				pairlist/twinrange_chargegroup_filter.h \
				pairlist/twinrange_chargegroup_filter.tcc \
				pairlist/chargegroup_range_pairlist_algorithm.h \
				pairlist/chargegroup_range_pairlist_algorithm.tcc

lib_LTLIBRARIES	= 		libinteraction.la

libinteraction_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libinteraction_la_SOURCES =	interaction.cc

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
