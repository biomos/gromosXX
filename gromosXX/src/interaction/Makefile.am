INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		interaction.h \
				forcefield/forcefield.h \
				forcefield/forcefield.tcc \
				forcefield/parameter.h \
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
				pairlist/simple_pairlist.h \
				pairlist/simple_pairlist.tcc \
				pairlist/twin_range_pairlist.h \
				pairlist/twin_range_pairlist.tcc \
				pairlist/twin_range_pairlist_cg.h \
				pairlist/twin_range_pairlist_cg.tcc

lib_LTLIBRARIES	= 		libinteraction.la

libinteraction_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libinteraction_la_SOURCES =	interaction.cc

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
