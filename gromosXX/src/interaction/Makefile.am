INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		forcefield/forcefield.h \
				forcefield/forcefield.tcc \
				interaction/interaction.h \
				interaction/nonbonded_interaction.h \
				interaction/nonbonded_interaction.tcc \
				interaction/harmonic_bond_interaction.h \
				interaction/harmonic_bond_interaction.tcc \
				pairlist/simple_pairlist.h \
				pairlist/simple_pairlist.tcc

lib_LTLIBRARIES	= 		libinteraction.la

libinteraction_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libinteraction_la_SOURCES =	interaction.cc

bin_PROGRAMS =

check_PROGRAMS = 		interaction.t

interaction_t_SOURCES =		interaction.t.cc

interaction_t_LDADD = 		libinteraction.la \
				../simulation/libsimulation.la

interaction_t_LDFLAGS =		-lboost_unit_test_framework

interaction_t_CXXFLAGS =	-DBZ_DEBUG -Wno-unused

TESTS = interaction.t
