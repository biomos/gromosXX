INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/interaction -I$(top_srcdir)/src/interaction/nonbonded -I$(top_srcdir)/src/interaction/nonbonded/pairlist -I$(top_srcdir)/src/interaction/nonbonded/interaction -I$(top_srcdir)/src/interaction/bonded -I$(top_srcdir)/src/interaction/forcefield -I$(top_srcdir)/src/interaction/special


xxincludedir = 			$(includedir)/gromosXX/interaction

xxinclude_HEADERS =		interaction.h \
				interaction_types.h \
				forcefield/forcefield.h \
				forcefield/create_forcefield.h \
				bonded/create_bonded.h \
				bonded/quartic_bond_interaction.h \
				bonded/harmonic_bond_interaction.h \
				bonded/angle_interaction.h \
				bonded/improper_dihedral_interaction.h \
				bonded/dihedral_interaction.h \
				bonded/perturbed_quartic_bond_interaction.h \
				bonded/perturbed_harmonic_bond_interaction.h \
				bonded/perturbed_angle_interaction.h \
				bonded/perturbed_improper_dihedral_interaction.h \
				bonded/perturbed_dihedral_interaction.h \
				nonbonded/create_nonbonded.h \
				nonbonded/interaction_spec.h \
				nonbonded/interaction/nonbonded_parameter.h \
				nonbonded/interaction/nonbonded_term.h \
				nonbonded/interaction/nonbonded_term.cc \
				nonbonded/interaction/nonbonded_innerloop.h \
				nonbonded/interaction/nonbonded_innerloop.cc \
				nonbonded/interaction/nonbonded_outerloop.h \
				nonbonded/interaction/perturbed_nonbonded_term.h \
				nonbonded/interaction/perturbed_nonbonded_term.cc \
				nonbonded/interaction/perturbed_nonbonded_innerloop.h \
				nonbonded/interaction/perturbed_nonbonded_innerloop.cc \
				nonbonded/interaction/perturbed_nonbonded_outerloop.h \
				nonbonded/interaction/perturbed_nonbonded_pair.h \
				nonbonded/interaction/nonbonded_set.h \
				nonbonded/interaction/perturbed_nonbonded_set.h \
				nonbonded/interaction/nonbonded_interaction.h \
				nonbonded/interaction/storage.h \
				nonbonded/pairlist/pairlist.h \
				nonbonded/pairlist/pairlist_algorithm.h \
				nonbonded/pairlist/standard_pairlist_algorithm.h \
				special/create_special.h \
				special/position_restraint_interaction.h \
				special/jvalue_restraint_interaction.h \
				special/pscale.h


noinst_LTLIBRARIES =		libgrointeraction.la

libgrointeraction_la_SOURCES =	interaction.cc \
				forcefield/forcefield.cc \
				forcefield/create_forcefield.cc \
				nonbonded/create_nonbonded.cc \
				nonbonded/pairlist/pairlist.cc \
				nonbonded/pairlist/standard_pairlist_algorithm.cc \
				nonbonded/pairlist/standard_pairlist_algorithm_atomic.cc \
				nonbonded/interaction/nonbonded_interaction.cc \
				nonbonded/interaction/nonbonded_set.cc \
				nonbonded/interaction/nonbonded_outerloop.cc \
				nonbonded/interaction/perturbed_nonbonded_set.cc \
				nonbonded/interaction/perturbed_nonbonded_outerloop.cc \
				nonbonded/interaction/perturbed_nonbonded_pair.cc \
				bonded/create_bonded.cc \
				bonded/quartic_bond_interaction.cc \
				bonded/perturbed_quartic_bond_interaction.cc \
				bonded/harmonic_bond_interaction.cc \
				bonded/perturbed_harmonic_bond_interaction.cc \
				bonded/angle_interaction.cc \
				bonded/perturbed_angle_interaction.cc \
				bonded/improper_dihedral_interaction.cc \
				bonded/perturbed_improper_dihedral_interaction.cc \
				bonded/dihedral_interaction.cc \
				bonded/perturbed_dihedral_interaction.cc \
				special/create_special.cc \
				special/position_restraint_interaction.cc \
				special/jvalue_restraint_interaction.cc \
				special/pscale.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =





