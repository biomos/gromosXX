INCLUDES = -I$(top_srcdir)/src

xxincludedir = 			$(includedir)/gromosXX/interaction
xxffincludedir = 		$(includedir)/gromosXX/interaction/forcefield
xxbincludedir =			$(includedir)/gromosXX/interaction/bonded
xxnbincludedir =		$(includedir)/gromosXX/interaction/nonbonded
xxspincludedir =		$(includedir)/gromosXX/interaction/special

xxinclude_HEADERS =		interaction.h \
				interaction_types.h

xxffinclude_HEADERS =		forcefield/forcefield.h \
				forcefield/create_forcefield.h

xxbinclude_HEADERS =		bonded/create_bonded.h \
				bonded/quartic_bond_interaction.h \
				bonded/harmonic_bond_interaction.h \
				bonded/angle_interaction.h \
				bonded/improper_dihedral_interaction.h \
				bonded/dihedral_interaction.h \
				bonded/perturbed_quartic_bond_interaction.h \
				bonded/perturbed_harmonic_bond_interaction.h \
				bonded/perturbed_angle_interaction.h \
				bonded/perturbed_improper_dihedral_interaction.h \
				bonded/perturbed_dihedral_interaction.h

xxnbinclude_HEADERS =		nonbonded/create_nonbonded.h \
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
				nonbonded/pairlist/standard_pairlist_algorithm.h

xxspinclude_HEADERS =		special/create_special.h \
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
				nonbonded/pairlist/perturbed_standard_pairlist_algorithm.cc \
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





