INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		simulation.h \
				simulation/simulation.h \
				simulation/simulation.tcc \
				simulation/nonbonded.h \
				simulation/nonbonded.tcc \
				simulation/multibath.h \
				simulation/multibath.tcc \
				system/system.h \
				system/system.tcc \
				topology/bond.h \
				topology/bond.tcc \
				topology/angle.h \
				topology/angle.tcc \
                                topology/improper_dihedral.h \
                                topology/improper_dihedral.tcc \
				topology/dihedral.h \
				topology/dihedral.tcc \
				topology/compound.h \
				topology/compound.tcc \
				topology/solute.h \
				topology/solute.tcc \
				topology/solvent.h \
				topology/solvent.tcc \
				topology/topology.h \
				topology/topology.tcc \
				core/atom_iterator.h \
				core/atom_group_iterator.h \
				core/molecule_iterator.h \
				core/chargegroup_iterator.h \
				core/parameter.h

lib_LTLIBRARIES	= 		libsimulation.la

libsimulation_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libsimulation_la_SOURCES =	simulation.cc core.cc

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
