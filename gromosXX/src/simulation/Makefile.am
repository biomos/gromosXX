INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		simulation.h \
				core.h \
				perturbation.h \
				simulation/simulation.h \
				simulation/simulation.tcc \
				simulation/nonbonded.h \
				simulation/nonbonded.tcc \
				simulation/multibath.h \
				simulation/multibath.tcc \
				system/system.h \
				system/system.tcc \
				system/energy.h \
				system/energy.tcc \
				system/energy_average.h \
				system/energy_average.tcc \
				topology/bond.h \
				topology/angle.h \
                                topology/improper_dihedral.h \
				topology/dihedral.h \
				topology/compound.h \
				topology/compound.tcc \
				topology/solute.h \
				topology/solute.tcc \
				topology/solvent.h \
				topology/solvent.tcc \
				topology/topology.h \
				topology/topology.tcc \
				topology/perturbation_topology.h \
				topology/perturbation_topology.tcc \
				topology/perturbed_atom.h \
				topology/perturbed_atom.tcc \
				topology/perturbed_atompair.h \
				topology/perturbed_bond.h \
				topology/perturbed_angle.h \
				topology/perturbed_improper_dihedral.h \
				topology/perturbed_dihedral.h \
				topology/perturbed_solute.h \
				topology/perturbed_solute.tcc \
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
