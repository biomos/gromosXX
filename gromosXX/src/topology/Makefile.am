INCLUDES = -I$(top_srcdir)/src

include_HEADERS =	core/atomgroup_iterator.h	core/chargegroup_iterator.h \
			core/atom_iterator.h		core/compound.h \
			core/body_term.h 		core/molecule_iterator.h \
			core/core.h\
			solute.h			solvent.h \
			perturbed_solute.h \
			topology.h			perturbed_atom.h

lib_LTLIBRARIES	=	libtopology.la

libtopology_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libtopology_la_SOURCES =	topology.cc	topology_global.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

