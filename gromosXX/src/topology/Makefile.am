INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX

xxinclude_HEADERS =	core/atomgroup_iterator.h	core/chargegroup_iterator.h \
			core/atom_iterator.h		core/compound.h \
			core/body_term.h 		core/molecule_iterator.h \
			core/core.h\
			solute.h			solvent.h \
			perturbed_solute.h \
			topology.h			perturbed_atom.h

noinst_LTLIBRARIES =	libgrotopology.la

libgrotopology_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libgrotopology_la_SOURCES =	topology.cc	topology_global.cc

libgrotopology_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

