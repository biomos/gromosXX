AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)

xxincludedir = $(includedir)/md++/topology
xxcoreincludedir = $(includedir)/md++/topology/core

xxcoreinclude_HEADERS =	core/atomgroup_iterator.h\
	core/chargegroup_iterator.h\
	core/atom_iterator.h\
	core/compound.h\
	core/body_term.h\
	core/molecule_iterator.h\
	core/temperaturegroup_iterator.h\
	core/pressuregroup_iterator.h\
	core/core.h\
	core/qm_atom.h

xxinclude_HEADERS =	solute.h\
	solvent.h\
	perturbed_atom.h\
	perturbed_solute.h\
	topology.h\
	sd.h

noinst_LTLIBRARIES = libgrotopology.la

libgrotopology_la_SOURCES =	topology.cc\
	topology_global.cc

libgrotopology_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

