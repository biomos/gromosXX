INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		simulation.h \
				simulation/simulation.h \
				simulation/simulation.tcc \
				system/system.h \
				system/system.tcc \
				topology/bond.h \
				topology/bond.tcc \
				topology/soluteatom.h \
				topology/soluteatom.tcc \
				topology/solvent.h \
				topology/solvent.tcc \
				topology/topology.h \
				topology/topology.tcc

lib_LTLIBRARIES	= 		libsimulation.la

libsimulation_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libsimulation_la_SOURCES =	simulation.cc

bin_PROGRAMS =

check_PROGRAMS = 		simulation.t

simulation_t_SOURCES =		simulation.t.cc

simulation_t_LDADD = 		libsimulation.la

simulation_t_CXXFLAGS =		-DBZ_DEBUG

TESTS = simulation.t
