INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		simulation/simulation.h \
				simulation/simulation.tcc \
				system/system.h \
				system/system.tcc \
				topology/topology.h \
				topology/topology.tcc

lib_LTLIBRARIES	= 		libsimulation.la

libsimulation_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libsimulation_la_SOURCES =	simulation.cc

bin_PROGRAMS =

check_PROGRAMS = 		simulation.t

simulation_t_SOURCES =		simulation.t.cc

simulation_t_LDADD = 		libsimulation.la

simulation_t_LDFLAGS =		-lboost_unit_test_framework

simulation_t_CXXFLAGS =		-DBZ_DEBUG -Wno-unused-variable

TESTS = simulation.t
