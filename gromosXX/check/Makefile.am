INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		check.h

lib_LTLIBRARIES = 

bin_PROGRAMS =			

check_PROGRAMS = 		simple

simple_SOURCES = 		simple.cc

simple_LDADD =			../src/simulation/libsimulation.la \
				../src/algorithm/libalgorithm.la \
				../src/interaction/libinteraction.la \
				../src/math/libmath.la \
				../src/io/libio.la


simple_CXXFLAGS =		-DBZ_DEBUG

TESTS = simple

