INCLUDES = -I$(top_srcdir)/src

include_HEADERS =		create_forcefield.h

bin_PROGRAMS =			vtest

vtest_SOURCES =			vtest.cc
vtest_LDADD =			../src/simulation/libsimulation.la \
				../src/algorithm/libalgorithm.la \
				../src/interaction/libinteraction.la \
				../src/math/libmath.la \
				../src/io/libio.la
