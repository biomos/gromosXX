INCLUDES = -I$(top_srcdir)/src

include_HEADERS =		create_forcefield.h

bin_PROGRAMS =			md

md_SOURCES =			md.cc
md_LDADD =			../src/simulation/libsimulation.la \
				../src/algorithm/libalgorithm.la \
				../src/interaction/libinteraction.la \
				../src/math/libmath.la \
				../src/io/libio.la
