INCLUDES = -I$(top_srcdir)/src

include_HEADERS =		flexi.h

bin_PROGRAMS =			md flexi

md_SOURCES =			md.cc
md_LDADD =			../src/simulation/libsimulation.la \
				../src/algorithm/libalgorithm.la \
				../src/interaction/libinteraction.la \
				../src/math/libmath.la \
				../src/io/libio.la

flexi_SOURCES =			flexi.cc
flexi_LDADD =			../src/simulation/libsimulation.la \
				../src/algorithm/libalgorithm.la \
				../src/interaction/libinteraction.la \
				../src/math/libmath.la \
				../src/io/libio.la
