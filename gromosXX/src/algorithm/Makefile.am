INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		algorithm.h \
				integration/leap_frog.h \
				integration/leap_frog.tcc \
				integration/runge_kutta.h \
				integration/runge_kutta.tcc \
				constraint/shake.h \
				constraint/shake.tcc

lib_LTLIBRARIES	= 		libalgorithm.la

libalgorithm_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libalgorithm_la_SOURCES =	algorithm.cc

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
