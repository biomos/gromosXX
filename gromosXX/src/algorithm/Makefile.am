INCLUDES = -I$(top_srcdir)/src

include_HEADERS =	algorithm.h \
			algorithm/algorithm_sequence.h \
			integration/leap_frog.h \
			temperature/temperature_calculation.h \
			temperature/berendsen.h \
			create_md_sequence.h \
			constraints/shake.h \
			constraints/shake.tcc

lib_LTLIBRARIES	=	libalgorithm.la

libalgorithm_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libalgorithm_la_SOURCES =	algorithm.cc \
				create_md_sequence.cc \
				algorithm/algorithm_sequence.cc \
				integration/leap_frog.cc \
				temperature/temperature_calculation.cc \
				temperature/berendsen.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

