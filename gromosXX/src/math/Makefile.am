INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		gmath.h \
						boundary_implementation.h \
						boundary_implementation.tcc \
						periodicity.h \
						periodicity.tcc \
						volume.h

lib_LTLIBRARIES	= 		libmath.la

libmath_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libmath_la_SOURCES =	math.cc \
			volume.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

