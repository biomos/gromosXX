AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)

xxincludedir = $(includedir)/gromosXX/math
xxinclude_HEADERS = 		gmath.h \
				boundary_implementation.h \
				boundary_implementation.cc \
				periodicity.h \
				periodicity.cc \
				volume.h \
                                random.h

noinst_LTLIBRARIES = 		libgromath.la

libgromath_la_SOURCES =		math.cc \
				volume.cc \
                                random.cc

libgromath_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

