AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/md++/math
xxinclude_HEADERS = gmath.h\
	boundary_implementation.h\
	boundary_implementation.cc\
	periodicity.h\
	periodicity.cc\
	volume.h\
	random.h\
	transformation.h\
	boundary_checks.h\
    fft.h

noinst_LTLIBRARIES = libgromath.la

libgromath_la_SOURCES =	math.cc\
	volume.cc\
	random.cc\
	transformation.cc\
	boundary_checks.cc

libgromath_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

