INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX

xxinclude_HEADERS = 		gmath.h \
				boundary_implementation.h \
				boundary_implementation.tcc \
				periodicity.h \
				periodicity.tcc \
				volume.h

noinst_LTLIBRARIES = 		libgromath.la

libgromath_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libgromath_la_SOURCES =		math.cc \
				volume.cc

libgromath_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

