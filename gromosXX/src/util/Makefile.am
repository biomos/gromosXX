INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 	stdheader.h \
			error.h \
			debug.h	\
			timing.h \
			generate_velocities.h \
			parse_tcouple.h \
			parse_verbosity.h \
			prepare_virial.h

lib_LTLIBRARIES	= 	libutil.la

libutil_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libutil_la_SOURCES =	debug.cc \
			timing.cc \
			generate_velocities.cc \
			parse_tcouple.cc \
			parse_verbosity.cc \
			prepare_virial.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

