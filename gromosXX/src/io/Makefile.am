INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 			io.h \
							blockinput.h

lib_LTLIBRARIES	= 			libio.la

libio_la_LDFLAGS =	 		-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libio_la_SOURCES =			io.cc

bin_PROGRAMS =			

check_PROGRAMS = 			blockinput.t

blockinput_t_SOURCES =		blockinput.t.cc

blockinput_t_CXXFLAGS =     -DBZ_DEBUG

TESTS = blockinput.t
