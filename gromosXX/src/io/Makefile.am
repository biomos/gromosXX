INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		blockinput.h

lib_LTLIBRARIES	= 		libio.la

libio_la_LDFLAGS =	 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libio_la_SOURCES =		

bin_PROGRAMS =			

check_PROGRAMS = 		blockinput.t

io_t_SOURCES =		blockinput.t.cc

io_t_CXXFLAGS =               -DBZ_DEBUG

TESTS = io.t
