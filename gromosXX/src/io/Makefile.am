INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		io.h \
				blockinput.h \
				blockinput.tcc \
				GInStream.h \
				topology/InTopology.h \
				topology/InTopology.tcc

lib_LTLIBRARIES	= 		libio.la

libio_la_LDFLAGS =	 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libio_la_SOURCES =		io.cc \
				GInStream.cc

bin_PROGRAMS =			

check_PROGRAMS = 		blockinput.t

blockinput_t_SOURCES =		blockinput.t.cc

blockinput_t_CXXFLAGS =     	-DBZ_DEBUG

blockinput_t_LDFLAGS =		-lboost_unit_test_framework

TESTS = 			blockinput.t
