INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		gmath.h

lib_LTLIBRARIES	= 		libmath.la

libmath_la_LDFLAGS =	 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libmath_la_SOURCES =		math.cc

bin_PROGRAMS =			

check_PROGRAMS = 		math.t

math_t_SOURCES =		math.t.cc

math_t_CXXFLAGS =		-DBZ_DEBUG

TESTS = math.t

