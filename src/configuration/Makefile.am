AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)

xxincludedir = $(includedir)/md++/configuration
xxinclude_HEADERS =	configuration_global.h\
	configuration.h\
	energy.h\
	average.h\
	state_properties.h\
	mesh.h\
	kspace.h\
	influence_function.h

noinst_LTLIBRARIES = libgroconfiguration.la

libgroconfiguration_la_SOURCES = configuration_global.cc\
	configuration.cc\
	energy.cc\
	average.cc\
	state_properties.cc	\
	mesh.cc\
	kspace.cc\
	influence_function.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

