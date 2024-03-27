AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)

xxincludedir = $(includedir)/md++/simulation
xxinclude_HEADERS =	simulation.h\
	multibath.h\
	mpiControl.h\
	parameter.h

noinst_LTLIBRARIES = libgrosimulation.la

libgrosimulation_la_SOURCES =	simulation.cc\
	mpiControl.cpp\
	multibath.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

