INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX

xxinclude_HEADERS =		configuration_global.h \
				configuration.h \
				energy.h \
				energy_average.h \
				state_properties.h


lib_LTLIBRARIES	=		libgroconfiguration.la

libgroconfiguration_la_LDFLAGS = -version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libgroconfiguration_la_SOURCES = configuration_global.cc	\
				configuration.cc \
				energy.cc \
				energy_average.cc \
				state_properties.cc	


bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

