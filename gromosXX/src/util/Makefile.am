AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX/util
xxinclude_HEADERS = error.h\
	debug.h\
	timing.h\
	usage.h\
	generate_velocities.h\
	parse_tcouple.h\
	parse_verbosity.h\
	prepare_virial.h\
	create_simulation.h\
	replica_data.h\
	replica_exchange.h\
	template_split.h\
	virtual_atom.h\
	virtual_grain.h\
	coding.h\
	feature.h\
	le_coordinate.h\
	umbrella.h\
	umbrella_weight.h

noinst_LTLIBRARIES = libgroutil.la

libgroutil_la_SOURCES =	debug.cc\
	timing.cc\
	usage.cc\
	generate_velocities.cc\
	parse_tcouple.cc\
	parse_verbosity.cc\
	prepare_virial.cc\
	replica_exchange.cc\
	replica_exchange_master.cc\
	replica_exchange_slave.cc\
	replica_exchange_control.cc\
	create_simulation.cc\
	virtual_atom.cc\
	virtual_grain.cc\
	coding.cc\
	feature.cc\
	le_coordinate.cc\
	umbrella.cc\
	umbrella_weight.cc

libgroutil_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

