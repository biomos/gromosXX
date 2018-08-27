AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/util -I$(top_srcdir)/src/util/replicaExchange

xxincludedir = $(includedir)/md++/util
xxinclude_HEADERS = error.h\
	bs_coordinate.h\
	bs_potentials.h\
	bs_subspace.h\
	bs_umbrella.h\
	bs_vector.h\
	cycle_thread.h\
	debug.h\
	timing.h\
	usage.h\
	generate_velocities.h\
	parse_tcouple.h\
	parse_verbosity.h\
	prepare_virial.h\
	create_simulation.h\
	template_split.h\
	virtual_atom.h\
	coding.h\
	feature.h\
	le_coordinate.h\
	umbrella.h\
	umbrella_weight.h\
	thread.h\
	system_call.h\
	undirected_graph.h


#include replica_exchange Files	
xxutilrepexdir = $(includedir)/md++/util/replicaExchange
xxutilrepex_HEADERS = replicaExchange/repex_mpi.h\
        replicaExchange/replica_data.h\
	replicaExchange/replica.h\
	replicaExchange/replica_exchange_base.h\
	replicaExchange/replica_exchange_master.h\
	replicaExchange/replica_exchange_slave.h
	
noinst_LTLIBRARIES = libgroutil.la

libgroutil_la_SOURCES =	cycle_thread.cc\
	bs_coordinate.cc\
	bs_potentials.cc\
	bs_subspace.cc\
	bs_umbrella.cc\
	bs_vector.cc\
	debug.cc\
	timing.cc\
	usage.cc\
	generate_velocities.cc\
	parse_tcouple.cc\
	parse_verbosity.cc\
	prepare_virial.cc\
	create_simulation.cc\
	virtual_atom.cc\
	coding.cc\
	feature.cc\
	le_coordinate.cc\
	umbrella.cc\
	umbrella_weight.cc\
	thread.cc\
	system_call.cc\
	undirected_graph.cc\
	replicaExchange/replica.cc\
	replicaExchange/replica_exchange_base.cc\
	replicaExchange/replica_exchange_master.cc\
	replicaExchange/replica_exchange_slave.cc

libgroutil_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

