INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX/io

xxinclude_HEADERS = 	argument.h \
			blockinput.h \
			instream.h \
			message.h \
			print_block.h \
			read_input.h \
			read_special.h \
			parameter/in_parameter.h \
			ifp.h \
			forcefield/in_ifp.h

xxiotopodir   = $(includedir)/gromosXX/io/topology
xxiotopo_HEADERS =      topology/in_jvalue.h \
                        topology/in_posres.h \
                        topology/in_topology.h \
                        topology/in_perturbation.h \
                        topology/in_distanceres.h \
                        topology/in_dihrest.h \
                        topology/in_friction.h

xxioconfdir   = $(includedir)/gromosXX/io/configuration
xxioconf_HEADERS =      configuration/in_configuration.h \
                        configuration/out_configuration.h \
                        configuration/inframe.h

xxioparamdir  = $(includedir)/gromosXX/io/parameter
xxioparam_HEADERS =     parameter/in_parameter.h

xxioffdir     = $(includedir)/gromosXX/io/forcefield
xxioff_HEADERS =        forcefield/in_ifp.h

noinst_LTLIBRARIES =	libgroio.la

libgroio_la_SOURCES =	io.cc \
			argument.cc \
			blockinput.cc \
			instream.cc \
			message.cc \
			print_block.cc \
			read_input.cc \
			read_special.cc \
			parameter/in_parameter.cc \
			topology/in_topology.cc \
			topology/in_perturbation.cc \
			configuration/inframe.cc \
			configuration/in_configuration.cc \
			configuration/out_configuration.cc \
			topology/in_posres.cc \
			topology/in_distanceres.cc \
			topology/in_dihrest.cc \
			topology/in_jvalue.cc \
                        topology/in_friction.cc \
			forcefield/in_ifp.cc

libgroio_la_LIBADD =

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

