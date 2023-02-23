AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)

xxincludedir = $(includedir)/md++/io

xxinclude_HEADERS =	argument.h\
	blockinput.h\
	instream.h\
	message.h\
	print_block.h\
	read_input.h\
	read_special.h\
	gzstream.h \
	ifp.h

xxiotopodir = $(includedir)/md++/io/topology
xxiotopo_HEADERS = topology/in_jvalue.h\
	topology/in_posres.h\
	topology/in_topology.h\
	topology/in_perturbation.h\
	topology/in_distanceres.h\
	topology/in_angrest.h\
	topology/in_dihrest.h\
	topology/in_friction.h\
	topology/in_xray.h\
	topology/in_gamd.h\
	topology/in_order.h\
	topology/in_leus.h\
	topology/in_bsleus.h\
	topology/in_qmmm.h\
	topology/in_symrest.h\
	topology/in_rdc.h

xxioconfdir = $(includedir)/md++/io/configuration
xxioconf_HEADERS = configuration/in_configuration.h\
  configuration/out_configuration.h\
  configuration/inframe.h

xxioparamdir = $(includedir)/md++/io/parameter
xxioparam_HEADERS = parameter/in_parameter.h\
  parameter/check_parameter.h


noinst_LTLIBRARIES =libgroio.la libgroiocheck.la

libgroiocheck_la_SOURCES = parameter/check_parameter.cc 
libgroiocheck_la_CXXFLAGS = @CXX_NDEBUG_FLAGS@

libgroio_la_SOURCES =	io.cc\
	argument.cc\
	blockinput.cc\
	instream.cc\
	message.cc\
	print_block.cc\
	read_input.cc\
	read_special.cc\
	gzstream.cc\
	parameter/in_parameter.cc\
	topology/in_topology.cc\
	topology/in_perturbation.cc\
	configuration/inframe.cc\
	configuration/in_configuration.cc\
	configuration/out_configuration.cc\
	topology/in_posres.cc\
	topology/in_distanceres.cc\
	topology/in_angrest.cc\
	topology/in_dihrest.cc\
	topology/in_jvalue.cc\
	topology/in_friction.cc\
	topology/in_xray.cc\
	topology/in_gamd.cc\
	topology/in_order.cc\
	topology/in_leus.cc\
	topology/in_bsleus.cc\
	topology/in_qmmm.cc\
	topology/in_symrest.cc\
	topology/in_rdc.cc

libgroio_la_LIBADD = libgroiocheck.la

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

