INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX

xxinclude_HEADERS = 	argument.h \
			blockinput.h \
			instream.h \
			message.h \
			print_block.h \
			read_input.h \
			read_special.h \
			parameter/in_parameter.h \
			topology/in_topology.h \
			topology/in_perturbation.h \
			configuration/inframe.h \
			configuration/in_configuration.h \
			configuration/out_configuration.h \
			topology/in_posres.h

lib_LTLIBRARIES	=	libgroio.la

libgroio_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

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
			topology/in_posres.cc

libgroio_la_LIBADD =	../util/libgroutil.la

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

