INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 	argument.h \
			blockinput.h \
			instream.h \
			message.h \
			read_input.h \
			parameter/in_parameter.h \
			topology/in_topology.h \
			topology/in_perturbation.h \
			configuration/inframe.h \
			configuration/in_configuration.h \
			configuration/out_configuration.h

lib_LTLIBRARIES	=			libio.la

libio_la_LDFLAGS =	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libio_la_SOURCES =	io.cc \
			argument.cc \
			blockinput.cc \
			instream.cc \
			message.cc \
			read_input.cc \
			parameter/in_parameter.cc \
			topology/in_topology.cc \
			topology/in_perturbation.cc \
			configuration/inframe.cc \
			configuration/in_configuration.cc \
			configuration/out_configuration.cc

libio_la_LIBADD =	../util/libutil.la

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

