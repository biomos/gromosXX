INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		io.h \
				blockinput.h \
				blockinput.tcc \
				message.h \
				GInStream.h \
				topology/InTopology.h \
				topology/InTopology.tcc \
				trajectory/InTrajectory.h \
				trajectory/InTrajectory.tcc \
				trajectory/OutTrajectory.h \
				trajectory/OutTrajectory.tcc \
				input/InInput.h \
				input/InInput.tcc \
				argument.h

lib_LTLIBRARIES	= 		libio.la

libio_la_LDFLAGS =	 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libio_la_SOURCES =		io.cc \
				message.cc \
				GInStream.cc \
				argument.cc

bin_PROGRAMS =			

check_PROGRAMS = 		blockinput.t \
				message.t

blockinput_t_SOURCES =		blockinput.t.cc

blockinput_t_CXXFLAGS =     	-DBZ_DEBUG

blockinput_t_LDADD =		libio.la \
				../simulation/libsimulation.la \
				../interaction/libinteraction.la

message_t_SOURCES = 		message.t.cc

message_t_CXXFLAGS =		-DBZ_DEBUG

message_t_LDADD =		libio.la

TESTS = 			blockinput.t \
				message.t
