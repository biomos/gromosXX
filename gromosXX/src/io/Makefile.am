INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		io.h \
				blockinput.h \
				blockinput.tcc \
				message.h \
				GInStream.h \
				topology/InTopology.h \
				topology/InTopology.tcc \
				topology/InPerturbationTopology.h \
				topology/InPerturbationTopology.tcc \
				trajectory/InTrajectory.h \
				trajectory/InTrajectory.tcc \
				trajectory/OutTrajectory.h \
				trajectory/OutTrajectory.tcc \
				trajectory/OutG96Trajectory.h \
				trajectory/OutG96Trajectory.tcc \
				trajectory/InFlexibleConstraints.h \
				trajectory/InFlexibleConstraints.tcc \
				input/InInput.h \
				input/InInput.tcc \
				argument.h \
				print_block.h \
				print_block.tcc

lib_LTLIBRARIES	= 		libio.la

libio_la_LDFLAGS =	 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libio_la_SOURCES =		io.cc \
				message.cc \
				GInStream.cc \
				argument.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

