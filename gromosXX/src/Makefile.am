AUTOMAKE_OPTIONS = foreign

SUBDIRS = simulation util math topology configuration algorithm interaction io .

lib_LTLIBRARIES = 	libgroxx.la

libgroxx_la_LDFLAGS =	-no-undefined -version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libgroxx_la_SOURCES =	

libgroxx_la_LIBADD =	algorithm/libgroalgorithm.la \
			configuration/libgroconfiguration.la \
			interaction/libgrointeraction.la \
			math/libgromath.la \
			simulation/libgrosimulation.la \
			topology/libgrotopology.la \
			util/libgroutil.la \
			io/libgroio.la \
			-lblitz


