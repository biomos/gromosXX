SUBDIRS = util simulation math topology configuration algorithm interaction io check .

noinst_HEADERS = 	stdheader.h

lib_LTLIBRARIES = 	libgroxx.la

libgroxx_la_LDFLAGS =	-no-undefined -version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0 $(GSL_LDFLAGS)

libgroxx_la_SOURCES =	

libgroxx_la_LIBADD =	algorithm/libgroalgorithm.la \
			configuration/libgroconfiguration.la \
			interaction/libgrointeraction.la \
			math/libgromath.la \
			simulation/libgrosimulation.la \
			topology/libgrotopology.la \
			util/libgroutil.la \
			io/libgroio.la \
			$(GSL_LIB)
