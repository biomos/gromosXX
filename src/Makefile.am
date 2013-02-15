AM_CXXFLAGS = @MY_CXXFLAGS@

if CUDA

SUBDIRS = util simulation math topology configuration algorithm interaction io check cuda .

else

SUBDIRS = util simulation math topology configuration algorithm interaction io check .

endif

noinst_HEADERS = stdheader.h

lib_LTLIBRARIES = libmdpp.la

libmdpp_la_LDFLAGS =	-no-undefined -release $(MD_VERSION)

libmdpp_la_SOURCES =	

if CUDA

libmdpp_la_LIBADD = algorithm/libgroalgorithm.la\
	configuration/libgroconfiguration.la\
	interaction/libgrointeraction.la\
	math/libgromath.la\
	simulation/libgrosimulation.la\
	topology/libgrotopology.la\
	util/libgroutil.la\
	io/libgroio.la\
	cuda/libgrocuda.la\
	$(GSL_LIBS)\
	$(CLIPPER_LIBS)

else

libmdpp_la_LIBADD = algorithm/libgroalgorithm.la\
	configuration/libgroconfiguration.la\
	interaction/libgrointeraction.la\
	math/libgromath.la\
	simulation/libgrosimulation.la\
	topology/libgrotopology.la\
	util/libgroutil.la\
	io/libgroio.la\
	$(GSL_LIBS)\
	$(CLIPPER_LIBS)

endif