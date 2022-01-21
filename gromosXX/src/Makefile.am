AM_CXXFLAGS = @MY_CXXFLAGS@

if WITH_CUDA
  CUKERNEL = cukernel
endif

SUBDIRS = $(CUKERNEL) util simulation math topology configuration algorithm interaction io check replicaExchange .

noinst_HEADERS = stdheader.h

lib_LTLIBRARIES = libmdpp.la

libmdpp_la_LDFLAGS =	-no-undefined -release $(MD_VERSION)

libmdpp_la_SOURCES =	

libmdpp_la_LIBADD = algorithm/libgroalgorithm.la\
	configuration/libgroconfiguration.la\
	interaction/libgrointeraction.la\
	math/libgromath.la\
	simulation/libgrosimulation.la\
	topology/libgrotopology.la\
    replicaExchange/libgrorepex.la\
	util/libgroutil.la\
	io/libgroio.la\
	$(GSL_LIBS)\
	$(CLIPPER_LIBS)
