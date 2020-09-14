AM_CXXFLAGS = @MY_CXXFLAGS@

if WITH_CUDA
  CUKERNEL = cukernel
endif

SUBDIRS = $(CUKERNEL) util simulation math topology configuration algorithm interaction io check .

noinst_HEADERS = stdheader.h

lib_LTLIBRARIES = libmdpp.la

libmdpp_la_LDFLAGS =	-no-undefined -release $(MD_VERSION)

libmdpp_la_SOURCES =	

libmdpp_la_LIBADD = cukernel/libgrocukernel.a\
	algorithm/libgroalgorithm.la\
	configuration/libgroconfiguration.la\
	interaction/libgrointeraction.la\
	math/libgromath.la\
	simulation/libgrosimulation.la\
	topology/libgrotopology.la\
	util/libgroutil.la\
	io/libgroio.la\
	$(GSL_LIBS)\
	$(CLIPPER_LIBS)

##libmdpp_la_LDADD = cukernel/libgrocukernel.la
