AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src

gdatadir = $(datadir)/../contrib
dist_gdata_DATA = schnet_v1.py schnet_v2.py

if WITH_CUDA
  AM_CPPFLAGS += -I$(top_srcdir)/src/cukernel
endif

include_HEADERS =

bin_PROGRAMS = rng_gsl tabulate_spc

rng_gsl_SOURCES = rng_gsl.cc
tabulate_spc_SOURCES = tabulate_spc.cc

LDADD = ../src/libmdpp.la \
	${GSL_LIBS}
if WITH_CUDA
  LDADD += ../src/cukernel/libgrocukernel.a
endif



# EXTRA_DIST =			BUILD_NUMBER
