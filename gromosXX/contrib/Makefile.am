AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src


if WITH_CUDA
  AM_CPPFLAGS += -I$(top_srcdir)/src/cukernel
endif

include_HEADERS =

bin_PROGRAMS = rng_gsl tabulate_spc

rng_gsl_SOURCES = rng_gsl.cc
tabulate_spc_SOURCES = tabulate_spc.cc
rng_gsl_LDFLAGS = ${GSL_LIBS}

if WITH_CUDA
  CUKERNEL = ../src/cukernel/libgrocukernel.a
else
  CUKERNEL = ../src/cukernel/libgrocukernel.la
endif
LDADD = $(CUKERNEL) ../src/libmdpp.la

# EXTRA_DIST =			BUILD_NUMBER

