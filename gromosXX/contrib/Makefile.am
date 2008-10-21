AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src

include_HEADERS =

bin_PROGRAMS =			rep_rewrite \
                                rng_gsl\
                                tabulate_spc

rep_rewrite_SOURCES =		rep_rewrite.cc
rng_gsl_SOURCES =               rng_gsl.cc
tabulate_spc_SOURCES = tabulate_spc.cc

LDADD = ../src/libgroxx.la $(extra_lib)

# EXTRA_DIST =			BUILD_NUMBER



