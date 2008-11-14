AM_CXXFLAGS = @MY_CXXFLAGS@
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/interaction -I$(top_srcdir)/src/interaction/nonbonded -I$(top_srcdir)/src/interaction/nonbonded/pairlist -I$(top_srcdir)/src/interaction/nonbonded/interaction -I$(top_srcdir)/src/interaction/bonded -I$(top_srcdir)/src/interaction/forcefield -I$(top_srcdir)/src/interaction/special -I${top_srcdir}/src/algorithm

check_PROGRAMS =  aladip.t\
	aladip_atomic.t\
	aladip_special.t\
	c16_cg.t\
	lambdas.t

aladip_t_SOURCES = aladip.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc

aladip_t_LDADD = ../libgroxx.la $(extra_lib)

aladip_special_t_SOURCES = aladip_special.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc


aladip_special_t_LDADD = ../libgroxx.la $(extra_lib)

aladip_atomic_t_SOURCES = aladip_atomic.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc

aladip_atomic_t_LDADD = ../libgroxx.la $(extra_lib)

c16_cg_t_SOURCES = c16_cg.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc

c16_cg_t_LDADD = ../libgroxx.la $(extra_lib)

lambdas_t_SOURCES = lambdas.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc

lambdas_t_LDADD = ../libgroxx.la $(extra_lib)

TESTS = aladip.t aladip_atomic.t aladip_special.t c16_cg.t lambdas.t

EXTRA_DIST = data/aladip.topo\
	data/aladip.pttopo\
	data/aladip.conf\
	data/aladip.in\
	data/aladip_special.in\
	data/c16_cg.in\
        data/c16_cg.topo\
        data/c16_cg.conf\
	data/aladip.distrest\
	data/aladip.dihrest\
	data/aladip_atomic.in\
	data/lambdas.off.in\
	data/lambdas.on.in\
	data/lambdas.pttopo

