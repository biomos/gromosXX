AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/interaction -I$(top_srcdir)/src/interaction/nonbonded -I$(top_srcdir)/src/interaction/nonbonded/pairlist -I$(top_srcdir)/src/interaction/nonbonded/interaction -I$(top_srcdir)/src/interaction/bonded -I$(top_srcdir)/src/interaction/forcefield -I$(top_srcdir)/src/interaction/special -I${top_srcdir}/src/algorithm

check_PROGRAMS =  aladip.t\
	aladip_unperturbed.t\
	aladip_atomic.t\
	aladip_special.t\
	aladip_ls.t\
	c16_cg.t\
	lambdas.t

aladip_t_SOURCES = aladip.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc
aladip_t_LDADD = ../libmdpp.la

aladip_unperturbed_t_SOURCES = aladip_unperturbed.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc 
aladip_unperturbed_t_LDADD = ../libmdpp.la

aladip_special_t_SOURCES = aladip_special.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc
aladip_special_t_LDADD = ../libmdpp.la

aladip_atomic_t_SOURCES = aladip_atomic.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc
aladip_atomic_t_LDADD = ../libmdpp.la

aladip_ls_t_SOURCES = aladip_ls.t.cc\
  check.h\
  check_forcefield.h\
  check_forcefield.cc\
  check_state.h\
  check_state.cc
aladip_ls_t_LDADD = ../libmdpp.la

c16_cg_t_SOURCES = c16_cg.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc

c16_cg_t_LDADD = ../libmdpp.la

lambdas_t_SOURCES = lambdas.t.cc\
	check.h\
	check_forcefield.h\
	check_forcefield.cc\
	check_state.h\
	check_state.cc

lambdas_t_LDADD = ../libmdpp.la
TESTS = aladip.t aladip_unperturbed.t aladip_atomic.t aladip_ls.t aladip_special.t c16_cg.t lambdas.t

EXTRA_DIST = data/aladip.topo\
	data/aladip.pttopo\
	data/aladip.conf\
	data/aladip.in\
	data/aladip_unperturbed.in\
	data/aladip_special.in\
	data/aladip.distrest\
	data/aladip.dihrest\
	data/aladip.order\
	data/aladip.xrs\
  data/aladip.lud\
	data/aladip.led\
	data/aladip_atomic.in\
	data/aladip_ewald.in\
	data/aladip_ls.in\
	data/lambdas.off.in\
	data/lambdas.on.in\
	data/lambdas.pttopo\
	data/cg16.conf\
	data/cg16.in\
	data/cg16.topo\
	data/scaling.in\
	data/scaling.lambdadep.in\
	data/scaling.lambdadep.pttopo\
	data/scaling.pttopo
