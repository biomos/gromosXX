INCLUDES = -I$(top_srcdir)/src

check_PROGRAMS = 	rasn.t \
			aladip.t

rasn_t_SOURCES = rasn.t.cc check.h check_forcefield.h check_forcefield.cc

rasn_t_LDADD = ../libgroxx.la

aladip_t_SOURCES = aladip.t.cc check.h check_forcefield.h check_forcefield.cc \
			check_state.h check_state.cc

aladip_t_LDADD = ../libgroxx.la

TESTS = aladip.t

EXTRA_DIST =	data/rasn.topo \
		data/rasn.conf \
		data/rasn.pttopo \
		data/rasn.in \
		data/aladip.topo \
		data/aladip.pttopo \
		data/aladip.conf \
		data/aladip.in
