INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/interaction -I$(top_srcdir)/src/interaction/nonbonded -I$(top_srcdir)/src/interaction/nonbonded/pairlist -I$(top_srcdir)/src/interaction/nonbonded/interaction -I$(top_srcdir)/src/interaction/bonded -I$(top_srcdir)/src/interaction/forcefield -I$(top_srcdir)/src/interaction/special

xxincludedir = 			$(includedir)/gromosXX/algorithm

xxinclude_HEADERS =		algorithm.h \
				create_md_sequence.h \
			  	algorithm/algorithm_sequence.h \
				integration/energy_calculation.h \
				integration/leap_frog.h \
				integration/stochastic.h \
				integration/slow_growth.h \
				integration/steepest_descent.h \
				integration/analyze.h \
			 	temperature/temperature_calculation.h \
				temperature/thermostat.h \
				temperature/berendsen_thermostat.h \
				temperature/nosehoover_thermostat.h \
			  	pressure/pressure_calculation.h \
			  	pressure/berendsen_barostat.h \
			  	constraints/create_constraints.h \
				constraints/position_constraints.h \
				constraints/shake.h \
				constraints/perturbed_shake.h \
				constraints/lincs.h \
				constraints/flexible_constraint.h \
				constraints/perturbed_flexible_constraint.h \
				constraints/remove_com_motion.h \
				constraints/rottrans.h

noinst_LTLIBRARIES	=	libgroalgorithm.la

libgroalgorithm_la_SOURCES =	algorithm.cc \
				create_md_sequence.cc \
				algorithm/algorithm_sequence.cc \
				integration/leap_frog.cc \
				integration/stochastic.cc \
				temperature/temperature_calculation.cc \
				temperature/thermostat.cc \
				temperature/berendsen_thermostat.cc \
				temperature/nosehoover_thermostat.cc \
				pressure/pressure_calculation.cc \
				pressure/berendsen_barostat.cc \
				constraints/create_constraints.cc \
				constraints/position_constraints.cc \
				constraints/remove_com_motion.cc \
				constraints/shake.cc \
				constraints/perturbed_shake.cc \
				constraints/flexible_constraint.cc \
				constraints/perturbed_flexible_constraint.cc \
				constraints/lincs.cc \
				integration/energy_calculation.cc \
				integration/slow_growth.cc \
				integration/steepest_descent.cc \
				integration/analyze.cc \
				constraints/rottrans.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

