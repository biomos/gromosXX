INCLUDES = -I$(top_srcdir)/src

xxincludedir = 			$(includedir)/gromosXX/algorithm

xxinclude_HEADERS =		algorithm.h \
				create_md_sequence.h \
			  	algorithm/algorithm_sequence.h \
				integration/leap_frog.h \
				integration/slow_growth.h \
				integration/steepest_descent.h \
			 	temperature/temperature_calculation.h \
				temperature/berendsen_thermostat.h \
			  	pressure/pressure_calculation.h \
			  	pressure/berendsen_barostat.h \
			  	constraints/create_constraints.h \
				constraints/shake.h \
				constraints/shake.cc \
				constraints/perturbed_shake.h \
				constraints/perturbed_shake.cc \
				constraints/lincs.h \
				constraints/lincs.cc \
				constraints/flexible_constraint.h \
				constraints/flexible_constraint.cc \
				constraints/perturbed_flexible_constraint.h \
				constraints/perturbed_flexible_constraint.cc \
				constraints/remove_com_motion.h \
				constraints/rottrans.h \
				constraints/rottrans.cc

noinst_LTLIBRARIES	=	libgroalgorithm.la

libgroalgorithm_la_SOURCES =	algorithm.cc \
				create_md_sequence.cc \
				algorithm/algorithm_sequence.cc \
				integration/leap_frog.cc \
				temperature/temperature_calculation.cc \
				temperature/berendsen_thermostat.cc \
				pressure/pressure_calculation.cc \
				pressure/berendsen_barostat.cc \
				constraints/create_constraints.cc \
				constraints/remove_com_motion.cc \
				integration/slow_growth.cc \
				integration/steepest_descent.cc

bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

