INCLUDES = -I$(top_srcdir)/src

xxincludedir = $(includedir)/gromosXX

xxinclude_HEADERS =	algorithm.h \
			algorithm/algorithm_sequence.h \
			integration/leap_frog.h \
			temperature/temperature_calculation.h \
			temperature/berendsen_thermostat.h \
			pressure/pressure_calculation.h \
			pressure/berendsen_barostat.h \
			create_md_sequence.h \
			constraints/create_constraints.h \
			constraints/shake.h \
			constraints/shake.tcc \
			constraints/perturbed_shake.h \
			constraints/perturbed_shake.tcc \
			constraints/lincs.h \
			constraints/lincs.tcc \
			constraints/flexible_constraint.h \
			constraints/flexible_constraint.tcc \
			constraints/perturbed_flexible_constraint.h \
			constraints/perturbed_flexible_constraint.tcc \
			constraints/remove_com_motion.h \
			integration/slow_growth.h \
			integration/steepest_descent.h

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

