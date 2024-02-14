AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/src/interaction -I$(top_srcdir)/src/interaction/nonbonded -I$(top_srcdir)/src/interaction/nonbonded/pairlist -I$(top_srcdir)/src/interaction/nonbonded/interaction -I$(top_srcdir)/src/interaction/bonded -I$(top_srcdir)/src/interaction/forcefield -I$(top_srcdir)/src/interaction/special

if WITH_CUDA
  AM_CPPFLAGS += -I$(top_srcdir)/src/cukernel
endif

xxincludedir = $(includedir)/md++/algorithm
xxinclude_HEADERS =	algorithm.h\
	create_md_sequence.h

xxalgalgdir = $(includedir)/md++/algorithm/algorithm
xxalgalg_HEADERS = algorithm/algorithm_sequence.h

xxalgintegrationdir = $(includedir)/md++/algorithm/integration
xxalgintegration_HEADERS = integration/energy_calculation.h\
	integration/leap_frog.h\
	integration/scaled_leap_frog.h\
	integration/monte_carlo.h\
	integration/stochastic.h\
	integration/slow_growth.h\
	integration/steepest_descent.h\
	integration/conjugate_gradient.h\
	integration/analyze.h\
	integration/lattice_shift.h\
	integration/multigradient.h\
    integration/eds.h\
	integration/gamd.h

xxalgtempdir = $(includedir)/md++/algorithm/temperature
xxalgtemp_HEADERS = temperature/temperature_calculation.h\
	temperature/thermostat.h\
	temperature/berendsen_thermostat.h\
	temperature/nosehoover_thermostat.h

xxalgpressudir = $(includedir)/md++/algorithm/pressure
xxalgpressu_HEADERS = pressure/pressure_calculation.h\
	pressure/berendsen_barostat.h

xxalgconstdir = $(includedir)/md++/algorithm/constraints
xxalgconst_HEADERS = constraints/create_constraints.h\
	constraints/position_constraints.h\
	constraints/shake.h\
	constraints/perturbed_shake.h\
	constraints/angle_constraint.cc\
        constraints/perturbed_angle_constraint.cc\
	constraints/dihedral_constraint.cc\
	constraints/perturbed_dihedral_constraint.cc\
	constraints/lincs.h\
	constraints/flexible_constraint.h\
	constraints/perturbed_flexible_constraint.h\
	constraints/remove_com_motion.h\
	constraints/rottrans.h\
	constraints/settle.h\
        constraints/m_shake.h\
        constraints/gpu_shake.h\
        constraints/gpu_shake_thread.h

xxalgvirtualatdir = $(includedir)/md++/algorithm/virtualatoms
xxalgvirtualat_HEADERS = virtualatoms/prepare_virtualatoms.h\
	virtualatoms/propagate_forces.h

noinst_LTLIBRARIES=	libgroalgorithm.la

libgroalgorithm_la_SOURCES = algorithm.cc\
	create_md_sequence.cc\
	algorithm/algorithm_sequence.cc\
	integration/leap_frog.cc\
	integration/scaled_leap_frog.cc\
	integration/monte_carlo.cc\
	integration/stochastic.cc\
	temperature/temperature_calculation.cc\
	temperature/thermostat.cc\
	temperature/berendsen_thermostat.cc\
	temperature/nosehoover_thermostat.cc\
	pressure/pressure_calculation.cc\
	pressure/berendsen_barostat.cc\
	constraints/create_constraints.cc\
	constraints/position_constraints.cc\
	constraints/remove_com_motion.cc\
	constraints/shake.cc\
	constraints/perturbed_shake.cc\
	constraints/flexible_constraint.cc\
	constraints/perturbed_flexible_constraint.cc\
	constraints/lincs.cc\
	constraints/settle.cc\
	constraints/m_shake.cc\
	constraints/gpu_shake.cc\
	constraints/gpu_shake_thread.cc\
	integration/energy_calculation.cc\
	integration/slow_growth.cc\
	integration/steepest_descent.cc\
	integration/conjugate_gradient.cc\
	integration/analyze.cc\
	constraints/rottrans.cc\
	integration/lattice_shift.cc\
	integration/multigradient.cc\
    integration/eds.cc\
	integration/gamd.cc\
	virtualatoms/prepare_virtualatoms.cc\
	virtualatoms/propagate_forces.cc


bin_PROGRAMS =			

check_PROGRAMS =

TESTS =

