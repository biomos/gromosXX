INCLUDES = -I$(top_srcdir)/src

include_HEADERS = 		algorithm.h \
				algorithm/md_global.h \
				algorithm/md_spec.h \
				algorithm/interaction_spec.h \
				algorithm/md_base.h \
				algorithm/md_base.tcc \
				algorithm/md.h \
				algorithm/md.tcc \
				algorithm/perturbation_md.h \
				algorithm/perturbation_md.tcc \
				integration/leap_frog.h \
				integration/leap_frog.tcc \
				integration/runge_kutta.h \
				integration/runge_kutta.tcc \
				constraint/shake.h \
				constraint/shake.tcc \
				constraint/perturbed_shake.h \
				constraint/perturbed_shake.tcc \
				constraint/flexible_constraint.h \
				constraint/flexible_constraint.tcc \
				constraint/perturbed_flexible_constraint.h \
				constraint/perturbed_flexible_constraint.tcc \
				temperature/berendsen.h \
				temperature/berendsen.tcc \
				pressure/berendsen.h \
				pressure/berendsen.tcc

lib_LTLIBRARIES	= 		libalgorithm.la

libalgorithm_la_LDFLAGS = 	-version-info $(MD_MAJOR_VERSION):$(MD_MINOR_VERSION):0

libalgorithm_la_SOURCES =	algorithm.cc algorithm/md_global.cc \
	algorithm/grid_scaled_molvir_perturbed_md.cc 	algorithm/grid_scaled_atvir_perturbed_md.cc \
	algorithm/grid_scaled_molvir_md.cc 	   	algorithm/grid_scaled_atvir_md.cc \
	algorithm/grid_scaled_perturbed_md.cc 	   	algorithm/grid_scaled_md.cc \
	algorithm/grid_molvir_perturbed_md.cc 	   	algorithm/grid_atvir_perturbed_md.cc \
	algorithm/grid_molvir_md.cc 		   	algorithm/grid_atvir_md.cc \
	algorithm/grid_perturbed_md.cc 		   	algorithm/grid_md.cc \
	algorithm/scaled_molvir_perturbed_md.cc 	algorithm/scaled_atvir_perturbed_md.cc \
	algorithm/scaled_molvir_md.cc 		  	algorithm/scaled_atvir_md.cc \
	algorithm/scaled_perturbed_md.cc 		algorithm/scaled_md.cc \
	algorithm/molvir_perturbed_md.cc 		algorithm/atvir_perturbed_md.cc \
	algorithm/molvir_md.cc 			   	algorithm/atvir_md.cc \
	algorithm/perturbed_md.cc 		   	algorithm/md.cc 


bin_PROGRAMS =

check_PROGRAMS =

TESTS =
