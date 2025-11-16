# This is the programs automake template.
AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src

EXTRA_DIST =   0index_contrib.doxy

bin_PROGRAMS = \
	aggregates\
	rottemp\
	ener_rplc\
	espmap\
	lb_top\
	m_ener\
	prep_bb\
	cg_potential\
	countbb\
	dipole_pot\
	stacking\
	shake\
	model\
	gca_ener\
	renumber\
	distmat\
	iwdcf\
	bb2tex\
	crossdih_top\
	saltbridge\
	bilayer_gel\
	make_rdc_spec\
	fit_rdc\
	add_atom\
	top_to_build\
	helix_axis\
	gromacs2gromos\
	diff_top\
	grid_dens\
	pdb2seq\
	nscatt\
	structal\
	cgLJpot\
	hvap\
	depsi2native\
	rep_rewrite\
	local_comp\
	prep_hybrid\
	sasa_new\
	ener_test\
	force\
	mut_top\
    mod_top\
    rmsd_cmp\
    eds_mult_all\
    stat_ana\
    sphericity\
	expand_ua\
    fit_ener\
    fit_ener_traj

noinst_LTLIBRARIES = libvmdplugin.la

aggregates_SOURCES = aggregates.cc
rottemp_SOURCES = rottemp.cc
ener_rplc_SOURCES = ener_rplc.cc
espmap_SOURCES = espmap.cc
lb_top_SOURCES = lb_top.cc
m_ener_SOURCES = m_ener.cc
prep_bb_SOURCES = prep_bb.cc
cg_potential_SOURCES = cg_potential.cc
countbb_SOURCES = countbb.cc
dipole_pot_SOURCES = dipole_pot.cc
stacking_SOURCES = stacking.cc
shake_SOURCES = shake.cc
model_SOURCES = model.cc
gca_ener_SOURCES = gca_ener.cc
renumber_SOURCES = renumber.cc
distmat_SOURCES = distmat.cc
iwdcf_SOURCES = iwdcf.cc
bb2tex_SOURCES = bb2tex.cc
crossdih_top_SOURCES = crossdih_top.cc
saltbridge_SOURCES = saltbridge.cc
bilayer_gel_SOURCES = bilayer_gel.cc
make_rdc_spec_SOURCES = make_rdc_spec.cc
fit_rdc_SOURCES = fit_rdc.cc
add_atom_SOURCES = add_atom.cc
top_to_build_SOURCES = top_to_build.cc
helix_axis_SOURCES = helix_axis.cc
gromacs2gromos_SOURCES = gromacs2gromos.cc
diff_top_SOURCES = diff_top.cc
grid_dens_SOURCES = grid_dens.cc
pdb2seq_SOURCES = pdb2seq.cc
nscatt_SOURCES = nscatt.cc
structal_SOURCES = structal.cc
cgLJpot_SOURCES = cgLJpot.cc
hvap_SOURCES = hvap.cc
depsi2native_SOURCES = depsi2native.cc
rep_rewrite_SOURCES = rep_rewrite.cc
local_comp_SOURCES = local_comp.cc
prep_hybrid_SOURCES = prep_hybrid.cc
sasa_new_SOURCES = sasa_new.cc
ener_test_SOURCES = ener_test.cc
force_SOURCES = force.cc
mut_top_SOURCES = mut_top.cc
mod_top_SOURCES = mod_top.cc
rmsd_cmp_SOURCES = rmsd_cmp.cc
eds_mult_all_SOURCES = eds_mult_all.cc
stat_ana_SOURCES = stat_ana.cc
sphericity_SOURCES = sphericity.cc
expand_ua_SOURCES = expand_ua.cc
fit_ener_SOURCES = fit_ener.cc
fit_ener_traj_SOURCES = fit_ener_traj.cc

libvmdplugin_la_SOURCES = vmdplugin.cc

LDADD = $(top_builddir)/src/libgromos.la

shake_LDADD = $(MDPP_LIB) $(top_builddir)/src/libgromos.la
shake_LDFLAGS = $(MDPP_LDFLAGS)

gromacs2gromos_LDADD = $(GMX_LIB) $(top_builddir)/src/libgromos.la
gromacs2gromos_LDFLAGS = $(GMX_LDFLAGS) $(top_builddir)/src/libgromos.la

gincludedir = $(includedir)/gromos++
