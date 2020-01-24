/**
 * @file bfactor.cc
 * bfactorbusiness
 */
#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

// special interactions
#include "../../../interaction/interaction_types.h"
#include "../../../util/template_split.h"
#include "../../../util/debug.h"

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "../../../interaction/special/xray/dens.h"
#include "../../../interaction/special/xray/sf.h"
#include "../../../interaction/special/xray/bfactor.h"
#endif

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#ifdef HAVE_CLIPPER
// this struct is used to pass the C++ objects to the C funciton as a void*
struct fit_param {
  topology::Topology * topo;
        configuration::Configuration * conf;
        simulation::Simulation * sim;
        const clipper::Cell * cell;
        clipper::Atom_list * atoms;
        clipper::HKL_data<clipper::data32::F_phi> * fphi_calc;
        clipper::HKL_data<clipper::data32::F_phi> * fphi;
        clipper::HKL_data<clipper::data32::F_phi> * fphi_obs;
        clipper::Xmap<clipper::ftype32> * rho_calc;
        clipper::FFTmap_p1 * D_k;
        clipper::Xmap<clipper::ftype32> * d_r;
};

double bfactor_residual(const gsl_vector * B, void * param);
void bfactor_deriv(const gsl_vector * B, void * param,
        gsl_vector * gradient);
void bfactor_residual_and_deriv(const gsl_vector * B, void * param,
        double * energy, gsl_vector * gradient);

void interaction::xray::fit_bfactor(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        const clipper::Cell & cell,
        clipper::Atom_list & atoms,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_calc,
        clipper::HKL_data<clipper::data32::F_phi> & fphi,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_obs,
        clipper::Xmap<clipper::ftype32> & rho_calc,
        clipper::FFTmap_p1 & D_k,
        clipper::Xmap<clipper::ftype32> & d_r) {
  fit_param fp;
  fp.topo = &topo;
  fp.conf = &conf;
  fp.sim = &sim;
  fp.cell = &cell;
  fp.atoms = &atoms;
  fp.fphi_calc = &fphi_calc;
  fp.fphi = &fphi;
  fp.fphi_obs = &fphi_obs;
  fp.rho_calc = &rho_calc;
  fp.D_k = &D_k;
  fp.d_r = &d_r;

  const double to_ang = sim.param().xrayrest.to_angstrom;
  const double sqpi2 = math::Pi * math::Pi * 8.0;

  unsigned int num_group = sim.param().xrayrest.bfactor.groups.size();
  gsl_multimin_function_fdf fit_func;
  fit_func.n = num_group;
  fit_func.f = &bfactor_residual;
  fit_func.df = &bfactor_deriv;
  fit_func.fdf = &bfactor_residual_and_deriv;
  fit_func.params = (void *) &fp;

  gsl_vector * B = gsl_vector_alloc(num_group);
  for(unsigned int i = 0; i < num_group; ++i) {
    gsl_vector_set(B, i, atoms[*(sim.param().xrayrest.bfactor.groups[i].begin())].u_iso() * sqpi2 / (to_ang * to_ang));
  }

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
       gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, num_group);

  gsl_multimin_fdfminimizer_set(s, &fit_func, B, sim.param().xrayrest.bfactor.min, sim.param().xrayrest.bfactor.terminate_gradient);
  int iter = 0;
  int status;
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;

    status = gsl_multimin_test_gradient(s->gradient, sim.param().xrayrest.bfactor.terminate_gradient);
  } while (status == GSL_CONTINUE && iter < sim.param().xrayrest.bfactor.terminate_iterations);

  std::cout << "R after " << iter << " iterations: " << conf.special().xray.R_inst << std::endl;

  for(unsigned int i = 0; i < num_group; ++i) {
    double B_i = gsl_vector_get(s->x, i);
    if (B_i < sim.param().xrayrest.bfactor.min) B_i = sim.param().xrayrest.bfactor.min;
    if (B_i > sim.param().xrayrest.bfactor.max) B_i = sim.param().xrayrest.bfactor.max;

    std::set<unsigned int>::const_iterator it = sim.param().xrayrest.bfactor.groups[i].begin(),
            to = sim.param().xrayrest.bfactor.groups[i].end();
    for (; it != to; ++it) {
      atoms[*it].set_u_iso(B_i * to_ang * to_ang / sqpi2);
      conf.special().xray_bfoc[*it].b_factor = B_i;
    }
  }
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(B);
}

double overall_bfactor_residual(const gsl_vector * B, void * param);
void overall_bfactor_deriv(const gsl_vector * B, void * param, gsl_vector * gradient);
void overall_bfactor_residual_and_deriv(const gsl_vector * B, void * param,
  double * energy, gsl_vector * gradient);

void interaction::xray::fit_overall_bfactor(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        const clipper::Cell & cell,
        clipper::Atom_list & atoms,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_calc,
        clipper::HKL_data<clipper::data32::F_phi> & fphi,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_obs,
        clipper::Xmap<clipper::ftype32> & rho_calc,
        clipper::FFTmap_p1 & D_k,
        clipper::Xmap<clipper::ftype32> & d_r) {
  fit_param fp;
  fp.topo = &topo;
  fp.conf = &conf;
  fp.sim = &sim;
  fp.cell = &cell;
  fp.atoms = &atoms;
  fp.fphi_calc = &fphi_calc;
  fp.fphi = &fphi;
  fp.fphi_obs = &fphi_obs;
  fp.rho_calc = &rho_calc;
  fp.D_k = &D_k;
  fp.d_r = &d_r;

  gsl_multimin_function_fdf fit_overall_func;
  fit_overall_func.n = 1;
  fit_overall_func.f = &overall_bfactor_residual;
  fit_overall_func.df = &overall_bfactor_deriv;
  fit_overall_func.fdf = &overall_bfactor_residual_and_deriv;
  fit_overall_func.params = (void *) &fp;

  double & overall_B = conf.special().xray.B_overall;
  gsl_vector * B = gsl_vector_alloc(1);
  gsl_vector_set(B, 0, overall_B);

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, 1);

  gsl_multimin_fdfminimizer_set(s, &fit_overall_func, B, 1E-4, 1E-4);
  int iter = 0;
  int status;
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;

    overall_B = gsl_vector_get(s->x, 0);

    status = gsl_multimin_test_gradient(s->gradient, sim.param().xrayrest.overall_bfactor.terminate_gradient);
  } while (status == GSL_CONTINUE && iter < sim.param().xrayrest.overall_bfactor.terminate_iterations);

  overall_B = gsl_vector_get(s->x, 0);
  std::cout << "R after " << iter << " iterations: " << conf.special().xray.R_inst << std::endl;
  std::cout << "overall B factor: " << overall_B << std::endl;
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(B);
}

using namespace interaction::xray;

double bfactor_residual(const gsl_vector * B, void * param) {
  fit_param & fp = *static_cast<fit_param*>(param);
  topology::Topology & topo = *fp.topo;
  configuration::Configuration & conf = *fp.conf;
  simulation::Simulation & sim = *fp.sim;
  clipper::Atom_list & atoms = *fp.atoms;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_calc = *fp.fphi_calc;
  clipper::HKL_data<clipper::data32::F_phi> & fphi = *fp.fphi;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_obs = *fp.fphi_obs;
  clipper::FFTmap_p1 & D_k = *fp.D_k;
  clipper::Xmap<clipper::ftype32> & rho_calc = *fp.rho_calc;

  const double to_ang = sim.param().xrayrest.to_angstrom;
  const double sqpi2 = math::Pi * math::Pi * 8.0;

  // adjust the B factors of the atoms
  for(unsigned int i = 0; i < sim.param().xrayrest.bfactor.groups.size(); ++i) {
    std::set<unsigned int>::const_iterator it = sim.param().xrayrest.bfactor.groups[i].begin(),
            to = sim.param().xrayrest.bfactor.groups[i].end();
    for(; it != to; ++it)
      atoms[*it].set_u_iso(gsl_vector_get(B, i) * to_ang * to_ang / sqpi2);
  }
  calculate_electron_density(rho_calc, atoms);
  // calculate structure factors and scale them
  rho_calc.fft_to(fphi_calc);
  scale_sf(topo, conf, sim, fphi_calc, fphi, fphi_obs);

  // get the energy
  double energy;
  calculate_energy_sf(sim, fphi, topo.xray_restraints(),
            conf.special().xray_rest,
            sim.param().xrayrest.xrayrest,
            conf.special().xray.k_inst, conf.special().xray.k_avg,
            D_k, 1.0, energy);
  return energy;
}

double overall_bfactor_residual(const gsl_vector * B, void * param) {
  fit_param & fp = *static_cast<fit_param*>(param);
  topology::Topology & topo = *fp.topo;
  configuration::Configuration & conf = *fp.conf;
  simulation::Simulation & sim = *fp.sim;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_calc = *fp.fphi_calc;
  clipper::HKL_data<clipper::data32::F_phi> & fphi = *fp.fphi;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_obs = *fp.fphi_obs;
  clipper::FFTmap_p1 & D_k = *fp.D_k;

  conf.special().xray.B_overall = gsl_vector_get(B, 0);

  scale_sf(topo, conf, sim, fphi_calc, fphi, fphi_obs);

  // get the energy
  double energy;
  calculate_energy_sf(sim, fphi, topo.xray_restraints(),
            conf.special().xray_rest,
            sim.param().xrayrest.xrayrest,
            conf.special().xray.k_inst, conf.special().xray.k_avg,
            D_k, 1.0, energy);
  return energy;
}

void bfactor_residual_and_deriv(const gsl_vector * B, void * param,
  double * energy, gsl_vector * gradient) {
  fit_param & fp = *static_cast<fit_param*>(param);
  topology::Topology & topo = *fp.topo;
  simulation::Simulation & sim = *fp.sim;
  clipper::Atom_list & atoms = *fp.atoms;
  clipper::FFTmap_p1 & D_k = *fp.D_k;
  clipper::Xmap<clipper::ftype32> & d_r = *fp.d_r;

  const double to_ang = sim.param().xrayrest.to_angstrom;

  // get the energy by just calling the residual
  *energy = bfactor_residual(B, param);
  // get the gradients of the bfactors
  math::VArray fdummy(topo.num_atoms());
  math::SArray b_deriv(topo.num_atoms());
  calculate_force_sf(true, D_k, d_r, atoms, fdummy, b_deriv, to_ang);
  for(unsigned int i = 0; i < sim.param().xrayrest.bfactor.groups.size(); ++i) {
    double grad_sum = 0.0;
    std::set<unsigned int>::const_iterator it = sim.param().xrayrest.bfactor.groups[i].begin(),
            to = sim.param().xrayrest.bfactor.groups[i].end();
    for(; it != to; ++it) {
      grad_sum += b_deriv(*it);
    }
    gsl_vector_set(gradient, i, grad_sum);
  }
}

void overall_bfactor_residual_and_deriv(const gsl_vector * B, void * param, double * energy, gsl_vector * gradient) {
  fit_param & fp = *static_cast<fit_param*>(param);
  topology::Topology & topo = *fp.topo;
  simulation::Simulation & sim = *fp.sim;
  configuration::Configuration & conf = *fp.conf;
  clipper::HKL_data<clipper::data32::F_phi> & fphi = *fp.fphi;
  std::vector<topology::xray_restraint_struct> & refl = topo.xray_restraints();
  std::vector<configuration::Configuration::special_struct::xray_struct> & refl_curr = conf.special().xray_rest;

  double & overall_B = conf.special().xray.B_overall;

  // get the energy by just calling the residual
  *energy = overall_bfactor_residual(B, param);
  // get the gradients of the bfactors

  // calculate normalisation factor to get rid of resolution dependence.
  double sum_xray_normalisation_factor = 0.0;
  const double invresolsq = 1.0 / (sim.param().xrayrest.resolution *
        sim.param().xrayrest.resolution *
        sim.param().xrayrest.to_angstrom *
        sim.param().xrayrest.to_angstrom);
  for (unsigned int i = 0; i < refl.size(); i++) {
    const topology::xray_restraint_struct & xrs = refl[i];
    const clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;

    double inv_var = 1.0;
    if (xrs.stddev_sf > math::epsilon) inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);
    sum_xray_normalisation_factor += inv_var * xrs.sf * xrs.sf;
  }
  double xray_normalisation_factor = 1.0;
  if (sum_xray_normalisation_factor > math::epsilon)
    xray_normalisation_factor = 1.0 / sum_xray_normalisation_factor;

  const double & k_inst = conf.special().xray.k_inst;
  double dterm = 0;
  const double force_constant = sim.param().xrayrest.force_constant;

  for (unsigned int i = 0; i < refl.size(); i++) {
    const topology::xray_restraint_struct & xrs = refl[i];
    const clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;
    const double fobs = xrs.sf;
    double inv_var = 1.0;
    if (xrs.stddev_sf > math::epsilon) inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);
    const double fcalc = refl_curr[i].sf_curr;
    const double term = fobs - k_inst * fcalc;
    DEBUG(8, "\tterm: " << term << " inv_var: " << inv_var);
    // calculate derivatives of target function
    dterm += force_constant * xray_normalisation_factor * inv_var *
            term * k_inst * fcalc * (-overall_B) * invresolsq;
  }

  gsl_vector_set(gradient, 0, dterm);
}

void bfactor_deriv(const gsl_vector * B, void * param,
  gsl_vector * gradient) {
  double dummy;
  bfactor_residual_and_deriv(B, param, &dummy, gradient);
}
void overall_bfactor_deriv(const gsl_vector * B, void * param,
  gsl_vector * gradient) {
  double dummy;
  overall_bfactor_residual_and_deriv(B, param, &dummy, gradient);
}
#endif
