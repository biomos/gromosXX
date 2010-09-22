/**
 * @file bfactor.cc
 * bfactorbusiness
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

// special interactions
#include <interaction/interaction_types.h>
#include <util/template_split.h>
#include <util/debug.h>

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <interaction/special/xray/dens.h>
#include <interaction/special/xray/sf.h>
#include <interaction/special/xray/bfactor.h>
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
  fp.fphi = &fphi;
  fp.fphi_obs = &fphi_obs;
  fp.rho_calc = &rho_calc;
  fp.D_k = &D_k;
  fp.d_r = &d_r;

  const double to_ang = sim.param().xrayrest.to_angstrom;
  const double sqpi2 = math::Pi * math::Pi * 8.0;

  gsl_multimin_function_fdf fit_func;
  fit_func.n = atoms.size();
  fit_func.f = &bfactor_residual;
  fit_func.df = &bfactor_deriv;
  fit_func.fdf = &bfactor_residual_and_deriv;
  fit_func.params = (void *) &fp;

  gsl_vector * B = gsl_vector_alloc(atoms.size());
  for(unsigned int i = 0; i < atoms.size(); ++i) {
    gsl_vector_set(B, i, atoms[i].u_iso() * sqpi2 / (to_ang * to_ang));
  }

  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
       gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, atoms.size());

  gsl_multimin_fdfminimizer_set(s, &fit_func, B, 0.01, 1e-4);
  unsigned int iter = 0;
  int status;
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;

    status = gsl_multimin_test_gradient(s->gradient, sim.param().xrayrest.bfactor.terminate_gradient);
  } while (status == GSL_CONTINUE && iter < sim.param().xrayrest.bfactor.terminate_iterations);

  std::cout << "R after " << iter << " iterations: " << conf.special().xray.R_inst << std::endl;

  for(unsigned int i = 0; i < atoms.size(); ++i) {
    double B_i = gsl_vector_get(s->x, i);
    if (B_i < sim.param().xrayrest.bfactor.min) B_i = sim.param().xrayrest.bfactor.min;
    if (B_i > sim.param().xrayrest.bfactor.max) B_i = sim.param().xrayrest.bfactor.max;
    atoms[i].set_u_iso(B_i * to_ang * to_ang / sqpi2);
    conf.special().xray_bfoc[i].b_factor = B_i;
  }
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
  clipper::HKL_data<clipper::data32::F_phi> & fphi = *fp.fphi;
  clipper::HKL_data<clipper::data32::F_phi> & fphi_obs = *fp.fphi_obs;
  clipper::FFTmap_p1 & D_k = *fp.D_k;
  clipper::Xmap<clipper::ftype32> & rho_calc = *fp.rho_calc;

  const double to_ang = sim.param().xrayrest.to_angstrom;
  const double sqpi2 = math::Pi * math::Pi * 8.0;

  // adjust the B factors of the atoms
  for(unsigned int i = 0; i < atoms.size(); ++i) {
    atoms[i].set_u_iso(gsl_vector_get(B, i) * to_ang * to_ang / sqpi2);
  }

  calculate_electron_density(rho_calc, atoms);
  // calculate structure factors and scale them
  rho_calc.fft_to(fphi);
  scale_sf(topo, conf, sim, fphi, fphi_obs);

  // get the energy
  double energy;
  calculate_energy_sf(topo.xray_restraints(),
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
  calculate_force_sf(D_k, d_r, atoms, fdummy, b_deriv, to_ang);
  for(unsigned int i = 0; i < b_deriv.size(); ++i) {
    gsl_vector_set(gradient, i, b_deriv(i));
  }
}

void bfactor_deriv(const gsl_vector * B, void * param,
  gsl_vector * gradient) {
  double dummy;
  bfactor_residual_and_deriv(B, param, &dummy, gradient);
}
#endif
