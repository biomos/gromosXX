/**
 * @file sf.cc
 * structure factor business
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
#include "../../../interaction/special/xray/sf.h"

#include "dens.h"
#endif

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

#ifdef HAVE_CLIPPER
void interaction::xray::scale_sf(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const clipper::HKL_data<clipper::data32::F_phi> & fphi_calc,
        clipper::HKL_data<clipper::data32::F_phi> & fphi,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_obs) {

  // scale here with overall B
  if (sim.param().xrayrest.overall_bfactor.B_overall_switcher == simulation::B_overall_on){
    /*clipper::HKL_data<clipper::data32::F_phi> * fphi_calc;
    interaction::xray::calculate_electron_density(conf.special().xray_conf.rho_calc, conf.special().xray_conf.atoms);
    // calculate structure factors and scale them
    conf.special().xray_conf.rho_calc.fft_to(fphi);*/
    clipper::HKL_data<clipper::data32::F_phi>::HKL_reference_index ix = fphi.first();
    clipper::HKL_data<clipper::data32::F_phi>::HKL_reference_index ix_calc = fphi_calc.first();
    for(; !ix.last(); ix.next(), ix_calc.next()) {
      std::complex<clipper::ftype32> exp_term(exp(-conf.special().xray.B_overall * ix.invresolsq()), 0.0f);
      fphi[ix] = exp_term * std::complex<clipper::ftype32>(fphi_calc[ix_calc]);
    }
  } else {
    fphi = fphi_calc;
  }
  // sqr_calc:       sum of w*Fcalc^2
  // obs:            sum of Fobs
  // obs_calc:       sum of w*Fobs*Fcalc
  // obs_k_calc:     sum of |Fobs-k*Fcalc|
  // zero all the sums
  double sqr_calc = 0.0, obs = 0.0, obs_free = 0.0, obs_calc = 0.0, obs_k_calc = 0.0,
         sqr_calcavg = 0.0, obs_calcavg = 0.0, obs_k_calcavg = 0.0,
         obs_k_calc_free = 0.0, obs_k_calcavg_free = 0.0, obs_calc_free = 0.0, obs_calcavg_free = 0.0,
         sqr_calc_free = 0.0, sqr_calcavg_free = 0.0;
  // Number of reflections
  const unsigned int num_xray_rest = topo.xray_restraints().size();
  const unsigned int num_xray_rfree = topo.xray_rfree().size();
  // e-term for time-average
  const double eterm = exp(-sim.time_step_size() / sim.param().xrayrest.tau);

  const double invresolsq = 1.0 / (sim.param().xrayrest.resolution *
          sim.param().xrayrest.resolution *
          sim.param().xrayrest.to_angstrom *
          sim.param().xrayrest.to_angstrom);

  // loop over structure factors
  unsigned int j = 0;
  for (unsigned int i = 0; i < num_xray_rest; i++, j++) {
    // filter calculated structure factors: save phases and amplitudes
    const clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
    conf.special().xray_rest[j].sf_curr = fabs(fphi[hkl].f());
    conf.special().xray_rest[j].phase_curr = fphi[hkl].phi();
    DEBUG(15,"HKL:" << hkl.h() << "," << hkl.k() << "," << hkl.l());
    DEBUG(15,"\tSF: " << conf.special().xray_rest[j].sf_curr);

    // reset the averages at the beginning if requested
    if (!sim.param().xrayrest.readavg && sim.steps() == 0) {
      conf.special().xray_rest[j].sf_av = conf.special().xray_rest[j].sf_curr;
      conf.special().xray_rest[j].phase_av = conf.special().xray_rest[j].phase_curr;
    }

    // calculate averages
    conf.special().xray_rest[j].sf_av = fabs((1.0 - eterm) * conf.special().xray_rest[j].sf_curr + eterm * conf.special().xray_rest[j].sf_av);
    conf.special().xray_rest[j].phase_av = (1.0 - eterm) * conf.special().xray_rest[j].phase_curr + eterm * conf.special().xray_rest[j].phase_av;

    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;

    // calculate the inverse of the variance (weight factor)
    double inv_var = 1.0;
    if (topo.xray_restraints()[i].stddev_sf > math::epsilon)
      inv_var = 1.0 / (topo.xray_restraints()[i].stddev_sf * topo.xray_restraints()[i].stddev_sf);

    // calc sums
    obs_calc += inv_var * conf.special().xray_rest[j].sf_curr * topo.xray_restraints()[i].sf;
    obs_calcavg += inv_var * conf.special().xray_rest[j].sf_av * topo.xray_restraints()[i].sf;
    sqr_calc += inv_var * conf.special().xray_rest[j].sf_curr * conf.special().xray_rest[i].sf_curr;
    obs += topo.xray_restraints()[i].sf;
    sqr_calcavg += inv_var * conf.special().xray_rest[j].sf_av * conf.special().xray_rest[j].sf_av;
  }
  // loop over structure factors in R free set
  for (unsigned int i = 0; i < num_xray_rfree; i++, j++) {
    // filter calculated structure factors: save phases and amplitudes for R free HKLs
    const clipper::HKL hkl(topo.xray_rfree()[i].h, topo.xray_rfree()[i].k, topo.xray_rfree()[i].l);
    conf.special().xray_rest[j].sf_curr = fabs(fphi[hkl].f());
    conf.special().xray_rest[j].phase_curr = fphi[hkl].phi();
    DEBUG(15,"HKL:" << hkl.h() << "," << hkl.k() << "," << hkl.l());
    DEBUG(15,"\tSF: " << conf.special().xray_rest[j].sf_curr);

    // reset the averages at the beginning if requested
    if (!sim.param().xrayrest.readavg && sim.steps() == 0) {
      conf.special().xray_rest[j].sf_av = conf.special().xray_rest[j].sf_curr;
      conf.special().xray_rest[j].phase_av = conf.special().xray_rest[j].phase_curr;
    }

    // calculate averages
    conf.special().xray_rest[j].sf_av = fabs((1.0 - eterm) * conf.special().xray_rest[j].sf_curr + eterm * conf.special().xray_rest[j].sf_av);
    conf.special().xray_rest[j].phase_av = (1.0 - eterm) * conf.special().xray_rest[j].phase_curr + eterm * conf.special().xray_rest[j].phase_av;

    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;
    // calculate the inverse of the variance (weight factor)
    double inv_var = 1.0;
    if (topo.xray_rfree()[i].stddev_sf > math::epsilon)
      inv_var = 1.0 / (topo.xray_rfree()[i].stddev_sf * topo.xray_rfree()[i].stddev_sf);

    // calc sums
    obs_free += topo.xray_rfree()[i].sf;
    obs_calc_free += inv_var * conf.special().xray_rest[j].sf_curr * topo.xray_rfree()[i].sf;
    obs_calcavg_free += inv_var * conf.special().xray_rest[j].sf_av * topo.xray_rfree()[i].sf;
    sqr_calc_free += inv_var * conf.special().xray_rest[j].sf_curr * conf.special().xray_rest[j].sf_curr;
    sqr_calcavg_free += inv_var * conf.special().xray_rest[j].sf_av * conf.special().xray_rest[j].sf_av;
  }
  // check for possible resolution problems
#ifdef HAVE_ISNAN
  if (std::isnan(sqr_calc)){
    io::messages.add("Structure factors were NaN. This can be due to numerical problems. "
                     "Try to slighlty increase the resolution.", "X-Ray Restraints", io::message::error);
    return;
  }
#endif

  // calculate the scaling constants for inst and avg.
  double & k_inst = conf.special().xray.k_inst;
  k_inst = obs_calc / sqr_calc;
  double & k_avg = conf.special().xray.k_avg;
  k_avg = obs_calcavg / sqr_calcavg;
  double & k_free_inst = conf.special().xray.k_free_inst;
  if (num_xray_rfree)
    k_free_inst = obs_calc_free / sqr_calc_free;
  else
    k_free_inst = 0.0;
  double & k_free_avg = conf.special().xray.k_free_avg;
  if (num_xray_rfree)
    k_free_avg = obs_calcavg_free / sqr_calcavg_free;
  else
    k_free_avg = 0.0;
  DEBUG(10, "k_inst value: " << k_inst);
  DEBUG(10, "k_avg  value: " << k_avg);
  DEBUG(10, "k_free_inst value: " << k_free_inst);
  DEBUG(10, "k_free_avg  value: " << k_free_avg);

  // calculate sums needed for R factors
  // and "observed" structure factors
  j = 0;
  for (unsigned int i = 0; i < num_xray_rest; i++, j++) {
    const topology::xray_restraint_struct & xrs = topo.xray_restraints()[i];
    const clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;

    double inv_var = 1.0;
    if (xrs.stddev_sf > math::epsilon)
      inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);

    obs_k_calc += fabs(xrs.sf - k_inst * conf.special().xray_rest[j].sf_curr);
    obs_k_calcavg += fabs(xrs.sf - k_avg * conf.special().xray_rest[j].sf_av);

    // save Fobs and PhiCalc for 2Fobs-kFcalc maps. This will be corrected
    // for symmetry in the FFT step.
    if (sim.param().xrayrest.xrayrest == simulation::xrayrest_inst)
      fphi_obs.set_data(hkl, 
              clipper::data32::F_phi(2.0 * xrs.sf - k_inst * conf.special().xray_rest[j].sf_curr,
              conf.special().xray_rest[j].phase_curr));
    else
      fphi_obs.set_data(hkl, 
              clipper::data32::F_phi(2.0 * xrs.sf - k_avg * conf.special().xray_rest[j].sf_av,
              conf.special().xray_rest[j].phase_av));
  }
  // and for R free
  for (unsigned int i = 0; i < num_xray_rfree; i++, j++) {
    const topology::xray_restraint_struct & xrs = topo.xray_rfree()[i];
    const clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;

    
    //calculate the inverse of the variance (weight factor)
    double inv_var = 1.0;
    if (xrs.stddev_sf > math::epsilon)
      inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);

    obs_k_calc_free += fabs(xrs.sf - k_inst * conf.special().xray_rest[j].sf_curr);
    obs_k_calcavg_free += fabs(xrs.sf - k_avg * conf.special().xray_rest[j].sf_av);
  }

  // calculate R factors: R_inst and R_avg
  double & R_inst = conf.special().xray.R_inst;
  R_inst = obs_k_calc / obs;
  double & R_avg = conf.special().xray.R_avg;
  R_avg = obs_k_calcavg / obs;
  double & R_free_inst = conf.special().xray.R_free_inst;
  if (num_xray_rfree)
    R_free_inst = obs_k_calc_free / obs_free;
  else
    R_free_inst = 0.0;
  double & R_free_avg = conf.special().xray.R_free_avg;
  if (num_xray_rfree)
    R_free_avg = obs_k_calcavg_free / obs_free;
  else
    R_free_avg = 0.0;
  DEBUG(10, "R_inst value: " << std::setw(15) << std::setprecision(8) << R_inst);
  DEBUG(10, "R_avg  value: " << std::setw(15) << std::setprecision(8) << R_avg);
  DEBUG(10, "R_free_inst value: " << std::setw(15) << std::setprecision(8) << R_free_inst);
  DEBUG(10, "R_free_avg  value: " << std::setw(15) << std::setprecision(8) << R_free_avg);
}

void interaction::xray::calculate_energy_sf(
        const simulation::Simulation & sim,
        const clipper::HKL_data<clipper::data32::F_phi> & fphi,
        const std::vector<topology::xray_restraint_struct> & refl,
        const std::vector<configuration::Configuration::special_struct::xray_struct> & refl_curr,
        simulation::xrayrest_enum averaging,
        const double k_inst, const double k_avg,
        clipper::FFTmap_p1 & D_k,
        const double force_constant,
        double & energy) {
  const double invresolsq = 1.0 / (sim.param().xrayrest.resolution *
          sim.param().xrayrest.resolution *
          sim.param().xrayrest.to_angstrom *
          sim.param().xrayrest.to_angstrom);
  
  // calculate normalisation factor to get rid of resolution dependence.
  double sum_xray_normalisation_factor = 0.0;
  for (unsigned int i = 0; i < refl.size(); i++) {
    const topology::xray_restraint_struct & xrs = refl[i];
    const clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;

    double inv_var = 1.0;
    if (xrs.stddev_sf > math::epsilon) inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);
    if (averaging == simulation::xrayrest_biq)
      sum_xray_normalisation_factor += inv_var * xrs.sf * xrs.sf * xrs.sf * xrs.sf;
    else
      sum_xray_normalisation_factor += inv_var * xrs.sf * xrs.sf;
  }
  double xray_normalisation_factor = 1.0;
  if (sum_xray_normalisation_factor > math::epsilon) 
    xray_normalisation_factor = 1.0 / sum_xray_normalisation_factor;

  // zero the reciprocal space difference map
  D_k.reset();

  double energy_sum = 0.0;

  // loop over retraints and calculate energy and difference map
  for (unsigned int i = 0; i < refl.size(); i++) {
    const topology::xray_restraint_struct & xrs = refl[i];
    const clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // skip them for sums if they are out of the requested resolution range
    if (invresolsq < fphi.invresolsq(fphi.hkl_info().index_of(hkl)))
      continue;
    // SWITCH FOR DIFFERENT METHODS
    switch (averaging) {
      case simulation::xrayrest_inst :
      {
        // INSTANTANEOUS
        // calculate energy-sum
        const double fobs = xrs.sf;
        double inv_var = 1.0;
        if (xrs.stddev_sf > math::epsilon) inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);
        const double fcalc = refl_curr[i].sf_curr;
        const double term = fobs - k_inst * fcalc;
        DEBUG(8, "\tterm: " << term << " inv_var: " << inv_var);
        energy_sum += xray_normalisation_factor * inv_var * term * term;
        // calculate derivatives of target function
        const double dterm = xray_normalisation_factor * inv_var * (k_inst * fcalc - fobs) * k_inst;
        // Here, I tried to apply symmetry operations for non P1 spacegroups
        // but this had the effect the forces were not in agreement with
        // the finite difference result anymore. So we just safe the relection
        // given in the reflection list and not all symmetric copies. It's
        // up to the user to decide whether he should provide also the symmetric
        // copies for the refinement.
        D_k.set_hkl(hkl, clipper::data32::F_phi(force_constant * dterm, refl_curr[i].phase_curr));
        break;
      }
      case simulation::xrayrest_avg :
      {
        // TIMEAVERAGED
        // calculate energy-sum
        const double fobs = xrs.sf;
        double inv_var = 1.0;
        if (xrs.stddev_sf > math::epsilon) inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);
        const double fcalc = refl_curr[i].sf_av;
        const double term = fobs - k_avg * fcalc;
        energy_sum += xray_normalisation_factor * inv_var * term * term;
        // calculate derivatives of target function
        // here we omit the 1-exp(-dt/tau) term.
        const double dterm = xray_normalisation_factor * inv_var * (k_avg * fcalc - fobs) * k_avg;
        D_k.set_hkl(hkl, clipper::data32::F_phi(force_constant * dterm, refl_curr[i].phase_curr));
        break;
      }
      case simulation::xrayrest_biq :
      {
        // BIQUADRATIC TIME-AVERAGED/INSTANTANEOUS
        // calculate energy-sum
        const double fobs = xrs.sf;
        double inv_var = 1.0;
        if (xrs.stddev_sf > math::epsilon) inv_var = 1.0 / (xrs.stddev_sf * xrs.stddev_sf);
        const double finst = refl_curr[i].sf_curr;
        const double favg = refl_curr[i].sf_av;
        const double inst_term = fobs - k_inst * finst;
        const double av_term = fobs - k_avg * favg;
        energy_sum += xray_normalisation_factor * inv_var * (inst_term * inst_term)*(av_term * av_term);
        // calculate derivatives of target function
        // here we omit the 1-exp(-dt/tau) term.
        const double dterm = xray_normalisation_factor * inv_var * ((k_inst * finst - fobs)*(av_term * av_term) * k_inst
                + (k_avg * favg - fobs)*(inst_term * inst_term) * k_avg);
        D_k.set_hkl(hkl, clipper::data32::F_phi(force_constant * dterm, refl_curr[i].phase_curr));
        break;
      }
      default: break;
    }
  }

  // finally calculate the energy
  energy = 0.5 * force_constant * energy_sum;
  DEBUG(6, "energy: " << energy);
}

void interaction::xray::calculate_force_sf(bool update, clipper::FFTmap_p1 & D_k,
        clipper::Xmap<clipper::ftype32> & d_r,
        const clipper::Atom_list & atoms,
        math::VArray & force,
        math::SArray & b_deriv,
        double to_ang) {
  const double sqpi2 = math::Pi * math::Pi * 8.0;

  if (update) {
    // these are just shortcuts to avoid many calls to the same functions
    const clipper::Spacegroup & spgr = d_r.spacegroup();

    // calculate the inverse symmetry operations of the spacegroup
    std::vector<clipper::Isymop> isymop;
    isymop.resize(spgr.num_symops());
    for (int j = 0; j < spgr.num_symops(); j++) {
      isymop[j] = clipper::Isymop(spgr.symop(j), d_r.grid_sampling());
    }
    const double volume = d_r.cell().volume();
    // convert from Angstrom to nm and add very annyoing scaling constants
    // to make the force volume AND resolution independent.
    const double scale = to_ang / 2.0 * volume / (d_r.grid_sampling().size());
    // perform FFT of the difference map
    D_k.fft_h_to_x(scale);
    // loop over the (symmetry corrected map - even though this doesn't matter).
    for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = d_r.first(); !ix.last(); ix.next()) {
      // set initial data value
      const clipper::Coord_grid & coord = ix.coord();
      d_r[ix] = D_k.real_data(coord);
      // loop over symmetric copies of the grid point and add the data from these points
      for (int j = 1; j < spgr.num_symops(); j++) {
        d_r[ix] += D_k.real_data(coord.transform(isymop[j]).unit(D_k.grid_real()));
      }
      // correct for points mapped on themselves
      d_r[ix] /= d_r.multiplicity(coord);
    }
  }

  // 3.5 is hardcoded atom radius for grid sampling.
  const double radius = 3.5;
  const clipper::Grid_sampling & grid = d_r.grid_sampling();
  // determine the range of the atomic electron density gradient on the gird
  clipper::Grid_range gd(d_r.cell(), grid, radius);

  // loop over the atoms - has to be int and not unsigned due to
  // stupid OpenMP rules
  const int atoms_size = atoms.size();
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < atoms_size; i++) {
    if (!atoms[i].is_null()) {
      math::Vec gradient(0.0, 0.0, 0.0);
      double Uiso_deriv = 0.0;
      clipper::AtomShapeFn sf(atoms[i].coord_orth(), atoms[i].element(),
              atoms[i].u_iso(), atoms[i].occupancy());

      // specify the derivatives we are interested in.
      sf.agarwal_params().resize(4);
      sf.agarwal_params()[0] = clipper::AtomShapeFn::X;
      sf.agarwal_params()[1] = clipper::AtomShapeFn::Y;
      sf.agarwal_params()[2] = clipper::AtomShapeFn::Z;
      sf.agarwal_params()[3] = clipper::AtomShapeFn::Uiso;
      // determine grid-ranges of this atom
      clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac(d_r.cell());
      clipper::Coord_grid g0 = uvw.coord_grid(grid) + gd.min();
      clipper::Coord_grid g1 = uvw.coord_grid(grid) + gd.max();
      clipper::Xmap<clipper::ftype64>::Map_reference_coord i0, iu, iv, iw;
      i0 = clipper::Xmap<clipper::ftype64>::Map_reference_coord(d_r, g0);
      std::vector<clipper::ftype> rho_grad(4, 0.0f);
      clipper::ftype temp_rho = 0.0f;
      // loop over grid and convolve with the atomic density gradient
      for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
          for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
            // get gradient of the atomic electron density
            sf.rho_grad(iw.coord_orth(), temp_rho, rho_grad);

            // convolve it with difference map
            const double d_r_iw = d_r[iw];
            gradient(0) += d_r_iw * rho_grad[0];
            gradient(1) += d_r_iw * rho_grad[1];
            gradient(2) += d_r_iw * rho_grad[2];

            Uiso_deriv += d_r_iw * rho_grad[3];
          }
        }
      } // loop over map
      // add to force
      force(i) -= gradient;
      b_deriv(i) = Uiso_deriv * 2.0 / to_ang * sqpi2 / to_ang / to_ang / to_ang;
      DEBUG(10, "grad(" << i << "): " << math::v2s(gradient));
    } // if atom not null
  } // for atoms
}

#endif

