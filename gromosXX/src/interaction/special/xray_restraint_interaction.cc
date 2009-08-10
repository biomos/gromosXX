/**
 * @file xray_restraint_interaction.cc
 * template methods of Xray_Restraint_Interaction
 */
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

// special interactions
#include <interaction/interaction_types.h>
#include <util/umbrella_weight.h>
#include <interaction/special/xray_restraint_interaction.h>

#include <math/periodicity.h>
#include <util/template_split.h>
#include <util/debug.h>

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special
//#define HAVE_CLIPPER
interaction::Xray_Restraint_Interaction::Xray_Restraint_Interaction() : Interaction("XrayRestraint") {
}

interaction::Xray_Restraint_Interaction::~Xray_Restraint_Interaction() {
}

#ifdef HAVE_CLIPPER
/**
 * calculate the electron density from the model structure
 * @param[out] rho_calc the calculated electron density on a map
 * @param[in] atoms the list of the atoms
 */
void calculate_electron_density(clipper::Xmap<clipper::ftype32> & rho_calc,
        const clipper::Atom_list & atoms) {
  // this code is basically copied from the clipper library but parallelised
  // some hardcoded settings
  const double radius = 2.5;

  // zero the map
  rho_calc = 0.0;
  // create the range (size of atom)
  const clipper::Cell & cell = rho_calc.cell();
  const clipper::Grid_sampling & grid = rho_calc.grid_sampling();
  clipper::Grid_range gd(cell, grid, radius);

  const int atoms_size = atoms.size();

  // loop over atoms
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < atoms_size; i++) {
    if (!atoms[i].is_null()) {
      clipper::AtomShapeFn sf(atoms[i].coord_orth(), atoms[i].element(),
              atoms[i].u_iso(), atoms[i].occupancy());
      // determine grad range of atom
      clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac(cell);
      clipper::Coord_grid g0 = uvw.coord_grid(grid) + gd.min();
      clipper::Coord_grid g1 = uvw.coord_grid(grid) + gd.max();

      // loop over atom's grid
      clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
      i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(rho_calc, g0);
      for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
        for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
          for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
            // calculate the electron density and assign it to the gird point
            const double density = sf.rho(iw.coord_orth());
#ifdef OMP
#pragma omp critical
#endif
            rho_calc[iw] += density;
          }
        }
      } // loop over grid
    }
  } // loop over atoms
  // loop over the grid again and correct the multiplicity
  for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = rho_calc.first();
          !ix.last(); ix.next())
    rho_calc[ix] *= rho_calc.multiplicity(ix.coord());
}

/**
 * calculates the force from a reciprocal space difference map
 * @param[inout] D_k the reciprocal space difference map
 * @param[out] d_r the real space difference map
 * @param[in] atoms the list containing the atoms
 * @param[out] the force vector
 */
void calculate_force_sf(clipper::FFTmap_p1 & D_k,
        clipper::Xmap<clipper::ftype32> & d_r,
        const clipper::Atom_list & atoms,
        math::VArray & force) {
  // these are just shortcuts to avoid many calls to the same functions
  const clipper::Spacegroup & spgr = d_r.spacegroup();

  // calculate the inverse symmetry operations of the spacegroup
  std::vector<clipper::Isymop> isymop;
  isymop.resize(spgr.num_symops());
  for(int j = 0; j < spgr.num_symops(); j++) {
    isymop[j] = clipper::Isymop(spgr.symop(j), d_r.grid_sampling());
  }
  const double volume = d_r.cell().volume();
  // convert from Angstrom to nm and add very annyoing scaling constants
  // to make the force volume AND resolution independent.
  const double scale = 10.0 / 2.0 * volume / (d_r.grid_sampling().size());
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
      clipper::AtomShapeFn sf(atoms[i].coord_orth(), atoms[i].element(),
              atoms[i].u_iso(), atoms[i].occupancy());

      // specify the derivatives we are interested in.
      sf.agarwal_params().resize(3);
      sf.agarwal_params()[0] = clipper::AtomShapeFn::X;
      sf.agarwal_params()[1] = clipper::AtomShapeFn::Y;
      sf.agarwal_params()[2] = clipper::AtomShapeFn::Z;
      // determine grid-ranges of this atom
      clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac(d_r.cell());
      clipper::Coord_grid g0 = uvw.coord_grid(grid) + gd.min();
      clipper::Coord_grid g1 = uvw.coord_grid(grid) + gd.max();
      clipper::Xmap<clipper::ftype64>::Map_reference_coord i0, iu, iv, iw;
      i0 = clipper::Xmap<clipper::ftype64>::Map_reference_coord(d_r, g0);
      std::vector<clipper::ftype> rho_grad(3, 0.0f);
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
          }
        }
      } // loop over map
      // add to force
      force(i) -= gradient;
      DEBUG(10, "grad(" << i << "): " << math::v2s(gradient));
    } // if atom not null
  } // for atoms
}

/**
 * calculate the energy for structure factor restraining
 * @param[in] refl the observed relefections
 * @param[in] refl_curr the calculated reflections
 * @param[in] averaging the averaging mode of the restraining
 * @param[in] k_inst the inst. scaling factor
 * @param[in] k_avg the avg. scaling factor
 * @param[out] D_k the difference map for gradients
 * @param[in] force_constant the force constant
 * @param[out] energy the energy obtained
 */
void calculate_energy_sf(const std::vector<topology::xray_restraint_struct> & refl,
        const std::vector<configuration::Configuration::special_struct::xray_struct> & refl_curr,
        simulation::xrayrest_enum averaging,
        const double k_inst, const double k_avg,
        clipper::FFTmap_p1 & D_k, 
        const double force_constant,
        double & energy) {
  // zero the reciprocal space difference map
  D_k.reset();

  double energy_sum = 0.0;
  // loop over retraints and calculate energy and difference map
  for (unsigned int i = 0; i < refl.size(); i++) {
    const topology::xray_restraint_struct & xrs = refl[i];
    clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // SWITCH FOR DIFFERENT METHODS
    switch (averaging) {
      case simulation::xrayrest_inst :
      {
        // INSTANTANEOUS
        // calculate energy-sum
        const double fobs = xrs.sf;
        const double fcalc = refl_curr[i].sf_curr;
        const double term = fobs - k_inst * fcalc;
        energy_sum += term * term;
        // calculate derivatives of target function
        const double dterm = (k_inst * fcalc - fobs) * k_inst;
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
        const double fcalc = refl_curr[i].sf_av;
        const double term = fobs - k_avg * fcalc;
        energy_sum += term * term;
        // calculate derivatives of target function
        // here we omit the 1-exp(-dt/tau) term.
        const double dterm = (k_avg * fcalc - fobs) * k_avg;
        D_k.set_hkl(hkl, clipper::data32::F_phi(force_constant * dterm, refl_curr[i].phase_curr));
        break;
      }
      case simulation::xrayrest_biq :
      {
        // BIQUADRATIC TIME-AVERAGED/INSTANTANEOUS
        // calculate energy-sum
        const double fobs = xrs.sf;
        const double finst = refl_curr[i].sf_curr;
        const double favg = refl_curr[i].sf_av;
        const double inst_term = fobs - k_inst * finst;
        const double av_term = fobs - k_avg * favg;
        energy_sum += (inst_term * inst_term)*(av_term * av_term);
        // calculate derivatives of target function
        // here we omit the 1-exp(-dt/tau) term.
        const double dterm = (k_inst * finst - fobs)*(av_term * av_term) * k_inst
                + (k_avg * favg - fobs)*(inst_term * inst_term) * k_avg;
        D_k.set_hkl(hkl, clipper::data32::F_phi(force_constant * dterm, refl_curr[i].phase_curr));
        break;
      }
      default: break;
    }
  }

  // finally calculate the energy
  energy = 0.5 * force_constant * energy_sum;
}

/**
 * calculate energy of the electron density restraining
 * @param[in] atoms the list of the atoms
 * @param[in] rho_obs the "observed" electron density
 * @param[in] force_constant the force constant
 * @param[out] energy the energy obtained
 * @param[out] force the forces are added to this vector
 */
void calculate_energy_rho(const clipper::Atom_list & atoms,
        clipper::Xmap<clipper::ftype32> & rho_obs,
        const clipper::Xmap<clipper::ftype32> & rho_calc,
        const double force_constant,
        double & energy,
        math::VArray & force) {
  const double radius = 3.5;
  energy = 0.0;

  // create the range (size of atom)
  const clipper::Cell & cell = rho_calc.cell();
  const clipper::Grid_sampling & grid = rho_calc.grid_sampling();
  clipper::Grid_range gd(cell, grid, radius);

  const int atoms_size = atoms.size();
  const double volume = cell.volume();
  // convert from Angstrom to nm and add very annyoing scaling constants
  // to make the force volume AND resolution independent.
  const double scale = volume / grid.size();

  // energy
  for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = rho_obs.first(),
          ix_c = rho_calc.first(); !ix.last(); ix.next(), ix_c.next()) {
    const double term = rho_obs[ix] - rho_calc[ix_c];
    energy += term * term;
  }
  energy *= 0.5 * force_constant * scale;
  // loop over atoms
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < atoms_size; i++) {
    if (!atoms[i].is_null()) {
      math::Vec gradient(0.0, 0.0, 0.0);
      clipper::AtomShapeFn sf(atoms[i].coord_orth(), atoms[i].element(),
              atoms[i].u_iso(), atoms[i].occupancy());
      sf.agarwal_params().resize(3);
      sf.agarwal_params()[0] = clipper::AtomShapeFn::X;
      sf.agarwal_params()[1] = clipper::AtomShapeFn::Y;
      sf.agarwal_params()[2] = clipper::AtomShapeFn::Z;

      // determine grad range of atom
      clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac(cell);
      clipper::Coord_grid g0 = uvw.coord_grid(grid) + gd.min();
      clipper::Coord_grid g1 = uvw.coord_grid(grid) + gd.max();

      clipper::ftype rho;
      std::vector<clipper::ftype> rho_grad(3, 0.0f);

      // loop over atom's grid
      clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
      i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(rho_obs, g0);
      clipper::Xmap<clipper::ftype32>::Map_reference_coord i0_c, iu_c, iv_c, iw_c;
      i0_c = clipper::Xmap<clipper::ftype32>::Map_reference_coord(rho_calc, g0);
      for (iu = i0, iu_c = i0_c; iu.coord().u() <= g1.u(); iu.next_u(), iu_c.next_u()) {
        assert(iu.coord().u() == iu_c.coord().u());
        for (iv = iu, iv_c = iu_c; iv.coord().v() <= g1.v(); iv.next_v(), iv_c.next_v()) {
          assert(iv.coord().v() == iv_c.coord().v());
          for (iw = iv, iw_c = iv_c; iw.coord().w() <= g1.w(); iw.next_w(), iw_c.next_w()) {
            assert(iw.coord().w() == iw_c.coord().w());
            // calculate electron density and gradient of it.
            sf.rho_grad(iw.coord_orth(), rho, rho_grad);
            const float term = rho_obs[iw] - rho_calc[iw_c];
            gradient(0) += -term * rho_grad[0];
            gradient(1) += -term * rho_grad[1];
            gradient(2) += -term * rho_grad[2];
          }
        }
      } // loop over grid
      // Angstrom -> nm
      gradient *= 10.0 * force_constant * scale;
      force(i) -= gradient;
    }
  } // loop over atoms
}

#endif
/**
 * calculate xray restraint interactions
 */
int interaction::Xray_Restraint_Interaction
::calculate_interactions(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
#ifdef HAVE_CLIPPER
  m_timer.start();
  // get number of atoms in simulation
  const int atoms_size = topo.num_atoms();
  // update clipper atomvec: convert the position to Angstrom
  math::Periodicity<math::triclinic> periodicity(conf.current().box);
  for (int i = 0; i < atoms_size; i++) {
    math::Vec in_box = conf.current().pos(i);
    periodicity.put_into_positive_box(in_box);
    in_box *= 10;
    atoms[i].set_coord_orth(clipper::Coord_orth(in_box(0),
            in_box(1), in_box(2)));
  }
  // Calculate structure factors
  m_timer.start("structure factor");
  calculate_electron_density(rho_calc, atoms);
  // FFT the electron density to obtain the structure factors
  rho_calc.fft_to(fphi);
  m_timer.stop("structure factor");

  // sqr_calc:       sum of squared Fcalc
  // obs:            sum of Fobs
  // calc:           sum of Fcalc
  // obs_calc:       sum of Fobs*Fcalc
  // obs_calcavg:    sum of Fobs*Fcalc(averaged)
  // obs_k_calcavg:  sum of Fobs-k_avg*Fcalc(averaged)
  // obs_k_calc:     sum of Fobs-k*Fcalc
  // sqr_calcavg:    sum of squared time-averaged Fcalc
  // calcavg:        sum of time-averaged Fcalc
   m_timer.start("scaling");
  // zero all the sums
  double sqr_calc = 0.0, obs = 0.0, obs_free = 0.0, calc = 0.0, obs_calc = 0.0, obs_k_calc = 0.0,
         sqr_calcavg = 0.0, calcavg = 0.0, obs_calcavg = 0.0, obs_k_calcavg = 0.0,
         obs_k_calc_free = 0.0, obs_k_calcavg_free = 0.0, obs_calc_free = 0.0, obs_calcavg_free = 0.0,
         sqr_calc_free = 0.0, sqr_calcavg_free = 0.0;
  // Number of reflections
  const unsigned int num_xray_rest = topo.xray_restraints().size();
  const unsigned int num_xray_rfree = topo.xray_rfree().size();
  // e-term for time-average
  const double eterm = exp(-sim.time_step_size() / sim.param().xrayrest.tau);

  // loop over structure factors
  unsigned int j = 0;
  for (unsigned int i = 0; i < num_xray_rest; i++, j++) {
    // filter calculated structure factors: save phases and amplitudes
    clipper::HKL hkl(topo.xray_restraints()[i].h, topo.xray_restraints()[i].k, topo.xray_restraints()[i].l);
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

    // calc sums
    obs_calc += conf.special().xray_rest[j].sf_curr * topo.xray_restraints()[i].sf;
    obs_calcavg += conf.special().xray_rest[j].sf_av * topo.xray_restraints()[i].sf;
    sqr_calc += conf.special().xray_rest[j].sf_curr * conf.special().xray_rest[i].sf_curr;
    obs += topo.xray_restraints()[i].sf;
    calc += conf.special().xray_rest[j].sf_curr;
    sqr_calcavg += conf.special().xray_rest[j].sf_av * conf.special().xray_rest[j].sf_av;
    calcavg += conf.special().xray_rest[j].sf_av;
  }
  // loop over structure factors in R free set
  for (unsigned int i = 0; i < num_xray_rfree; i++, j++) {
    // filter calculated structure factors: save phases and amplitudes for R free HKLs
    clipper::HKL hkl(topo.xray_rfree()[i].h, topo.xray_rfree()[i].k, topo.xray_rfree()[i].l);
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

    // calc sums
    obs_free += topo.xray_rfree()[i].sf;
    obs_calc_free += conf.special().xray_rest[j].sf_curr * topo.xray_rfree()[i].sf;
    obs_calcavg_free += conf.special().xray_rest[j].sf_av * topo.xray_rfree()[i].sf;
    sqr_calc_free += conf.special().xray_rest[j].sf_curr * conf.special().xray_rest[i].sf_curr;
    sqr_calcavg_free += conf.special().xray_rest[j].sf_av * conf.special().xray_rest[j].sf_av;
  }
  // check for possible resolution problems
#ifdef HAVE_ISNAN
  if (std::isnan(calc)){
    io::messages.add("Structure factors were NaN. This can be due to numerical problems. "
                     "Try to slighlty increase the resolution.", "X-Ray Restraints", io::message::error);
    return 1;
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
    obs_k_calc += fabs(xrs.sf - k_inst * conf.special().xray_rest[j].sf_curr);
    obs_k_calcavg += fabs(xrs.sf - k_avg * conf.special().xray_rest[j].sf_av);
    
    clipper::HKL hkl(xrs.h, xrs.k, xrs.l);
    // save Fobs and PhiCalc for density maps. This will be corrected
    // for symmetry in the FFT step.
    if (sim.param().xrayrest.xrayrest == simulation::xrayrest_inst)
      fphi_obs.set_data(hkl, clipper::data32::F_phi(xrs.sf / k_inst, conf.special().xray_rest[j].phase_curr));
    else
      fphi_obs.set_data(hkl, clipper::data32::F_phi(xrs.sf / k_avg, conf.special().xray_rest[j].phase_av));
  }
  // and for R free
  for (unsigned int i = 0; i < num_xray_rfree; i++, j++) {
    const topology::xray_restraint_struct & xrs = topo.xray_rfree()[i];
    obs_k_calc_free += fabs(xrs.sf - k_free_inst * conf.special().xray_rest[j].sf_curr);
    obs_k_calcavg_free += fabs(xrs.sf - k_free_avg * conf.special().xray_rest[j].sf_av);
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
  m_timer.stop("scaling");

  if (sim.param().xrayrest.local_elevation) {
    ///////////////////////////////////////////////
    // LOCAL ELEVATION
    ///////////////////////////////////////////////
    m_timer.start("obs. electron density");
    rho_obs.fft_from(fphi_obs);
    m_timer.stop("obs. electron density");
  } else if (sim.param().xrayrest.mode != simulation::xrayrest_mode_electron_density) {
    ///////////////////////////////////////////////
    // STRUCTURE FACTOR RESTRAINING
    ///////////////////////////////////////////////
    m_timer.start("energy");
    calculate_energy_sf(topo.xray_restraints(), conf.special().xray_rest,
            sim.param().xrayrest.xrayrest, k_inst, k_avg,
            D_k, sim.param().xrayrest.force_constant, conf.current().energies.xray_total);
    m_timer.stop("energy");

    // start to calculate the forces
    m_timer.start("force");
    calculate_force_sf(D_k, d_r, atoms, conf.current().force);
    m_timer.stop("force");
  } else {
    ///////////////////////////////////////////////
    // ELECTRON DENSITY RESTRAINING
    ///////////////////////////////////////////////
    m_timer.start("energy & force");
    rho_obs.fft_from(fphi_obs);
    calculate_energy_rho(atoms, rho_obs, rho_calc, sim.param().xrayrest.force_constant,
            conf.current().energies.xray_total, conf.current().force);
    m_timer.stop("energy & force");
  }
  DEBUG(10, "energy: " << conf.current().energies.xray_total);

  m_timer.stop();

  // write xmap to external file
  if (sim.param().xrayrest.writexmap != 0 && sim.steps() % sim.param().xrayrest.writexmap == 0) {
    if (sim.param().xrayrest.mode != simulation::xrayrest_mode_electron_density)
      rho_obs.fft_from(fphi_obs);
    clipper::CCP4MAPfile mapfile;
    std::ostringstream file_name, asu_file_name;
    file_name << "density_frame_" << std::setw(int(log10(sim.param().step.number_of_steps)))
            << std::setfill('0') << sim.steps() << ".ccp4";
    if (sim.param().xrayrest.writedensity == 1 || sim.param().xrayrest.writedensity == 3) {
      mapfile.open_write(file_name.str());
      mapfile.export_xmap(rho_obs);
      mapfile.close_write();
    }
    // Non Cristallographic Map
    clipper::NXmap<double> asu(rho_obs.grid_asu(), rho_obs.operator_orth_grid());
    asu_file_name << "density_asu_frame_" << std::setw(int(log10(sim.param().step.number_of_steps)))
            << std::setfill('0') << sim.steps() << ".ccp4";
    if (sim.param().xrayrest.writedensity == 2 || sim.param().xrayrest.writedensity == 3) {
      mapfile.open_write(asu_file_name.str());
      mapfile.export_nxmap(asu);
      mapfile.close_write();
    }
  }
#endif

  return 0;
}

int interaction::Xray_Restraint_Interaction::init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {
#ifdef HAVE_CLIPPER
  DEBUG(15, "Xray_Restraint_Interaction: init")
  const double sqpi2 = (math::Pi * math::Pi * 8.0);

  // Redirect clipper errors
  clipper::Message message;
  message.set_stream(os);

  // Construct clipper objects
  clipper::Spacegroup spacegr;
  try {
    clipper::Spgr_descr spgrinit(clipper::String(sim.param().xrayrest.spacegroup), clipper::Spgr_descr::HM);
    spacegr.init(spgrinit);
  } catch (const clipper::Message_fatal & msg) {
    io::messages.add("Xray_restraint_interaction", msg.text(), io::message::error);
    return 1;
  }

  math::Box & box = conf.current().box;
  const double a = math::abs(box(0));
  const double b = math::abs(box(1));
  const double c = math::abs(box(2));
  const double alpha = acos(math::costest(dot(box(1), box(2)) / (b * c)));
  const double beta = acos(math::costest(dot(box(0), box(2)) / (a * c)));
  const double gamma = acos(math::costest(dot(box(0), box(1)) / (a * b)));
  clipper::Cell_descr cellinit(a * 10.0, b * 10.0, c * 10.0,
          alpha, beta, gamma);
  clipper::Cell cell(cellinit);

  clipper::Resolution reso(sim.param().xrayrest.resolution * 10.0);

    // create a grid and a crystalographic map
  const double shannon_rate = 1.5;
  const clipper::Grid_sampling grid_rho_calc(spacegr, cell, reso, shannon_rate);
  rho_calc.init(spacegr, cell, grid_rho_calc);
  const clipper::Grid_sampling grid_rho_obs(spacegr, cell, reso, shannon_rate);
  rho_obs.init(spacegr, cell, grid_rho_obs);

  hkls.init(spacegr, cell, reso, true);
  fphi.init(hkls, hkls.cell());
  fphi_obs.init(hkls, hkls.cell());

  // The difference map has to be a P 1 map in order to get agreement with
  // the finite difference results. However, the reasons for this are
  // not 100% clear. In principle (from theory) a spacegroup depdendent
  // Xmap should do the job.
  clipper::Spgr_descr spgrinit(clipper::String("P 1"), clipper::Spgr_descr::HM);
  clipper::Spacegroup p1_spacegr;
  p1_spacegr.init(spgrinit);
  // create a grid and a P1 FFT map. Here we can use the FFTmap_p1 which was
  // designed for fast P1.
  // 1.5 is shannon-rate for oversampled FFT
  const clipper::Grid_sampling fftgrid(p1_spacegr, cell, reso, 1.5);
  D_k.init(fftgrid);
  // we still need an Xmap for convenient looping over the data in order
  // to do the convolution.
  const clipper::Grid_sampling grid_d_r(spacegr, cell, reso, 1.5);
  d_r.init(spacegr, cell, grid_d_r);

  // Fill clipper atom-vector
  std::vector<clipper::Atom> atomvec;
  // Fill solute
  for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
    clipper::Atom atm;
    assert(i < topo.xray_occupancies().size());
    atm.set_occupancy(topo.xray_occupancies()[i]);
    atm.set_coord_orth(clipper::Coord_orth(0.0, 0.0, 0.0));
    assert(i < topo.xray_b_factors().size());
    atm.set_u_iso(topo.xray_b_factors()[i] * 100.0 / sqpi2);
    assert(i < topo.xray_elements().size());
    atm.set_element(topo.xray_elements()[i]);
    atomvec.push_back(atm);
  }
  // Fill solvent
  for (unsigned int i = 0; i < topo.num_solvent_atoms(); i++) {
    clipper::Atom atm;
    unsigned int index = i % topo.solvent(0).num_atoms();
    assert(index < topo.xray_solv_occupancies().size());
    atm.set_occupancy(topo.xray_solv_occupancies()[index]);
    atm.set_coord_orth(clipper::Coord_orth(0.0, 0.0, 0.0));
    assert(index < topo.xray_solv_b_factors().size());
    atm.set_u_iso(topo.xray_solv_b_factors()[index] * 100 / sqpi2);
    assert(index < topo.xray_solvelements().size());
    atm.set_element(topo.xray_solvelements()[index]);
    atomvec.push_back(atm);
  }
  atoms = clipper::Atom_list(atomvec);

  conf.special().xray_rest.resize(topo.xray_restraints().size() + topo.xray_rfree().size());

  // Scale Fobs (needed for constant force-constant) -> scaled to sfscale
  const double sfscale = 100.0;
  double maxsf = 0.0;
  // Get max structure factor
  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    if (maxsf < topo.xray_restraints()[i].sf)
      maxsf = topo.xray_restraints()[i].sf;
  }
  for (unsigned int i = 0; i < topo.xray_rfree().size(); i++) {
    if (maxsf < topo.xray_rfree()[i].sf)
      maxsf = topo.xray_rfree()[i].sf;
  }
  const double scalefactor = maxsf / sfscale;
  // scale
  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    topo.xray_restraints()[i].sf /= scalefactor;
  }
  for (unsigned int i = 0; i < topo.xray_rfree().size(); i++) {
    topo.xray_rfree()[i].sf /= scalefactor;
  }

  if (sim.param().xrayrest.local_elevation) {
    std::vector<topology::xray_umbrella_weight_struct>::const_iterator it =
            topo.xray_umbrella_weights().begin(), to =
            topo.xray_umbrella_weights().end();
    for (; it != to; ++it) {
      std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
              umb_to = conf.special().umbrellas.end();
      bool found = false;
      for (; umb_it != umb_to; ++umb_it) {
        if (umb_it->id == it->id) {
          found = true;
          umb_it->umbrella_weight_factory = new interaction::Electron_Density_Umbrella_Weight_Factory(it->atoms,
                  it->threshold, it->cutoff, conf, atoms, rho_calc, rho_obs);
        }

      }
      if (!found) {
        std::ostringstream msg;
        msg << "Cannot find umbrella " << it->id;
        io::messages.add(msg.str(), "Xray_Restraint_Interaction", io::message::error);
        return 1;
      }
    } // for weights
  } // if local elev

  if (!quiet) {
    os.precision(2);
    os << "\nXRAYREST\n";
    os << "Restraint type              : ";
    switch (sim.param().xrayrest.xrayrest) {
      case simulation::xrayrest_off :
        os << "No xray restraining";
        break;
      case simulation::xrayrest_inst :
        os << "Instantaneous restraining";
        break;
      case simulation::xrayrest_avg :
        os << "Time-averaged restraining";
        break;
      case simulation::xrayrest_biq :
        os << "Biquadratic time-averaged/instantaneous restraining";
        break;
    }
    if (sim.param().xrayrest.mode == simulation::xrayrest_mode_electron_density) {
      os << " on the electron density";
      if (sim.param().xrayrest.local_elevation)
        os << " using local elevation";
      os << ".\n";
    } else {
      os << " on the structure factors.\n";
    }

    if (sim.param().xrayrest.readavg)
      os << "\treading xray averages from file\n";

    os << "Restraint force-constant    : " << sim.param().xrayrest.force_constant << std::endl;
    os << "Spacegroup                  : " << sim.param().xrayrest.spacegroup << std::endl;
    os.precision(4);
    os << "Resolution                  : " << sim.param().xrayrest.resolution << std::endl;
    os << "Num experimental reflections: " << topo.xray_restraints().size() << std::endl;
    os << "Num R-free reflections      : " << topo.xray_rfree().size() << std::endl;
    os << "Num expected reflections    : " << hkls.num_reflections() << std::endl;
    os << "Writeing electron density   : " << sim.param().xrayrest.writexmap << std::endl << std::endl;
    if (sim.param().xrayrest.local_elevation) {
      os << "The following local elevation umbrellas are weighted by the electron density:" << std::endl;
      std::vector<topology::xray_umbrella_weight_struct>::const_iterator it =
              topo.xray_umbrella_weights().begin(), to =
              topo.xray_umbrella_weights().end();
      for (; it != to; ++it) {
        os.precision(5);
        os << "  - " << std::setw(5) << it->id << ": Threshold " << std::setw(15) << it->threshold << " atoms: ";
        for (std::vector<unsigned int>::const_iterator a_it = it->atoms.begin(),
                a_to = it->atoms.end(); a_it != a_to; ++a_it) {
          os << std::setw(5) << (*a_it) + 1;
        }
        os << std::endl;
      }
    }
    os << "END\n\n";
  }
  // Check if too low resolution
  double expnorm = 0.0, calcnorm = 0.0, tempnorm = 0.0;
  for (unsigned int i = 0; i < topo.xray_restraints().size(); i++) {
    // calc max. experimental-index norm
    tempnorm = sqrt(double((topo.xray_restraints()[i].h * topo.xray_restraints()[i].h)+
            (topo.xray_restraints()[i].k * topo.xray_restraints()[i].k)+
            (topo.xray_restraints()[i].l * topo.xray_restraints()[i].l)));
    if (tempnorm > expnorm)
      expnorm = tempnorm;
  }
  for (unsigned int i = 0; i < topo.xray_rfree().size(); i++) {
    // calc max. experimental-index norm
    tempnorm = sqrt(double((topo.xray_rfree()[i].h * topo.xray_rfree()[i].h)+
            (topo.xray_rfree()[i].k * topo.xray_rfree()[i].k)+
            (topo.xray_rfree()[i].l * topo.xray_rfree()[i].l)));
    if (tempnorm > expnorm)
      expnorm = tempnorm;
  }
  for (int i = 0; i < hkls.num_reflections(); i++) {
    // calc max. calculation-index norm
    tempnorm = sqrt(double((hkls.hkl_of(i).h() * hkls.hkl_of(i).h())+
            (hkls.hkl_of(i).k() * hkls.hkl_of(i).k())+
            (hkls.hkl_of(i).l() * hkls.hkl_of(i).l())));
    if (tempnorm > calcnorm)
      calcnorm = tempnorm;
  }

  if (expnorm > calcnorm) {
    io::messages.add("Xray_restraint_interaction", "Too little reflections. Set higher resolution!", io::message::error);
    return 1;
  }
#endif
  return 0;
}

void interaction::Electron_Density_Umbrella_Weight::increment_weight() {
#ifdef HAVE_CLIPPER
  // convert to angstrom
  const float cutoff2 = cutoff * cutoff * 100.0;
  // create the range (size of atom)
  const clipper::Cell & cell = rho_calc.cell();
  const clipper::Grid_sampling & grid = rho_calc.grid_sampling();
  clipper::Grid_range gd(cell, grid, cutoff * 10.0);

  const int atoms_size = variable_atoms.size();
  const double volume = cell.volume();

  // loop over atoms to find covered grid points
  std::set<int> grid_points;

  for (int i = 0; i < atoms_size; i++) {
    DEBUG(1, "Atom index: " << variable_atoms[i]);
    const clipper::Coord_orth & atom_pos = atoms[variable_atoms[i]].coord_orth();

    // determine grad range of atom
    const clipper::Coord_frac uvw = atom_pos.coord_frac(cell);
    const clipper::Coord_grid g0 = uvw.coord_grid(grid) + gd.min();
    const clipper::Coord_grid g1 = uvw.coord_grid(grid) + gd.max();

    // loop over atom's grid
    clipper::Xmap<clipper::ftype32>::Map_reference_coord i0, iu, iv, iw;
    i0 = clipper::Xmap<clipper::ftype32>::Map_reference_coord(rho_obs, g0);
    for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
      for (iv = iu; iv.coord().v() <= g1.v(); iv.next_v()) {
        for (iw = iv; iw.coord().w() <= g1.w(); iw.next_w()) {
          // check if this cell is within the cutoff
          if ((iw.coord_orth() - atom_pos).lengthsq() < cutoff2) {
            grid_points.insert(iw.index());
          } // if within cutoff
        }
      }
    } // loop over grid
  } // loop over atoms
  // now we have a list of grid points without duplicates. Calculate density
  // differences (weight)
 
  const double scale = volume / grid.size();
  const double threshold_scaled = threshold; // / sqrt(scale);
  DEBUG(1, "threshold scaled: " << threshold_scaled);
  double diff = 0.0;
   // loop over grid points
  for (std::set<int>::const_iterator it = grid_points.begin(), to = grid_points.end();
          it != to; ++it) {
    const double obs = rho_obs.get_data(*it);
    const double calc = rho_calc.get_data(*it);
    DEBUG(1, "rho_obs: " << obs << " rho_calc: " << calc << " diff: " << fabs(obs - calc));
    // this is the implementation of the flat bottom potential
    const double rhoobs_m_thres = obs - threshold_scaled;
    const double rhoobs_p_thres = obs + threshold_scaled;

    if (calc < rhoobs_m_thres) {
      const double term = rhoobs_m_thres - calc;
      diff += term * term;
    } else if (calc > rhoobs_p_thres) {
      const double term = rhoobs_p_thres - calc;
      diff += term * term;
    }
  }
  // correct for volume
  weight += diff * scale;
#endif
}

void interaction::Electron_Density_Umbrella_Weight::write(std::ostream& os) const {
  os.precision(8);
  os.flags(std::ios::scientific);
  os << std::setw(15) << weight;
}

util::Umbrella_Weight * interaction::Electron_Density_Umbrella_Weight_Factory::get_instance() {
#ifdef HAVE_CLIPPER
  util::Umbrella_Weight * instance =
          new interaction::Electron_Density_Umbrella_Weight(variable_atoms,
            threshold, cutoff, conf, atoms, rho_calc, rho_obs);
  // make sure the weight is increased as soon as the instance is created
  instance->increment_weight();
  instances.push_back(instance);
  return instance;
#endif
}

