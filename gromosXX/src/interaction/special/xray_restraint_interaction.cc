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
#include <set>
#include <vector>

#include "xray_restraint_interaction.h"

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
 * @param[out] force the force vector
 * @param[in] to_ang conversion factor for unit length
 */
void calculate_force_sf(clipper::FFTmap_p1 & D_k,
        clipper::Xmap<clipper::ftype32> & d_r,
        const clipper::Atom_list & atoms,
        math::VArray & force,
        double to_ang) {
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
 * fit two electron densities on top of each other. This is done using a
 * linear regression such that
 * @f[ \rho_1 = \alpha + \beta\rho_2 @f]
 * 
 * @param[in] rho1 the first electron density
 * @param[in] rho2 the second electron density
 * @param[in] points the set of grid point to consider
 * @param[out] slope slope @f$\beta@f$ of the linear regression
 * @param[out] intercept intercept @f$\alpha@f$  of the linrar regression
 */
void fit_rho(
        const clipper::Xmap<clipper::ftype32> & rho1,
        const clipper::Xmap<clipper::ftype32> & rho2,
        std::set<int> & points,
        double & slope, double & intercept) {
  DEBUG(6, "fitting electron densities");
  // if the set is empty just fill it with all grid points
  if (points.empty()) {
    DEBUG(10, "\tfitting whole grid.");
    for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = rho1.first();
          !ix.last(); ix.next()) {
      points.insert(ix.index());
    }
  }

  double sum_xy = 0.0, sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0;
  const double n = points.size();
  for(std::set<int>::const_iterator it = points.begin(), to = points.end(); it != to; ++it) {
    // get the data
    const double & x = rho2.get_data(*it);
    const double & y = rho1.get_data(*it);
    // calculate the suns
    sum_xy += x*y;
    sum_x += x;
    sum_y += y;
    sum_xx += x*x;
  }
  const double mean_x = sum_x / n;
  DEBUG(9, "mean rho2: " << mean_x);
  const double mean_y = sum_y / n;
  DEBUG(9, "mean rho1: " << mean_y);
  slope = (sum_xy - sum_x*mean_y) / (sum_xx - sum_x*mean_x);
  intercept = mean_y - slope * mean_x;
  DEBUG(9, "slope = " << slope << " intercept = " << intercept);
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
 * @parma[in] to_ang converison factor for length unit
 */
void calculate_energy_rho(const clipper::Atom_list & atoms,
        clipper::Xmap<clipper::ftype32> & rho_obs,
        const clipper::Xmap<clipper::ftype32> & rho_calc,
        const double force_constant,
        double & energy,
        math::VArray & force,
        double to_ang
) {
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

  // the fitting parameters
  double a, b;
  std::set<int> points; // empty -> whole grid
  fit_rho(rho_calc, rho_obs, points, a, b);

  // energy
  for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = rho_obs.first(),
          ix_c = rho_calc.first(); !ix.last(); ix.next(), ix_c.next()) {
    DEBUG(1, "rho_obs: " << rho_obs[ix] << " rho_calc: " << rho_calc[ix_c]);
    const double term = a*rho_obs[ix]+b - rho_calc[ix_c];
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
            const float term = a*rho_obs[iw]+b - rho_calc[iw_c];
            gradient(0) += -term * rho_grad[0];
            gradient(1) += -term * rho_grad[1];
            gradient(2) += -term * rho_grad[2];
          }
        }
      } // loop over grid
      // Angstrom -> nm
      gradient *= to_ang * force_constant * scale;
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
  const double & to_ang = sim.param().xrayrest.to_angstrom;
  m_timer.start();
  // get number of atoms in simulation
  const int atoms_size = topo.num_atoms();
  // update clipper atomvec: convert the position to Angstrom
  math::Periodicity<math::triclinic> periodicity(conf.current().box);
  for (int i = 0; i < atoms_size; i++) {
    math::Vec in_box = conf.current().pos(i);
    periodicity.put_into_positive_box(in_box);
    in_box *= to_ang;
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
    sqr_calc_free += conf.special().xray_rest[j].sf_curr * conf.special().xray_rest[j].sf_curr;
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
      fphi_obs.set_data(hkl, clipper::data32::F_phi(xrs.sf, conf.special().xray_rest[j].phase_curr));
    else
      fphi_obs.set_data(hkl, clipper::data32::F_phi(xrs.sf, conf.special().xray_rest[j].phase_av));
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
    calculate_force_sf(D_k, d_r, atoms, conf.current().force, to_ang);
    m_timer.stop("force");
  } else {
    ///////////////////////////////////////////////
    // ELECTRON DENSITY RESTRAINING
    ///////////////////////////////////////////////
    m_timer.start("energy & force");
    rho_obs.fft_from(fphi_obs);
    calculate_energy_rho(atoms, rho_obs, rho_calc, sim.param().xrayrest.force_constant,
            conf.current().energies.xray_total, conf.current().force, to_ang);
    m_timer.stop("energy & force");
  }
  DEBUG(10, "energy: " << conf.current().energies.xray_total);

  if (sim.param().xrayrest.ncsrest != simulation::xray_ncsrest_off) {
    DEBUG(6, "NCS restraints");
    m_timer.start("NCS restraints");
    std::vector<unsigned int>::const_iterator
            atom_it = topo.xray_ncs_restraints().begin(),
            atom_to = topo.xray_ncs_restraints().end();

    const clipper::Cell & cell = rho_calc.cell();

    for (; atom_it != atom_to; ++atom_it) {
      const unsigned int atom_p = *atom_it - topo.xray_asu()[0];

      switch (sim.param().xrayrest.ncsrest) {
        case simulation::xray_ncsrest_off : break;
        case simulation::xray_ncsrest_ind :
        {
          DEBUG(6, "NCS individual atoms");
          // loop over images
          for (unsigned int i = 1; i < topo.xray_asu().size(); ++i) {
            const unsigned int atom_img = topo.xray_asu()[i] + atom_p;
            // optain the image position
            const clipper::Coord_orth pos_img_ang(ncs_spacegroup.symop(i).rtop_orth(cell) * atoms[*atom_it].coord_orth());
            math::Vec pos_img(pos_img_ang.x(), pos_img_ang.y(), pos_img_ang.z());
            pos_img /= to_ang;
            DEBUG(8, "pos    : " << math::v2s(conf.current().pos(atom_img)));
            DEBUG(8, "refpos : " << math::v2s(pos_img));

            // obtain the distance between the image and the real atom
            math::Vec r;
            periodicity.nearest_image(conf.current().pos(atom_img), pos_img, r);
            DEBUG(9, "r      : " << math::v2s(r));
            const math::Vec & f = (-sim.param().xrayrest.ncs_force_constant) * r;
            DEBUG(7, "f      : " << math::v2s(f));

            conf.current().force(atom_img) += f;
            const double V = 0.5 * sim.param().xrayrest.ncs_force_constant * abs2(r);
            conf.current().energies.xray_total += V;
            DEBUG(7, "energy : " << V);
          } // loop over images
          break;
        }
        case simulation::xray_ncsrest_avg:
        {
          // loop over images
          DEBUG(6, "NCS averaged atoms");
          math::Vec avg(0.0, 0.0, 0.0);
          const math::Vec & atom_pos = conf.current().pos(*atom_it);
          DEBUG(9, "atom_pos: " << math::v2s(atom_pos));
          for (unsigned int i = 1; i < topo.xray_asu().size(); ++i) {
            const unsigned int atom_img = topo.xray_asu()[i] + atom_p;
            // optain the image position transformed to the fist ASU position
            const clipper::Coord_orth pos_img_ang(ncs_spacegroup.symop(i).inverse().rtop_orth(cell) * atoms[atom_img].coord_orth());
            math::Vec pos_img(pos_img_ang.x(), pos_img_ang.y(), pos_img_ang.z());
            pos_img /= to_ang;
            DEBUG(9, "pos_img: " << math::v2s(pos_img));

            math::Vec shift;
            periodicity.nearest_image(atom_pos, pos_img, shift);
            DEBUG(9, "shift: " << math::v2s(shift));
            avg += atom_pos + shift;
          }

          const double inv_nsym_m_one = 1.0 / (topo.xray_asu().size() - 1);

          avg *= inv_nsym_m_one;
          DEBUG(8, "avg     : " << math::v2s(avg));

          // obtain the distance between the image and the real atom
          const math::Vec r(avg - atom_pos);
          DEBUG(8, "r       : " << math::v2s(r));
          const double V = 0.5 * sim.param().xrayrest.ncs_force_constant * abs2(r);
          DEBUG(7, "V       : " << V);
          conf.current().energies.xray_total += V;

          for (unsigned int i = 1; i < topo.xray_asu().size(); ++i) {
            const unsigned int atom_img = topo.xray_asu()[i] + atom_p;
            const clipper::Mat33<double> & Sinv = ncs_spacegroup.symop(i).inverse().rtop_orth(cell).rot();


            // create the transformation matrix.
            const clipper::Coord_orth Sinv_ex_ang(Sinv * clipper::Vec3<double>(1.0, 0.0, 0.0));
            const clipper::Coord_orth Sinv_ey_ang(Sinv * clipper::Vec3<double>(0.0, 1.0, 0.0));
            const clipper::Coord_orth Sinv_ez_ang(Sinv * clipper::Vec3<double>(0.0, 0.0, 1.0));

            math::Matrix trans(
                    math::Vec(Sinv_ex_ang.x(), Sinv_ex_ang.y(), Sinv_ex_ang.z()),
                    math::Vec(Sinv_ey_ang.x(), Sinv_ey_ang.y(), Sinv_ey_ang.z()),
                    math::Vec(Sinv_ez_ang.x(), Sinv_ez_ang.y(), Sinv_ez_ang.z()));
            trans /= to_ang;
            DEBUG(8, "trans  :\n\t" << math::m2s(trans));

            const math::Vec & f = math::product(trans, (sim.param().xrayrest.ncs_force_constant*inv_nsym_m_one) * r);
            DEBUG(7, "f      : " << math::v2s(f));
            conf.current().force(atom_img) += f;
          }
          break;
        }
      } // method switch
    } // loop over restrained atoms
    m_timer.stop("NCS restraints");
  }

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

  // increase the local elevation thresholds here if required
  if (sim.steps() > 0) {
    std::vector<topology::xray_umbrella_weight_struct>::iterator it =
            topo.xray_umbrella_weights().begin(), to =
            topo.xray_umbrella_weights().end();
    for (; it != to; ++it) {
      if (!it->threshold_freeze) {
        it->threshold += it->threshold_growth_rate;
        if (it->threshold < 0.0) {
          it->threshold = 0.0;
          it->threshold_freeze = true;
          std::ostringstream msg;
          msg << "Threshold of umbrella " << it->id << " was zero and thus frozen.";
          io::messages.add(msg.str(), "Xray_Restraint_Interaction", io::message::warning);
        }
      }
    } // for umbrellas
  } // steps?

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

  const double & to_ang = sim.param().xrayrest.to_angstrom;

  math::Box & box = conf.current().box;
  const double a = math::abs(box(0));
  const double b = math::abs(box(1));
  const double c = math::abs(box(2));
  const double alpha = acos(math::costest(dot(box(1), box(2)) / (b * c)));
  const double beta = acos(math::costest(dot(box(0), box(2)) / (a * c)));
  const double gamma = acos(math::costest(dot(box(0), box(1)) / (a * b)));
  clipper::Cell_descr cellinit(a * to_ang, b * to_ang, c * to_ang,
          alpha, beta, gamma);
  clipper::Cell cell(cellinit);

  clipper::Resolution reso(sim.param().xrayrest.resolution * to_ang);

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
    atm.set_u_iso(topo.xray_b_factors()[i] * to_ang * to_ang / sqpi2);
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
    atm.set_u_iso(topo.xray_solv_b_factors()[index] * to_ang * to_ang / sqpi2);
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
    std::vector<topology::xray_umbrella_weight_struct>::iterator it =
            topo.xray_umbrella_weights().begin(), to =
            topo.xray_umbrella_weights().end();
    for (; it != to; ++it) {
      std::vector<util::Umbrella>::iterator umb_it = conf.special().umbrellas.begin(),
              umb_to = conf.special().umbrellas.end();
      bool found = false;
      for (; umb_it != umb_to; ++umb_it) {
        if (umb_it->id == it->id) {
          found = true;
          umb_it->umbrella_weight_factory = new interaction::Electron_Density_Umbrella_Weight_Factory(*it, conf, atoms, rho_calc, rho_obs, to_ang);
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

  if (sim.param().xrayrest.ncsrest != simulation::xray_ncsrest_off) {
    try {
      clipper::Spgr_descr spgrinit(clipper::String(sim.param().xrayrest.ncs_spacegroup), clipper::Spgr_descr::HM);
      ncs_spacegroup.init(spgrinit);
    } catch (const clipper::Message_fatal & msg) {
      io::messages.add("Xray_Restraint_Interaction", "NCS spacegroup: " + msg.text(), io::message::error);
      return 1;
    }
  } // if NCS rest

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

    if (sim.param().xrayrest.ncsrest != simulation::xray_ncsrest_off) {
      os << std::endl;
      switch (sim.param().xrayrest.ncsrest) {
        case simulation::xray_ncsrest_off : break;
        case simulation::xray_ncsrest_ind:
          os << "NCS restraints on individual atom positions" << std::endl;
          break;
        case simulation::xray_ncsrest_avg:
          os << "NCS restraints on averaged atom positions" << std::endl;
          break;
      }
      os << " - force constant : " << sim.param().xrayrest.ncs_force_constant << std::endl;
      os << " - Spacegroup     : " << sim.param().xrayrest.ncs_spacegroup << std::endl;
      os << " - " << topo.xray_asu().size() << " asymmetric units. First atoms: " << std::endl;
      for (unsigned int i = 0; i < topo.xray_asu().size(); ++i) {
        os << "       - " << topo.xray_asu()[i]+1 << std::endl;
      }
      os << " - the restraint atoms: " << std::endl;
      for(unsigned int i = 0; i < topo.xray_ncs_restraints().size(); ++i) {
        os << std::setw(6) << (topo.xray_ncs_restraints()[i]+1);
        if ((i+1) % 10 == 0) os << std::endl;
      }
      os << std::endl;
    }
    os << "END" << std::endl;
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
  const float cutoff2 = param.cutoff * param.cutoff * to_ang * to_ang;
  // create the range (size of atom)
  const clipper::Cell & cell = rho_calc.cell();
  const clipper::Grid_sampling & grid = rho_calc.grid_sampling();
  clipper::Grid_range gd(cell, grid, param.cutoff * to_ang);

  const int atoms_size = param.atoms.size();

  // loop over atoms to find covered grid points
  std::set<int> grid_points;

  for (int i = 0; i < atoms_size; i++) {
    DEBUG(10, "Atom index: " << param.atoms[i]);
    const clipper::Coord_orth & atom_pos = atoms[param.atoms[i]].coord_orth();

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

   // the fitting parameters
  double a, b;
  fit_rho(rho_calc, rho_obs, grid_points, a, b);
  double sum_obs_m_calc = 0.0;
  double sum_obs_p_calc = 0.0;
  DEBUG(10, "old weight: " << weight);
  // loop over grid points
  for (std::set<int>::const_iterator it = grid_points.begin(), to = grid_points.end();
          it != to; ++it) {
    // bring obs on the same scale as calc
    const double obs = a*rho_obs.get_data(*it)+b;
    DEBUG(12, "obs: " << obs);
    const double calc = rho_calc.get_data(*it);
    DEBUG(12, "calc: " << calc);

    // for R-real
    sum_obs_m_calc += fabs(obs - calc);
    sum_obs_p_calc += fabs(obs + calc);
  }
  const double r_real = sum_obs_m_calc / sum_obs_p_calc;
  DEBUG(10, "R real space: " << r_real);

  if (r_real > param.threshold) {
    const double term = r_real - param.threshold;
    weight += 0.5 * term * term;
  } else {
    param.threshold_freeze = true;
  }

  DEBUG(10, "new weight: " << weight);
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
          new interaction::Electron_Density_Umbrella_Weight(param, conf, atoms, rho_calc, rho_obs, to_ang);
  // make sure the weight is increased as soon as the instance is created
  instance->increment_weight();
  instances.push_back(instance);
  return instance;
#endif
  return NULL;
}

