/**
 * @file dens.cc
 * electron density business
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
#endif

#ifdef OMP
#include <omp.h>
#endif

#ifdef HAVE_CLIPPER
void interaction::xray::calculate_electron_density(clipper::Xmap<clipper::ftype32> & rho_calc,
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

void interaction::xray::fit_rho(
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

void interaction::xray::calculate_energy_rho(const clipper::Atom_list & atoms,
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
