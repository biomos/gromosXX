/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file dihedral_constraint.cc
 * contains the dihedral constraint iteration
 * as its a template, the file is included in shake.cc
 */

/**
 * Dihedral Constraints
 *
 * see: 
 * J Phys Chem B. 2006 Apr 27;110(16):8488-98 (Appendix)
 * and J. Chem. Phys. 152, 024109 (2020)
 *
 */
#ifndef DIHEDRAL_CONSTRAINT_CC
#define DIHEDRAL_CONSTRAINT_CC
template <math::boundary_enum B, math::virial_enum V>
int algorithm::Shake::dih_constr_iteration(
    topology::Topology const &topo,
    configuration::Configuration &conf,
    simulation::Simulation const &sim,
    bool &convergence,
    std::vector<bool> &skip_now,
    std::vector<bool> &skip_next,
    std::vector<topology::dihedral_restraint_struct> const &dihedral_restraints,
    math::Periodicity<B> const &periodicity)
{
  const double tolerance = sim.param().dihrest.tolerance;

  convergence = true;

  math::VArray &pos = conf.current().pos;
  math::VArray &ref = conf.old().pos;

  std::vector<topology::dihedral_restraint_struct>::const_iterator
      it = dihedral_restraints.begin(),
      to = dihedral_restraints.end();

  for (; it != to; ++it)
  {

    // check whether we can skip this constraint
    if (skip_now[it->i] && skip_now[it->j] &&
        skip_now[it->k] && skip_now[it->l])
      continue;

    // calculate dihedral angle
    DEBUG(9, "dihedral angle " << it->i << "-" << it->j << "-" << it->k << "-" << it->l);

    math::Vec r12, r32, r34;
    periodicity.nearest_image(pos(it->i), pos(it->j), r12);
    periodicity.nearest_image(pos(it->k), pos(it->j), r32);
    periodicity.nearest_image(pos(it->k), pos(it->l), r34);

    // eq numbers from J. Chem. Phys. 152, 024109 (2020); doi: 10.1063/1.5124923

    // eq 3
    const math::Vec r52 = math::cross(r12, r32);
    const double d52 = math::abs(r52);
    // eq 4
    const math::Vec r63 = math::cross(r32, r34);
    const double d63 = math::abs(r63);
    // eq 11
    const int sign_phi = (math::dot(r12, r63) >= 0.0) ? 1 : -1;
    // eq 35
    const double cos_phi = math::dot(r52, r63) / (d52 * d63);

    double phi= 0.0;
    // cos_phi can be >1 or <-1 because of precision limits
    if (cos_phi > 1)
      phi = 0.0;
    else if (cos_phi < -1)
      phi = math::Pi;
    else
      phi = sign_phi * acos(cos_phi);

    const double phi0 = it->phi;

    while (phi < phi0 - math::Pi)
      phi += 2 * math::Pi;
    while (phi > phi0 + math::Pi)
      phi -= 2 * math::Pi;

    // decide if constraint is fulfilled using phi
    // and not cos(phi)

    const double diff = phi - phi0;
    DEBUG(8, "phi=" << std::setprecision(7) << 180 * phi / math::Pi
                    << "\tphi0=" << 180 * phi0 / math::Pi
                    << std::setprecision(7) << "\tdiff=" << 180 * diff / math::Pi);

    if (fabs(diff) >= tolerance)
    {
      // we have to shake
      DEBUG(10, "shaking");

      math::Vec ref12, ref32, ref34;
      periodicity.nearest_image(ref(it->i), ref(it->j), ref12);
      periodicity.nearest_image(ref(it->k), ref(it->j), ref32);
      periodicity.nearest_image(ref(it->k), ref(it->l), ref34);

      // eq 3
      const math::Vec ref52 = math::cross(ref12, ref32);
      const double dref52 = math::abs(ref52);
      // eq 4
      const math::Vec ref63 = math::cross(ref32, ref34);
      const double dref63 = math::abs(ref63);

      const double dref32 = math::abs(ref32);

      const math::Vec a1 = dref32 / (dref52 * dref52) * ref52;

      const math::Vec a2 =
          (math::dot(ref12, ref32) / (dref32 * dref32) - 1) *
              dref32 / (dref52 * dref52) * ref52 +
          math::dot(ref34, ref32) / (dref32 * dref32) *
              dref32 / (dref63 * dref63) * ref63;

      const math::Vec a3 = -((math::dot(ref34, ref32) / (dref32 * dref32) - 1) *
                                 dref32 / (dref63 * dref63) * ref63 +
                             math::dot(ref12, ref32) / (dref32 * dref32) *
                                 dref32 / (dref52 * dref52) * ref52);

      const math::Vec a4 = -dref32 / (dref63 * dref63) * ref63;

      const double m1 = topo.mass(it->i);
      const double m2 = topo.mass(it->j);
      const double m3 = topo.mass(it->k);
      const double m4 = topo.mass(it->l);

      // eq 58
      const math::Vec b123 =
          math::cross(r12, a3 / m3 - a2 / m2) -
          math::cross(r32, a1 / m1 - a2 / m2);

      // eq 59, 60
      const math::Vec b234 =
          math::cross(r32, a3 / m3 - a4 / m4) -
          math::cross(r34, a3 / m3 - a2 / m2);

      const double c3 = d52 * d63;
      const double c4 = d52 * math::dot(r63, b234) / d63 + d63 * math::dot(r52, b123) / d52;

      //////////////////////////////////////////////////

      double l_dt2 = 0;
      double phi0_abs = fabs(it->phi);
      if (phi0_abs <= math::Pi / 4 || phi0_abs > 3 * math::Pi / 4)
      { // sine case
        math::Vec c5 = math::cross(r52, r63);
        math::Vec c6 = math::cross(r52, b234) - math::cross(r63, b123);
        double d32 = math::abs(r32);

        math::Vec am32 = (a3 / m3 - a2 / m2);

        double sinphi0 = sin(it->phi);
        double sin_nom = 0.0, sin_denom = 0.0;

        sin_nom = sinphi0 * d32 * c3 - math::dot(c5, r32);

        double T1 = c3 / d32 * math::dot(r32, am32) + c4 * d32;
        double T2 = math::dot(c5, am32) + math::dot(c6, r32);

        sin_denom = sinphi0 * T1 - T2;

        // denominator might be zero
        const double epsilon = math::abs(r12) * 0.0000000001;
        if (fabs(sin_denom) < epsilon)
        {
          sin_denom = epsilon;
        }

        l_dt2 = sin_nom / sin_denom;
      }
      else if (phi0_abs > math::Pi / 4 && phi0_abs <= 3 * math::Pi / 4)
      { // cosine case
        // eq 61,62
        const double c1234 =
            math::dot(
                math::cross(r12, r32),
                math::cross(r32, r34));

        const double d1234 =
            math::dot(
                math::cross(r12, r32),
                b234) +
            math::dot(
                math::cross(r32, r34),
                b123);

        double nominator = 0.0, denominator = 0.0;
        nominator = cos(it->phi) * c3 - c1234;
        denominator = cos(it->phi) * c4 - d1234;

        // denominator might be zero
        const double epsilon = math::abs(r12) * 0.0000000001;
        if (fabs(denominator) < epsilon)
        {
          denominator = epsilon;
        }

        l_dt2 = nominator / denominator;
      }

      //////////////////////////////////////////////////

      // eq. 48, 50
      pos(it->i) -= l_dt2 * a1 / m1;
      pos(it->j) -= l_dt2 * a2 / m2;
      pos(it->k) -= l_dt2 * a3 / m3;
      pos(it->l) -= l_dt2 * a4 / m4;

      if (V == math::atomic_virial)
      {
        io::messages.add("atomic virial not implemented. copy from dihedral angle interaction!",
                         "Dihedral Constraints",
                         io::message::error);
      }

      convergence = false;

      // consider atoms in the next step
      skip_next[it->i] = false;
      skip_next[it->j] = false;
      skip_next[it->k] = false;
      skip_next[it->l] = false;

    } // we have to shake
  }   // constraints

  return 0;
}
#endif //DIHEDRAL_CONSTRAINT_CC
