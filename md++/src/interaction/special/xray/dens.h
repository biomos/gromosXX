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
 * @file dens.h
 * All kind of electron density routines
 */

#ifndef DENS_H
#define	DENS_H

namespace interaction {
  namespace xray {
    /**
     * calculate the electron density from the model structure
     * @param[out] rho_calc the calculated electron density on a map
     * @param[in] atoms the list of the atoms
     */
    void calculate_electron_density(clipper::Xmap<clipper::ftype32> & rho_calc,
            const clipper::Atom_list & atoms);

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
            double & slope, double & intercept);

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
            double to_ang);
  }
}
#endif	/* DENS_H */

