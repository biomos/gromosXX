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

/* 
 * File:   cremer_pople.h
 * Author: dahahn
 *
 * Created on December 14, 2015, 3:15 PM
 */

#ifndef CREMER_POPLE_H
#define	CREMER_POPLE_H

#include <vector>

#include "../gmath/Vec.h"

#ifdef	__cplusplus
extern "C" {
#endif

  
#define Delta(x,y) ((x)==(y)?5./6.:-1./6.)
#define SQRT3_2 .86602540378443864676

namespace util {
  namespace cremerpople {
    // WEIGHTS FOR SUMMATION OVER i
    // weights for R' and R'' definitions, w_1
    const double puckW[] = {0., SQRT3_2, SQRT3_2, 0., -SQRT3_2, -SQRT3_2};
    const double puckV[] = {1., 0.5, -0.5, -1., -0.5, 0.5};
    // weights for puckering coordinate definitions, w_2
    const double puckBIGW[] = {0., SQRT3_2, -SQRT3_2, 0., SQRT3_2, -SQRT3_2};
    const double puckBIGV[] = {1., -0.5, -0.5, 1., -0.5, -0.5};
    const double puckSIGN[] = {1., -1., 1., -1., 1., -1.}; // (-1)^(j-1) term in summation

    void calcCOM(const std::vector<gmath::Vec> &R, gmath::Vec &Rcm);
    void shiftToCOM(std::vector<gmath::Vec> &R, const gmath::Vec &Rcm);

    void calcPlane(const std::vector<gmath::Vec> &R, gmath::Vec &Rp, gmath::Vec &Rdp, gmath::Vec &RpxRdp);

    double calcZeta(const std::vector<gmath::Vec> &R, const gmath::Vec &RpxRdp, std::vector<double> &zeta);
    void calcGradZeta(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &gizj, unsigned int i, unsigned int j);

    void calcGradQ(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &giQ, unsigned int index);

    double calcPhi(const std::vector<double> &zeta);
    void calcGradPhi(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &giphi, unsigned int index);

    double calcTheta(const std::vector<double> &zeta);
    void calcGradTheta(const std::vector<gmath::Vec> &R, const gmath::Vec &Rp, const gmath::Vec &Rdp, const gmath::Vec &RpxRdp, const std::vector<double> &zeta, gmath::Vec &githeta, unsigned int index);

  }
}

#ifdef	__cplusplus
}
#endif

#endif	/* CREMER_POPLE_H */

