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
 * @file bfactor.h
 * all kind of bfactor fitting routines
 */


#ifndef BFACTOR_H
#define	BFACTOR_H

namespace interaction {
  namespace xray {
    void fit_bfactor(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            const clipper::Cell & cell,
            clipper::Atom_list & atoms,
            clipper::HKL_data<clipper::data32::F_phi> & fphi_calc,
            clipper::HKL_data<clipper::data32::F_phi> & fphi,
            clipper::HKL_data<clipper::data32::F_phi> & fphi_obs,
            clipper::Xmap<clipper::ftype32> & rho_calc,
            clipper::FFTmap_p1 & D_k,
            clipper::Xmap<clipper::ftype32> & d_r);

    void fit_overall_bfactor(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation & sim,
        const clipper::Cell & cell,
        clipper::Atom_list & atoms,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_calc,
        clipper::HKL_data<clipper::data32::F_phi> & fphi,
        clipper::HKL_data<clipper::data32::F_phi> & fphi_obs,
        clipper::Xmap<clipper::ftype32> & rho_calc,
        clipper::FFTmap_p1 & D_k,
        clipper::Xmap<clipper::ftype32> & d_r);
  }
}

#endif	/* BFACTOR_H */

