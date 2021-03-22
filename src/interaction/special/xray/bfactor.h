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

