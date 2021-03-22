/**
 * @file sf.h
 * all kind of structure factor routines
 */

#ifndef SF_H
#define	SF_H

namespace interaction {
  namespace xray {
    /**
     * Scale the structure factor amplitudes, store them in the configuration
     * @param[in] topo the topology
     * @param[inout] conf the configuration
     * @param[in] sim the simulation
     * @param[in] fphi_calc the unscaled structure factors
     * @param[out] fphi the scaled structure factors
     * @param[out] fphi_obs the observed structure factors
     */
    void scale_sf(const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            const clipper::HKL_data<clipper::data32::F_phi> & fphi_calc,
            clipper::HKL_data<clipper::data32::F_phi> & fphi,
            clipper::HKL_data<clipper::data32::F_phi> & fphi_obs
            );

    /**
     * calculate the energy for structure factor restraining
     * @param[in] the simulation
     * @param[in] the reflection list for hkl info
     * @param[in] refl the observed relefections
     * @param[in] refl_curr the calculated reflections
     * @param[in] averaging the averaging mode of the restraining
     * @param[in] k_inst the inst. scaling factor
     * @param[in] k_avg the avg. scaling factor
     * @param[out] D_k the difference map for gradients
     * @param[in] force_constant the force constant
     * @param[out] energy the energy obtained
     */
    void calculate_energy_sf(const simulation::Simulation & sim,
            const clipper::HKL_data<clipper::data32::F_phi> & fphi,
            const std::vector<topology::xray_restraint_struct> & refl,
            const std::vector<configuration::Configuration::special_struct::xray_struct> & refl_curr,
            simulation::xrayrest_enum averaging,
            const double k_inst, const double k_avg,
            clipper::FFTmap_p1 & D_k,
            const double force_constant,
            double & energy);

    /**
     * calculates the force from a reciprocal space difference map
     * @param[in] update this is set to true if the SF were updated. It means
     * that the difference map has to be updated too.
     * @param[inout] D_k the reciprocal space difference map
     * @param[out] d_r the real space difference map
     * @param[in] atoms the list containing the atoms
     * @param[out] force the force vector
     * @param[in] to_ang conversion factor for unit length
     */
    void calculate_force_sf(bool update,
            clipper::FFTmap_p1 & D_k,
            clipper::Xmap<clipper::ftype32> & d_r,
            const clipper::Atom_list & atoms,
            math::VArray & force,
            math::SArray & b_deriv,
            double to_ang);
  }
}

#endif	/* SF_H */

