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
 * @file energy.h
 * storage of the energies.
 */

#ifndef INCLUDED_ENERGY_H
#define INCLUDED_ENERGY_H

namespace configuration
{
  /**
   * @class Energy
   * storage of energies.
   */
  class Energy
  {
  public:
    /**
     * Constructor
     */
    Energy();
    /**
     * Destructor
     */
    ~Energy(){}
    
    /**
     * total energy of the system
     */
    double total;
    /**
     * total kinetic energy
     */
    double kinetic_total;
    /**
     * total potential energy
     */
    double potential_total;

    /**
     * total bond stretching energy
     */
    double bond_total;
    /**
     * total bond angle bending energy
     */
    double angle_total;
    /**
     * total improper dihedral torsional energy
     */
    double improper_total;
    /**
     * total dihedral torsional energy
     */
    double dihedral_total;
    /**
     * total crossdihedral torsional energy
     */
    double crossdihedral_total;
    /**
     * total potential energy of the "bonded" interaction terms
     */
    double bonded_total;
    /**
     * total potential energy of the "non-bonded" interaction terms
     */
    double nonbonded_total;
    /**
     * total energy of the Lennard-Jones interaction
     */
    double lj_total;
    /**
     * total energy of the coulomb reaction-field interaction
     */
    double crf_total;
    /**
     * total energy of the coulomb LS interaction
     */
    double ls_total;
    /**
     * total pairwise LS interaction energy
     */
    double ls_pair_total;
    /**
     * total energy of the coulomb LS real space interaction
     */
    double ls_realspace_total;
    /**
     * total energy of the coulomb LS k space interaction
     */
    double ls_kspace_total;
    /**
     * total electrostatic self energy (LS)
     */
    double ls_self_total;
    /**
     * total electrostatic self energy for constant volume simulations (LS)
     */
    double ls_self_total_nvt;
    /**
     * total surface energy
     */
    double ls_surface_total;
    /**
     * total A term energy
     */
    double ls_a_term_total;
    /**
     * total A term energy for constant volume simulations
     */
    double ls_a_term_total_nvt;
    /**
     * total QM energy
     */
    double qm_total;
    /**
     * shifting energy
     */
    std::vector<std::vector<double> > shift_extra_orig;
    double shift_extra_orig_total;
    std::vector<std::vector<double> > shift_extra_phys;
    double shift_extra_phys_total;
    /**
     * total energy of the "special" interactions
     */
    double special_total;
    /**
     * total energy of the position restraint interaction
     */
    double posrest_total;
    /** 
     * total energy of the distance restraint interaction
     */
    double distanceres_total;
    /**
     *  total energy of the distancefield restraint interaction
     */
    double disfieldres_total;
    /**
     *  total energy of the angle restraint interaction
     */
    double angrest_total;
    /**
     *  total energy of the dihedral restraint interaction
     */
    double dihrest_total;
    /**
     * total energy of the J-value restraint interaction
     */
    double jvalue_total;
    /**
     * total energy of the X-ray restaint interaction
     */
    double xray_total;
    /**
     * total energy of the local elevation interaction
     */
    double leus_total;
    /**
     * total energy of the bsleus interaction
     */
    double bsleus_total;
    /**
     * total energy of the order parameter restraint interaction
     */
    double oparam_total;
    /**
     * total energy of the rdc restraint interaction
     */
    double rdc_total;
    /**
     * symmetry restraints energy
     */
    double symrest_total;
    /**
     * total energy (=0.0) of the (distance) constraint interaction(s).
     */
    double constraints_total;
    
    /**
     * total energy of the dipole-dipole interaction (self energy)
     */
    double self_total;

    /**
     * total energy of sasa interaction
     */
    double sasa_total;
    /**
     * total energy of volume term
     */
    double sasa_volume_total;
    /**
     * energy difference of validation NN model
     */
    double nn_valid;

    /** ANITA
    * total A_lj for each lambda
    */
    std::vector<double> A_lj_total;

    /** ANITA
    * total B_lj for each lambda
    */
    std::vector<double> B_lj_total;

    /** ANITA
    * total A_crf for each lambda
    */
    std::vector<double> A_crf_total;

    /** ANITA
    * total B_crf for each lambda
    */
    std::vector<double> B_crf_total;

    /** ANITA
    * kinetic energy for each lambda
    */
//    double A_kinetic;
//    double B_kinetic;
    std::vector<double> AB_kinetic;

    /** ANITA
    * bond stretching energy for each lambda
    */
    std::vector<double> AB_bond;

    /** ANITA
    * bond angle energy for each lambda
    */
    std::vector<double> AB_angle;

    /** ANITA
    * improper dihedral energy for each lambda
    */
    std::vector<double> AB_improper;

    /** ANITA
    * dihedral energy for each lambda
    */
    double A_dihedral;
    double B_dihedral;
//    std::vector<double> AB_dihedral;

    /** Betty
     * special interaction energies for each lambda
     */
    std::vector<double> AB_disres;
    std::vector<double> AB_angres;
    std::vector<double> AB_dihres;
    std::vector<double> AB_disfld;

    /**
     * total energy of the dipole-dipole interaction (self energy)
     */

    /**
     * energy of the reference state in eds
     */
    double eds_vr;

    /**
     * energy of the gromos reference state in eds with shifting
     */
    double eds_vr_shift_orig;
    double eds_vr_shift_phys;
    
    /**
    *
    */
    double eds_vmix;

    /**
    *
    */
    double eds_emax;

    /**
    *
    */
    double eds_emin;

    /**
    *
    */
    double eds_globmin;

    /**
    *
    */
    double eds_globminfluc;

    /**
     * nonbonded energy of the endstates in eds
     */
    std::vector<double> eds_vi;

    /**
     * shifting energy of V_i states
     */
    std::vector<double> eds_vi_shift_extra_orig;
    std::vector<double> eds_vi_shift_extra_phys;

    /**
    * offset energy of the endstates in eds
    */
    std::vector<double> eds_eir;
    
    /**
     * special energy of the endstates in eds
     */
    std::vector<double> eds_vi_special;

    // ORIOL_GAMD
    /**
     * dihedral energie for each accel group (dihedral + improper + crossdihedral)
     */
    std::vector<double> gamd_dihedral_total;
    /**
     * total potential energie for each accel group without dihedral
     */
    std::vector<double> gamd_potential_total;
    /**
     * boosting potential used per accel group
     */
    std::vector<double> gamd_DV;
    /**
     * Energy threshold for the dihedral acceleration per accel group
     */
    std::vector<double> gamd_ED;
    /**
     * Energy threshold for the potential energy acceleration per accel group
     */
    std::vector<double> gamd_ET;
    /**
     * Force constant for the dihedral acceleration
     */
    std::vector<double> gamd_KD;
    /**
     * Force constant for the potential energy acceleration
     */
    std::vector<double> gamd_KT;
    /**
     * total boosting potential used
     */
    double gamd_DV_total;
           
    // this should be size of bath
    /**
     * kinetic energy term
     */
    std::vector<double> kinetic_energy;
    /**
     * molecular translational kinetic energy term
     */
    std::vector<double> com_kinetic_energy;
    /**
     * molecular internal and rotational kinetic energy term
     */
    std::vector<double> ir_kinetic_energy;

    /**
     * bond stretching energy term
     */
    std::vector<double> bond_energy;
    /**
     * bond angle energy term
     */
    std::vector<double> angle_energy;
    /**
     * improper dihedral energy term
     */
    std::vector<double> improper_energy;
    /**
     * dihedral angle energy term
     */
    std::vector<double> dihedral_energy;
    /**
     * crossdihedral angle energy term
     */
    std::vector<double> crossdihedral_energy;

    /**
     * lennard-jones interaction energy term
     */
    std::vector<std::vector<double> > lj_energy;
    /**
     * coulomb reaction field energy term
     */
    std::vector<std::vector<double> > crf_energy;
    
     /**
     * coulomb LS realspace energy term
     */
    std::vector<std::vector<double> > ls_real_energy;
     /**
     * coulomb LS kspace energy term
     */
    std::vector<std::vector<double> > ls_k_energy;

    /**
     * position restraint energy term
     */
    std::vector<double> posrest_energy;
    /**
     * distance restraint energy term
     */
    std::vector<double> distanceres_energy;
    /**
     * distancefield restraint energy term
     */
    std::vector<double> disfieldres_energy;
    /**
     * angle restraint energy term
     */
    std::vector<double> angrest_energy;
    /**
     * dihedral restraint energy term
     */
    std::vector<double> dihrest_energy;
    /**
     * jvalue restraint energy term
     */
    std::vector<double> jvalue_energy;
    /**
     * RDC restraint energy term
     */
    std::vector<double> rdc_energy;  
    /**
     * (distance) constraints energy term
     * (has to be 0.0 always)
     */
    std::vector<double> constraints_energy;
    /**
     * self energy term (polarisation)
     */
    std::vector<double> self_energy;
    /**
     * entropy estimation
     * dH/dl * dH
     */
    double entropy_term;
    /**
     * energy group names
     */
    std::vector<std::string> group_name;

    /**
     * sasa interaction energy term
     */
    std::vector<double> sasa_energy;
    /**
     * volume (sasa) interaction energy term
     */
    std::vector<double> sasa_volume_energy;

    /** ANITA
     * A_lj energies for [lam][groupi][groupj]    
     */
    std::vector<std::vector<std::vector<double> > > A_lj_energy;

    /** ANITA
     * B_lj energies for [lam][groupi][groupj]    
     */
    std::vector<std::vector<std::vector<double> > > B_lj_energy;

    /** ANITA
     * A_crf energies for [lam][groupi][groupj]    
     */
    std::vector<std::vector<std::vector<double> > > A_crf_energy;

    /** ANITA
     * B_crf energies for [lam][groupi][groupj]    
     */
    std::vector<std::vector<std::vector<double> > > B_crf_energy;

    /**
     * reset the energy terms to zero
     */
    void zero(bool potential = true, bool kinetic = true);
    /**
     * resize the arrays for energy groups and temperature baths
     */
    void resize(unsigned int energy_groups, unsigned int multi_baths = 0,
                unsigned int nr_lambdas = 0); //ANITA
    /**
     * calculate the totals of the individual energy terms.
     */
    int  calculate_totals();

    /**
     * set energy warning
     */
    void ewarn(double ew) { m_ewarn = ew; }
    
    /**
     * gets the energy by an index
     */
    double get_energy_by_index(const unsigned int & index);
    
    static const unsigned int MAX_ENERGY_INDEX = 45; //44

  private:
    double m_ewarn;
    
  };

} // configuration

#endif
