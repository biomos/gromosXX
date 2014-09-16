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
     * total QM energy
     */
    double qm_total;

    /**
     * total energy of the dipole-dipole interaction (self energy)
     */

    /**
     * energy of the reference state in eds
     */
    double eds_vr;
    
    /**
     * nonbonded energy of the endstates in eds
     */
    std::vector<double> eds_vi;
    
    /**
     * special energy of the endstates in eds
     */
    std::vector<double> eds_vi_special;
           
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
    
    /**
     * reset the energy terms to zero
     */
    void zero(bool potential = true, bool kinetic = true);
    /**
     * resize the arrays for energy groups and temperature baths
     */
    void resize(unsigned int energy_groups, unsigned int multi_baths = 0);
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
    
    static const unsigned int MAX_ENERGY_INDEX = 38;

  private:
    double m_ewarn;
    
  };

} // configuration

#endif
