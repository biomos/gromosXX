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
    Energy(){}

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
     * total potential energy of the "bonded" interaction terms
     */
    double bonded_total;
    /**
     * total potential energy of the "non-bonded" interaction terms
     */
    double nonbonded_total;
    /**
     * total energy of the lennard-jones interaction
     */
    double lj_total;
    /**
     * total energy of the coulomb reaction-field interaction
     */
    double crf_total;
    /**
     * total energy of the "special" interactions
     */
    double special_total;
    /**
     * total energy of the position restraint interaction
     */
    double posrest_total;
    /** total energy of the distance restraint interaction
     */
    double distrest_total;
    /**
     * total energy of the J-value restraint interaction
     */
    double jvalue_total;
    /**
     * total energy (=0.0) of the (distance) constraint interaction(s).
     */
    double constraints_total;
    /**
     * total energy of an external interaction
     */
    double external_total;
    
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
     * lennard-jones interaction energy term
     */
    std::vector<std::vector<double> > lj_energy;
    /**
     * coulomb reaction field energy term
     */
    std::vector<std::vector<double> > crf_energy;

    /**
     * position restraint energy term
     */
    std::vector<double> posrest_energy;
    /**
     * distance restraint energy term
     */
    std::vector<double> distrest_energy;
    /**
     * jvalue restraint energy term
     */
    std::vector<double> jvalue_energy;
    /**
     * (distance) constraints energy term
     * (has to be 0.0 always)
     */
    std::vector<double> constraints_energy;

    /**
     * energy group names
     */
    std::vector<std::string> group_name;
    
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
    
  };

} // configuration

#endif
