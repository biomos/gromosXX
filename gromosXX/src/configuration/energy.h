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
    double total;
    double kinetic_total;
    double potential_total;
    
    double bond_total;
    double angle_total;
    double improper_total;
    double dihedral_total;
    double bonded_total;
    double nonbonded_total;
    double lj_total;
    double crf_total;
    double special_total;
    double posrest_total;

    // this should be size of bath
    std::vector<double> kinetic_energy;
    std::vector<double> com_kinetic_energy;
    std::vector<double> ir_kinetic_energy;

    std::vector<double> bond_energy;
    std::vector<double> angle_energy;
    std::vector<double> improper_energy;
    std::vector<double> dihedral_energy;

    std::vector<std::vector<double> > lj_energy;
    std::vector<std::vector<double> > crf_energy;

    std::vector<double> posrest_energy;
    
    std::vector<std::string> group_name;
    
    void zero(bool potential = true, bool kinetic = true);
    void resize(size_t const energy_groups, size_t const multi_baths = 0);
    void calculate_totals();
    
  };

} // configuration

#endif
