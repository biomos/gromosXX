/**
 * @file system.h
 * the system class
 */

#ifndef INCLUDED_SYSTEM_H
#define INCLUDED_SYSTEM_H

namespace simulation
{
  /**
   * @class System
   * holds the state information of
   * the simulated system.
   * @sa simulation::simulation
   * @sa simulation::topology
   */
  template<math::boundary_enum b>
  class System
  {
  public:
    /**
     * @struct index_struct
     * struct to hold the box indices for
     * the atoms.
     * @TODO would allow for Phil's method to gather
     * -- not implemented --
     */
    struct index_struct
    {
      /**
       * Constructor.
       */
      index_struct() : k(0), l(0), m(0) {};
      /**
       * box index k.
       */
      int k;
      /**
       * box index l.
       */
      int l;
      /**
       * box index m.
       */
      int m;
    };
    
    /**
     * Constructor
     */
    explicit System();
    /**
     * set the number of atoms in the system.
     */
    void resize(size_t s);
    /**
     * const position accessor
     */
    math::VArray const &pos()const;
    /**
     * position accessor
     */
    math::VArray &pos();
    /**
     * old position
     */
    math::VArray const &old_pos()const;
    /**
     * const velocity accessor
     */
    math::VArray const &vel()const;
    /**
     * velocity accessor
     */
    math::VArray &vel();
    /**
     * const old velocity
     */
    math::VArray const &old_vel()const;
    /**
     * old velocity
     */
    math::VArray & old_vel();
    /**
     * force accessor
     */
    math::VArray const &force()const; 
    /**
     * force accessor
     */
    math::VArray &force();
    /**
     * old force
     */
    math::VArray const &old_force()const;

    /**
     * const constraint force
     */
    math::VArray const &constraint_force()const;

    /**
     * constraint force
     */
    math::VArray &constraint_force();
    
    /**
     * exchange positions
     */
    void exchange_pos();
    /**
     * exchange velocities
     */
    void exchange_vel();
    /**
     * exchage forces
     */
    void exchange_force();
    /**
     * const periodicity accessor
     */
    math::Periodicity<b> const & periodicity()const;
    /**
     * periodicity accessor
     */
    math::Periodicity<b> & periodicity();
    /**
     * box index accessor
     */
    std::vector<index_struct> & box_indices();
    /**
     * box index i accessor
     */
    index_struct & box_index(size_t i);
    
    /**
     * const rel_mol_com_pos
     * relative molecular com position
     */
    math::VArray const & rel_mol_com_pos()const;
    
    /**
     * rel_mol_com_pos
     * relative molecular com position
     */
    math::VArray & rel_mol_com_pos();
     /**
     * const virial accessor.
     */
    math::Matrix const & virial()const;

    /**
     * virial accessor.
     */
    math::Matrix & virial();

    /**
     * const molecular kinetic energy accessor.
     */
    math::Matrix const & molecular_kinetic_energy()const;
    
    /**
     * molecular kinetic energy accessor.
     */
    math::Matrix & molecular_kinetic_energy();

    /**
     * pressure accessor.
     */
    math::Matrix & pressure();
    
    /**
     * const pressure accessor.
     */
    math::Matrix const & pressure()const;

    /**
     * const energy arrays
     */
    simulation::Energy const & energies()const;

    /**
     * energy arrays
     */
    simulation::Energy & energies();
    
    /**
     * const energy averages.
     */
    simulation::Energy_Average const & energy_averages()const;
    
    /**
     * energy averages.
     */
    simulation::Energy_Average & energy_averages();

    /**
     * lambda derivative of the energies.
     */
    simulation::Energy & lambda_energies();

    /**
     * lambda derivative averages.
     */
    simulation::Energy_Average & lambda_derivative_averages();
    /**
     * const lambda derivative averages.
     */
    simulation::Energy_Average const & lambda_derivative_averages()const;

    /**
     * generate initial velocities.
     */
    void generate_velocities(double const temp, math::SArray const &mass,
			     unsigned int const seed);

    /**
     * check state
     */
    int check_state()const;

    /**
     * calculate the center of mass and the
     * translational kinetic energy of a group of
     * atoms.
     * @param start begin of a group of atoms.
     * @param end of a group of atoms.
     * @param mass the masses of (all) atoms.
     * @param com_pos returns the center of mass.
     * @param com_e_kin returns the tranlational kinetic energy tensor.
     * 
     * @TODO the gathering of the molecule is hardcoded in here.
     * Maybe this should be changed to a generic implementation.
     * Gathering is done in respect to the previous atom. An idea would
     * be to gather as default with respect to the previous atom but
     * letting the user override this (GATHER block).
     * This does not yield the same answer as Phils approach for all cases
     * but maybe for the practical ones???
     */
    void center_of_mass(Atom_Iterator start, Atom_Iterator end,
			math::SArray const &mass, 
			math::Vec &com_pos, math::Matrix &com_e_kin);
    
    void molecular_translational_ekin(Atom_Iterator start, Atom_Iterator end,
				      math::SArray const &mass, 
				      math::Vec &mol_v, double &com_e_kin,
				      double &e_kin,
				      int mean = 0);
    
  protected:
    /**
     * position 1
     */
    math::VArray m_position1;
    /**
     * position 2
     */
    math::VArray m_position2;
    /**
     * velocity 1
     */
    math::VArray m_velocity1;
    /**
     * velocity 2
     */
    math::VArray m_velocity2;
    /**
     * force 1
     */
    math::VArray m_force1;
    /**
     * force 2
     */
    math::VArray m_force2;
    /**
     * position
     */
    math::VArray *m_pos;
    /**
     * old position
     */
    math::VArray *m_old_pos;
    /**
     * velocity
     */
    math::VArray *m_vel;
    /**
     * old velocity
     */
    math::VArray *m_old_vel;
    /**
     * force
     */
    math::VArray *m_force;
    /**
     * old force
     */
    math::VArray *m_old_force;
    /**
     * constraint force
     */
    math::VArray m_constraint_force;
    
    /**
     * the box.
     */
    math::Box m_box;
    /**
     * the periodicity.
     * hard coded to any periodicity...
     */
    math::Periodicity<math::any> m_periodicity;
    /**
     * the box indices for every atom.
     */
    std::vector<index_struct> m_box_index;

    /*
     * rel_mol_com_pos
     * relative molecular com position
     */
    math::VArray m_rel_mol_com_pos;
    
    /**
     * the virial.
     */
    math::Matrix m_virial;
    
    /**
     * the molecular kinetic energy.
     */
    math::Matrix m_molecular_kinetic_energy;

    /**
     * the pressure.
     */
    math::Matrix m_pressure;
    
    /**
     * the energies.
     */
    Energy m_energy;

    /**
     * the average energies and fluctuations.
     */
    Energy_Average m_energy_average;
    
    /**
     * the lambda derivative of the energy.
     */
    Energy m_lambda_energy;

    /**
     * the lambda derivative averages.
     */
    Energy_Average m_lambda_derivative_average;

  }; // System
  
} // simulation

// template and inline definitions
#include "system.tcc"

#endif

  
