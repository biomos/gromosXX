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
     */
    struct index_struct
    {
      index_struct() : k(0), l(0), m(0) {};
      int k;
      int l;
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
    math::VArray &force();
    /**
     * old force
     */
    math::VArray const &old_force()const;
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
     * check state
     */
    int check_state()const;
    
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

  }; // System
  
} // simulation

// template and inline definitions
#include "system.tcc"

#endif

  
