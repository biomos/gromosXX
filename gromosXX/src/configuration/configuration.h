/**
 * @file configuration.h
 * the configuration class
 */

#ifndef INCLUDED_CONFIGURATION_H
#define INCLUDED_CONFIGURATION_H

namespace configuration
{
  /**
   * @class Configuration
   * holds the state information of
   * the simulated system.
   */
  class Configuration
  {
  public:
    
    /**
     * Constructor
     */
    explicit Configuration();
    
    struct state_struct
    {
      /**
       * position
       */
      math::VArray pos;
      /**
       * velocity
       */
      math::VArray vel;
      /**
       * force
       */
      math::VArray force; 

      /**
       * the box.
       */
      math::Box box;
      
      /**
       * virial tensor.
       */
      math::Matrix virial_tensor;

      /**
       * molecular kinetic energy tensor.
       */
      math::Matrix kinetic_energy_tensor;

      /**
       * pressure tensor.
       */
      math::Matrix pressure_tensor;

      /**
       * energy arrays
       */
      configuration::Energy energies;

      /**
       * energy averages.
       */
      configuration::Energy_Average energy_averages;
    
      /**
       * perturbed energy derivatives.
       */
      configuration::Energy perturbed_energy_derivatives;

      /**
       * perturbed energy derivative averages.
       */
      configuration::Energy_Average perturbed_energy_derivative_averages;
      
      /**
       * resize the arrays for s atoms.
       */
      void resize(size_t s);

    };

    /**
     * @struct special_struct
     * special information storage.
     */
    struct special_struct
    {
      /**
       * const rel_mol_com_pos
       * relative molecular com position
       */
      math::VArray rel_mol_com_pos;

      /**
       * the dihedral angle minima for monitoring
       */
      std::vector<double> dihedral_angle_minimum;

      /**
       * flexible constraints velocity.
       */
      std::vector<double> flexible_vel;
      
      /**
       * flexible constraints kinetic energy.
       */
      std::vector<double> flexible_ekin;
      
    };
    
    /**
     * get the current state
     */
    state_struct & current() { return *m_current; }
    /**
     * get the old state
     */
    state_struct & old() { return *m_old; }
    /**
     * get the current state (const)
     */
    state_struct const & current()const { return *m_current; }
    /**
     * get the old state (const)
     */
    state_struct const & old()const { return *m_old; }

    /**
     * exchange the old and the current state
     */
    void exchange_state() { state_struct * dummy = m_current; m_current = m_old; m_old = dummy; }
              
    /**
     * special information accessor.
     */
    special_struct & special() { return m_special; }

    /**
     * special information accessor (const).
     */
    special_struct const & special()const { return m_special; }
    
    /**
     * set the number of atoms in the system.
     */
    void resize(size_t s);
    
    /**
     * boundary type.
     */
    math::boundary_enum boundary_type;
    
  protected:

    /**
     * current state information.
     */
    state_struct * m_current;
    /**
     * state information from last step.
     */
    state_struct * m_old;
    /**
     * state information.
     */
    state_struct m_state1;
    /**
     * state information.
     */
    state_struct m_state2;

    /**
     * special information
     */
    special_struct m_special;
    

  }; // Configuration
  
} // configuration

#endif

  
