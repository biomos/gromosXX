/**
 * @file configuration.h
 * the configuration class
 */

#ifndef INCLUDED_CONFIGURATION_H
#define INCLUDED_CONFIGURATION_H

// headers needed for configuration
#include <gromosXX/configuration/energy.h>
#include <gromosXX/configuration/average.h>

namespace topology
{
  class Topology;
}

namespace simulation
{
  class Parameter;
}

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
    
    /**
     * @struct state_struct
     * holds state information.
     */
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
       * averages.
       */
      configuration::Average averages;
      
      /**
       * perturbed energy derivatives.
       */
      configuration::Energy perturbed_energy_derivatives;

      /**
       * resize the position / velocities and force arrays
       */
      void resize(size_t num_atoms);

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
    
    /**
     * initialise.
     * should be called after a (initial) configuration has been read in or
     * constructed otherwise.
     * - resizes energy / perturbed energy arrays
     * - resizes averages
     * - resizes special data
     * - resizes state_struct current() and old()
     * - if coordinates have been read in it is usually necessary to gather them
     * to keep chargegroups together.
     */
    void initialise(topology::Topology & topo,
		    simulation::Parameter const & param,
		    bool gather = true);

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

  
