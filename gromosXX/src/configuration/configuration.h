/**
 * @file configuration.h
 * the configuration class
 */

#ifndef INCLUDED_CONFIGURATION_H
#define INCLUDED_CONFIGURATION_H

// headers needed for configuration
#include "energy.h"
#include "average.h"
#include "mesh.h"
#include "influence_function.h"
#include "kspace.h"
#include "../util/umbrella.h"

namespace topology {
  class Topology;
}

namespace simulation {
  class Parameter;
}

namespace configuration {

  /**
   * @class Configuration
   * holds the state information of
   * the simulated system.
   */
  class Configuration {
  public:

    /**
     * Constructor
     */
    Configuration();

    /**
     * copy constructor
     */
    Configuration(Configuration const & conf);
    /**
     * assignment
     */
    Configuration & operator=(Configuration const & conf);

    /**
     * @struct state_struct
     * holds state information.
     */
    struct state_struct {
      /**
       * position
       */
      math::VArray pos;
      /** 
       * charge-on-spring (distance vector between cos and real atom)
       */
      math::VArray posV;
      /**
       * velocity
       */
      math::VArray vel;
      /**
       * force
       */
      math::VArray force;

      /**
       * sasa and volume forces; only needed for testing force
       */
      //math::VArray fsasa;
      /**
       * sasa and volume forces; only needed for testing force
       */
      //math::VArray fvolume;
      
      /**
       * the constraint force
       */
      math::VArray constraint_force;
      /**
       * stochastic integrals (SD)
       */
      math::VArray stochastic_integral;
      /**
       * stochastic dynamics random number seed
       */
      std::string stochastic_seed;

      /**
       * sasa of atom i
       */
      std::vector<double> sasa_area;
      /**
       * the volume interior of a molecule
       */
      std::vector<double> sasa_vol;
      /**
       * total sasa
       */
      double sasa_tot;
      /**
       * total volume
       */
      double sasavol_tot;

      /**
       * the box.
       */
      math::Box box;

      /**
       * the Euler angles 
       */
      double phi; //yaw   (z   axis)
      double theta; //pitch (y'  axis)
      double psi; //roll  (x'' axis)

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
      void resize(unsigned int num_atoms);
    };

    /**
     * @struct special_struct
     * special information storage.
     */
    struct special_struct {
      /**
       * the dihedral angle minima for monitoring
       */
      std::vector<double> dihedral_angle_minimum;

      /**
       * the local elevation umbrellas
       */
      std::vector<util::Umbrella> umbrellas;
      //////////////////////////////////////////////////

      /**
       * @struct flexible_constraint
       * flexible constraint data
       */
      struct flexible_constraint_struct {
        /**
         * flexible constraints velocity.
         */
        std::vector<double> flexible_vel;
        /**
         * flexible constraints kinetic energy.
         */
        std::vector<double> flexible_ekin;
        /**
         * flexible constraint lengths
         */
        std::vector<double> flex_len;

      } /** flexible constraint data */ flexible_constraint;

      //////////////////////////////////////////////////
      /**
       * j value average
       */
      std::vector<double> jvalue_av;
      /**
       * current j value
       */
      std::vector<double> jvalue_curr;
      /**
       * local elevation counter
       */
      std::vector<std::vector<double> > jvalue_epsilon;
      //////////////////////////////////////////////////

      /**
       * @struct xray_struct
       * struct that holds the average and current structure factors and their
       * phases.
       */
      struct xray_struct {
        /**
         * structure factor time averaged
         */
        double sf_av;
        /**
         * current structure factor
         */
        double sf_curr;
        /**
         * average phase of structure factor
         */
        double phase_av;
        /**
         * current pahse of structure factor
         */
        double phase_curr;
      };
      
      /**
       * xray restraint averages
       */
      std::vector<xray_struct> xray_rest;

            /**
       * @struct xray_rvalue_struct
       * holds some of the X-ray restraints configuration data
       */
      struct xray_rvalue_struct {
        /**
         * the instantaneous structure factor scaling constant
         */
        double k_inst;
        /**
         * the instantaneous R-value
         */
        double R_inst;
        /**
         * the time-averaged structure factor scaling constant
         */
        double k_avg;
        /**
         * the time-averaged R-value
         */
        double R_avg;
      } /** xray informaton */ xray;

      /**
       * distance restraint average
       */
      std::vector<double> distanceres_av;

      /**
       * @struct pscale_struct
       * stores periodic scaling information
       */
      struct pscale_struct {
        /**
         * maps J-value restraints to dihedral angles.
         */
        std::vector<int> JtoDihedral;
        /**
         * stores original dihedral angle potential force constants
         */
        std::vector<double> KDIH;
        /**
         * stores original J-value restraint force constants.
         */
        std::vector<double> KJ;
        /**
         * scale time (no scaling if 0.0)
         */
        std::vector<double> t;
        /**
         * scaling right now?
         */
        std::vector<int> scaling;

      } /** periodic scaling information */ pscale;

      //////////////////////////////////////////////////

      /**
       * roto-translational constraints
       * @struct rottrans_constr_struct
       */
      struct rottrans_constr_struct {
        /**
         * inverse theta: translational part (diagonal)
         */
        math::Matrix theta_inv_trans;
        /**
         * inverse theta: rotational part
         */
        math::Matrix theta_inv_rot;
        /**
         * reference positions
         */
        math::VArray pos;
      } /** roto-translational constraints information */ rottrans_constr;

      struct ramd_struct {
        /**
         * directionality of the force
         */
        math::Vec force_direction;
        /**
         * old center of mass of RAMD atoms
         */
        math::Vec old_com;
        /**
         * total mass of the RAMD atoms
         */
        double total_mass;
        /**
         * time averaged distance that we have travelled
         */
        double ta_average;

      } /** ramd information */ ramd;

      struct eds_struct {
        /**
         * (longrange) force storage (perturbed part of multiple states)
         */
        std::vector<math::VArray> force_endstates;
        /**
         * virial tensor (perturbed part of multiple states)
         */
        std::vector<math::Matrix> virial_tensor_endstates;

      } /** enveloping distribution sampling information */ eds;

      /**
       * lattice shifts
       */
      math::VArray lattice_shifts;

      /**
       * there was a SHAKE failure
       */
      bool shake_failure_occurred;

      /**
       * position restraints reference positions
       */
      math::VArray reference_positions;

      /**
       * position restraints bfactors
       */
      math::SArray bfactors;
    }; // special

    /**
     * @struct lattice_sum_struct
     * lattice sum information
     */
    struct lattice_sum_struct {

      /**
       * constructor
       */
      lattice_sum_struct() : charge_density(NULL), potential(NULL),
      a2_tilde(0.0), a2_tilde_derivative(0.0), squared_charge(NULL) {
      }
      /**
       * destructor
       *
      ~lattice_sum_struct() {
        if (charge_density!=NULL) delete charge_density;
        if (potential!=NULL) delete potential;
      }*/
      /**
       * a vector holding the k-space elements
       * (vectors, absolute values, fourier coefficients)
       */
      std::vector<KSpace_Element> kspace;

      /**
       * influence function (real)
       */
      Influence_Function influence_function;
      /**
       * charge density
       */
      Mesh* charge_density;
      /**
       * electrostatic potential
       */
      Mesh* potential;

      /**
       * electric field
       */
      struct electric_field_mesh {

        electric_field_mesh() : x(NULL), y(NULL), z(NULL) {
        }
        /*
        ~electric_field_mesh() {
          if (x!=NULL) delete x;
          if (y!=NULL) delete y;
          if (z!=NULL) delete z;
        }*/
        Mesh* x;
        Mesh* y;
        Mesh* z;
      } electric_field;

      /**
       * methodology dependent A2 term
       */
      double a2_tilde;

      /**
       * derivative of the methodology dependent A2 term
       */
      math::SymmetricMatrix a2_tilde_derivative;

      /**
       * indices of domain decomposition
       */
      std::vector<int> domain;
      /**
       * the mesh with the squared charges
       */
      Mesh* squared_charge;

      /**
       * init lattice sum
       */
      void init(topology::Topology const & topo, simulation::Simulation & sim);
    };
    //////////////////////////////////////////////////////////////////////
    // accessors
    //////////////////////////////////////////////////////////////////////

    /**
     * get the current state
     */
    state_struct & current() {
      return *m_current;
    }

    /**
     * get the old state
     */
    state_struct & old() {
      return *m_old;
    }

    /**
     * get the current state (const)
     */
    state_struct const & current()const {
      return *m_current;
    }

    /**
     * get the old state (const)
     */
    state_struct const & old()const {
      return *m_old;
    }

    /**
     * exchange the old and the current state
     */
    void exchange_state() {
      state_struct * dummy = m_current;
      m_current = m_old;
      m_old = dummy;
    }

    /**
     * special information accessor.
     */
    special_struct & special() {
      return m_special;
    }

    /**
     * special information accessor (const).
     */
    special_struct const & special()const {
      return m_special;
    }

    /**
     * lattice sum information accessor.
     */
    lattice_sum_struct & lattice_sum() {
      return m_lattice_sum;
    }

    /**
     * lattice sum information accessor (const).
     */
    lattice_sum_struct const & lattice_sum()const {
      return m_lattice_sum;
    }
    /**
     * set the number of atoms in the system.
     */
    void resize(unsigned int s);

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
    void init(topology::Topology const & topo,
            simulation::Parameter & param,
            bool gather = true);

    //////////////////////////////////////////////////////////////////////
    // data
    //////////////////////////////////////////////////////////////////////

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

    /**
     * lattice sum information
     */
    lattice_sum_struct m_lattice_sum;


  }; // Configuration

} // configuration

#endif


