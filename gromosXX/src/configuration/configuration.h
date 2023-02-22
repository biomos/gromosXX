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
#include "../util/bs_umbrella.h"

// Additional Clipper Headers
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#endif

namespace topology {
  class Topology;
  struct four_body_term_struct;
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
       * the constraint force
       */
      math::VArray constraint_force;
      /**
       * conjugate gradient search direction
       */
      math::VArray cgrad;
      /**
       * stochastic integrals (SD)
       */
      math::VArray stochastic_integral;
      /**
       * stochastic dynamics random number seed
       */
      std::string stochastic_seed;

      /**
       * sasa of each atom
       */
      std::vector<double> sasa_area;
      /**
       * volume of each buried atom
       */
      std::vector<double> sasa_buriedvol;
      /**
       * switching function for volume term
       */
      std::vector<double> gvol;
      /**
       * derivative of switching function for volume term
       */
      std::vector<double> dgvol;
      /**
       * total sasa of entire solute molecule
       */
      double sasa_tot;
      /**
       * total volume of buried atoms of solute molecule
       */
      double sasa_buriedvol_tot;

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
       * @struct dihangle_trans
       * dihedral angle transition monitoring data
       */
      struct dihangle_trans_struct {
      /**
       * the residues for monitoring: number
       */
      std::vector<int> resid;
      /**
       * the atoms for monitoring
       */
      std::vector<int> i;
      /**
       * the atoms for monitoring
       */
      std::vector<int> j;
      /**
       * the atoms for monitoring
       */
      std::vector<int> k;
      /**
       * the atoms for monitoring
       */
      std::vector<int> l;
      /**
       * the old dihedral angle minima
       */
      std::vector<double> old_minimum;
      /**
       * the dihedral angle minima for monitoring
       */
      std::vector<double> dihedral_angle_minimum;
      } /** dihedral angle transition data */ dihangle_trans;
      ///////////////////////////////////////////////////////////////
      /**
       * the local elevation umbrellas
       */
      std::vector<util::Umbrella> umbrellas;
      //////////////////////////////////////////////////
      /**
       * The Umbrella for the BS&LEUS algorithm. Contains all the B&S 
       * potentials.
       */
      util::BS_Umbrella bs_umbrella;
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


      /**
       * @struct rdc_struct
       * struct that holds the average and current rdc values
       * and other parameters that occur per rdc interaction.
       */
      struct rdc_struct{
        /**
         * RDC averages
         */
        std::vector<double> av;
        /**
         * current RDC values
         */
        std::vector<double> curr;
        /**
         * Magnetic field endpoints
         */
        math::VArray MFpoint;
        /**
         * Magnetic field endpoint velocities
         */
        math::VArray MFpointVel;
        /**
         * Magnetic field endpoint masses
         */
        std::vector<double> MFpointMass;
        /**
         * Alignment tensor components 
         */
        std::vector<double> Tensor;
        /**
         * Tensor component velocities
         */
        std::vector<double> TensorVel;      
        /**
         * Alignment tensor component masses
         */
        std::vector<double> TensorMass;
        /**
         * Shperical harmonics coefficients 
         */
        std::vector<double> clm;
        /**
         * Shperical harmonics coefficient masses
         */
        std::vector<double> clmMass;
        /**
         * Shperical harmonics coefficient velocities
         */
        std::vector<double> clmVel;      
        /**
         * Conversion of the frequency
         */
        double factorFreq;
        /**
         * Conversion of the gyromagnetic ratio
         */
        double factorGyr;
        /**
         * 'kinetic energy'
         */
        double Ekin;
        /**
         * stochastic integrals for the magnetic field vector representation
         */
        math::VArray stochastic_integral_mf;
        /**
         * stochastic integrals for the tensor representation
         */
        std::vector<double> stochastic_integral_t;
        /**
         * stochastic integrals for the spherical harmonic representation
         */
        std::vector<double> stochastic_integral_sh;
      }; // rdc related variables

      /**
       * rdc restraint averages and other related data
       */
      std::vector<rdc_struct> rdc;


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
         * current phase of structure factor
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
        double k_inst, k_free_inst;
        /**
         * the instantaneous R-value
         */
        double R_inst, R_free_inst;
        /**
         * the time-averaged structure factor scaling constant
         */
        double k_avg, k_free_avg;
        /**
         * the time-averaged R-value
         */
        double R_avg, R_free_avg;
        /**
         *overall B factor XROB
         */
        double B_overall;
      } /** xray information */ xray;
      /**
       * @struct xray_bfoc_struct
       * the atomic parameters for X-ray
       */
      struct xray_bfoc_struct {
        /**
         * isotropic B factor
         */
        double b_factor;
        /**
         * gradient of isotropic B factor
         */
        double b_factor_gradient;
        /**
         * occupancy of isotropic B factor
         */
        double occupancy;
      };

      /**
       * the atomic parameters
       */
      std::vector<xray_bfoc_struct> xray_bfoc;
      /**
       * @struct adde_struct
       * struct that holds the average and current potential energies
       * for reweighting when using adiabatic decoupling
       */
      struct adde_struct {
        /**
         * potential energy between heavy particles
         */
        double vhh;
        /**
         * mean potential of light paricles on heay particle
         */
        double evhl;
         /**
         *potential of light paricles on heay particle at time 0
         */
        double vhl0;
      } /** adde information */ adde;

      /**
       * @struct disres_struct
       * holds the distance restraints configuration data
       */
      struct disres_struct {
        /**
         * the distance
         */
        std::vector<double> d;
        /**
         * the energy
         */
        std::vector<double> energy;
        /**
         * the running average
         */
        std::vector<double> av;

      } /** disres informaton */ distanceres;
      disres_struct pertdistanceres;

      /**
       * @struct disfield_struct
       * holds the distance field restraints grid
       */
      struct disfield_struct {
	/**
	 * the number of gridpoints
	 */
	std::vector<int> ngrid;
	/**
	 * the distance for the gridpoints
	 */
	std::vector<double> distance;
	/**
	 * the actual distance
	 */
	double dist;
	/**
	 * energy
	 */
	double energy;
	/**
	 * energy derivative
	 */
	double energy_deriv;
      }  /** disfield information */ distancefield;
  

      /**
       * @struct angres_struct
       * holds the angle restraints configuration data
       */
      struct angres_struct {
        /**
         * the angle
         */
        std::vector<double> d;
        /**
         * the energy
         */
        std::vector<double> energy;

      } /** angres informaton */ angleres;
      angres_struct pertangleres;

      /**
       * @struct dihres_struct
       * holds the dihedral restraints configuration data
       */
      struct dihres_struct {
        /**
         * the dihedral
         */
        std::vector<double> d;
        /**
         * the energy
         */
        std::vector<double> energy;
      } /** dihres informaton */ dihedralres;
      dihres_struct pertdihedralres;
      
	
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

      // ORIOL_GAMD
      struct gamd_struct {
        /**
         *  force storage for the dihedral forces (copy of the dihedral forces)
         */
        std::vector<math::VArray> dihe_force;
        /**
         *  force storage by pairs of acceleration groups (copy of the forces)
         */
        std::vector<math::VArray> total_force;
        /**
         * virial tensor dihedral contribution by charge group
         */
        std::vector<math::Matrix> virial_tensor_dihe;
        /**
         * virial tensor by charge group pairs
         */
        std::vector<math::Matrix> virial_tensor;

      } /** enveloping distribution sampling information */ gamd;
      
      struct nemd_conf_struct {
        /**
         * Variable to accumulate the momentum
         */
        double Px;
        /**
         * Vector to accumulate the data per slab
         */
        std::vector<double> stored_data_per_bin;
        /**
         * Vector to accumulate velocities per atom
         */
        std::vector<double> vel_per_atom;
        /**
         * Vector to accumulate Dvx per atom
         */
        std::vector<double> dvx_per_atom;
        /**
         * Counter
         */
        unsigned int counter;
      } nemd_conf;

      struct xray_conf_struct {
        #ifdef HAVE_CLIPPER
            /**
             * the atoms
             */
            clipper::Atom_list atoms;
            /**
             * the atoms used for structure factor calculation
             */
            clipper::Atom_list atoms_sf;
            /**
             * the calculated electron density
             */
            clipper::Xmap<clipper::ftype32> rho_calc;
            /**
             * the observed electron density
             */
            clipper::Xmap<clipper::ftype32> rho_obs;
            /**
             * the HKLs (reflections)
             */
            clipper::HKL_info hkls;
            /**
             * structure factos built from electron density. Unmodified, unscaled.
             */
            clipper::HKL_data<clipper::data32::F_phi> fphi_calc;
            /**
             * the structure factors
             */
            clipper::HKL_data<clipper::data32::F_phi> fphi;
            /**
             * structure factos built from the observed amplitudes and the
             * calculated phases
             */
            clipper::HKL_data<clipper::data32::F_phi> fphi_obs;
            /**
             * the gradient map
             */
            clipper::FFTmap_p1 D_k;
            /**
             * the map for the gradient convolution
             */
            clipper::Xmap<clipper::ftype32> d_r;
            /**
             * spacegroup for symmetry restraints
             */
            clipper::Spacegroup sym_spacegroup;
        #endif

      } xray_conf;

      /**
       * @struct oparam_struct
       * holds the order parameter restraints configuration data
       */
      struct oparam_struct {
        /**
         * the averaged order parameter
         */
        std::vector<double> S2_avg;
        /**
         * the energy
         */
        std::vector<double> energy;
        /**
         * the running average of Q
         */
        std::vector<math::Matrix> Q_avg;
        /**
         * the running average of D
         */
        std::vector<double> D_avg;
        /**
         * the averaging window of Q
         */
        std::vector<std::list<math::Matrix> > Q_winavg;
        /**
         * the averaging window of D
         */
        std::vector<std::list<double> > D_winavg;
        
      } /** disres informaton */ orderparamres;

      /**
       * lattice shifts
       */
      math::VArray lattice_shifts;

      /**
       * there was a SHAKE failure
       */
      bool shake_failure_occurred;

      unsigned int change_on_slave;

      /**
       * position restraints reference positions
       */
      math::VArray reference_positions;

      /**
       * position restraints bfactors
       */
      math::SArray bfactors;
      
      /**
       * group-wise forces
       */
      std::vector<std::vector<math::VArray > > force_groups;
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
    /**
     * check configuration
     */
    bool check(topology::Topology const & topo, simulation::Simulation & sim);

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
    
    /**
     * check the positions for overlapping atoms.
     * @param topo topology
     * @param error zero if fine, non-zero if not
     */
    template<math::boundary_enum B> 
    void check_positions(topology::Topology const & topo, int & error) const;

    /** 
     * check dihedrals for atoms at the same position
     * @param dihedrals vector with (improper) dihedrals
     * @param error zero if fine, non-zero if not
     */
    template<math::boundary_enum B> 
    void check_dihedrals(std::vector<topology::four_body_term_struct> const & dihedrals, int & error) const;
    
    /**
     * check the positions for excluded atoms which are further away than
     * the inner cut-off. As this exclusions are ignored, it is worth a check.
     * See Gitlab issue #18.
     * Pascal, Feb 2017
     */
    template<math::boundary_enum B> 
    void check_excluded_positions(topology::Topology const & topo, 
                                  simulation::Simulation & sim);


  }; // Configuration
} // configuration

#endif


