/**
 * @file parameter.h
 * input parameters
 */

#ifndef INCLUDED_PARAMETER_H
#define INCLUDED_PARAMETER_H

namespace simulation
{
  /**
   * @enum constr_enum
   * constraints enumeration.
   */
  enum constr_enum{
    /**
     * no constraints
     */
    constr_off,
    /**
     * shake
     */
    constr_shake,
    /**
     * lincs
     */
    constr_lincs,
    /**
     * flexible shake
     */
    constr_flexshake,
    /**
     * settle
     */
    constr_settle
  };
  /**
   * @enum special_loop_enum
   * special solvent loop
   */
  enum special_loop_enum{
    /**
     * no special loop
     */
    special_loop_off = 0,
    /**
     * spc loop
     */
    special_loop_spc = 1,
    /**
     * special solvent loop (generic)
     */
    special_loop_generic = 2,
    /**
     * spc loop table
     */
    special_loop_spc_table = 3
  };
  
  /**
   * @enum special_loop_solvent_enum
   * holds the solvent used in a special loop
   */
  enum special_loop_solvent_enum {
    /**
     * not solvent specific. use topology
     */
    sls_topo,
    /**
     * special loop spc
     */
    sls_spc
  };
  
  /**
   * @enum special_loop_acceleration_enum
   * holds the acceleration method
   */
  enum special_loop_acceleration_enum {
    /**
     * no acceleration. use standard loops
     */
    sla_off,
    /**
     * generic solvent properties
     */
    sla_generic,
    /**
     * hardcoded parameters
     */
    sla_hardcode,
    /**
     * tabulated forces and energies
     */
    sla_table,
    /**
     * CUDA library acceleration
     */
    sla_cuda
  };
  
  /**
   * @enum jvalue_restr_enum
   * J-value restraints enumeration
   */
  enum jvalue_restr_enum{
    /**
     * no restraints
     */
    jvalue_restr_off = 0,
    /**
     * instantaneous restraints
     */
    jvalue_restr_inst = 1,
    /**
     * instantaneous restraints, weighted
     */
    jvalue_restr_inst_weighted = 2,
    /**
     * time-averaged restraints
     */
    jvalue_restr_av = -1,
    /**
     * time-averaged restraints, weighted
     */
    jvalue_restr_av_weighted = -2,
    /**
     * biquadratic (time averaged & instantaneous) restraints, weighted
     */
    jvalue_restr_biq_weighted = -3
  };

  /**
   * @enum dihedral_restr_enum
   * Dihedral restraints enumeration
   */
  enum dihedral_restr_enum{
    /**
     * no restraints
     */
    dihedral_restr_off = 0,
    /**
     * instantaneous restraints
     */
    dihedral_restr_inst = 1,
    /**
     * instantaneous restraints, weighted
     */
    dihedral_restr_inst_weighted = 2,
    /**
     * dihedral constraints
     */
    dihedral_constr = 3
  };

  /**
   * @enum integrate_enum
   * integration method
   */
  enum integrate_enum{
    /**
     * off
     */
    integrate_off = 0,
    /**
     * leap-frog
     */
    integrate_leap_frog = 1
  };

  /**
   * @enum interaction_func_enum
   * which interaction function to use.
   * if a interaction function is added, it has to be added
   * to the innerloop_template.h file as well,
   * and (of course) to the switch in the Nonbonded_Innerloop
   * and Perturbed_Nonbonded_Innerloop
   */
  enum interaction_func_enum{
    /** lj_crf_function */ lj_crf_func,
    /** lj_ls_function */ lj_ls_func,
    /** pol_lj_crf_function */ pol_lj_crf_func,
    /** cgrain_function */ cgrain_func
  };
  
  /**
   * @enum electrostatic_method_enum
   * which electrostatic method to use
   */
  enum electrostatic_method_enum {
    /** reaction field */ el_reaction_field,
    /** Ewald */ el_ewald,
    /** particle-particle-particle mesh (P3M) */ el_p3m,
    /** smooth particle mesh Ewald */ el_spme
  };
  
  /**
   * @enum ls_a2_method_enum
   * how to calculate the A2 term for lattice sum self energy
   */
  enum ls_a2_method_enum {
    /** a2 and a2_tilde set to zero */ ls_a2_zero,
    /** a2_tilde exact, a2=a2_tilde */ ls_a2t_exact,
    /** a2 numerical, a2_tilde = a2 */ ls_a2_numerical,
    /** a2_tilde exact (ewald or mesh+coords), a2 numerical */ ls_a2t_exact_a2_numerical,
    /** a2_tilde averaged from mesh only, a2 numerical */ ls_a2t_ave_a2_numerical
  };
  
  /**
   * @enum efield_site_enum
   * determines on which site the electric field is calculated
   */
  enum efield_site_enum {
    /**
     * electric field at the atom
     */
    ef_atom = 0, 
    /**
     * electric field at the carge-on-spring
     */
    ef_cos = 1
  };
  
  /**
   * @enum randomgenerator_enum
   * determines which random number generator is used
   */
  enum randomgenerator_enum {
    /**
     * g96 algorithm
     */
    random_g96 = 0, 
    /**
     * GSL library
     */
    random_gsl = 1
  };

  /**
   * @enum posrest_enum
   * position restraints enumeration
   */
  enum posrest_enum{
    /**
     * no restraints
     */
    posrest_off = 0,
    /**
     * normal harmonic position restraints
     */
    posrest_on = 1,
    /**
     * position restraints with force constant weighted by B-factor
     */
    posrest_bfactor = 2,
    /**
     * position constraints
     */
    posrest_const = 3,
  };

  /**
   * @enum xrayrest_enum
   * xray restraints enumeration
   */
  enum xrayrest_enum {
    /**
     * no restraints
     */
    xrayrest_off = 0,
    /**
     *instantaneous xray restraints
     */
    xrayrest_inst = 1,
    /**
     *timeaveraged xray restraints
     */
    xrayrest_avg = 2,
    /**
     *biquadratic instantaneous/timeaveraged xray restraints
     */
    xrayrest_biq = 3
  };

  /**
   * @enum xrayrestmode_enum
   * xray restraints mode enumeration
   */
  enum xrayrestmode_enum {
    /**
     * restrain structure factors
     */
    xrayrest_mode_structure_factor = 0,
    /**
     * restrain electron density
     */
    xrayrest_mode_electron_density = 1
  };

  /**
   * @enum eds_enum
   * eds functional form enumeration
   */
  enum eds_enum {
    /**
     * single s parameter i.e.
     * @f$ V_R = - \left(\beta s \right)^{-1} \ln \sum_i e^{-\beta s \left(V_i-E_i^R\right)} @f$
     */
    single_s = 1,
    /** 
     * pairwise s parameters i.e.
     * @f$ V_R = - \beta ^{-1} \ln \left\{
       \left[
       \sum_{i,j pairs}
       \left(
       e^{-\beta s_{ij} \left(V_i-E_i^R\right)} + e^{-\beta s_{ij} \left(V_j-E_j^R\right)}
       \right)^{1/s_{ij}}
       \right]
       \frac{1}{N-1}
       \right\}
       @f$
     */
    multi_s = 2,
    /**
     * pairwise s parameters using only (N-1) pairs
     *
     * @f$ V_R = - \beta ^{-1} \ln \left\{
       \left[
       \sum_{spec. pairs}
       \left(
       e^{-\beta s_{ij} \left(V_i-E_i^R\right)} + e^{-\beta s_{ij} \left(V_j-E_j^R\right)}
       \right)^{1/s_{ij}}
       \right]
       \frac{N}{(N-1)*2}
       \right\}
       @f$
     */
    pair_s = 3,
  };

  /**
   * @enum interaction_lambda_enum
   * used to refer to interaction with their own lambda dependence
   */
  enum interaction_lambda_enum{
    /**
     * bond interaction
     */
    bond_lambda = 0,
    /**
     * angle interaction
     */
    angle_lambda = 1,
    /**
     * dihedral interaction
     */
    dihedral_lambda = 2,
    /**
     * improper interaction
     */
    improper_lambda = 3,
    /** 
     * Van der Waals interaction
     */
    lj_lambda = 4,
    /**
     * Van der Waals softness value
     */
    lj_softness_lambda = 5,
    /**
     * Coulomb-reaction field interaction
     */
    crf_lambda = 6,
    /**
     * Coulomb-reaction field softness value
     */
    crf_softness_lambda = 7,
    /**
     * position restraint interaction
     */
    disres_lambda = 8,
    /**
     * dihedral restraint interaction
     */
    dihres_lambda = 9,
    /**
     * mass-scaling value
     */
    mass_lambda = 10,
    /**
     * one extra interaction for looping
     */
    last_interaction_lambda=11
  };

  /**
   * @enum localelev_enum
   * local elevation
   */
  enum localelev_enum {
    /**
     * don't use local elevation
     */
    localelev_off = 0,
    /**
     * use local elevation
     */
    localelev_on = 1
  };

  /**
   * @class Parameter
   * input parameters.
   */
  class Parameter
  {
  public:
    Parameter() : title("GromosXX") {}
    
    /**
     * title of the simulation (from the input file)
     */
    std::string title;
    
    /**
     * @struct system_struct
     * system block
     */
    struct system_struct
    {
      /**
       * Constructor.
       * Default values:
       * - npm 1 (1 solute)
       * - nsm 0 (no solvent)
       */
      system_struct() : npm(1), nsm(0) {}
      
      /**
       * Number of protein molecules
       */
      int npm;
      /** 
       * Number of solvent molecules
       */
      int nsm;
    } /** the system paramters */ system;
    
    /**
     * @struct minimise_struct
     * minimise block
     */
    struct minimise_struct
    {
      /**
       * Constructor.
       * Default values:
       * - ntem 0      (no energy minimisation)
       * - ncyc 0      (unused, conjugate gradient not implemented)
       * - dele 0.0001 (minimal energy difference)
       * - dx0  0.1    (initial step size)
       * - dxm  0.5    (maximum step size)
       * - nmin 1      (at least 1 step)
       * - flim 0.0    (no force limitation)
       */
      minimise_struct() : ntem(0), ncyc(0), dele(0.0001),
			  dx0(0.1), dxm(0.5), nmin(1), flim(0.0)
      {}
      /**
       * minimisation method.
       */
      int ntem;
      /**
       * cycle numbers.
       */
      int ncyc;
      /**
       * minimum energy criterion.
       */
      double dele;
      /**
       * starting step size.
       */
      double dx0;
      /**
       * maximum step size.
       */
      double dxm;
      /**
       * minimum number of steps.
       */
      int nmin;
      /**
       * force limit.
       */
      double flim;

    } /** energy minimisation parameters */ minimise;

    /**
     * @struct start_struct
     * start block
     */
    struct start_struct
    {
      /**
       * Constructor.
       * Default values:
       * - shake_pos                false  (no initial SHAKE of positions)
       * - shake_vel                false  (no initial SHAKE of velocities)
       * - remove_com_translation   false  (no initial removal of COM translation)
       * - remove_com_rotation      false  (no initial removal of COM rotation)
       * - generate_velocities      false  (no generation of initial velocities)
       * - read_nosehoover_chains   true   (read them from configuration)
       * - read_nosehoover_barostat true   (read them from configuration)
       * - read_rottrans            true   (read initial setting of positions
       *                                    and orientations for rot-trans constraints)
       * - read_lattice_shifts      true   (read initial lattice shifts)
       * - ig                       0      (random number seed)
       * - tempi                    0.0    (temperature to generate initial velocities)
       */
      start_struct() : shake_pos(false), shake_vel(false), 
                       remove_com_translation(false), remove_com_rotation(false),
		       generate_velocities(false), ig(0), tempi(0.0),
                       read_nosehoover_chains(true), read_nosehoover_barostat(true),
                       read_rottrans(true), read_lattice_shifts(true) {}
      
      /**
       * shake initial positions
       */
      bool shake_pos;
      /**
       * shake initial velocities
       */
      bool shake_vel;
      /**
       * COM translation removal.
       */
      bool remove_com_translation;
      /**
       * COM rotation removal.
       */
      bool remove_com_rotation;
      /**
       * generate velocities.
       */
      bool generate_velocities;
      /**
       * Random number seed
       */
      unsigned int ig;
      /**
       * Initial temperature
       */
      double tempi;
      /**
       * Read Nose-Hoover Chain variables from configuration or reset them
       */
      bool read_nosehoover_chains;
      /**
       * Read Nose-Hoover Chain barostat variables from configuration or reset them
       */
      bool read_nosehoover_barostat;
      /**
       * Read initial setting of positions and orientations for roto-translational
       * constraints from configuration or reset them
       */
      bool read_rottrans;
      /**
       * Read initial lattice shifts from configuration
       */
      bool read_lattice_shifts;
    } /** startup parameters */ start;

    /**
     * @struct step_struct
     * step block
     */
    struct step_struct
    {
      /**
       * Constructor
       * Default values:
       * - number_of_steps  0
       * - t0               0.0 (initial time)
       * - dt               0.0 (time step)
       */
      step_struct() : number_of_steps(0), t0(0.0), dt(0.0) {}
      
      /**
       * Number of steps
       */
      int number_of_steps;
      /**
       * initial time
       */
      double t0;
      /**
       * time step
       */
      double dt;
    } /** (time) step paramters */ step;

    /**
     * @struct boundary_struct
     * BOUNDARY block
     */
    struct boundary_struct
    {
      /**
       * Constructor
       * Default values:
       * - boundary math::vacuum
       * - dof_to_subtract 0
       */
      boundary_struct() : boundary(math::vacuum), dof_to_subtract(0) {}
      
      /**
       * NTB switch
       */
      math::boundary_enum boundary;
      /**
       * Number of degrees of freedom subtracted for temperature
       * NDFMIN switch
       */
      int dof_to_subtract;
    } /** boundary parameters */ boundary;

    /**
     * @struct multibath_struct
     * multibath block
     */
    struct multibath_struct
    {
      /**
       * Constructor
       * Default values:
       * - couple false (no temperature coupling)
       * - found multibath false
       * - found tcouple false
       * - nosehoover 0 (weak coupling)
       */
      multibath_struct() : couple(false), found_multibath(false), found_tcouple(false), nosehoover(0) {}
      
      /**
       * do temperature coupling?
       */
      bool couple;
      /**
       * ready made multibath
       */
      Multibath multibath;
      /**
       * tcouple struct
       * is translated to the multibath before the 
       * configuration / topology is read in.
       */
      struct tcouple_struct
      {
	/**
	 * Constructor
	 * Default values:
	 * - ntt     0    (no temperature coupling)
	 * - temp0 300    (300 K)
	 * - tau     0.1  (coupling time)
	 */
	tcouple_struct(){
	  ntt[0] = ntt[1] = ntt[2] = 0;
	  temp0[0] = temp0[1] = temp0[2] = 300.0;
	  tau[0] = tau[1] = tau[2] = 0.1;
	}
	
	/**
	 * ntt array
	 */
	int ntt[3];
	/**
	 * temp0
	 */
	double temp0[3];
	/**
	 * tau
	 */
	double tau[3];
      } /** TCOUPLE paramter */ tcouple;
      
      /**
       * have multibath
       */
      bool found_multibath;
      /**
       * have tcouple
       */
      bool found_tcouple;
      /**
       * Nose-Hoover?
       *  0 : berendsen
       *  1 : Nose-Hoover
       * >1 : Nose-Hoover-Chains
       */
      int nosehoover;

    } /** temperature coupling parameters */ multibath;
    
    /**
     * @struct pcouple_struct
     * PCOUPLE block
     */
    struct pcouple_struct
    {
      /**
       * Constructor
       * Default values:
       * - scale pcouple_off   (no pressure coupling)
       * - calculate false     (no pressure calculation)
       * - virial no_virial    (no virial calculation)
       * - pres0 diag(0.06102) (1atm in Gromos96 units)
       * - tau 0.5
       * - compressibility 0.000751
       */
      pcouple_struct()
      {
	scale=math::pcouple_off;
	calculate=false;
	virial=math::no_virial;
	pres0 = 0.0;
	pres0(0,0) = pres0(1,1) = pres0(2,2) = 0.06102;
	tau = 0.5;
	compressibility = 0.000751;
      }
      /**
       * calculate pressure?
       */
      bool calculate;
      /**
       * scale pressure?
       */
      math::pressure_scale_enum scale;
      /**
       * virial type
       */
      math::virial_enum virial;
      /**
       * reference pressure
       */
      math::Matrix pres0;
      /**
       * pressure coupling relaxation time
       */
      double tau;
      /**
       * isothermal compressibility
       */
      double compressibility;
    } /** pressure coupling parameters */ pcouple;

    /**
     * @struct centreofmass_struct
     * CENTREOFMASS block
     */
    struct centreofmass_struct
    {
      /**
       * Constructor
       * Default values:
       * - skip_step 0         (number of steps to skip between removal of com motion)
       * - remove_rot false    (remove center of mass rotation)
       * - remove_trans false  (remove center of mass translation)
       */
      centreofmass_struct() : skip_step(0), remove_rot(false), remove_trans(false) {}
      
      /**
       * NSCM parameter
       */
      int skip_step;
      /**
       * remove angular momentum.
       */
      bool remove_rot;
      /**
       * remove translational momentum.
       */
      bool remove_trans;
      
    } /** centre of mass motion removal parameters */ centreofmass;

    /**
     * @struct print_struct
     * PRINT block
     */
    struct print_struct
    {
      /**
       * Constructor
       * Default values:
       * - stepblock 0               (no printing of energies)
       * - centreofmass 0            (no printing of centre of mass information)
       * - monitor_dihedrals false   (do not monitor dihedral angle transitions)
       */
      print_struct() : stepblock(0), centreofmass(0), monitor_dihedrals(false) {}
      
      /**
       * print stepblock
       */
      int stepblock;
      /**
       * print centre of mass
       */
      int centreofmass;
      /**
       * dihedral angle transitions
       */
      bool monitor_dihedrals;
    } /** output parameters */ print;

    /**
     * @struct write_struct
     * WRITE block
     */
    struct write_struct
    {
      /**
       * Constructor
       * Default values:
       * - position 0        (write position trajectory)
       * - velocity 0        (write velocity trajectory)
       * - force    0        (write force trajectory)
       * - energy   0        (write energy trajectory)
       * - free_energy 0     (write energy lambda derivative trajectory)
       * - block_average 0   (write block averaged energy trajectories)
       * - position_solute_only false (write solute and solvent)
       * - velocity_solute_only false (write solute and solvent)
       * - force_solute_only false (write solute and solvent)
       * - energy_index 0    (don't write minimum energy trajectory)
       */
      write_struct() : position(0), velocity(0), force(0), energy(0), free_energy(0), 
		       block_average(0), position_solute_only(false),
                       velocity_solute_only(false), force_solute_only(false),
                       energy_index(0) {}
      
      /**
       * position.
       */
      int position;
      /**
       * velocity.
       */
      int velocity;
      /**
       * force
       */
      int force;
      /**
       * energy
       */
      int energy;
      /**
       * free energy.
       */
      int free_energy;
      /**
       * block averages.
       */
      int block_average;
      /**
       * write solute only for position trajectory
       */
      bool position_solute_only;
      /**
       * write solute only for velocity trajectory
       */
      bool velocity_solute_only;
      /**
       * write solute only for force trajectory
       */
      bool force_solute_only;
      /**
       * index of the energy array taken to write minimum energy
       * trajectory
       */
      int energy_index;
      
    } /** write out paramters (trajectories) */ write;

    /**
     * @struct constraint_struct
     * SHAKE block
     */
    struct constraint_struct
    {
      /**
       * Constructor
       * Default values:
       * - ntc = 1
       */
      constraint_struct() : ntc(1) {}
      
      /**
       * NTC parameter (off=1, hydrogens=2, all=3, specified=4)
       * specified shakes everything in the constraint block in the topology.
       * hydrogens or all add the bonds containing hydrogens or all bonds to
       * the constraint block and shake those.
       */
      int ntc;
      /**
       * @struct constr_param_struct
       * constraint parameter for
       * solute and solvent.
       */
      struct constr_param_struct
      {
	/**
	 * Constructor
	 * Default values:
	 * - algorithm constr_off
	 * - shake_tolerance 0.0001
	 * - lincs_order 4
	 * - flexshake_readin false
	 * - flexshake_mode 0
	 */
	constr_param_struct()
	  : algorithm(constr_off),
	    shake_tolerance(0.0001),
	    lincs_order(4),
	    flexshake_readin(false),
	    flexshake_mode(0)
	{}
	
	/**
	 * constraint algorithm to use.
	 */
	constr_enum algorithm;
	/**
	 * SHAKE tolerance
	 */
	double shake_tolerance;
	/**
	 * LINCS order.
	 */
	int lincs_order;
	/**
	 * read flexible constraint information
	 * from configuration file.
	 */
	bool flexshake_readin;
	/**
	 * mode of flexible constraints.
	 * - 0: kinetic and potential energy, approximate (standard)
	 * - 1: only potential energy, approximate
	 * - 2: kinetic and potential energy, exact (but not conservative)
	 * - 3: only potential energy, exact (conservative ???)
	 */
	int flexshake_mode;
      };
      /**
       * parameter for solute.
       */
      constr_param_struct solute;
      /**
       * parameter for solvent.
       */
      constr_param_struct solvent;
      
    } /** Constraint method parameters */ constraint;

    /**
     * @struct force_struct
     * FORCE block
     */
    struct force_struct
    {
      /**
       * Constructor
       * Default values:
       * - bond 1
       * - angle 1
       * - improper 1
       * - dihedral 1
       * - nonbonded 1
       * - energy_group empty
       * - special_loop -1
       * - interaction function lj_crf_func
       * - external interaction 0
       */
      force_struct() : bond(1), angle(1), improper(1),
		       dihedral(1), nonbonded_vdw(1),
		       nonbonded_crf(1), special_loop(special_loop_off),
		       interaction_function(lj_crf_func),
		       external_interaction(0)
      {}
      
      /**
       * bonds?
       */
      int bond;
      /**
       * angles?
       */
      int angle;
      /**
       * improper?
       */
      int improper;
      /**
       * dihedral?
       */
      int dihedral;
      /**
       * nonbonded van der Waals?
       */
      int nonbonded_vdw;
      /**
       * nonbonded Coulomb and reaction field?
       */
      int nonbonded_crf;
      /**
       * Energy groups
       */
      std::vector<unsigned int> energy_group;
      /**
       * special loops
       */
      int special_loop;
      /**
       * nonbonded interaction function
       */
      interaction_func_enum interaction_function;
      /**
       * add an external interaction
       */
      int external_interaction;
      
    } /** Force(field) parameters */ force;

    /**
     * @struct plist_struct
     * PLIST block
     */
    struct plist_struct
    {
      /**
       * Constructor
       * Default values:
       * - grid 0
       * - skip_step 5
       * - cutoff_short 0.8
       * - cutoff_long 1.4
       * - grid_size 0.4
       * - atomic_cutoff false
       */
      plist_struct() : grid(0), skip_step(5), cutoff_short(0.8),
		       cutoff_long(1.4), grid_size(0.4),
		       atomic_cutoff(false), print(false) {}
      
      /**
       * algorithm.
       */
      int grid;
      /**
       * skip step
       */
      int skip_step;
      /** 
       * short range cutoff
       */
      double cutoff_short;
      /**
       * long range cutoff
       */
      double cutoff_long;
      /**
       * grid size
       */
      double grid_size;
      /**
       * atomic cutoff
       */
      bool atomic_cutoff;
      
      /**
       * print the pairlist
       */
      bool print;
      
    } /** Pairlist method parameters */ pairlist;

    /**
     * @struct nonbonded_struct
     * NONBONDED block
     */
    struct nonbonded_struct
    {
      /**
       * Constructor
       * Default values:
       * - rf_epsilon 66 (spc water)
       * - rf_kappa    0
       * - rf_cutoff 1.4
       * - rf_excluded true (new standard)
       * - epsilon     1
       */
      nonbonded_struct() :      
        method(el_reaction_field),
        lserf(false),
        rf_kappa(0.0),
        rf_cutoff(1.4),
        rf_epsilon(66.0),
        ls_charge_shape(-1),
        ls_charge_shape_width(0.27),
        ls_calculate_a2(ls_a2_zero),
        ls_a2_tolerance(0.0001),
        ls_epsilon(66.0),
        ewald_max_k_x(10),
        ewald_max_k_y(10),
        ewald_max_k_z(10),
        ewald_kspace_cutoff(0.8),
        p3m_grid_points_x(10),
        p3m_grid_points_y(10),
        p3m_grid_points_z(10),
        p3m_charge_assignment(0),
        p3m_finite_differences_operator(0),
        p3m_mesh_alias(0),
        spme_bspline(0),
        accuracy_evaluation(0),
        influence_function_rms_force_error(0.001),
        influence_function_read(false),
        influence_function_write(false),
        lj_correction(false),
        lj_solvent_density(997.0),
        rf_excluded(true),
        epsilon(1.0) {}

      /**
       * method for longrange electrostatics
       */
      electrostatic_method_enum method;
      /**
       * using lattice-sum emulated reaction field (LSERF)?
       */
      bool lserf;
      /**
       * inverse debye scrining length
       */
      double rf_kappa;
      /**
       * reaction field cutoff
       */
      double rf_cutoff;
      /**
       *reaction field permittivity
       */
      double rf_epsilon;
      /**
       * lattice-sum charge shaping function
       */
      int ls_charge_shape;
      /** 
       * lattice-sum charge shaping function width
       */
      double ls_charge_shape_width;
      /**
       * lattice-sum A2 calculation
       */
      ls_a2_method_enum ls_calculate_a2;
      /**
       * lattice-sum A2 tolaerance
       */
      double ls_a2_tolerance;
      /**
       * lattice-sum external permittivity
       */
      double ls_epsilon;
      /**
       * ewald maximum k compontents
       */
      int ewald_max_k_x;
      int ewald_max_k_y;
      int ewald_max_k_z;
      /**
       * ewald k-space cut-off
       */
      double ewald_kspace_cutoff;
      /**
       * p3m number of grid points
       */
      int p3m_grid_points_x;
      int p3m_grid_points_y;
      int p3m_grid_points_z;
      /**
       * p3m charge assignment function
       */
      int p3m_charge_assignment;
      /**
       * p3m order of mesh finite-difference operator
       */
      int p3m_finite_differences_operator;
      /**
       * p3m number of mesh alias vectors consideres
       */
      int p3m_mesh_alias;
      /**
       * oder of SPME B-spline function
       */
      int spme_bspline;
      /**
       * accuracy evaluation
       */
      int accuracy_evaluation;
      /**
       * rms force error
       */
      double influence_function_rms_force_error;
      /**
       * read influence function from file
       */
      bool influence_function_read;
      /**
       * write influence function to conf
       */
      bool influence_function_write;
      /**
       * longrange LJ correction
       */
      bool lj_correction;
      /**
       * solvent density
       */
      double lj_solvent_density;
      /**
       * RF contribution of excluded atoms
       */
      bool rf_excluded;
      /**
       * Epsilon 1 within the cutoff.
       * in GROMOS this is hardcoded to be 1;
       * we do so in In_Parameter
       */
      double epsilon;
    } nonbonded;
    /**
     * @struct posrest_struct
     * POSITONRES block
     */
    struct posrest_struct
    {
      /** 
       * Constructor
       * Default values:
       * - posrest 0 (no position restraints)
       * - read true
       * - force_constant 10000
       */
      posrest_struct() : posrest(posrest_off), read(true), force_constant(1E4),
                         scale_reference_positions(false) {}
      
      /**
       * posrest
       */
      posrest_enum posrest;
      /**
       * read from specification file (otherwise from startup file)
       */
      bool read;
      /**
       * CPOR
       */
      double force_constant;
      /**
       * scale reference positions.
       */
      bool scale_reference_positions;
    } /** Position restraint parameters */ posrest;



    /**
     * @struct xrayrest_struct
     * XRAYRES block (for xrayrest and force_constant)
     *
     * spacegroup, cell-values, resolution, bfactor filled by in_xray.cc
     */
    struct xrayrest_struct {

      /**
       * Constructor
       * Default values:
       * - xrayrest 0 (no xray restraints)
       * - mode 0 (structure factor)
       * - force_constant 10000
       * - write 0
       * - writexmap 0
       * - tau 0
       * - spacegroup "P 1"
       * - resolution 1.0
       * - readavg 0
       */
      xrayrest_struct() : xrayrest(xrayrest_off), mode(xrayrest_mode_structure_factor), force_constant(1E4), tau(0),
      write(0), writedensity(0), writexmap(0), spacegroup("P 1"),
      resolution(1.0), readavg(0) {
      }

      /**
       * xrayrest
       */
      xrayrest_enum xrayrest;
      /**
       * restraining mode
       */
      xrayrestmode_enum mode;
      /**
       * CXR
       */
      double force_constant;
      /**
       * CXTAU
       */
      double tau;
      /**
       * NTWXR
       */
      unsigned int write;
      /**
       * NTWDE
       */
      unsigned int writedensity;
      /**
       * NTXMAP
       */
      unsigned int writexmap;
      /**
       * spacegroup
       */
      std::string spacegroup;
      /**
       * scattering resolution
       */
      double resolution;
      /**
       * decision-boolean for reading averages or not
       */
      bool readavg;
    } /** Xray restraint parameters */ xrayrest;

    /**
     * @struct distanceres_struct
     * DISTANCERES block
     */
    struct distanceres_struct
    {
      /**
       * Constructor
       * Default values:
       * - distanceres 0 (no distance restraints)
       * - K 0
       * - r_linear 0
       * - tau 10
       * - read 0
       */

      distanceres_struct()
	: distanceres(0),
	  K(0),
	  r_linear(0),
	  tau(10),
	  read(0)
      {
      }
      
      /** 
       * distance restraints on/off
       */
      int distanceres;
      
      /**
       * force constant K
       */
      double K;
      
      /**
       * distance where potential gets linear
       */
      double r_linear;
      
      /**
       * memory time for time averaging 
       */
      double tau;
      
      /**
       * read on/off (not supported)
       */
      bool read;
      
    }/** Distance restraints parameters */ distanceres;

    /**
     * @struct dihrest_struct
     * DIHREST block
     */
    struct dihrest_struct
    {
      /**
       * Constructor
       * Default values:
       * - dihrest 0 (no dihedral restraints)
       * - K 0
       */
      dihrest_struct()
	: dihrest(dihedral_restr_off),
	  K(0.0),
	  phi_lin(0.0) {}
      
      /** 
       * dihedral restraints
       * method:
       * - 0: off
       * - 1: uniform K
       * - 2: K * Ki (weight by Ki in dihedral restraint file)
       * - 3: constraints
       */
      dihedral_restr_enum dihrest;
      /**
       * force constant K
       */
      double K;
      /**
       * deviation larger phi_lin leads to linear potential
       */
      double phi_lin;
      
    }/** dihedral restraint parameters */ dihrest;

    /**
     * @struct perturb_struct
     * PERTURB block
     */
    struct perturb_struct
    {
      /**
       * Constructor
       * Default values:
       * - perturbation false
       * - lambda 0
       * - lambda_exponent 1
       * - dlamt 0
       * - scaling false
       * - scaled_only false
       * - soft_lj 0.0
       * - soft_crf 0.0
       */
      perturb_struct() : perturbation(false), read_initial(false), 
                         lambda(0), lambda_exponent(1),
			 dlamt(0), scaling(false), scaled_only(false),
			 soft_vdw(0.0), soft_crf(0.0) {}
      
      /**
       * perturbation?
       */
      bool perturbation;
      /**
       * read initial lambda from configuration
       */
      bool read_initial;
      /**
       * lambda
       */
      double lambda;
      /**
       * lambda exponent
       */
      int lambda_exponent;
      /**
       * change of lambda per time unit (you call it picosecond)
       */
      double dlamt;
      /**
       * scaling?
       */
      bool scaling;
      /**
       * perturb only scaled interactions.
       */
      bool scaled_only;
      /**
       * soft van der Waals interaction
       */
      double soft_vdw;
      /**
       * soft crf interaction
       */
      double soft_crf;
      
    } /** Perturbation parameters */ perturbation;

    /**
     * @struct jvalue_struct
     * j-value restraint parameters.
     */
    struct jvalue_struct
    {
      /**
       * Constructor
       * Default values:
       * - mode restr_off
       * - le 0
       * - tau 0
       * - ngrid 1
       * - K 1.0
       * - delta 0.0
       * - read_av false
       * - write 0
       */
      jvalue_struct()
	: mode(jvalue_restr_off),
	  le(0),
	  tau(0.0),
	  ngrid(1),
	  K(1.0),
	  delta(0.0),
	  read_av(false),
          write(0)
      {
      }
      /**
       * restraining mode.
       */
      jvalue_restr_enum mode;
      /**
       * local elevation restraining
       */
      int le;
      /**
       * coupling time.
       */
      double tau;
      /**
       * number of local elevation grid points
       */
      int ngrid;
      /**
       * force constant
       * (multiplied by individual restraint weighting)
       */
      double K;
      /**
       * no elevation of potential if J is whitin delta to J0
       */
      double delta;
      /**
       * read averages.
       */
      bool read_av;
      /**
       * write averages and LE grid to special trajectory avery n-th step
       */
      unsigned int write;
    } /** jvalue-parameters */ jvalue;
    
    /**
     * @struct pscale_struct
     * periodic scaling parameters.
     */
    struct pscale_struct
    {
      /**
       * Constructor
       * Default values:
       * - jrest false
       * - KDIH 1.0
       * - KJ 1.0
       * - T 1.0
       * - difference 1.0
       * - ratio 1.0
       * - read_data false
       */
      pscale_struct() : jrest(false), KDIH(1.0), KJ(1.0), T(1.0), difference(1.0), ratio(1.0), read_data(false)
      {
      }
      
      /**
       * do J-Value restraints dependent periodic scaling?
       */
      bool jrest;
      /**
       * maximum scaling factor for the dihedral interaction force constant
       */
      double KDIH;
      /**
       * maximum scaling factor for the J-Value interaction force constant
       */
      double KJ;
      /**
       * periodicity of the consine scaling function.
       */
      double T;
      /**
       * difference in J that starts a scaling.
       */
      double difference;
      /**
       * ration between non-scaled and scaled time
       */
      double ratio;
      /**
       * read data for continuous runs
       */
      bool read_data;
      
    } /** pscale parameters */ pscale;

    /**
     * @struct rottrans_struct
     * rotational translational constraints
     */
    struct rottrans_struct
    {
      /**
       * Constructor
       * Default values:
       * - rottrans false
       * - last 0
       */
      rottrans_struct() : rottrans(false), last(0)
      {
      }
      /**
       * apply rotational translational constraints?
       */
      bool rottrans;
      /**
       * last atom to be roto-translationally constrained
       */
      int last;
    } /** rottrans parameters */ rottrans;

    /**
     * @struct replica_struct
     * information to do replica-exchange simulations
     */
    struct replica_struct
    {
      /**
       * Constructor
       * Default values:
       * - num_T 0
       * - num_l 0
       * - temperature <empty>
       * - scale (false)
       * - lambda <empty>
       * - dt <empty>
       * - trials 0
       * - equilibrate 0
       * - slave_runs 0
       * - write 0
       */
      replica_struct() : num_T(0), num_l(0), scale(false), trials(0),
			 equilibrate(0), slave_runs(0), write(0)
      {
      }
      /**
       * number of replicas with different temperature
       */
      int num_T;
      /**
       * number of replicas with different lambdas
       */
      int num_l;
      /**
       * tempeartures
       */
      std::vector<double> temperature;
      /**
       * scale
       */
      bool scale;
      /**
       * lambdas
       */
      std::vector<double> lambda;
      /**
       * time step to use when running at corresponding lambda
       */
      std::vector<double> dt;
      /**
       * trial moves
       */
      int trials;
      /**
       * equilibrate: no switching for the first N trials
       */
      int equilibrate;
      /**
       * runs per slave
       */
      int slave_runs;
      /**
       * write replicas (to master trajectory) every n steps
       */
      int write;
      
    } /** replica exchange parameters */ replica;

    /**
     * @struct cgrain_struct
     * coarse grain potential
     */
    struct cgrain_struct
    {
      /**
       * Constructor
       * Default values:
       * - level (0)
       * - EPS (20.0)
       */
      cgrain_struct()
	: level(0), EPS(20.0)
      {
      }
      /**
       * do coarse graining
       * - 0 atomistic (off)
       * - 1 coarse-graining (on)
       * - 2 multi-graining (mixed)
       *
       * multigraining requires a set of topology,
       * configuration and input parameters for an
       * atomistic and for a coarse-grained simulation.
       * The coarse-grained files are indicated by
       * cg_topo, cg_conf and cg_input.
       * the coarse-grained topology should contain an
       * additional block VIRTUALGRAINS that specifies
       * the virtual atoms from "real" atomistic atoms.
       */
      int level;
      /**
       * EPS for the CG coulomb
       */
      double EPS;
    } /** coarse grain parameters */ cgrain;

    /**
     * @struct multicell_struct
     * multiple unit cell simulations
     */
    struct multicell_struct
    {
      /**
       * Constructor
       * Default values:
       * - multicell false
       * - x 1
       * - y 1
       * - z 1
       */
      multicell_struct() : multicell(false), x(1), y(1), z(1)
      {
      }
      /**
       * do multicell
       */
      bool multicell;
      /**
       * multiples in x direction
       */
      int x;
      /**
       * multiples in y direction
       */
      int y;
      /**
       * multiples in z direction
       */
      int z;
    } /** multicell parameter */ multicell;

    /**
     * @struct analyze_struct
     * re-analyze trajectory
     */
    struct analyze_struct
    {
      /**
       * Constructor
       * Default values:
       * - analyze(false)
       * - copy_pos(false)
       * - trajectory("")
       */
      analyze_struct() : analyze(false), copy_pos(false), trajectory("")
      {
      }
      /** 
       * re-analyze trajectory
       */
      bool analyze;
      /**
       * copy position (to old position)
       */
      bool copy_pos;
      /**
       * trajectory filename
       */
      std::string trajectory;
      
    } /** analyze parameter */ analyze;

    /**
     * @struct integrate_struct
     * select integration routine
     */
    struct integrate_struct
    {
      /**
       * Constructor
       * Default values:
       * - method(integrate_leap_frog)
       */
      integrate_struct() : method(integrate_leap_frog)
      {
      }
      /** 
       * select integration method
       */
      integrate_enum method;

    } /** integration parameter */ integrate;

    /** 
     * @struct lambdas_struct
     * individual lambdas
     */
    struct lambdas_struct
    {
      /**
       * constructor
       * Default values:
       * - individual_lambdas(false)
       * - a(empty)
       * - b(empty)
       * - c(empty)
       * - d(empty)
       * - e(empty)
       */
      lambdas_struct() : individual_lambdas(false), 
			 a(last_interaction_lambda),
			 b(last_interaction_lambda),
			 c(last_interaction_lambda),
			 d(last_interaction_lambda),
			 e(last_interaction_lambda)
      {
      }
      /**
       * use individual values for lambda
       */
      bool individual_lambdas;
      /**
       * polynomial coefficient a
       */
      std::vector< std::vector< std::vector < double > > > a;
       /**
       * polynomial coefficient b
       */
      std::vector< std::vector< std::vector < double > > > b;
       /**
       * polynomial coefficient c
       */
      std::vector< std::vector< std::vector < double > > > c;
       /**
       * polynomial coefficient d
       */
      std::vector< std::vector< std::vector < double > > > d;
        /**
       * polynomial coefficient e
       */
      std::vector< std::vector< std::vector < double > > > e;

    } /** lambdas struct */ lambdas;
 
    struct stochastic_struct
    {
      /**
       * Constructor
       * Default values:
       * - SD(0)
       * - NTFR(0)
       * - NSFR(0)
       * - NBREF(0)
       * - RCUTF(0.0)
       * - CFRIC(0.0)
       * - TEMP(0.0)
       * - generate_integral(false)
       */
      stochastic_struct() : sd(0), ntfr(0), nsfr(0), nbref(0), rcutf(0.0),
			    cfric(0.0), temp(0.0), generate_integral(false)
      {
      }
      /**
       * do stochastic dynamics
       */
      int sd;
      /**
       * calculate friction coefficients?
       * - 0: set gamma to 0.0
       * - 1: set gamma to CFRIC
       * - 2: set gamma to CFRIC * gamma(0), gamma(0) read from file (not implemented)
       * - 3: calculate gamma
       */
      int ntfr;
      /**
       * recalculate gamma every nsfr steps
       */
      int nsfr;
      /**
       * number of neighbour atoms within RCUTF distance to be considered buried
       */
      int nbref;
      /**
       * inter atom distance to consider when calculating gamma
       */
      double rcutf;
      /**
       * global weighting for gamma
       */
      double cfric;
      /**
       * temperature of the stochastic bath
       */
      double temp;
      /**
       * initially generate stochastic integral
       */
      bool generate_integral;
      
    } /** stochastic dynamics */ stochastic;

    /**
     * @struct ewarn_struct
     * warn about too high energy terms
     */
    struct ewarn_struct
    {
      /**
       * Constructor
       * default values
       * - ewarn 1e99
       */
      ewarn_struct() : limit(1E99)
      {
      }
      /**
       * maximum allowed energy (per term)
       */
      double limit;
    } /** ewarn */ ewarn;

    /**
     * @struct multistep_struct
     * multiple time stepping
     */
    struct multistep_struct
    {
      /**
       * constructor
       */
      multistep_struct() : steps(1), boost(false)
      {
      }
      /**
       * number of steps to boost the
       * nonbonded (and external; multigraining ;-)
       * terms
       */
      int steps;
      /**
       * use boost method (impulse)
       */
      bool boost;
      
    } /** multistep */ multistep;

    /**
     * @struct montecarlo_struct
     * monte-carlo simulation
     */
    struct montecarlo_struct
    {
      /**
       * Constructor
       */
      montecarlo_struct() : mc(0), steps(0), dlambda(0)
      {
      }
      /**
       * chemical monte-carlo
       */
      int mc;
      /**
       * number of md steps between mc trials
       */
      int steps;
      /**
       * value of dlambda in MC move 
       */
      double dlambda;
    } /** chemical monte-carlo */ montecarlo;

    
    /**
     * @struct ramd_struct
     * random acceleration simulation
     */
    struct ramd_struct
    {
      /**
       * Constructor
       */
      ramd_struct() : fc(0.0), steps(0), r_min(0.0), every(0), 
		      do_ta(false), tau(0.0), ta_min(0.0)
      {
      }
      
      /**
       * force constant
       */
      double fc;
      /**
       * number of md steps between change of direction
       */
      int steps;
      /**
       * minimum distance to travel 
       */
      double r_min;
      /**
       * print the ramd force and com position every steps
       */
      int every;
      /**
       * logical to turn time-averaged thingies on
       */
      bool do_ta;
      /**
       * averaging time for time averaged distances
       */
      double tau;
      /**
       * minimum distance to travel on the time averaged scale
       */
      double ta_min;
      /**
       * atoms to apply force to
       */
      std::set<unsigned int> atom;
      
    } /** ramd */ ramd;
    
    /**
     * @struct polarise_struct
     * polarisation simulation
     */
    struct polarise_struct {
      /**
       * Constructor
       * Default values:
       * - cos, no charge-on-spring polarisation
       * - minfield: 2.5
       * - damp, no damping of polarisability
       * - output cospositions every 'write'th block to special trajectory
       */
      polarise_struct() : cos(false), minfield(2.5), efield_site(ef_atom),
                          damp(false), write(0)
      {}
      /**
       * use charge-on-spring polarisation
       */
      bool cos;
      /** 
       * minfield
       */
      double minfield;
      /**
       * site to calculate the electric field
       */
      efield_site_enum efield_site;
      /**
       * use damping
       */
      bool damp;
      /**
       * write cos positions every write-th step to special trajectory
       */
      int write;
    } /** polarise */ polarise;
    
    /**
     * @struct rng_struct
     * random number generator settings
     */
    struct rng_struct {
      /**
       * Constructor
       * Default values:
       * - g96
       */
      rng_struct() : rng(random_g96), gsl_rng(-1) {}
      /**
       * random number generator
       */
      randomgenerator_enum rng;
      /** 
       * GSL random number generator
       * use the rng_gsl contrib program to find out which values are supported.
       */
      int gsl_rng;
    } /** random number generator */ rng;    
    
    /**
     * @struct eds_struct
     * parameters for enveloping distribution sampling (eds)
     */
    struct eds_struct{
      /** 
       * Constructor:
       * Default values:
       * - eds: no eds sampling
       * - form: single_s
       */
      eds_struct() : eds(false), form(single_s), numstates(0) {}
      /**
       * do enveloping distribution sampling using the Hamiltonian:
       * 
       */
      bool eds;
      /**
       * functional form of eds Hamiltonian
       */
      eds_enum form;
      /**
       * number of eds states
       */
      unsigned int numstates;
      /**
       * smoothness parameter(s) @f$ s@f$ of @f$ s_{ij}@f$ used in reference state Hamiltonian.
       */
      std::vector<double> s;
      /**
       * vector of indices of specified pairs (for form = pair_s)
       */
      struct state_pair{
        int i,j;
      };
      std::vector<state_pair> pairs;
      /**
       * energy offsets @f$E_i^R@f$ of states
       */
      std::vector<double> eir;
    } /** enveloping distribution sampling*/ eds;
    
    /**
     * @struct innerloop_struct
     * Constructor:
     * Default values:
     * - solvent: from topology
     * - method: off
     */
    struct innerloop_struct {
      /**
       * constructor
       */
      innerloop_struct() : solvent(sls_topo), method(sla_off), cuda_device(0) {}
      /**
       * the solvent
       */
      special_loop_solvent_enum solvent;
      /**
       * the acceleration method
       */
      special_loop_acceleration_enum method;
      /**
       * the GPU CUDA device number to use
       */
      unsigned int cuda_device;
    } /** special inner loops */ innerloop;

    /**
     * @struct localelev_struct
     * Constructor:
     * Default values:
     * - localelev: off
     * - read: false
     * - umbrellas: empty
     */
    struct localelev_struct {
      /**
       * constructor
       */
      localelev_struct() : localelev(localelev_off), read(false) {}
      /**
       * use the method or not
       */
      localelev_enum localelev;
      /**
       * read from external file
       */
      bool read;
      /**
       * ids of the umbrella potentials to apply
       * true if building up
       */
      std::map<int, bool> umbrellas;
    } localelev;
  };
}

#endif
