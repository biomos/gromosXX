/**
 * @file parameter.h
 * input parameters
 */

#ifdef XXMPI
    #include <mpi.h>
#endif

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
    constr_settle,
    /**
     * m_shake
     */
    constr_m_shake,
    /**
     * gpu_settle
     */
    constr_gpu_settle,
    /**
     * gpu_shake
     */
    constr_gpu_shake

  };
   /**
    * @enum replica_exchange_interruptor_enum
    * mode switcher
    * of replica exchange
    */
  enum replica_exchange_interruptor_enum{
      /**
       * no replica exchange
       */
      replica_exchange_off = 0,
      /**
       * replica exchange with force constant
       */
      replica_exchange_force = 1,
      /**
       * replica exchange with resolution
       */
      replica_exchange_resolution = 2
  };
     /**
    * @enum rep_ex_energy_interruptor_enum
    * mode switcher
    * of the energy
    */
  enum rep_ex_energy_interruptor_enum{
      /**
       * total energy
       */
      energy_tot = 0,
      /**
       * physical energy
       */
      energy_phys = 1,
      /**
       * special energy
       */
      energy_special = 2
  };
   /**
    * @enum B_overall_enum
    * interruptor of the overall B factor
    */
  enum B_overall_enum{
      /**
       * B_overall on
       */
      B_overall_off = 0,
      /**
       * B_overall off
       */
      B_overall_on = 1,
  };
  /**
   * @enum dfunct_enum
   * dfunct enum, useful for sampling reaction transition states
   */
  enum dfunct_enum {
    /**
     * dfunct off
     */
    dfunct_off, 
    /**
     * restrain substitution type geometry
     */
    dfunct_substitution,
    /**
     * restrain cycloaddition type geometry
     */
    dfunct_cycloaddition
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

#ifdef HAVE_HOOMD
  /**
   * @enum hoomd
   * HOOMD code settings enumeration
   */
  enum hoomd {
    unknown = -1,
    cpu = 0,  // One CPU
	gpus = 1  // All GPUs
  };
#endif

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
    jvalue_restr_biq_weighted = -3,
    /**
     * time-averaged restraints with no force scaling
     */
    jvalue_restr_tar_no_scaling = 0,
    /**
     * time-averaged restraints with force scaling
     */
    jvalue_restr_tar_scaling = 1,
    /**
     * biquadratic with equal weights of the two terms
     */
    jvalue_restr_biq_equal_weight = 0,
    /**
     * biquadratic with exponential weight of the average term
     */
    jvalue_restr_biq_exp_weight = 1,
    /**
     * biquadratic with zero weight of the average term
     */
    jvalue_restr_biq_zero_weight = 2
  };

    /**
   * @enum oparam_restr_enum
   * order parameter restraints enumeration
   */
  enum oparam_restr_enum {
    /**
     * no restraints
     */
    oparam_restr_off = 0,
    /**
     * time-averaged restraints
     */
    oparam_restr_av = -1,
    /**
     * time-averaged restraints, weighted
     */
    oparam_restr_av_weighted = -2,
    /**
     * time-averaged restraints (window)
     */
    oparam_restr_winav = 1,
    /**
     * time-averaged restraints (window), weighted
     */
    oparam_restr_winav_weighted = 2
  };

   /**
   * @enum rdc_restr_enum
   * RDC restraints enumeration
   */
  enum rdc_restr_enum{
    /**
     * no restraints
     */
    rdc_restr_off = 0,
    /**
     * instantaneous restraints
     */
    rdc_restr_inst = 1,
    /**
     * instantaneous restraints, weighted
     */
    rdc_restr_inst_weighted = 2,
    /**
     * time-averaged restraints
     */
    rdc_restr_av = -1,
    /**
     * time-averaged restraints, weighted
     */
    rdc_restr_av_weighted = -2,
    /**
     * biquadratic (time averaged & instantaneous) restraints
     */
    rdc_restr_biq = -3,
    /**
     * biquadratic (time averaged & instantaneous) restraints, weighted
     */
    rdc_restr_biq_weighted = -4
  };

  /**
   * @enum rdc_mode_enum
   * Method of updating RDC magnetic field vectors enumeration
   */

  enum rdc_mode_enum {
      /**
       * Energy minimisation
       */
      rdc_em = 0,
      /**
       * Stochastic dynamics
       */
      rdc_sd = 1,
      /**
       * Molecular dynamics
       */
      rdc_md = 2
  };

  /**
   * @enum rdc_type_enum
   * Type of magnetic field representation
   */

  enum rdc_type_enum {
      /**
       * Magnetic field vectors
       */
      rdc_mf = 0,
      /**
       * Alignment tensor
       */
      rdc_t = 1,
      /**
       * Spherical harmonics
       */
      rdc_sh = 2
  };

  /**
 * @enum angle_restr_enum
 * Angle restraints enumeration
 */
  enum angle_restr_enum{
      /**
       * no restraints
       */
      angle_restr_off = 0,
      /**
       * instantaneous restraints
       */
      angle_restr_inst = 1,
      /**
       * instantaneous restraints, weighted
       */
      angle_restr_inst_weighted = 2,
      /**
       * angle constraints
       */
      angle_constr = 3
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
    /** lj_function */ lj_func,
    /** pol_lj_crf_function */ pol_lj_crf_func,
    /** pol_off_lj_crf_function */ pol_off_lj_crf_func,
    /** cgrain_function (MARTINI)*/ cgrain_func,
    /** cgrain_function (GROMOS) */ cggromos_func,
    /** default */ default_func
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
   * @enum charge_type_enum
   * use standard MM charges or special charge of QM buffer atoms
   */
  enum charge_type_enum {
    /**
     * standard MM charge
     */
    mm_charge = 0, 
    /**
     * special charge calculation for QM buffer atoms
     */
    qm_buffer_charge = 1
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
   * @enum xray_symrest_enum
   * x-ray symmetry restraints enum
   */
  enum xray_symrest_enum {
    /**
     * no symmetry restraints
     */
    xray_symrest_off = 0,
    /**
     * use symmetry restraints on the individual atoms
     */
    xray_symrest_ind = 1,
    /*
     * use symmetry contraints
     */
    xray_symrest_constr = 2,
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
    /**
    * A-EDS using Emax and Emin
    */
    aeds = 4,
    /**
    * A-EDS using Emax and Emin, search for EiR
    */
    aeds_search_eir = 5,
    /**
    * A-EDS using Emax and Emin, search for Emax and Emin
    */
    aeds_search_emax_emin = 6,
    /**
    * A-EDS using Emax and Emin, search for Eir, Emax an Emin
    */
    aeds_search_all = 7,
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
     * distance restraint interaction
     */
    disres_lambda = 8,
    /**
     * distancefield restraint interaction
     */
    disfield_lambda = 9,
    /**
     * dihedral restraint interaction
     */
    dihres_lambda = 10,
    /**
     * mass-scaling value
     */
    mass_lambda = 11,
    /**
     * angle restraint interaction
     */
    angres_lambda = 12,

    /**
     * one extra interaction for looping
     */
    last_interaction_lambda=13
  };

  /**
   * @enum nemd_enum
   * non-equilibrium molecular dynamics
   */
  enum nemd_enum {
    /**
     * don't use nemd
     */
    nemd_off = 0,
    /**
     * use nemd
     */
    nemd_on = 1
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
   * @enum bsleus_enum
   * do B&S-LEUS or not
   */
  enum bsleus_enum {
    /**
     * Don't use the B&S-LEUS method
     */
    bsleus_off = 0,
    /**
     * Use B&S-LEUS
     */
    bsleus_on = 1
  };

  /**
   * @enum electric_enum
   * electric field
   */
  enum electric_enum {
    /**
     * don't apply electric field
     */
    electric_off = 0,
    /**
     * apply electric field
     */
    electric_on = 1
  };

  /**
   * @enum qmmm_enum
   * do QM/MM, or not
   */
  enum qmmm_enum {
    /**
     * disable QM/MM
     */
    qmmm_off = 0,
    /**
     * enable QM/MM - mechanical embedding
     * QM charges are used for nonbonded QM-MM interaction.
     * QM atoms do not see any MM atoms.
     */
    qmmm_mechanical = 1,
    /**
     * enable QM/MM - electrostatic embedding
     * Nonbonded QM-MM interaction is modelled on QM level using
     * MM atoms as point charges. Only LJ interactions are
     * calculated clasically.
     */
    qmmm_electrostatic = 2,
    /**
     * enable QM/MM - polarisable embedding
     * Nonbonded QM-MM interaction is modelled on QM level using
     * MM atoms and their charge-on-spring as point charges.
     * Self-consistent field iteration is performed every step.
     * LJ interactions are calculated clasically.
     */
    qmmm_polarisable = 3
  };

  /**
   * @enum qm_ch_enum
   * update charges of QM atoms from the QM calculation
   * (only in mechanical embedding)
   */
  enum qm_ch_enum {
    /**
     * use constant charges from the topology
     */
    qm_ch_constant = 0,
    /**
     * update charges every step
     */
    qm_ch_dynamic = 1
  };

  /**
   * @enum qm_lj_enum
   * apply LJ between QM atoms
   */
  enum qm_lj_enum {
  /**
     * don't apply LJ dispersion between QM atoms
     */
    qm_lj_off = 0,
    /**
     * apply LJ dispersion between QM atoms
     */
    qm_lj_on = 1
  };

  /**
   * @enum qm_constr_enum
   * keep distance constraints within QM zone
   */
  enum qm_constr_enum {
  /**
     * remove constraints in QM zone
     */
    qm_constr_off = 0,
    /**
     * keep constraints in QM zone
     */
    qm_constr_on = 1
  };

  /**
   * @enum qmmm_software_enum
   * which QM software to use
   */
  enum qm_software_enum {
    /**
     * use MNDO
     */
    qm_mndo = 0,
    /**
     * use Turbomole
     */
    qm_turbomole = 1,
    /**
     * use DFTB
     */
    qm_dftb = 2,
    /**
     * use MOPAC
     */
    qm_mopac = 3,
    /**
     * use Gaussian
     */
    qm_gaussian = 4,
    /*
     * use Schnetpack NN
     */
    qm_nn = 5,
    /**
     * use Orca
     */
    qm_orca = 6,
    /**
     * use XTB
     */
    qm_xtb = 7
  };

  /**
   * @enum qmmm_nn_device_enum
   * which device to run NN on
   */
  enum qm_nn_device_enum {
    /**
     * Try CUDA, otherwise CPU
     */
    nn_device_auto = 0,
    /**
     * use CUDA
     */
    nn_device_cuda = 1,
    /**
     * use CPU
     */
    nn_device_cpu = 2
  };

  /**
   * @enum qmmm_nn_model_type_enum
   * specify if the model was trained on both the QM+buffer region and the buffer region
   * or on the difference of QM+buffer and buffer region
   */
  enum qmmm_nn_model_type_enum {
    /**
     * BuRNN model - trained on the difference between the QM+buffer region and the buffer region
     */
    nn_model_type_burnn = 0,
    /**
     * Standard model - trained on the QM+buffer region and the buffer region
     */
    nn_model_type_standard = 1
  };

  /**
   * @enum qmmm_nn_learning_type_enum
   * which device to run NN on
   */
  enum qmmm_nn_learning_type_enum {
    /**
     * 1: model was learned on all atoms (QMZONE + BUFFERZONE)
     */
    nn_learning_type_all = 1,
    /**
     * 2: model was learned by assigning energies only to the QMZONE atoms
     */
    nn_learning_type_qmonly = 2
  };

  /**
   * @class Parameter
   * input parameters.
   */
  class Parameter
  {
  public:
    Parameter() : title(GROMOSXX) {
        develop.develop = false;
    }

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
       * - ncyc 0      (number of steps before resetting of conjugate-gradient search direction)
       * - dele 0.0001 (minimal energy difference)
       * -             (conjugate-gradient - RMS force threshold)
       * - dx0  0.1    (initial step size)
       * -             (conjugate-gradient - initial and minimum step size)
       * - dxm  0.5    (maximum step size)
       * - nmin 1      (at least 1 step)
       * - flim 0.0    (no force limitation)
       * - cgim 3      (conjugate-gradient - maximum number of interpolations)
       * - cgic 1e-3   (conjugate-gradient - displacement threshold after interpolation)
       */
      minimise_struct() : ntem(0), ncyc(0), dele(0.0001),
			  dx0(0.1), dxm(0.5), nmin(1), flim(0.0), cgim(3), cgic(1e-3)
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
      /**
       * maximum interpolations.
       */
      int cgim;
      /**
       * interpolation displacement threshold.
       */
      double cgic;

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
       * - algorithm 0 (weak coupling)
       */
      multibath_struct() : couple(false), found_multibath(false), found_tcouple(false), algorithm(0) {}

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
       * algorithm?
       *  0 : berendsen
       *  1 : Nose-Hoover
       * >1 : Nose-Hoover-Chains
       */
      int algorithm;

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
        x_semi = y_semi = z_semi = 1;
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
      /**
       * semianisotropic couplings
       */
      int x_semi, y_semi, z_semi;
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
      print_struct() : stepblock(0), centreofmass(0), monitor_dihedrals(0) {}

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
      int monitor_dihedrals;
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
    /**
     * Number of GPUs
     */
    unsigned int number_gpus;
    /**
     * GPU IDs
     */
    std::vector<int> gpu_device_number;
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
       * - crossdihedral 1
       * - nonbonded 1
       * - energy_group empty
       * - special_loop -1
       * - interaction function lj_crf_func
       * - force_groups false
       */
      force_struct() : bond(1), angle(1), improper(1),
		       dihedral(1), crossdihedral(1), nonbonded_vdw(1),
		       nonbonded_crf(1), special_loop(special_loop_off),
		       interaction_function(lj_crf_func),
		       force_groups(false)
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
       * crossdihedral?
       */
      int crossdihedral;
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
       * use energy groups also for forces
       */
      bool force_groups;

    } /** Force(field) parameters */ force;

#ifdef HAVE_HOOMD
	/**
	 * @strict hoomd_struct
	 * HOOMD block
	 */
	struct hoomd_struct
	{
		hoomd_struct() : processor(unknown) {}
	    ~hoomd_struct() {}
		enum hoomd processor;
	} /** Hoomd parameter */ hoomd;
#endif

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
       * - selfterm_excluded_atoms 1 (selfterm considered)
       * - rf_excluded true (new standard)
       * - epsilon     1
       */
      nonbonded_struct() :
        method(el_reaction_field),
        lserf(false),
        rf_kappa(0.0),
        rf_cutoff(1.4),
        rf_epsilon(66.0),
        selfterm_excluded_atoms(1),
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
       * self term of the excluded atoms for reaction field
       */
      int selfterm_excluded_atoms;
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
       * - local_elevatin false
       * - force_constant 10000
       * - write 0
       * - writexmap 0
       * - tau 0
       * - spacegroup "P 1"
       * - resolution 1.0
       * - readavg 0
       * - to_angstrom 10.0
       * - symrest (off)
       * - sym_force_constant 0.0
       * - sym_spacegroup "P 1"
       */
      xrayrest_struct() : xrayrest(xrayrest_off),
      mode(xrayrest_mode_structure_factor),
      local_elevation(false),
      force_constant(1E4),
      tau(0),
      write(0),
      writedensity(0),
      writexmap(0),
      spacegroup("P 1"),
      resolution(1.0),
      readavg(0),
      to_angstrom(10.0),
      symrest(xray_symrest_off),
      sym_force_constant(0.0),
      sym_spacegroup("P 1") {
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
       * local elevation
       */
      bool local_elevation;
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
      /**
       * converison factor for length unit to angstrom
       */
      double to_angstrom;
      /**
       * do symmetry restraints?
       */
      xray_symrest_enum symrest;
      /**
       * force constant for symmetry restraints
       */
      double sym_force_constant;
      /**
       * symmetry spacegroup
       */
      std::string sym_spacegroup;
/**
       * @struct bfactor_struct
       * B factor settings structure
       * default value:
       * - step 0: don't fit B factors
       * - terminate_iterations 10
       * - terminate_gradient 0.1
       * - max 1.0
       * - min 0.001
       */
      struct bfactor_struct {
        bfactor_struct() : init(true),
        step(0),
        terminate_iterations(10),
        terminate_gradient(0.1),
        max(1.0),
        min(0.001) {
        }
        /**
         * init them from topology
         */
        bool init;
        /**
         * fit B factors every step
         */
        int step;
        /**
         * terminate after the number of iterations
         */
        int terminate_iterations;
        /**
         * terminate after |gradient| &lt; value
         */
        double terminate_gradient;
        /**
         * maximum B-factor
         */
        double max;
        /**
         * minimum B-factor
         */
        double min;
        /**
         * B factor groups
         */
        std::vector<std::set<unsigned int> > groups;
      } /** B factor settings */ bfactor;
     /**
     * @struct replica_exchange_parameters_struct
     * default values:
      * -no replica exchange
      * -minimum value set to 0.0
      * -maximum value set to 0.0
     */
      struct replica_exchange_parameters_struct
    {
        replica_exchange_parameters_struct() :
        switcher(replica_exchange_off),
        lambda_dependant_min(0.0),
        lambda_dependant_max(0.0),
        energy_switcher(energy_tot)
        {
        }
        /**
         *replica exchange switcher
         */
        replica_exchange_interruptor_enum switcher;
        /**
         * minimal force constant
         */
        double lambda_dependant_min;
        /**
         * maximal force constant
         */
        double lambda_dependant_max;
        /**
         *switcher for the energy
         */
        rep_ex_energy_interruptor_enum energy_switcher;
     }replica_exchange_parameters;
      /**
     * @struct overall_bfactor_struct
     * default values: 0
     */
     struct overall_bfactor_struct
    {
        overall_bfactor_struct() :
        B_overall_switcher(B_overall_off), init(0.0)
        {
        }
        /**
         * overall B factor switching
         */
        B_overall_enum B_overall_switcher;
        /**
         * initial value (read from specification file)
         */
        double init;
        /**
         * terminate after the number of iterations
         */
        int terminate_iterations;
        /**
         * terminate after |gradient| &lt; value
         */
        double terminate_gradient;
     }overall_bfactor;
     /**
     * @struct structure_factor_calculation_struct
     * default values:
     */
      struct structure_factor_calculation_struct
    {
        structure_factor_calculation_struct() :
        atom_move_tolerance(0.0),
        steps_nb_constant(1)
        {
        }
        /**
         * tolerance for atom move
         */
        double atom_move_tolerance;
        /**
         * every how many steps
         */
        unsigned int steps_nb_constant;
     } structure_factor_calculation;
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
       * - write 0
       */

      distanceres_struct()
	: distanceres(0),
	  K(0),
	  r_linear(0),
	  tau(10),
	  read(0),
	  virial(0),
	  forcescale(0),
          write(0)
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

      /**
       * compute virial contribution
       */
      unsigned int virial;

      /**
       * force scaling according to equation 8.17
       */
      unsigned int forcescale;

      /**
       * write on/off
       */
      unsigned int write;

    }/** Distance restraints parameters */ distanceres;
     /**
     * @struct distancefield_struct
     * DISTANCEFIELD block
     */
    struct distancefield_struct
    {
      /**
       * Constructor
       * Default values:
       * - distancefield 0 (no distance restraints)
       * - grid 0
       * - proteinoffset 0
       * - proteincutoff 0
       * - smooth 0
       * - r_l 0
       * - write 0
       * - update 1
       */

      distancefield_struct()
	: distancefield(0),
	  grid(0),
	  proteinoffset(0),
	  proteincutoff(0),
	  smooth(0),
          r_l(0),
          write(0),
          printgrid(false),
	  update(1),
          protect(0)
      {
      }

      /**
       * distance field restraints on/off
       */
      int distancefield;
      /**
       * grid length
       */
      double grid;
      /**
       * protein offset as penalty
       */
      double proteinoffset;
      /**
       * cutoff to determine gridpoints within the protein
       */
      double proteincutoff;
      /**
       * flag to smoothen the forces
       */
      int smooth;
      /**
       * distance where potential gets linear
       */
      double r_l;
      /**
       * write on/off
       */
      unsigned int write;
      /**
       * print the final grid to file
       */
      bool printgrid;
      /**
       * update frequency
       */
      int update;
      /**
       * radius around the zero-grid-point which is protected
       * from being flagged as protein
       */
      double protect;
    }/** Distancefield restraints parameters */ distancefield;

    struct angrest_struct
            {
        /**
         * Constructor
         * Default values:
         * - angrest 0 (no angle restraints)
         * - K 0
         */
        angrest_struct()
        : angrest(angle_restr_off),
        K(0.0),
        virial(0),
        write(0) {}

        /**
         * angle restraints
         * method:
         * - 0: off
         * - 1: uniform K
         * - 2: K * Ki (weight by Ki in angle restraint file)
         * - 3: constraints
         */
        angle_restr_enum angrest;

        /**
         * force constant K
         */
        double K;

        /**
        * compute virial contribution
        */
        unsigned int virial;

        /**
         * write on/off
         */
        unsigned int write;

        /**
         * tolerance
         */
        double tolerance;

    }/** angle restraint parameters */ angrest;


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
        phi_lin(0.0),
        virial(0),
        write(0) {}

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

        /**
        * compute virial contribution
        */
        unsigned int virial;

        /**
         * write on/off
         */
        unsigned int write;

        /**
         * tolerance
         */
        double tolerance;

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
			 soft_vdw(0.0), soft_crf(0.0), perturbed_par(false) {}

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
      /**
       * will be set to true if any perturbed parameter
       * is read in read_special  or read_topology
       */
      bool perturbed_par;

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
          tarfscale(0),
          biqweight(0),
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
       * force scaling in time-averaged restraining
       */
      unsigned int tarfscale;
      /**
       * weighting of the two terms in biquadratic restraining
       */
      unsigned int biqweight;
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
       * - temperature \<empty\>
       * - scale (false)
       * - lambda \<empty\>
       * - dt \<empty\>
       * - trials 0
       * - equilibrate 0
       * - cont 0
       */
      replica_struct() : retl(false), num_T(0), num_l(0), scale(false), trials(0),
                         equilibrate(0), cont(0)
      {
      }
      /**
       * Shall a Replica_exchange Temperature or lambdaDep simulation be exectued
       */
      bool retl;
      /**
       * number of replicas with different temperature
       */
      int num_T;
      /**
       * number of replicas with different lambdas
       */
      int num_l;
      /**
       * temperatures
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
       * do continuation run
       */
      int cont;

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
       * - EPSM (20.0)
       */
      cgrain_struct()
	: level(0), EPS(20.0), EPSM(20.0)
      {
      }
      /**
       * do coarse graining
       * - 0 atomistic (off)
       * - 1 coarse-graining using MARTINI model (on)
       * - 2 coarse-graining using GROMOS model (on)
       * - 3 mixed-graining using GROMOS model (on)
       */
      int level;
      /**
       * EPS for the pure CG coulomb
       */
      double EPS;
      /**
       * EPS for the mixed FG-CG coulomb
       */
      double EPSM;
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
       * - no_constraints (false)
       * - stride (1)
       * - trajectory("")
       */
      analyze_struct() : analyze(false), copy_pos(false), trajectory(""),
                         stride(1), no_constraints(false)
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
      /**
       * stride
       */
      int stride;
      /**
       * do not apply any constraints (also not on solvent)
       */
      bool no_constraints;

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


    /**
     * ANITA
     * @struct precalclam_struct
     * pre-calculate energies for other lambda values
     */
    struct precalclam_struct
    {
      /**
       * constructor
       * default values:
       * - nr_lambdas (0)
       * - min_lam (0.0)
       * - max_lam (1.0)
       */
      precalclam_struct() : nr_lambdas(0),
                          min_lam(0.0),
                          max_lam(1.0)
      {
      }
      /**
       * calculate nr_lambdas extra lambda points
       */
       unsigned int nr_lambdas;
      /**
       * starting from lambda
       */
       double min_lam;
      /**
       * up to lambda
       */
       double max_lam;

    } /** precalculate lambdas struct */ precalclam;
    // END ANITA

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
       * nonbonded terms
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
      polarise_struct() : cos(0), minfield(37.3), efield_site(ef_atom),
                          damp(false), write(0)
      {}
      /**
       * use charge-on-spring polarisation
       */
      int cos;
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
       * - gsl
       */
      rng_struct() : rng(random_gsl), gsl_rng(-1) {}
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
      eds_struct() : eds(false), soft_vdw(0.0), soft_crf(0.0), form(single_s), numstates(0) {}
      /**
       * do enveloping distribution sampling using the Hamiltonian:
       */
      unsigned int eds;
      /**
       * soft core van der Waals interactions
       */
      double soft_vdw;
      /**
       * soft core electrostatic interactions
       */
      double soft_crf;
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
       * position information first: start position; second: current position of coord_ID
       */
      std::pair<int, int> pos_info;
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
      /**
      * parameter emax for aeds
      */
      double emax;
      /**
      * parameter emin for aeds
      */
      double emin;
      /**
      * do we want to init an aeds parameter search?
      */
      bool initaedssearch;
      /**
      * current maximum transition energy within a state round-trip
      */
      double searchemax;
      /**
      * how many emaxes did we already find?
      */
      unsigned int emaxcounts;
      /**
      * ln of exponential energy differences between the states and the reference state
      */
      std::vector<double> lnexpde;
      /**
      * free energy differences between the states and the reference state
      */
      std::vector<double> statefren;
      /**
      * states that were already visited within a state round-trip
      */
      std::vector<bool> visitedstates;
      /**
      * how many times did we visit a state?
      */
      std::vector<unsigned int> visitcounts;
      /**
      * state of the last simulation step
      */
      unsigned int oldstate;
      /**
      * average energy of an end-state
      */
      std::vector<double> avgenergy;
      /**
      * average energy including offset of an end-state
      */
      std::vector<double> eiravgenergy;
      /**
      * helper variable for calculation of running standard deviation of the end-state energies
      */
      std::vector<double> bigs;
      /**
      * running standard deviation of the end-state energies
      */
      std::vector<double> stdevenergy;
      /**
      * which kind of bmax is given in the input parameters?
      */
      unsigned int bmaxtype;
      /**
      * the maximum energy barrier parameter
      */
      double setbmax;
      /**
      * do we want to accelerate over the minimum average energy of the end-states?
      */
      bool fullemin;
      /**
      * half-life of the offset parameters at the beginning of the run
      */
      unsigned int asteps;
      /**
      * half-life of the offset parameters at the beginning of the run
      */
      unsigned int bsteps;
    } /** enveloping distribution sampling*/ eds;

 struct reeds_struct : public replica_struct
    {
      /**
       * Constructor
       * Default values:
       * - num_s 0
       * - num_eoff 0
       * - temperature \<empty\>
       * - scale (false)
       * - svals \<empty\>
       * - dt \<empty\>
       * - trials 0
       * - equilibrate 0
       * - cont 0
       */
      reeds_struct() : reeds(0),
                       num_states(0),  num_s(0), num_eoff(0),
                       trials(0), equilibrate(0),
                       cont(0), eds_stat_out(true), periodic(false) {}
      /**
       * Check what kind of reed run.f this is
       **/
      int reeds;
      /**
       * num_states
       */
      int num_states;

      /**
       * number of replicas with different s in REEDS these are the smoothing values
       */
      int num_s;
      /**
       * * number of energy offsets param-vectors with different offsets for each state in REEDS
       * * length of one param-vector == NUMSTATES
       */
      int num_eoff;
      /**
       * temperatures
       */
      double temperature;
      /**
       * lambdas: contains all smoothness parameter of RE_EDS system
       */
      std::vector<double> svals;
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
       * do continuation run
       */
      int cont;
      /**
       * write output to stat_file (repdat)
       **/
      bool eds_stat_out;
      /**
       * periodic boundary?
       **/
      bool periodic;
       /**
       * for RE-EDS Sim many eds parameters have to be accessible for
       * energy calculation.
       */
      std::vector<eds_struct> eds_para;

    } /** replica exchange parameters */ reeds;

    /**
     * @struct sasa
     * parameters for calculating the sasa and volume term
     */
    struct sasa_struct {
      /**
       * Constructor: disable SASA
       */
      sasa_struct() : switch_sasa(false), switch_volume(false), p_12(0.0), p_13(0.0), p_1x(0.0),
      sigma_v(0.0), r_solv(0.0), max_cut(0.0), min_cut(0.0) {}
      /**
       * SASA switch parameter
       */
      bool switch_sasa;
      /**
       * volume switch parameter
       */
      bool switch_volume;
      /**
       * p_ij for one bond interaction, first neighbours
       */
      double p_12;
      /**
       * p_ij for three and four bond interaction, second and third neighbours
       */
      double p_13;
      /**
       * pij for > 1,4
       */
      double p_1x;
      /**
       * sigma value for volume term, will be deleted and stored an topology?
       */
      double sigma_v;
      /**
       * solvent radius
       */
      double r_solv;
      /**
       * parameters for switching function, upper cutoff
       */
      double max_cut;
      /**
       * parameters for switching function, lower cutoff
       */
      double min_cut;
      /**
       * difference between upper and lower cutoff
       */
      double cut_diff;

    } /** sasa */ sasa;

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
      innerloop_struct() : solvent(sls_topo), method(sla_off), number_gpus(0)  {}
      /**
       * the solvent
       */
      special_loop_solvent_enum solvent;
      /**
       * the acceleration method
       */
      special_loop_acceleration_enum method;
      /**
       * The number of GPUs used for CUDA
       */
      unsigned int number_gpus;
      /**
       * Which device number should be used for CUDA
       */
      std::vector<int> gpu_device_number;
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
      localelev_struct() : localelev(localelev_off), read(false), write(0) {}
      /**
       * use the method or not
       */
      localelev_enum localelev;
      /**
       * read from external file
       */
      bool read;
      /**
       * write to special traj
       */
      int write;
      /**
       * ids of the umbrella potentials to apply
       * true if building up
       */
      std::map<int, bool> umbrellas;
    } localelev;

    /**
     * @struct bsleus_struct
     * Constructor:
     * Default values:
     * - bslues: off
     */
    struct bsleus_struct {
      /**
       * constructor
       */
      bsleus_struct() : bsleus(bsleus_off), building(0), write(0), transition_conf(false) {}
      /**
       * Use B&S-LEUS or not
       */
      bsleus_enum bsleus;
      /**
       * Are we building the potential?
       */
      bool building;
      /**
       * How the memory is incremented (k_LE)
       */
      //double forceConstantIncrement;
      /**
       * Do we write the bsleus potential?
       */
      int write;
      /**
       * Is this just the configuration along the transition path, which doesn't
       * need a velocity?
       */
      bool transition_conf;
    } bsleus;

    /**
     * @struct electric_struct
     * Constructor:
     * Default values:
     * - field, dipole, current: off
     */
    struct electric_struct {
      /**
       * constructor
       */
      electric_struct() : electric(electric_off), Ef_x(0.0), Ef_y(0.0), Ef_z(0.0),
                          dipole(false), dip_groups(0), dip_write(0),
                          current(false), cur_groups(0), cur_write(0) {}
      /**
       * use the method or not
       */
      electric_enum electric;
       /**
       * external field components x, y and z
       */
      double Ef_x, Ef_y, Ef_z;
      /**
       * Dipole calc and write to special traj
       */
      bool dipole;
      /**
       * Calculate the box dipole considering the groups
       * 0 : solute
       * 1 : solvent
       * 2 : all
       */
      unsigned int dip_groups;
      /**
       * Write dipole to special traj every dip_write steps
       */
      unsigned int dip_write;
      /**
       * Current calc and write to special traj
       */
      bool current;
      /**
       * Number of current groups
       */
      unsigned int cur_groups;
      /**
       * Vector including the last atom of each group
       */
      std::vector<unsigned int> current_group;
      /**
       * Write current to special traj every cur_write steps
       */
      unsigned int cur_write;

    } electric;

    struct nemd_struct {
      /**
       * constructor
       * Default values
       * nemd (nemd_off): do not performd nemd
       * method
       * 0: velocity exchange
       * For now only velocity exchange is implemented
       */
      nemd_struct() : nemd(nemd_off),  property(0), method(0),
                      slabnum(1), pertfrq(1),
                      ampbath(0), stdyaft(0), write(0) {}
      /**
       * use the method or not
       */
      nemd_enum nemd;
      /**
       * define method
       */
      unsigned int property;
       /**
       * define method
       */
      unsigned int method;
      /**
       * number of slabs used in grid based methods (discretized in z-direction)
       */
      unsigned int slabnum;
      /**
       * pertfrq (frequency of applied perturbation)
       */
      unsigned int pertfrq;
      /**
       * amplitude of external bath (for possible implementation of weak-coupling method)
       * or amplitude of perturbation (for the PPM method)
       */
      double ampbath;
      /**
       * after this point a steady state is assumed and data is accumulated
       */
      unsigned int stdyaft;
      /**
       * write every nth timesteps (write the velocities and flux)
       */
      unsigned int write;

    } nemd;

    struct multigradient_struct {
      /**
       * @struct multigradient_struct
       * constructor
       * default values: disable multigradient
       */
      multigradient_struct() : multigradient(false),
      print_graph(true), print_curve(false)
      {
      }
      /**
       * enable the gradients
       */
      bool multigradient;
      /**
       * print the graphs to the output file
       */
      bool print_graph;
      /**
       * print the curve data to the output file
       */
      bool print_curve;
      /**
       * the variables to affect
       */
      std::vector<std::string> variable;
      /**
       * the functional form
       */
      std::vector<int> functional_form;
      /**
       * the control points
       */
      std::vector<std::vector<std::pair<double, double> > > control_points;
    } multigradient;

    /**
     * @struct addecouple_struct
     * Constructor:
     * Default values:
     * - adgr = 0
     */
    struct addecouple_struct {
      /**
       * constructor: no addiabatic decoupling
       */
      addecouple_struct() : adgr(0), tmf(0), write(0) {
      }
      /**
       * number of addiabatic decoupling groups
       */
      unsigned int adgr;

      struct adc_struct {

        adc_struct(int adstart, int adend, double sm, double sv,
                double st, int tir, int eg, int tg)
        : adstart(adstart), adend(adend), sm(sm), sv(sv), st(st),
        tir(tir), eg(eg), tg(tg) {
        }
        /**
         * first atom of addiabatic decoupling groups
         */
        int adstart;
        /**
         * last atom of addiabatic decoupling groups
         */
        int adend;
        /**
         * scaling factor mass
         */
        double sm;
        /**
         * scaling factor potential energy function
         */
        double sv;
        /**
         * scaling factor temperature
         */
        double st;
        /**
         * which temperature bath to scale
         */
        int tir;
        /**
         * energy group of decoupled group
         */
        int eg;
        /**
         * temperature group of decoupled group
         */
        unsigned int tg;
      };
      std::vector<adc_struct> m_adc_index;

      void add_adc(adc_struct adc) {
        m_adc_index.push_back(adc);
      }

      void add_adc(int adstart, int adend, double sm, double sv,
              double st, int tir, int eg, int tg) {
        m_adc_index.push_back(adc_struct(adstart, adend, sm, sv, st, tir, eg, tg));
      }

      int check_index_adc(int atom_number)const {
        for (unsigned int i = 0; i < m_adc_index.size(); i++) {
          if (atom_number >= m_adc_index[i].adstart && atom_number <= m_adc_index[i].adend)
            return i;
        }
        return -1;
      }


      int check_index_adc(int atom_number) {
        for (unsigned int i = 0; i < m_adc_index.size(); i++) {
          if (atom_number >= m_adc_index[i].adstart && atom_number <= m_adc_index[i].adend)
            return i;
        }
        return -1;
      }

      std::vector<adc_struct> & adc_index() {
        return m_adc_index;
      }

      std::vector<adc_struct> const & adc_index()const {
        return m_adc_index;
      }
      /**
       * tau mean field
       */
      double tmf;
      /**
       * printing of the special trajectory
       */
      int write;
    } /** addecouple */ addecouple;

    /**
     * @struct orderparamres_struct
     * ORDERPARAMRES block
     */
    struct orderparamrest_struct
    {
      /**
       * Constructor
       * Default values:
       * - orderparamrest 0 (no order parameter restraints)
       * - K 0.0
       * - tau 10.0
       * - read false
       * - write 0
       */

      orderparamrest_struct()
	: orderparamrest(oparam_restr_off),
	  K(0.0),
	  tau(10.0),
          update_step(1),
	  read(false),
          write(0)
      {
      }

      /**
       * order parameter restraing method
       */
      oparam_restr_enum orderparamrest;

      /**
       * force constant K
       */
      double K;
      /**
       * memory time for time averaging
       */
      double tau;
      /**
       * update order parameter step
       */
      unsigned int update_step;

      /**
       * read on/off
       */
      bool read;
      /**
       * write on/off, every n-th step
       */
      unsigned int write;
    }/** order parameter restraints parameters */ orderparamrest;

     /**
     * @struct rdc_struct
     * RDC restraint parameters.
     */
    struct rdc_struct
    {
      /**
       * Constructor
       * Default values:
       * - mode restr_off
       * - read_av false
       * - type rdc_mf
       * - read_mfv false
       * - method rdc_em
       * - emgradient 0.0
       * - emstepsize 0.0
       * - emmaxiter 0
       * - sdfric 0.0
       * - temp 0.0
       * - delta 0.0
       * - K 1.0
       * - tau 0
       * - tAVfactor 0
       * - biqfactor 0
       * - write 0
       */

      rdc_struct()
	: mode(rdc_restr_off),
          read_av(false),
          type(rdc_mf),
          read_align(false),
          method(rdc_em),
          emgradient(0.0),
          emstepsize(0.0),
          emmaxiter(0),
          sdfric(0.0),
          temp(0.0),
          delta(0.0),
          K(1.0),
          tau(0),
          tAVfactor(0),
          biqfactor(0),
          write(0)
      {
      }
      /**
       * restraining mode.
       */
      rdc_restr_enum mode;
      /**
       * read averages
       */
      bool read_av;
      /**
       * type of magnetic field vector representation
       */
      rdc_type_enum type;
      /**
       * read magnetic field vectors
       */
      bool read_align;
      /**
       * method of updating the magnetic field vectors
       */
      rdc_mode_enum method;
      /**
       * EM: stop if gradient is below emgradient
       */
      double emgradient;
      /**
       * EM: start with emstepsize
       */
      double emstepsize;
      /**
       * EM: stop after emmaxiter are reached
       */
      unsigned int emmaxiter;
      /**
       * SD: friction coefficient
       */
      double sdfric;
      /**
       * reference temperature for SD and for initial velocities
       */
      double temp;
      /**
       * half the width of the flat bottom potential
       */
      double delta;
      /**
       * force constant
       * (multiplied by individual restraint weighting)
       */
      double K;
      /**
       * coupling time.
       */
      double tau;
      /**
       * choose if (1 - e^(...)) factor is omitted in time AV
       */
      unsigned tAVfactor;
      /**
       * choose factor by which biquad term is multiplied (1, (1 - e^(...)), 0)
       */
      unsigned biqfactor;
       /**
       * write output to special trajectory every n-th step
       */
      unsigned int write;
    } /** RDC-parameters */ rdc;



    struct qmmm_struct {
      /**
       * Constructor
       * Default values:
       * - cutoff 0.0 (no cutoff)
       * - cap_length 0.109 (capping atom distance)
       * - mm_scale -1.0 (no scaling)
       * - qmmm qmmm_off (no QMMM)
       * - qm_lj qm_lj_off (no LJ dispersion within QM zone)
       * - qm_constraint qm_constr_off (no constraints in QM zone)
       * - software qm_mndo (MNDO)
       * - write 0 (no writing)
       * - atomic_cutoff false (using charge-group based cutoff)
       * - use_qm_buffer false (not using buffer zone)
       */
      qmmm_struct() :
                      cutoff(0.0)
                    , cap_length(0.109)
                    , mm_scale(-1.0)
                    , qmmm(qmmm_off)
                    , qm_ch(qm_ch_constant)
                    , qm_constraint(qm_constr_off)
                    , software(qm_mndo)
                    , write(0)
                    , atomic_cutoff(false)
                    , use_qm_buffer(false)
                    , dynamic_buffer_charges(false) {}
      /**
       *
       * Common QMMM parameters
       *
       */
      /**
       * cutoff to determine atoms included in QM calculation as point charges.
       */
      double cutoff;
      /**
       * Capping atom bond length
       */
      double cap_length;
      /**
       * scaling factor for the MM charges in the QM/MM interaction
       */
      double mm_scale;
      /**
       * QM-MM embedding scheme or disable
       */
      qmmm_enum qmmm;
      /**
       * QM-MM embedding scheme or disable
       */
      qm_ch_enum qm_ch;
      /**
       * apply LJ interaction in QM zone or not
       */
      qm_lj_enum qm_lj;
      /**
       * keep constraints in QM zone and QM-MM link
       */
      qm_constr_enum qm_constraint;
      /**
       * the QM software to use
       */
      qm_software_enum software;
      /**
       * write QM/MM related stuff to special trajectory
       */
      unsigned write; // What can be written here?
      /**
       * type of cutoff (atomic or chargegroup-based)
       */
      bool atomic_cutoff;
      /**
       * if QM buffer zone is used
       */
      bool use_qm_buffer;
      /**
       * if dynamic charges are used with QM buffer zone
       */
      bool dynamic_buffer_charges;

      /**
       * QM zone parameters
       */
      struct qm_zone_struct {
      /**
       * Constructor
       * Default values:
       * - charge 0 (neutral)
       * - spin_mult 1 (no unpaired electrons)
       */
      qm_zone_struct() :
                      charge(0)
                    , spin_mult(1) {}
        /**
         * net charge
         */
        int charge;
        /**
         * spin multiplicity
         */
        int spin_mult;
      } qm_zone;

      /**
       * QM buffer zone parameters
       */
      struct buffer_zone_struct : qm_zone_struct {
      /**
       * Constructor
       * Default values:
       * - cutoff 0.0 (no adaptive QM buffer)
       */
      buffer_zone_struct() :
                      cutoff(0.0) {}
        /**
         * Adaptive buffer zone cutoff
         */
        double cutoff;
      } buffer_zone;

      /**
       * QM program unspecific parameters
       */
      struct qm_param_struct{
      /**
       * Constructor
       * Default values:
       * - unit_factor_length 1.0
       * - unit_factor_energy 1.0
       * - unit_factor_force 1.0
       * - unit_factor_charge 1.0
       */
      qm_param_struct() :
                      unit_factor_length(1.0)
                    , unit_factor_energy(1.0)
                    , unit_factor_force(1.0)
                    , unit_factor_charge(1.0) {}
        /**
         * maps atomic number to elements name
         */
        std::map<unsigned, std::string> elements;
        /**
         * maps IAC numbers to atomic numbers; 
         */
        std::map<unsigned, unsigned> iac_elements;
        /**
         * path for the program binary
         */
        std::string binary;
        /**
         * path for the input file. Empty for a temporary file
         */
        std::string input_file;
        /**
         * path for the output file. Empty for a temporary file
         */
        std::string output_file;
        /**
         * header of the input file
         */
        std::string input_header;
        /**
         * factor to convert the QM length unit to the GROMOS one
         */
        double unit_factor_length;
        /**
         * factor to convert the QM energy unit to the GROMOS one
         */
        double unit_factor_energy;
        /**
         * factor to convert the QM energy unit to the GROMOS one
         */
        double unit_factor_force;
        /**
         * factor to convert the QM charge unit to the GROMOS one
         */
        double unit_factor_charge;
      };

      /**
       * MNDO specific parameters
       */
      struct mndo_param_struct : public qm_param_struct {
        /**
         * path for the gradient output file. Empty for a temporary file
         */
        std::string output_gradient_file;
        /**
         * path for the density matrix file. Empty for a temporary file
         */
        std::string density_matrix_file;
      } mndo;

      /**
       * Turbomole specific parameters
       */
      struct turbomole_param_struct : public qm_param_struct {
        /**
         * the tools to run in the working directory
         */
        std::vector<std::string> toolchain;
        /**
         * the working directory containing the control file
         */
        std::string working_directory;
        /**
         * the directory containing the turbomole binaries
         */
        std::string binary_directory;
        /**
         * the input file containing the atomic coordinates
         */
        std::string input_coordinate_file;
        /**
         * the input file containing the positions of the MM atoms
         */
        std::string input_mm_coordinate_file;
        /**
         * the output file containing the energy
         */
        std::string output_energy_file;
        /**
         * the output file containing the cartesian gradients
         */
        std::string output_gradient_file;
        /**
         * the output file containing the cartesian gradients of the MM atoms
         */
        std::string output_mm_gradient_file;
        /**
         * the output file containing the ESP charges of the QM atoms
         */
        std::string output_charge_file;
      } turbomole;

      /**
       * DFTB specific parameters
       */
      struct dftb_param_struct : public qm_param_struct {
        /**
         * the working directory containing the dftb_in.hsd
         */
        std::string working_directory;
        /**
         * path of the input geometry file
         */
        std::string input_coordinate_file;
        /**
         * path of the input MM charges geometry file
         */
        std::string input_mm_coordinate_file;
        /**
         * path for the stdout file
         */
        std::string stdout_file;
      } dftb;

      /**
       * MOPAC specific parameters
       */
      struct mopac_param_struct : public qm_param_struct {
        /**
         * Constructor
         * Default values:
         * - link_atom_mode 0
         */
        mopac_param_struct() :
                        link_atom_mode(0) {}
        /**
         * path for the output aux file. Empty for a temporary file
         */
        std::string output_aux_file;
        /**
         * path for the output arc file. Empty for a temporary file
         */
        std::string output_arc_file;
        /**
         * path for the stdout file. Empty for a temporary file
         */
        std::string stdout_file;
        /**
         * path for the output aux file. Empty for a temporary file
         */
        std::string output_dens_file;
        /**
         * path for the molin file. Empty for a temporary file
         */
        std::string molin_file;
        /**
         * link atom treatment mode
         */
        int link_atom_mode;
      } mopac;
      
      /**
       * Gaussian specific parameters
       */
     struct gaussian_param_struct : public qm_param_struct{
        /**
         * route section of the input file
         */
        std::string route_section;
        /**
         * total charge and spin multiplicity in the input file
         */
        std::string chsm;
      } gaussian;

     /**
       * ORCA specific parameters
       */
      struct orca_param_struct : public qm_param_struct { 
        /**
         * the input file containing the positions and element types of the QM atoms
         */
        std::string input_coordinate_file;
        /**
         * the input file containing the positions and charges of the MM atoms
         */
        std::string input_pointcharges_file;
        /**
         * the output file containing the cartesian gradients
         */
        std::string output_gradient_file;
        /**
         * the output file containing the cartesion gradients of the MM atoms
         */
        std::string output_mm_gradient_file;
      } orca; 

      /**
       * XTB specific parameters
       */
      struct xtb_param_struct : public qm_param_struct { 
        /**
         * the version of the XTB Hamiltonian used
         * options are 1 and 2 (for GFN1-xTB or GFN2-xTB, respectively)
         */
        unsigned int hamiltonian;
        /**
         * the verbosity level of XTB
         * options are 0, 1, and 2 corresponding to muted, minimal, or full verbosity, respectively
         */
        unsigned int verbosity;
        /**
         * the input file containing the positions and element types of the QM atoms
         */
        std::string input_coordinate_file;
        /**
         * the input file containing the positions and charges of the MM atoms
         */
        std::string input_pointcharges_file;
        /**
         * the output file containing the cartesian gradients
         */
        std::string output_gradient_file;
        /**
         * the output file containing the cartesion gradients of the MM atoms
         */
        std::string output_mm_gradient_file;
      } xtb; 

      /**
       * NN specific parameters
       */
      struct nn_param_struct : public qm_param_struct {
      /**
       * Constructor
       * Default values:
       * - model_path "" (empty string)
       * - val_model_path "" (empty string)
       * - model_type 0 (bool)
       * - device 0 (auto)
       */
      nn_param_struct() :
                      model_path()
                      , val_model_path() 
                      , val_thresh(0.0)
                      , val_steps(0)
                      , val_forceconstant(0.0)
                      , charge_model_path()
                      , charge_steps(0)
                      , model_type(nn_model_type_burnn) 
                      , learning_type(nn_learning_type_all)
                      , device(nn_device_auto) {}
        /**
         * Schnetpack model path
         */
        std::string model_path;
        /**
         * Schnetpack model path
         */
        std::string val_model_path;
        /**
         * Threshold of energy validation
         */
        double val_thresh;
        /**
         * Number of steps between validations
         */
        unsigned val_steps;
        /**
         * Force constant to enforce agreement between NN models
         */
        double val_forceconstant;
        /**
         * Schnetpack model path
         */
        std::string charge_model_path;
        /**
         * Number of steps between validations
         */
        unsigned charge_steps;
        /**
         * nn model type
         */
        qmmm_nn_model_type_enum model_type;
        /**
         * nn learning type
         */
        qmmm_nn_learning_type_enum learning_type;
        /**
         * Device to run model on
         */
        qm_nn_device_enum device;
      } nn;

    } qmmm;


    struct symrest_struct {
      /**
       * Constructor
       * - no symmetry restraints
       * - zero force constant
       */
      symrest_struct():
      symrest(xray_symrest_off),
      force_constant(0.0) {}
      /**
       * symmetry restraints
       */
      xray_symrest_enum symrest;
      /**
       * force constant
       */
      double force_constant;
      /**
       * symmetry operations
       */
      std::vector<std::pair<math::Matrix, math::Vec> > symmetry_operations;
    } /* symmetry restraints */symrest;

    /**
     * @struct amber_struct
     * AMBER block
     */
    struct amber_struct {
      /**
       * Constructor
       * - no AMBER topology
       * - electrostatic interaction scaling for 1,4-interactions 1.2
       */
      amber_struct():
      amber(false),
      coulomb_scaling(1.2) {}
      /**
       * amber topology
       */
      bool amber;
      /**
       * electrostatic interaction scaling for 1,4-interactions
       */
      double coulomb_scaling;
    } amber;

    struct dfunct_struct {

      dfunct_struct() : dfunct(dfunct_off), atom_i(0), atom_j(0), atom_k(0), atom_l(0), r_0(0.0), d(0), force(0.0) {}

      /**
       * dfunct enum 
       */
      dfunct_enum dfunct;
      /**
       * index of first atom involved in the potential
       */
      int atom_i, 
      /**
       * index of second atom involved in the potential
       */
      atom_j, 
      /**
       * index of third atom involved in the potential
       */
      atom_k, 
      /**
       * index of fourth atom involved in the potential
       */
      atom_l;
      /**
       * r_0 distance
       */
      double r_0;
      /**
       * addition or subtraction of distances r_ij and r_kl (can be scaled)
       */
      double d;
      /**
       * force constant of the bias
       */
      double force;
    } dfunct;
    
    /**
     A struct to mark parts of the code as "under development"
     */
    struct develop_struct {
    public:
      bool develop;
      std::string msg;
    } develop;

    /**
     set the development flag as true and specify the error message
     */
    void setDevelop(std::string s) {
      develop.develop = true;
      develop.msg = s;
    }

  };

}

#endif
