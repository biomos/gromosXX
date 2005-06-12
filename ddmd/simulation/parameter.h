/**
 * @file parameter.h
 * input parameters
 */

#ifndef INCLUDED_PARAMETER_H
#define INCLUDED_PARAMETER_H

#include "multibath.h"

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
  constr_flexshake
};
  
/**
 * @enum restr_enum
 * restraints enumeration
 */
enum restr_enum{
  /**
   * no restraints
   */
  restr_off = 0,
  /**
   * instantaneous restraints
   */
  restr_inst = 1,
  /**
   * time averaged restraints
   */
  restr_av = 2,
  /**
   * biquadratic (time averaged & instantaneous) restraints
   */
  restr_biq = 3,
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
  /** lj_crf_function */ lj_crf_func
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
     * - shake_pos           false  (no initial SHAKE of positions)
     * - shake_vel           false  (no initial SHAKE of velocities)
     * - remove_com          false  (no initial removal of COM motion)
     * - generate_velocities false  (no generation of initial velocities)
     * - ig                      0  (random number seed)
     * - tempi                 0.0  (temperature to generate initial velocities)
     */
    start_struct() : shake_pos(false), shake_vel(false), remove_com(false),
		     generate_velocities(false), ig(0), tempi(0.0) {}
      
    /**
     * shake initial positions
     */
    bool shake_pos;
    /**
     * shake initial velocities
     */
    bool shake_vel;
    /**
     * COM removal.
     */
    bool remove_com;
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
     * - boundary vacuum
     */
    boundary_struct() : boundary(vacuum) {}
      
    /**
     * NTB switch
     */
    boundary_enum boundary;
  } /** boundary parameters */ boundary;
 
  /**
   * @struct submolecules_struct
   * submolecules block
   */
  struct submolecules_struct
  {
    /**
     * Constructor
     * Default values:
     * - submolecules empty (if it stays empty, all atoms are considered to be in one submolecule)
     */
    submolecules_struct() {}
      
    /**
     * Vector containing the last atom of every molecule
     */
    std::vector<unsigned int> submolecules;
  } /** submolecule array */ submolecules;

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
      scale=pcouple_off;
      calculate=false;
      virial=no_virial;
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
    pressure_scale_enum scale;
    /**
     * virial type
     */
    virial_enum virial;
    /**
     * reference pressure
     */
    Matrix pres0;
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
     * - ndfmin 0            (number of degrees of freedom to substract)
     * - skip_step 0         (number of steps to skip between removal of com motion)
     * - remove_rot false    (remove center of mass rotation)
     * - remove_trans false  (remove center of mass translation)
     */
    centreofmass_struct() : ndfmin(0), skip_step(0), remove_rot(false), remove_trans(false) {}
      
    /**
     * Number of degrees of freedom to substract
     */
    int ndfmin;
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
     * - stepblock 0               (print out every step)
     * - centreofmass 0            (print centre of mass information every step)
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
     * - energy   0        (write energy trajectory)
     * - free_energy 0     (write energy lambda derivative trajectory)
     * - block_average 0   (write block averaged energy trajectories)
     * - solute_only false (write solute and solvent)
     */
    write_struct() : position(0), velocity(0), energy(0), free_energy(0), 
		     block_average(0), solute_only(false) {}
      
    /**
     * position.
     */
    int position;
    /**
     * velocity.
     */
    int velocity;
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
     * write solute only trajectory
     */
    bool solute_only;
      
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
     * - spc_loop -1
     */
    force_struct() : bond(1), angle(1), improper(1),
		     dihedral(1), nonbonded(1), spc_loop(-1),
		     interaction_function(lj_crf_func)
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
     * nonbonded?
     */
    int nonbonded;
    /**
     * Energy groups
     */
    std::vector<unsigned int> energy_group;
    /**
     * fast spc loops
     */
    int spc_loop;

    /**
     * nonbonded interaction function
     */
    interaction_func_enum interaction_function;

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
     * - grid false
     * - skip_step 5
     * - cutoff_short 0.8
     * - cutoff_long 1.4
     * - grid_size 0.4
     * - atomic_cutoff false
     */
    plist_struct() : grid(false), skip_step(5), cutoff_short(0.8),
		     cutoff_long(1.4), grid_size(0.4), atomic_cutoff(false) {}
      
    /**
     * algorithm.
     */
    bool grid;
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
      
  } /** Pairlist method parameters */ pairlist;

  /**
   * @struct longrange_struct
   * LONGRANGE block
   */
  struct longrange_struct
  {
    /**
     * Constructor
     * Default values:
     * - rf_epsilon 66 (spc water)
     * - rf_kappa    0
     * - rf_cutoff 1.4
     * - rf_excluded false (old standard)
     * - epsilon     1
     */
    longrange_struct() : rf_epsilon(66), rf_kappa(0), rf_cutoff(1.4),
			 rf_excluded(false), epsilon(1) {}
      
    /**
     * Reaction field epsilon
     */
    double rf_epsilon;
    /**
     * Reaction field Kappa
     */
    double rf_kappa;
    /**
     * Reaction field cutoff
     */
    double rf_cutoff;
    /**
     * include rf contributions from excluded atoms
     */
    bool rf_excluded;
    /**
     * Epsilon 1 within the cutoff.
     * in GROMOS this is hardcoded to be 1;
     * we do so in In_Parameter
     */
    double epsilon;
      
  } /** longrange treatment parameters */ longrange;

  /**
   * @struct posrest_struct
   * POSREST block
   */
  struct posrest_struct
  {
    /** 
     * Constructor
     * Default values:
     * - posrest 0 (no position restraints)
     * - nrdrx true
     * - force_constant 10000
     */
    posrest_struct() : posrest(0), nrdrx(true), force_constant(1E4) {}
      
    /**
     * posrest
     */
    int posrest;
    /**
     * NRDRX
     */
    bool nrdrx;
    /**
     * CHO
     */
    double force_constant;
  } /** Position restraint parameters */ posrest;

  /**
   * @struct distrest_struct
   * DISTREST block
   */

  struct distrest_struct
  {
    /**
     * Constructor
     * Default values:
     * - distrest 0 (no distance restraints)
     * - K 0
     * - r_linear 0
     * - tau 10
     * - read 0
     */

    distrest_struct()
      : distrest(0),
	K(0),
	r_linear(0),
	tau(10),
	read(0)
    {
    }
      
    /** 
     * distance restraints on/off
     */
    int distrest;
      
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
      
  }/** Distance restraints parameters */ distrest;
    

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
     */
    perturb_struct() : perturbation(false), lambda(0), lambda_exponent(1),
		       dlamt(0), scaling(false), scaled_only(false) {}
      
    /**
     * perturbation?
     */
    bool perturbation;
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
     * - tau 0
     * - read_av false
     */
    jvalue_struct()
      : mode(restr_off),
	tau(0.0),
	read_av(false)
    {
    }
    /**
     * restraining mode.
     */
    restr_enum mode;
    /**
     * coupling time.
     */
    double tau;
    /**
     * read averages.
     */
    bool read_av;
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
     */
    rottrans_struct() : rottrans(false)
    {
    }
    /**
     * apply rotational translational constraints?
     */
    bool rottrans;
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
     * - ID 0
     * - T 0
     * - scale 0
     */
    replica_struct() : number(0), ID(0), T(0.0), scale(0), trials(0)
    {
    }
    /**
     * num
     */
    int number;
    /**
     * job ID
     */
    int ID;
    /**
     * @deprecated
     * temperature
     */
    double T;
    /**
     * @deprecated
     * scale initial velocities to T?
     */
    int scale;
    /**
     * tempeartures
     */
    std::vector<double> temperature;
    /**
     * lambdas
     */
    std::vector<double> lambda;

    /**
     * trial moves
     */
    int trials;
      
  } /** replica exchange parameters */ replica;
    
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
     * - trajectory("")
     */
    analyze_struct() : analyze(false), trajectory("")
    {
    }
    /** 
     * re-analyze trajectory
     */
    bool analyze;
    /**
     * trajectory filename
     */
    std::string trajectory;
      
  } /** analyze parameter */ analyze;

};

#endif