/**
 * @file parameter.h
 * input parameters
 */

#ifndef INCLUDED_PARAMETER_H
#define INCLUDED_PARAMETER_H

namespace simulation
{
  /**
   * @class Parameter
   * input parameters.
   */
  class Parameter
  {
  public:
    /**
     * @struct system_struct
     * system block
     */
    struct system_struct
    {
      /**
       * Number of protein molecules
       */
      int npm;
      /** 
       * Number of solvent molecules
       */
      int nsm;
    } system;
    
    /**
     * @struct minimise_struct
     * minimise block
     */
    struct minimise_struct
    {
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
    } minimise;

    /**
     * @struct istart_struct
     * start block
     */
    struct start_struct
    {
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
    } start;

    /**
     * @struct step_struct
     * step block
     */
    struct step_struct
    {
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
    } step;

    /**
     * @struct boundary_struct
     * BOUNDARY block
     */
    struct boundary_struct
    {
      /**
       * NTB switch
       */
      math::boundary_enum boundary;
    } boundary;
 
    /**
     * @struct submolecules_struct
     * submolecules block
     */
    struct submolecules_struct
    {
      /**
       * Vector containing the last atom of every molecule
       */
      std::vector<size_t> submolecules;
    } submolecules;

    /**
     * multibath block
     */
    struct multibath_struct
    {
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
       */
      struct tcouple_struct
      {
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
      } tcouple;
      
      /**
       * have multibath
       */
      bool found_multibath;
      /**
       * have tcouple
       */
      bool found_tcouple;
    } multibath;
    
    
    /**
     * @struct pcouple_struct
     * PCOUPLE block
     */
    struct pcouple_struct
    {
      /**
       * default constructor
       */
      pcouple_struct()
      {
	scale=math::pcouple_off;
	calculate=false;
	virial=math::no_virial;
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
    } pcouple;

    /**
     * @struct centreofmass_struct
     * CENTREOFMASS block
     */
    struct centreofmass_struct
    {
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
      
    } centreofmass;

    /**
     * @struct print_struct
     * PRINT block
     */
    struct print_struct
    {
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
    } print;

    /**
     * @struct write_struct
     * WRITE block
     */
    struct write_struct
    {
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
    } write;

    /**
     * @struct shake_struct
     * SHAKE block
     */
    struct shake_struct
    {
      /**
       * NTC parameter (off=1, hydrogens=2, all=3, specified=4)
       * specified shakes everything in the constraint block in the topology.
       * hydrogens or all add the bonds containing hydrogens or all bonds to
       * the constraint block and shake those.
       */
      int ntc;
      /**
       * SHAKE tolerance
       */
      double tolerance;
    } shake;

    /**
     * @struct force_struct
     * FORCE block
     */
    struct force_struct
    {
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
      std::vector<size_t> energy_group;
    } force;

    /**
     * @struct plist_struct
     * PLIST block
     */
    struct plist_struct
    {
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
      
    } pairlist;
    /**
     * @struct longrange_struct
     * LONGRANGE block
     */
    struct longrange_struct
    {
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
      
      
    } longrange;

    /**
     * @struct posrest_struct
     * POSREST block
     */
    struct posrest_struct
    {
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
    } posrest;

    /**
     * @struct perturb_struct
     * PERTURB block
     */
    struct perturb_struct
    {
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
    } perturbation;
    
  };
  
}

#endif
