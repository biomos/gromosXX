/* 
 * File:   eds.h
 * Author: haniels
 *
 * Created on August 2, 2011, 1:41 PM
 */

#ifndef EDS_H
#define	EDS_H
namespace algorithm
{

  /**
   * @class EDS
   * implements EDS.
   */
  class EDS : public Algorithm
  {
  public:
    // stores current emin, emax, target_emax
    double emin, emax, target_emax;

    /**
     * Constructor.
     */
    EDS() : Algorithm("EDS"), conf2(NULL) {emin=0.; emax=0.; target_emax=0.;}
    
    void set_conf2(configuration::Configuration & conf) {
      conf2 = &conf;
    }
    
    // functions for Gaussian acceleration
      // calculate accelerated energy & scaling factor for the force (gaussian acceleration)
      void accelerate_E_F_gauss(double *E_a, double *f_der_x, double E, double E_offset/*=0*/);

      // calculate accelerated energy (gaussian acceleration)
      void accelerate_E_gauss(double *E_a, double E, double E_offset/*=0*/);
    

    // functions needed for gamma acceleration
      // sets emin, emax and target_emax & updates other gamma acceleration variables (e.g. target_frac, diff_emm, gamma_pow, ...)
      void set_EDS_params(double emin_value, double emax_value, double target_emax_value);

      // updates gamma acceleration variables (e.g. target_frac, diff_emm, gamma_pow, ...)
      void update_gamma_params();

      // transforms energy value into x variable for gamma acceleration
      double transform_E_2_x(double E);
  
      // calculate accelerated energy & scaling factor for the force (gamma acceleration)
      void accelerate_E_F_gamma(double *E_a, double *f_der_x, double E, double E_offset/*=0*/);

      // calculate accelerated energy (gamma acceleration)
      void accelerate_E_gamma(double *E_a, double E, double E_offset/*=0*/);


    /**
     * Destructor.
     */
    virtual ~EDS(){}
    
    /**
     * 
     */
    virtual int apply(topology::Topology &topo, 
		      configuration::Configuration &conf,
		      simulation::Simulation &sim);

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) 
    
    {
      if (!quiet)
	os << "\tEDS\nEND\n";
      return 0;
    };

   private:
     configuration::Configuration * conf2;


     // all vars needed for gamma acceleration
     double target_frac, diff_emm;
     double lin_frac, gamma_frac, f_gamma, f_der_gamma, gamma_c;
     int gamma_pow, gamma_pow_1;

     // all functions needed for gamma acceleration
      // finds the gamma power needed to fulfill the target_frac
      void find_gamma_pow();

      //calculates the linear (and the gamma term) scaling factor to get the target_frac
      // sets lin_frac, gamma_frac and gamma_pow based on target_frac
      void find_lin_scaling_frac();

      // calculates f(x) for a given gamma power
      // f(x) i the scaling factor for the energy (acceleration)
      // where x = (H - Emin) / (Emax - Emin)
      // gamma_pow has to be 0, 2, 4, 6...
      void f_k_pow(double x);

      // calculates f'(x) for a given gamma power
      // f'(x) i the scaling factor for the force (acceleration)
      // where x = (H - Emin) / (Emax - Emin)
      // gamma_pow has to be 0, 2, 4, 6...
      void f_der_k_pow(double x);
      
      // the 2 functions from above in one (optimized version)
      void f_f_der_k_pow(double x);

      // calculate f(x) and f'(x) needed to get H* and F*
      // f(x)
      double get_f_x_lin(double x);
      // f'(x)
      double get_f_der_x_lin(double x);

      // do both f(x) and f'(x) at the same time
      void get_f_f_der_x_lin(double x, double *f_x, double *f_der_x);

  // transition from f(x) to accelerated H
      void get_target_frac();
  
  };
   
} // algorithm

#endif



