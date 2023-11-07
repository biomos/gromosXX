/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>

#ifndef Haccel_H
#define Haccel_H
namespace algorithm
{
  class EnergyAcceleration{
    public:
      /** updates acceleration parameters/variables*/
      virtual void update_params() = 0;
      
      // calculate accelerated energy
      virtual void accelerate_E(double E, double * E_a) = 0;
      
      // calculate accelerated energy & scaling factor for the force
      virtual void accelerate_E_F(double E, double * E_a, double * f_k) = 0;
      
      /**
      * sets acceleration range (between `accel_emin` and `accel_emax`)
      * and `accel_target`, which defines the acceleration level
      * and updates the acceleration parameters/variables accordingly
      */
      virtual void set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value) = 0;
  };
  
  class PowerAcceleration:EnergyAcceleration{

    public:
      // current emin, emax & other varibles
      bool valid_params;
      double emin, emax, pow_frac, emax_a; // emax_a is the value of E at emax
      int pow2;
      // target acceleration vars
      double accel_emin, accel_emax, accel_target;

      /**
      * Constructor.
      */
      PowerAcceleration() {emin=0.; emax=0.; pow2=2; pow_frac=1.; valid_params=false;}

      /** updates power acceleration variables (e.g. target_frac, diff_emm, gamma_pow, ...)*/
      void update_params();

      /**
      * calculate accelerated energy
      * @param E energy to be accelerated
      * @param E_a pointer to a variable to assign the value of the accelerated energy
      */
      void accelerate_E(double E, double * E_a);

      /**
      * calculate accelerated energy & scaling factor for the force
      * @param E energy to be accelerated
      * @param E_a pointer to a variable to assign the value of the accelerated energy
      * @param f_k pointer to a variable to assign the value of the force scaling factor
      */
      void accelerate_E_F(double E, double * E_a, double * f_k);

      /**
      * sets acceleration range (between `accel_emin` and `accel_emax`)
      * and `accel_target`, which defines the acceleration level
      * @param accel_emax_value target emax value
      * @param accel_emin_value target emin value
      * @param accel_target_value target acceleration level (accelerated energy at emax)
      */
      void set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value);
    
      /**
      * Destructor.
      */
      virtual ~PowerAcceleration(){}

    private:
      // all vars needed for power acceleration
      double diff_emm;
      double lin_frac, pow_c;
      int pow2_1;
      // target acceleration var
      double target_frac;

      /**
      * transforms energy value into x variable for power acceleration
      * @param E
      */
      double transform_E_2_x(double E);

      /**
      * transforms f(x) value back to energy (for power acceleration)
      * @param E
      */
      double transform_fx_2_E(double f_x);

      /** 
      * calculate f(x) needed to get H*
      * @param x 
      * @return f(x)
      */
      double get_f_x_lin(double x);

      /**
      * calculate both f(x) and f'(x) at the same time
      * @param x 
      * @param f_x pointer to a double to assign the value of f(x)
      * @param f_der_x pointer to a double to assign the value of f'(x)
      */
      void get_f_f_der_x_lin(double x, double *f_x, double *f_der_x);

    /////// for target accel
      /**finds the power needed to fulfill the target_frac */
      void find_pow2();

      /**
      * calculates the linear (and the pow term) scaling factor to get the target_frac
      * sets lin_frac, pow_frac and pow2 based on target_frac
      */
      void find_lin_scaling_frac();

      /**get the target_frac based on emin, emax and target_accel*/
      void get_target_frac();

      /**updates the acceleration variables based on target acceleration*/
      void update_target_params();
  };


  class InverseGaussianAcceleration:EnergyAcceleration{

    public:
      double emin, emax;
      bool valid_params;
      // target acceleration vars
      double accel_emin, accel_emax, accel_target;

      /**
      * Constructor.
      */
      InverseGaussianAcceleration() {emin=0.; emax=0.; valid_params=false;}

      /**
      * Destructor.
      */
      virtual ~InverseGaussianAcceleration(){}

      /**
      * calculate accelerated energy
      * @param E energy to be accelerated
      * @param E_a pointer to a variable to assign the value of the accelerated energy
      */
      void accelerate_E(double E, double * E_a);

      /**
      * calculate accelerated energy & scaling factor for the force
      * @param E energy to be accelerated
      * @param E_a pointer to a variable to assign the value of the accelerated energy
      * @param f_k pointer to a variable to assign the value of the force scaling factor
      */
      void accelerate_E_F(double E, double * E_a, double * f_k);

      /** updates acceleration variables (e.g. diff_emm)*/
      void update_params();

      /**
      * sets acceleration range (between `accel_emin` and `accel_emax`)
      * and `accel_target`, which defines the acceleration level
      * @param accel_emax_value target emax value
      * @param accel_emin_value target emin value
      * @param accel_target_value target acceleration level (accelerated energy at emax)
      */
      void set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value);
    
    private:
      // all vars needed for invGauss acceleration
      double diff_emm;
      // target acceleration var
      double target_frac;

      void update_target_params();
      void get_emin();
  };



  class GaussianAcceleration:EnergyAcceleration{

    public:
      double emin, emax, k0;
      bool valid_params;
      // target acceleration vars
      double accel_emin, accel_emax, accel_target;

      /**
      * Constructor.
      */
      GaussianAcceleration() {emin=0.; emax=0.; k0=1., valid_params=false;flag_global_emin=false;}

      /**
      * Destructor.
      */
      virtual ~GaussianAcceleration(){}

      /**
      * calculate accelerated energy
      * @param E energy to be accelerated
      * @param E_a pointer to a variable to assign the value of the accelerated energy
      */
      void accelerate_E(double E, double * E_a);

      /**
      * calculate accelerated energy & scaling factor for the force
      * @param E energy to be accelerated
      * @param E_a pointer to a variable to assign the value of the accelerated energy
      * @param f_k pointer to a variable to assign the value of the force scaling factor
      */
      void accelerate_E_F(double E, double * E_a, double * f_k);

      /** updates acceleration variables (e.g. diff_emm)*/
      void update_params();

      /**
      * sets acceleration range (between `accel_emin` and `accel_emax`)
      * and `accel_target`, which defines the acceleration level
      * @param accel_emax_value target emax value
      * @param accel_emin_value target emin value
      * @param accel_target_value target acceleration level (accelerated energy at emax)
      */
      void set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value);
      
      // set global_emin value
      void set_global_emin(double global_emin_value);

      // unset global_emin value
      void unset_global_emin();

    private:
      // all vars needed for invGauss acceleration
      double diff_emm, global_emin;
      // target acceleration var
      double target_frac;
      bool flag_global_emin;
      
      void update_target_params();
      void get_emin();
  };

}
#endif
