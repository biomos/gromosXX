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

#include "../../stdheader.h"
#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../util/error.h"

#include "H_acceleration.h"


// PowerAcceleration
/** sets acceleration parameters/variables
if aparam_form==1
    params hold: {emax, emin, pow, pow_frac}
*/
void algorithm::PowerAcceleration::set_params(int aparam_form, std::vector<double> params){
    if (aparam_form==1){
        emax = params[0];
        emin = params[1];
        pow2 = params[2];
        assert((pow2 - int(pow2) == 0));
        pow_frac = params[3];
        update_params();
    }
}

/** gets acceleration parameters/variables
{emax, emin, pow, pow_frac}
*/
void algorithm::PowerAcceleration::get_params(std::vector<double> params){
    params.push_back(emax);
    params.push_back(emin);
    params.push_back(pow2);
    params.push_back(pow_frac);
}

// updates power acceleration variables (e.g. diff_emm, pow_c, ...)
void algorithm::PowerAcceleration::update_params(){
    if (emax > emin){
        valid_params = true;
        pow2_1 = pow2 + 1;
        pow_c = 1. / (2*pow2_1);
        diff_emm = emax - emin;
        lin_frac = 1. - pow_frac;
        double f_x = get_f_x_lin(1);
        emax_a = transform_fx_2_E(f_x);
    }
    else
        {valid_params = false;}
}

/**
 * transforms energy value into x variable for the acceleration
 * @param E
 */
inline double algorithm::PowerAcceleration::transform_E_2_x(double E){
    return (E - emin) / diff_emm;
}

/**
 * transforms f(x) value back to energy (for power acceleration)
 * @param E
 */
inline double algorithm::PowerAcceleration::transform_fx_2_E(double f_x){
    return emin + f_x * diff_emm;
}


/** 
 * calculate f_lin(x) needed to get H*
 * @param x 
 * @return f_lin(x)
 */
inline double algorithm::PowerAcceleration::get_f_x_lin(double x){
    double f_pow = pow(2.0, pow2) * pow(x - 0.5, pow2_1) / pow2_1 + pow_c;
    return pow_frac * f_pow + lin_frac * x;
}

/**
 * calculate both f_lin(x) and f_lin'(x) at the same time
 * @param x 
 * @param f_x pointer to a double to assign the value of f_lin(x)
 * @param f_der_x pointer to a double to assign the value of f_lin'(x)
 */
inline void algorithm::PowerAcceleration::get_f_f_der_x_lin(double x, double *f_x, double *f_der_x){
    double f_der_pow = pow(2 * x - 1, pow2);
    double f_pow = f_der_pow * (x - 0.5) / pow2_1 + pow_c;
    *f_x = pow_frac * f_pow + lin_frac * x;
    *f_der_x = pow_frac * f_der_pow + lin_frac;
}

/**
 * calculate accelerated energy (power acceleration)
 * @param E_a pointer to a variable to assign the value of the accelerated energy
 * @param E energy to be accelerated
 */
void algorithm::PowerAcceleration::accelerate_E(double E, double *E_a){
    *E_a = E;
    if (valid_params){
        if (E <= emin)
            *E_a = E;
        else if (E >= emax)
            *E_a = emax_a + E - emax;
        else{
            double x = transform_E_2_x(E);
            double f_x = get_f_x_lin(x);
            *E_a = transform_fx_2_E(f_x);
        }
    }
}

/**
 * calculate accelerated energy & scaling factor for the force (power acceleration)
 * @param E_a pointer to a variable to assign the value of the accelerated energy
 * @param f_der_x pointer to a variable to assign the value of the force scaling factor
 * @param E energy to be accelerated
 */
void algorithm::PowerAcceleration::accelerate_E_F(double E, double *E_a, double *f_der_x){
    *f_der_x = 1.;
    *E_a = E;
    if (valid_params){
        if (E <= emin)
            *E_a = E;
        else if (E >= emax)
            *E_a = emax_a + E - emax;
        else{
            double x = transform_E_2_x(E);
            double f_x;
            get_f_f_der_x_lin(x, &f_x, f_der_x);
            *E_a = transform_fx_2_E(f_x);
        }
    }
}


// finds the gamma power needed to fulfill the target_frac
inline void algorithm::PowerAcceleration::find_pow2(){
    double temp_pow = (1. / target_frac - 1) / 2;
    pow2 =  2*ceil(temp_pow);
    pow2_1 = pow2 + 1;
    pow_c = 1. / (2*pow2_1);
}

/**
 * calculates the linear (and the gamma term) scaling factor to get the target_frac
 * sets lin_frac, gamma_frac and gamma_pow based on target_frac
 */
inline void algorithm::PowerAcceleration::find_lin_scaling_frac(){
    double pow_frac_min = 1. / pow2_1;
    lin_frac = (target_frac - pow_frac_min) / (1. - pow_frac_min);
    pow_frac = 1. - lin_frac;
}

/**get the target_frac based on emin, emax and accel_target*/
inline void algorithm::PowerAcceleration::get_target_frac(){
    if (emax == emin)
      target_frac = 1.;
    else
      target_frac = accel_target / (emax - emin);
}

/**
* sets acceleration range (between `accel_emin` and `accel_emax`)
* and `accel_target`, which defines the acceleration level
* @param accel_emax_value target emax value
* @param accel_emin_value target emin value
* @param accel_target_value target acceleration level (accelerated energy at emax)
*/
void algorithm::PowerAcceleration::set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value){
    accel_emin = accel_emin_value;
    accel_emax = accel_emax_value;
    accel_target = accel_target_value;
    update_target_params();
}

/**updates the acceleration variables based on target acceleration*/
void algorithm::PowerAcceleration::update_target_params(){
    if (0 < accel_target && accel_target < accel_emax - accel_emin){
        valid_params=true;
        emin = accel_emin;
        emax = accel_emax;
        diff_emm = emax - emin;
        emax_a = accel_target + emin;
        get_target_frac();
        find_pow2();
        find_lin_scaling_frac();    
    }
    else{
        valid_params=false;
    }
}



// Inverse Gaussian Acceleration

/** sets acceleration parameters/variables
if aparam_form==1
    params hold: {emax, emin}
*/
void algorithm::InverseGaussianAcceleration::set_params(int aparam_form, std::vector<double> params){
    if (aparam_form==1){
        emax = params[0];
        emin = params[1];
        update_params();
    }
}

/** gets acceleration parameters/variables
{emax, emin}
*/
void algorithm::InverseGaussianAcceleration::get_params(std::vector<double> params){
    params.push_back(emax);
    params.push_back(emin);
}

// updates Inverse Gaussian Acceleration variables
void algorithm::InverseGaussianAcceleration::update_params(){
    if (emax > emin){
        valid_params = true;
        diff_emm = emax - emin;
    }
    else
        {valid_params = false;}
}


/**
* calculate accelerated energy & scaling factor for the force
* @param E energy to be accelerated
* @param E_a pointer to a variable to assign the value of the accelerated energy
* @param f_k pointer to a variable to assign the value of the force scaling factor
*/
void algorithm::InverseGaussianAcceleration::accelerate_E_F(double E, double * E_a, double * f_k){
    *E_a = E;
    *f_k = 1.;
    if (valid_params){
        // if (E <= emin)*E_a = E; // not needed
        if (E >= emax)
            *E_a = E - 0.5 * diff_emm;
        else if (E > emin){
            double demix = E - emin;
            double kfac = 1.0 / diff_emm;
            *f_k = 1.0 - kfac * demix;
            *E_a = E - 0.5 * kfac * demix * demix;
        }
    }
}

/**
* calculate accelerated energy
* @param E energy to be accelerated
* @param E_a pointer to a variable to assign the value of the accelerated energy
*/
void algorithm::InverseGaussianAcceleration::accelerate_E(double E, double * E_a){
    *E_a = E;
    if (valid_params){
        //if (E <= emin) *E_a = E; // this is not needed
        if (E >= emax)
            *E_a = E - 0.5 * diff_emm;
        else if (E > emin){
            double demix = E - emin;
            double kfac = 1.0 / diff_emm;
            *E_a = E - 0.5 * kfac * demix * demix;
        }
    }
}

/**
* sets acceleration range (between `accel_emin` and `accel_emax`)
* and `accel_target`, which defines the acceleration level
* @param accel_emax_value target emax value
* @param accel_emin_value target emin value
* @param accel_target_value target acceleration level (accelerated energy at emax)
*/
void algorithm::InverseGaussianAcceleration::set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value){
    accel_emin = accel_emin_value;
    accel_emax = accel_emax_value;
    accel_target = accel_target_value;
    update_target_params();
}

/**updates the acceleration variables based on target acceleration*/
void algorithm::InverseGaussianAcceleration::update_target_params(){
    if (0 < accel_target && accel_target < accel_emax - accel_emin){
        valid_params=true;
        emax = accel_emax;
        get_emin();
    }
    else{
        valid_params=false;
    }
}

// adjusts emin to target acceleration
void algorithm::InverseGaussianAcceleration::get_emin(){
    emin = 2 * (accel_target + accel_emin) - emax;
    if (emin < accel_emin){
        double temp = - emax*emax + 2*emax*accel_target + 2*emax*accel_emin - accel_emin*accel_emin;
        emin = temp / (2*accel_target);
    }
    diff_emm = emax - emin;
}
    


// Gaussian Acceleration

/** sets acceleration parameters/variables
if aparam_form==1
    params hold: {emax, k}
if aparam_form==2
    params hold: {emax, emin, k0}
*/
void algorithm::GaussianAcceleration::set_params(int aparam_form, std::vector<double> params){
    if (aparam_form==1){
        emax = params[0];
        k = params[1];
        valid_params = true;
    }
    if (aparam_form==2){
        emax = params[0];
        emin = params[1];
        k0 = params[2];
        update_params();
    }
}

/** gets acceleration parameters/variables
{emax, k}
*/
void algorithm::GaussianAcceleration::get_params(std::vector<double> params){
    params.push_back(emax);
    params.push_back(k);
}

// updates Gaussian Acceleration variables
void algorithm::GaussianAcceleration::update_params(){
    if (emax > emin){
        valid_params = true;
        diff_emm = emax - emin;
        k = k0 / diff_emm;
    }
    else
        {valid_params = false;}
}

/**
* calculate accelerated energy
* @param E energy to be accelerated
* @param E_a pointer to a variable to assign the value of the accelerated energy
*/
void algorithm::GaussianAcceleration::accelerate_E(double E, double * E_a){
    *E_a = E;
    if (valid_params){
        if (E < emax){
            double temp_dE = emax - E;
            *E_a = E + 0.5* k * temp_dE*temp_dE;
        }
    }
}

/**
* calculate accelerated energy & scaling factor for the force
* @param E energy to be accelerated
* @param E_a pointer to a variable to assign the value of the accelerated energy
* @param f_k pointer to a variable to assign the value of the force scaling factor
*/
void algorithm::GaussianAcceleration::accelerate_E_F(double E, double * E_a, double * f_k){
    *E_a = E;
    *f_k = 1.;
    if (valid_params){
        if (E < emax){
            double temp_dE = emax - E;
            *f_k = k  * temp_dE;
            *E_a = E + *f_k * 0.5 * temp_dE;
            *f_k = 1 - *f_k;
        }
    }
}

/**
* sets acceleration range (between `accel_emin` and `accel_emax`)
* and `accel_target`, which defines the acceleration level
* @param accel_emax_value target emax value
* @param accel_emin_value target emin value
* @param accel_target_value target acceleration level (accelerated energy at emax)
*/
void algorithm::GaussianAcceleration::set_target_acceleration(double accel_emax_value, double accel_emin_value, double accel_target_value){
    accel_emin = accel_emin_value;
    accel_emax = accel_emax_value;
    accel_target = accel_target_value;
    update_target_params();
}

/**updates the acceleration variables based on target acceleration*/
void algorithm::GaussianAcceleration::update_target_params(){
    if (0 < accel_target && accel_target < accel_emax - accel_emin){
        valid_params=true;
        emax = accel_emax;
        if (flag_global_emin)
            emin = global_emin;
        else
            emin = accel_emin;
        diff_emm = emax - emin;
        if (!flag_global_emin && accel_target >= diff_emm/2){
            k0 = 2*(diff_emm - accel_target) / diff_emm;
        }
        else{
            k0=1.; // reset k0 to 1 to ensure one always gets the same result
            double non_acell_d_emm = accel_emax - accel_emin; // normal difference
            double temp_fac = 2*(accel_target - non_acell_d_emm) / k0;
            double temp = emin * temp_fac + accel_emax*accel_emax - accel_emin*accel_emin;
            double temp_emax = temp / (temp_fac + 2*non_acell_d_emm);
            if (temp_emax > emax){
                emax = temp / (temp_fac + 2*non_acell_d_emm);
                diff_emm = emax - emin;
            }
            else{
                temp = 2*(non_acell_d_emm - accel_target) * diff_emm;
                k0 = temp / ((emax - accel_emin)*(emax - accel_emin));
            }
        }
        k = k0 / diff_emm;
    }
    else{
        valid_params=false;
    }
}

// set global_emin value
void algorithm::GaussianAcceleration::set_global_emin(double global_emin_value){
    flag_global_emin=true;
    global_emin = global_emin_value;
}

// unset global_emin value
void algorithm::GaussianAcceleration::unset_global_emin(){
    flag_global_emin=false;
}


// AccelerationContainer
void algorithm::AccelerationContainer::set_accel(std::string accel_name, int accel_type, int aparam_form, std::vector<double> params){
    switch (accel_type){
        case 1:{
            algorithm::InverseGaussianAcceleration temp_accel;
            if (aparam_form!=0){
                temp_accel.set_params(aparam_form, params);
            }
            IGAs.push_back(temp_accel);
            accel_map[accel_name] = &(IGAs.back());
            break;
        }
        case 2:{
            algorithm::GaussianAcceleration temp_accel;
            if (aparam_form!=0){
                temp_accel.set_params(aparam_form, params);
            }
            GAs.push_back(temp_accel);
            accel_map[accel_name] = &(GAs.back());
            break;
        }
        case 3:{
            algorithm::PowerAcceleration temp_accel;
            if (aparam_form!=0){
                temp_accel.set_params(aparam_form, params);
            }
            PAs.push_back(temp_accel);
            accel_map[accel_name] = &(PAs.back());
            break;
        }
    }
    if (aparam_form==0)
        (*accel_map[accel_name]).set_target_acceleration(params[0], params[1], params[2]);
}
