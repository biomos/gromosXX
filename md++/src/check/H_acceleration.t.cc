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

//#include <iostream>
//#include <iomanip>
//#include <map>
//#include <vector>

#include "../stdheader.h"
#include "check.h"

#include "../algorithm/integration/H_acceleration.cc"



void hard_coded_values(std::map<std::string, std::vector<double>> & m){
    m["PowerA_E"].assign({-25, 0, 11.6015625, 15, 18.3984375, 30, 55});
    m["PowerA_f"].assign({1, 1, 0.17968745, 0.125, 0.17968745, 1, 1});
    m["PowerA_E_90"].assign({-25.0, 0.0, 23.4375, 45, 66.5625, 90.0, 115.0});
    m["PowerA_f_90"].assign({1.0, 1.0, 0.8875, 0.85, 0.8875, 1.0, 1.0});
    
    m["InvGaussA_E"].assign({-56.25, -45.0, -36.25, -30.0, -26.25, -25.0});
    m["InvGaussA_f"].assign({0.5, 0.4, 0.3, 0.2, 0.1, 1.0});
    m["InvGaussA_E_90"].assign({-25.0, 0.0, 25.0, 50.0, 75.0, 90.0, 115.0});
    m["InvGaussA_f_90"].assign({1., 1., 1., 1., 1., 1.0});
    m["InvGaussA_Ef_90"].assign({87.5, 0.5});
    
    m["GaussA_E"].assign({126.25, 125.0, 126.25, 130.0, 136.25, 145.0});
    m["GaussA_f"].assign({-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5});
    m["GaussA_E_globalmin"].assign({162.5, 163.33333333333334, 165.83333333333334, 170.0, 175.83333333333331, 183.33333333333331});
    m["GaussA_f_globalmin"].assign({0.0, 0.06666666666666665, 0.1333333333333333, 0.19999999999999996, 0.2666666666666667, 0.33333333333333337});
    m["GaussA_E_90"].assign({-9.375, 10.0, 30.625, 52.5, 75.625, 100});
    m["GaussA_f_90"].assign({0.75, 0.8, 0.85, 0.9, 0.95, 1.0});
}

int main(){
    // set hard coded values to compare to
    std::map<std::string, std::vector<double>> ref_values;
    hard_coded_values(ref_values);
    double ref_value;
    double delta = 0.00001;
    int i, res=0, total=0;
    
    double E, E_a, E_a2, fk;
    double target_emax,target_emin, target_accel;
    double E_offset=0, E_scaling=1, E_step=25.;
    double E_offset_start=0, E_offset_end=201, E_offset_step=200;
    double E_scaling_start=1, E_scaling_end=3.1, E_scaling_step = 3;
    
    // ###################################### Power Acceleration test ######################################
    std::cout << "####### Power Acceleration test #######\n";
    algorithm::PowerAcceleration PA;
    
    std::cout << "\tcheck for target acceleration\n";
    target_emin = 0; target_emax=100;
    res = 0;
    for (target_accel=10; target_accel<100; target_accel+=10){
        PA.set_target_acceleration(target_emax, target_emin, target_accel);
        PA.accelerate_E_F(target_emin, &E_a, &fk);
        CHECK_APPROX_EQUAL(fk, 1, delta, res);
        PA.accelerate_E_F(target_emax, &E_a2, &fk);
        CHECK_APPROX_EQUAL(fk, 1, delta, res);
        CHECK_APPROX_EQUAL(E_a2 - E_a, target_accel, delta, res);
    }
    RESULT(res, total);
    
    target_emin = 0; target_emax=100; target_accel=90;
    std::cout << "\ttarget_emax: "<< target_emax<<
        "    target_emin: "<<target_emin <<
        "    target_accel: " << target_accel << "\n";
    for (E_scaling=E_scaling_start; E_scaling<E_scaling_end;E_scaling*=E_scaling_step) for (E_offset=E_offset_start; E_offset<E_offset_end; E_offset+=E_offset_step){
        E = (target_emin - E_step + E_offset) * E_scaling;
        res = 0;
        std::cout << "\t\tE_offset: "<< E_offset << "    scaling: " << E_scaling << "\n";
        PA.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        for (i=0;i<6;i++){
            PA.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["PowerA_E_90"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            CHECK_APPROX_EQUAL(fk, ref_values["PowerA_f_90"][i], delta, res);
            PA.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            E = E + E_step*E_scaling;
        }
        RESULT(res, total);
    }

    
    target_emin = 0; target_emax=100; target_accel=30;
    std::cout << "\ttarget_emax: "<< target_emax<<
        "    target_emin: "<<target_emin <<
        "    target_accel: " << target_accel << "\n";

    for (E_scaling=E_scaling_start; E_scaling<E_scaling_end;E_scaling*=E_scaling_step) for (E_offset=E_offset_start; E_offset<E_offset_end; E_offset+=E_offset_step){
        E = (target_emin - E_step + E_offset) * E_scaling;
        res = 0;
        std::cout << "\t\tE_offset: "<< E_offset << "    scaling: " << E_scaling << "\n";
        PA.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        for (i=0;i<6;i++){
            PA.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["PowerA_E"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            CHECK_APPROX_EQUAL(fk, ref_values["PowerA_f"][i], delta, res);
            PA.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            E = E + E_step*E_scaling;
        }
        RESULT(res, total);
    }
    
    // ###################################### Inverse Gaussian Acceleration test ######################################
    std::cout << "\n####### Inverse Gaussian Acceleration test #######\n";
    algorithm::InverseGaussianAcceleration IGA;

    std::cout << "\tcheck for target acceleration\n";
    target_emin = 0.; target_emax=100.;
    res = 0;
    for (target_accel=10.; target_accel<100; target_accel+=10){
        IGA.set_target_acceleration(target_emax, target_emin, target_accel);
        IGA.accelerate_E_F(target_emin, &E_a, &fk);
        IGA.accelerate_E_F(target_emax, &E_a2, &fk);
        CHECK_APPROX_EQUAL(double(E_a2 - E_a), target_accel, delta, res);
    }
    RESULT(res, total);

    target_emin = 0; target_emax=100; target_accel=90;
    std::cout << "\ttarget_emax: "<< target_emax<<
        "    target_emin: "<<target_emin <<
        "    target_accel: " << target_accel << "\n";

    for (E_scaling=E_scaling_start; E_scaling<E_scaling_end;E_scaling*=E_scaling_step) for (E_offset=E_offset_start; E_offset<E_offset_end; E_offset+=E_offset_step){
        E = (target_emin - E_step + E_offset) * E_scaling;
        res = 0;
        std::cout << "\t\tE_offset: "<< E_offset << "    scaling: " << E_scaling << "\n";
        IGA.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        for (i=0;i<6;i++){
            IGA.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["InvGaussA_E_90"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(fk, ref_values["InvGaussA_f_90"][i], delta, res);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            IGA.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            E = E + E_step*E_scaling;
        }
        E = (90+E_offset)*E_scaling;
        IGA.accelerate_E_F(E, &E_a, &fk);
        ref_value = (ref_values["InvGaussA_Ef_90"][0] + E_offset) * E_scaling;
        CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
        CHECK_APPROX_EQUAL(fk, ref_values["InvGaussA_Ef_90"][1], delta, res);
        
        RESULT(res, total);
    }
    
    
    target_emin = 0; target_emax=100; target_accel=20;
    std::cout << "\ttarget_emax: "<< target_emax<<
        "    target_emin: "<<target_emin <<
        "    target_accel: " << target_accel << "\n";

    for (E_scaling=E_scaling_start; E_scaling<E_scaling_end;E_scaling*=E_scaling_step) for (E_offset=E_offset_start; E_offset<E_offset_end; E_offset+=E_offset_step){
        E = (target_emin - E_step + E_offset) * E_scaling;
        res = 0;
        std::cout << "\t\tE_offset: "<< E_offset << "    scaling: " << E_scaling << "\n";
        IGA.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        for (i=0;i<6;i++){
            IGA.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["InvGaussA_E"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(fk, ref_values["InvGaussA_f"][i], delta, res);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            IGA.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            E = E + E_step*E_scaling;
        }
        RESULT(res, total);
    }

    // ###################################### Gaussian Acceleration test ######################################
    std::cout << "\n####### Gaussian Acceleration test #######\n";
    algorithm::GaussianAcceleration GA, GA_global_emin;
    std::cout << "\tcheck for target acceleration\n";
    target_emin = 0.; target_emax=100.;
    res = 0;
    for (target_accel=10.; target_accel<100; target_accel+=10){
        GA.set_target_acceleration(target_emax, target_emin, target_accel);
        GA.accelerate_E_F(target_emin, &E_a, &fk);
        GA.accelerate_E_F(target_emax, &E_a2, &fk);
        CHECK_APPROX_EQUAL(double(E_a2 - E_a), target_accel, delta, res);
    }
    RESULT(res, total);

    target_emin = 0; target_emax=100; target_accel=20;
    std::cout << "\ttarget_emax: "<< target_emax<<
        "    target_emin: "<<target_emin <<
        "    target_accel: " << target_accel << "\n";

    for (E_scaling=E_scaling_start; E_scaling<E_scaling_end;E_scaling*=E_scaling_step) for (E_offset=E_offset_start; E_offset<E_offset_end; E_offset+=E_offset_step){
        E = (target_emin - E_step + E_offset) * E_scaling;
        res = 0;
        std::cout << "\t\tE_offset: "<< E_offset << "    scaling: " << E_scaling << "\n";
        GA.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        GA_global_emin.set_global_emin(E);
        GA_global_emin.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        
        for (i=0;i<6;i++){
            GA.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["GaussA_E"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(fk, ref_values["GaussA_f"][i], delta, res);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            GA.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            
            GA_global_emin.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["GaussA_E_globalmin"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(fk, ref_values["GaussA_f_globalmin"][i], delta, res);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            GA_global_emin.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            
            E = E + E_step*E_scaling;
        }
        RESULT(res, total);
    }

    
    target_emin = 0; target_emax=100; target_accel=90;
    std::cout << "\ttarget_emax: "<< target_emax<<
        "    target_emin: "<<target_emin <<
        "    target_accel: " << target_accel << "\n";

    for (E_scaling=E_scaling_start; E_scaling<E_scaling_end;E_scaling*=E_scaling_step) for (E_offset=E_offset_start; E_offset<E_offset_end; E_offset+=E_offset_step){
        E = (target_emin - E_step + E_offset) * E_scaling;
        res = 0;
        std::cout << "\t\tE_offset: "<< E_offset << "    scaling: " << E_scaling << "\n";
        GA.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        GA_global_emin.set_global_emin(E);
        GA_global_emin.set_target_acceleration((target_emax+E_offset)*E_scaling, (target_emin+E_offset)*E_scaling, target_accel*E_scaling);
        
        for (i=0;i<6;i++){
            GA.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["GaussA_E_90"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(fk, ref_values["GaussA_f_90"][i], delta, res);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            GA.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            
            GA_global_emin.accelerate_E_F(E, &E_a, &fk);
            ref_value = (ref_values["GaussA_E_90"][i] + E_offset) * E_scaling;
            CHECK_APPROX_EQUAL(fk, ref_values["GaussA_f_90"][i], delta, res);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            GA_global_emin.accelerate_E(E, &E_a);
            CHECK_APPROX_EQUAL(E_a, ref_value, delta, res);
            
            E = E + E_step*E_scaling;
        }
        RESULT(res, total);
    }

    
    
    return total;
    
}

