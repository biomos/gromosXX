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

#ifndef INCLUDED_PHYSICS_H
#define INCLUDED_PHYSICS_H

#include <string>

namespace gmath {

  /**
   * The PhysConst class holds and handles the physical constants used in
   * GROMOS++.
   * 
   * GROMOS is basically unitless which means you can feed input with any
   * physical units an get the output with the corresponding units. Some values
   * may be derived by other (fundamental) constants which are defined  in the
   * PHYSICALCONSTANTS block of the molecular topology file. These are:
   * 
   *   - 1 / (4 pi epsilon_0) with epsilon_0 the permittivity of vacuum
   *   - the Planck constant
   *   - the speed of light
   *   - the Boltzmann constant
   * 
   * Therefore, no hard coded physical constants must exist. However, the class
   * is initialised using values from literature [1] in specified units and 
   * replaces these values if the variables are reset using the set_variable()
   * functions below (e.g. when reading the constants from a topology).
   * 
   * Every constant has a corresponding bool variable which indicates the state
   * of each physical variable (initialised from literature [1] or reset, e.g.
   * when read from topology). The accessors of each variable use these bool
   * variables to print a warning message on the screen in case the returned
   * physical constant was not read from the topology file or reset on purpose
   * in the program code. This makes it still possible to do some calculations
   * without a topology file (using the values from literature [1]) but warns
   * the user to avoid any confusion when doing so.
   *
   * 
   * References:
   *
   * [1] Mills, I. and Cvitas, N.K. and Homann, K. and Kuchitsu, K. "Quantities,
   *     Units and Symbols in Physical Chemsitry" Blackwell Science Ltd,
   *     <b>1993</b>, London.
   *
   * @class PhysConst
   * @author A.P. Eichenberger
   * @ingroup gmath
   */
  class PhysConst {

  public:

    // ====================================
    //    CONSTRUCTIR(S) & DESTUCTOR(S):
    // ====================================
    //
    /**
     * Cunstructor, inertialises all values with values from literature [1]
     */
    PhysConst();


    // ================
    //    ACCESSORS:
    // ================
    //
    // Numerical Constants:
    // --------------------
    /**
     * Accessor for the circle constant pi
     */
    double get_pi();
    /**
     * Accessor for the kilo variable
     */
    double get_kilo();
    /**
     * Accessor for the nano variable
     */
    double get_nano();
    /**
     * Accessor for the pico variable
     */
    double get_pico();
    /**
     * Accessor for the Euler number
     */
    double get_euler();
    /**
     * Accessor for the radian2degree converting constant
     */
    double get_radian2degree();
    /**
     * Accessor for the degree2radian converting constant
     */
    double get_degree2radian();
    /**
     * Accessor for the small number
     */
    double get_epsilon();
    //
    // Physical Constants:
    // -------------------
    /**
     * Accessor for the atomic mass unit
     */
    double get_atomic_mass_unit();
    /**
     * Accessor for the Avogadro number
     */
    double get_avogadro();
    /**
     * Accessor for the Boltzmann constant
     */
    double get_boltzmann();
    /**
     * Accessor for the Boltzmann constant without 
     * warning for cases without topology
     */
    double get_boltzmann_silent();
    /**
     * Accessor for the elementary charge constant
     */
    double get_elementary_charge();
    /**
     * Accessor for the dielectric permittivty of vacuum
     */
    double get_eps0();
    /**
     * Accessor for the four_pi_eps_i variable
     */
    double get_four_pi_eps_i();
    /**
     * Accessor for the Planck constant
     */
    double get_h();
    /**
     * Accessor for hbar
     */
    double get_hbar();
    /**
     * Accessor for the permeability of free space
     */
    double get_mu0();
    /**
     * Accessor for the speed of light
     */
    double get_speed_of_light();


    // =======================
    //    MEMBER FUNCTIONS:
    // =======================
    //
    /**
     * Sets the speed of light constant to the indicated value and turns off the
     * printWarning function when returning this variable (get_speed_of_light).
     * All dependant variables are recalculated again.
     */
    void set_speed_of_light(double value);
    /**
     * Sets the four_pi_epsilon constant to the indicated value and turns off the
     * printWarning function when returning this variable (get_four_pi_epsilon).
     * All dependant variables are recalculated again.
     */
    double set_four_pi_eps_i(double value);
    /**
     * Sets the boltz constant to the indicated value and turns off the
     * printWarning function when returning this variable (get_boltz).
     * All dependant variables are recalculated again.
     */
    double set_boltzmann(double value);
    /**
     * Sets the hbar constant to the indicated value and turns off the
     * printWarning function when returning this variable (get_hbar).
     * All dependant variables are recalculated again.
     */
    double set_hbar(double value);
    /**
     * Recalculates all physical constants which are dependant from the four
     * basic physical constants defined in the topology.
     * It sets also the bools to true if the variables are fully derived
     * by values e.g. from the topology and now inertialised (hard coded)
     * values were taken.
     */
    void calc();
    /**
     * Prints a warning message (used to warn the user that the physical
     * constant is not coming from the topology.
     * All dependant variables are recalculated again.
     */
    void printWarning(std::string name, double &value, bool &b);

    /**
     * Prints out all variables of the class PhysConst on the screen (std::cerr).
     */
    void printAll();

  
  private:

    // ==========================
    //    NUMERICAL CONSTANTS:
    // ==========================
    //
    /**
     * circle constant pi, calculated as 4 * atan(1.0)
     */
    double pi;
    /**
     * prefactor for kilo
     */
    double kilo;
    /**
     * prefactor for nano
     */
    double nano;
    /**
     * prefactor for pico
     */
    double pico;
    /**
     * Euler number
     */
    double euler;
    /**
     * factor to convert radian2 to degrees
     */
    double radian2degree;
    /**
     * factor toconvert degrees to radians
     */
    double degree2radian;
    /**
     * a small number, a factor for float comparison
     */
    double epsilon;

    
    // =========================
    //    PHYSICAL CONSTANTS:
    // =========================
    //
    /**
     * atomic mass unit (amu or u)\n
     * inertialised units: kg
     */
    double atomic_mass_unit;
    /**
     * to know if atomic_mass_unit was read from topology (true) ore not(false)
     */
    bool atomic_mass_unit_;
    /**
     * Avogadro number\n
     * inertialised units: 1/mol
     */
    double avogadro;
    /**
     * to know if avogadro was read from topology (true) ore not(false)
     */
    bool avogadro_;
    /**
     * the Boltzmann/gas constant\n
     * inertialised units: kJ/(mol K)
     */
    double boltzmann;
    /**
     * to know if boltzmann was read from topology (true) ore not(false)
     */
    bool boltzmann_;
    /**
     * elementary charge / C
     */
    double elementary_charge;
    /**
     * to know if elementary_charge was read from topology (true) ore not(false)
     */
    bool elementary_charge_;
    /**
     * dielectric permitivitty of vacuum / (mol e/kJ nm)
     */
    double eps0;
    /**
     * to know if eps0 was read from topology (true) ore not(false)
     */
    bool eps0_;
    /**
     * coulombic prefactor
     */
    double four_pi_eps_i;
    /**
     * to know if four_pi_eps_i was read from topology (true) ore not(false)
     */
    bool four_pi_eps_i_;
    /**
     * Plank constant (not devided by 2 pi) in gromos units
     */
    double h;
    /**
     * to know if h was read from topology (true) ore not(false)
     */
    bool h_;
    /**
     * Plank constant devided by 2 pi
     */
    double hbar;
    /**
     * to know if hbar was read from topology (true) ore not(false)
     */
    bool hbar_;
    /**
     * permeability of free space
     */
    double mu0;
    /**
     * to know if mu0 was read from topology (true) ore not(false)
     */
    bool mu0_;
    /**
     * speed of light / (nm/ps)
     */
    double speed_of_light;
    /**
     * to know if speed_of_light was read from topology (true) ore not(false)
     */
    bool speed_of_light_;
  };

  // an extern class variable to hold/handle the physical constants
  extern PhysConst physConst;
}

#endif	/* INCLUDED_PHYSICS_H */


