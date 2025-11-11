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

// utils_Energy.h

// Energy class: calculates all energies of a system

#ifndef INCLUDED_UTILS_ENERGY
#define INCLUDED_UTILS_ENERGY

#include <vector>
#include <set>

#include "SimplePairlist.h"
#include "../gmath/Vec.h"
#include "../gromos/Exception.h"

namespace gcore{
  class GromosForceField;
  class System;
}
namespace bound{
  class Boundary;
}
namespace gmath{
  class Vec;
}

namespace utils
{
  class AtomSpecifier;
  class PropertyContainer;
  class Property;
  class Energy;
  /**
   * Class Energy
   * Purpose: calculates the potential energy for properties or atoms
   *
   * Description:
   * The energy class can be used to calculate bonded potential energies 
   * for specified Properties, or non-bonded potential energies involving 
   * specific atoms. What you want to calculate is parsed as Property- and 
   * AtomSpecifiers. A call to the member function calc actually performs 
   * the calculation, after which a series of accessors are available.
   * <p>
   * For the non-bonded interactions, a chargegroup based cutoff is applied
   * with a Reaction field correction for the electrostatic interactions.<br>
   * The Vanderwaals interactions are calculated as
   * @f[ V^{LJ}=\left[\frac{C_{12}(i,j)}{r_{ij}^6}-C_6(i,j)\right] 
   *     \frac{1}{r_{ij}^6} @f]
   * And the electrostatic interactions from
   * @f[ V^{CRF}=\frac{q_iq_j}{4\pi\epsilon_0}\left[\frac{1}{r_{ij}} -
   *     \frac{\frac{1}{2}C_{rf}r_{ij}^2}{R_{rf}^3} - 
   *     \frac{(1-\frac{1}{2}C_{rf})}{R_rf}\right] @f]
   * <p>
   * The bonded interactions are calculated for every specified Property:<br>
   * bonds:
   * @f[ V^{bond}=\frac{1}{4}K_{b_n}\left[b_n^2 - b_{0_n}^2\right]^2@f]
   * angles:
   * @f[ V^{angle}=\frac{1}{2}K_{\theta_n}\left[\cos{\theta_n} - 
   *     \cos{\theta_{0_n}}\right]^2@f]
   * dihedrals:
   * @f[ V^{trig}=K_{\phi_n}\left[1+\cos(\delta_n)\cos(m_n\phi_n)\right]
   * @f]
   * impropers:
   * @f[ V^{har}=\frac{1}{2}K_{\xi_n}\left[\xi_n - \xi_{0_n}\right]^2
   * @f]
   *
   * @class Energy
   * @author C. Oostenbrink
   * @ingroup utils
   * @sa utils::PropertyContainer
   * @sa utils::AtomSpecifier
   */
  class Energy{
    gcore::System *d_sys;
    gcore::GromosForceField *d_gff;
    bound::Boundary *d_pbc;
    utils::AtomSpecifier *d_as;
    utils::PropertyContainer *d_pc;
    std::vector<double> d_cov, d_vdw_m, d_vdw_s, d_el_m, d_el_s;
    std::vector<std::vector<double> > d_covpar;
    utils::AtomSpecifier *d_soft;
    double d_lam, d_alj, d_ac, d_eps, d_kap, d_cut, d_p_vdw, d_p_el;
    std::vector<std::set<int> > d_ex, d_third;
    std::vector<utils::SimplePairlist> d_pl;
    std::vector<gmath::Vec> d_f_el_m, d_f_el_s;
    bool d_RFex;
    bool coulomb_scaling;
  public: 
    /**
     * Energy Constructor
     * @param sys We need to have the System information
     * @param gff The GromosForceField contains all force field parameters
     *            that are needed
     * @param pbc We also need to know about the boundary conditions that
     *            are to be applied
     */
    Energy(gcore::System &sys, gcore::GromosForceField &gff, 
           bound::Boundary &pbc);
    /**
     * Energy deconstructor
     */
    ~Energy(){}
   
    /**
     * The method setAtoms allows you to specify which atoms you want to 
     * calculate the non-bonded interactions for. Interactions of these 
     * atoms with ALL atoms in the system, are considered.
     * @param as An AtomSpecifier containing the atoms that you are 
     *           interested in
     */
    int setAtoms(utils::AtomSpecifier &as);

    /**
     * The method setProperties allows you to specify which covalent 
     * interactions you want to monitor
     * @param pc A PropertyContainer with all properties you are 
     *           interested in.
     */
    int setProperties(utils::PropertyContainer &pc);

    /**
     * The method setSoft allows you to specify some atoms that are 
     * considered soft in the calculation, together with the necessary 
     * parameters.
     *
     * In the current implementation a soft atom has the following 
     * Vanderwaals interaction:
     * @f[V^{LJ}_{soft}=\left[\frac{C_{12}(i,j)}{\alpha_{LJ}\lambda^2
     * C_{126}(i,j) + r_{ij}^6}-C_6(i,j)\right]\frac{1}{\alpha_{LJ}\lambda^2
     * C_{126}(i,j) + r_{ij}^6}@f]
     * with C126 = C12/C6 if C12!=0 && C6 !=0 and C126=0 otherwise
     * <p>
     * And for the Electrostatic interaction we use:
     * @f[ V^{CRF}_{soft}=\frac{q_iq_j}{4\pi\epsilon_0} 
     * \left[\frac{1}{\left[\alpha_C(i,j)\lambda^2 + r_{ij}^2\right]^{1/2}} - 
     * \frac{\frac{1}{2}C_{rf}r_{ij}^2}{\left[\alpha_C(i,j)\lambda^2 
     * + R_{rf}^2\right]^{3/2}} - \frac{(1-\frac{1}{2}C_{rf})}{R_{rF}}\right]
     * @f]
     * 
     * @param soft An AtomSpecifier to say which atoms are soft
     * @param lam  The value of lambda (@f$\lambda@f$)
     * @param alj  The soft Vdw parameter to use (@f$\alpha_{LJ}@f$)
     * @param ac   The soft Coulomb parameter to use (@f$\alpha_{C}@f$)
     */
    void setSoft(utils::AtomSpecifier &soft, double lam, double alj, double ac);

    /**
     * Method to set the Reaction field parameters that are required for 
     * the calculation of the electrostatic interaction energy.
     * 
     * These variables are needed to calculate the reaction field coefficient
     * @f[ C_{rf} = \frac{(2-2\epsilon_2)(1+\kappa R_{rf}) - 
     *     \epsilon_2(\kappa R_{rf})^2}{(1+2\epsilon_2)(1+\kappa R_{rf}) +
     *     \epsilon_2(\kappa R_{rf})^2} @f]
     *
     * @param eps The reaction field relative permitivity (@f$\epsilon_2@f$)
     * @param kap The inverse Debye screening length (@f$\kappa@f$)
     */
    void setRF(double eps, double kap);

    /**
     * Method to set the cutoff radius
     *
     * Only non-bonded interactions between atoms for which the centre of 
     * geometries of their respective charge groups lie within this cutoff
     * are calculated
     * @param cut The cutoff that is to be used
     */
    void setCutOff(double cut);
    
    /**
     * Method to set the pairlist type
     *
     * @param t ATOMIC or CHARGEGROUP
     */
    void setPairlistType(std::string t);

    /** 
     * Method to turn on the RF contribution for excluded atoms
     */
    void setRFexclusions(bool p);
    
    /** 
     * Method to turn on the RF contribution for excluded atoms
     */
    void setCoulombScaling(bool p);
    
    /**
     * Method to actually perform the calculations
     *
     * A call to this function will calculate all covalent interactions and
     * the non-bonded interactions for all specified properties and atoms
     */
    void calc();
    /**
     * Method to calculate the non-bonded interactions by first making 
     * the pairlist
     */
    void calcNb();
    /**
     * Method to calculate the electrostatic force on the atoms
     */
    void calcField();
    /**
     * Method to calculate the pairlists for all relevant atoms
     */
    void calcPairlist();
    /**
     * Method to calculate the reaction field contribution for the excluded
     * atoms and the self-interaction
     */
    void calcRFex();
    /**
     * Method to set the pairlist of particle i
     */
    void setPairlist(int i, SimplePairlist &as);
    /**
     * Method to calculate the non-bonded interactions, assuming that the
     * pairlist has already been made
     */
    void calcNb_interactions();
    
    /**
     * Method to calculate the covalent interactions
     */
    void calcCov();
    
    /**
     * Method to calculate the non-bonded interactions between two of the 
     * specified atoms
     * @param i,j The interaction between the i-th and the j-th atoms in the
     *            AtomSpecifier is calculated
     * @param vdw Is returned with the vdw energy between the atoms
     * @param el  Is returned with the electrostatic energy between the atoms
     */
    void calcPair(int i, int j, double &vdw, double &el);
    
    /**
     * Accessor, returns the total potential energy: all bonded and non-bonded
     * terms are added together
     */
    double tot() const;

    /**
     * Accessor, returns the total non-bonded energy
     */
    double nb() const;

    /**
     * Accessor, returns the total Vanderwaals energy
     */
    double vdw() const;

    /**
     * Accessor, returns the total electrostatic energy
     */
    double el() const;

    /**
     * Accessor, returns the total covalent energy
     */
    double cov() const;

    /**
     * Accessor, returns the total bonded energy
     */
    double dist() const;

    /**
     * Accessor, returns the total energy from angles
     */
    double angle() const;

    /**
     * Accessor, returns the total energy from improper dihedrals
     */
    double impdihed() const;

    /**
     * Accessor, returns the total energy from torsional dihedrals
     */
    double torsdihed() const;
    
    /** 
     * Accessor, returns the total Vanderwaals energy of the i-th atom 
     * in the AtomSpecifier
     * @param i The i-th atom in the AtomSpecifier
     */
    double vdw(unsigned int i) const;

    /**
     * Accessor, returns the total Vanderwaals energy with other solute atoms
     */
    double vdw_m() const;

    /**
     * Accessor, returns the total Vanderwaals energy with solvent atoms
     */
    double vdw_s() const;
    
    /**
     * Accessor, returns the total electrostatic energy with other solute atoms
     */
    double el_m() const;

    /**
     * Accessor, returns the total electrostatic energy with solvent atoms
     */
    double el_s() const;

    /**
     * Accessor, returns the total electrostatic energy of the i-th atom
     * in the AtomSpecifier
     * @param i The i-th atom in the AtomSpecifier
     */
    double el(unsigned int i) const;

    /**
     * Accessor, returns the potential energy corresponding to the i-th 
     * property in the PropertyContainer
     * @param i The i-th property in the PropertyContainer
     */
    double cov(unsigned int i) const;

    /**
     * Accessor, returns the total Vanderwaals interaction of the i-th atom
     * in the AtomSpecifier with other Solute atoms
     * @param i The i-th atom in the AtomSpecifier
     */
    double vdw_m(unsigned int i) const;

    /**
     * Accessor, returns the total Vanderwaals interaction of the i-th atom
     * in the AtomSpecifier with Solvent atoms
     * @param i The i-th atom in the AtomSpecifier
     */
    double vdw_s(unsigned int i) const;

    /**
     * Accessor, returns the total Electrostatic interaction of the i-the atom
     * in the AtomSpecifier with other Soluter atoms
     * @param i The i-th atom in the AtomSpecifier
     */
    double el_m(unsigned int i) const;

    /**
     * Accessor, returns the total Electrostratic interaction of the i-th atom
     * in the AtomSpecifier with Solvent atoms
     * @param i The i-th atom in the AtomSpecifier
     */
    double el_s(unsigned int i) const;
    /**
     * Accessor, returns the total electrostatic force of the i-th atom
     * i n the AtomSpecifier
     * @param i The i-th atom in the AtomSpecifier
     */
    gmath::Vec f_el(unsigned int i) const;
    /**
     * Accessor, returns the total electrostatic force of the i-th atom
     * i n the AtomSpecifier with Solute atoms
     * @param i The i-th atom in the AtomSpecifier
     */
    gmath::Vec f_el_m(unsigned int i) const;
    /**
     * Accessor, returns the total electrostatic force of the i-th atom
     * i n the AtomSpecifier with Solvent atoms
     * @param i The i-th atom in the AtomSpecifier
     */
    gmath::Vec f_el_s(unsigned int i) const;
    
    // Exception
    /**
     * Exception to be thrown if anything goes wrong
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string &what): gromos::Exception("Energy", what){}
    };

  protected:
    /**
     * A function to calculate the centre of geometry for the charge group
     * to which atom i belongs
     * @param i The i-th atom in the AtomSpecifier
     * @returns A gmath::Vec containing the coordinates of the centre of 
     *          geometry of the charge group
     */
    gmath::Vec calcChgrp(int i) const;
    /**
     * A function to determine the BondType of the Property that is specified.
     * If the bond is found in the topology, its type is returned, otherwise 
     * an exception is thrown.
     */
    int findBond(utils::Property &pp) const;
    /**
     * A function to determine the AngleType of the Property that is
     * specified. If the angle is found in the topology, its type is returned,
     * otherwise an exception is thrown.
     */
    int findAngle(utils::Property &pp) const;
    /**
     * A function to determine the DihedralType or ImproperType of the 
     * Property that is specified. If the Property matches an improper, 
     * -(type+1) is returned. If it is a Dihedral, its type is returned.
     * If it cannot be found at all, an exception is thrown.
     */
    std::vector<int> findDihedral(utils::Property &pp) const;
    /**
     * A function to determine the DihedralType of the
     * Property that is specified.
     * If it cannot be found at all, an exception is thrown.
     */
    std::vector<int> findCrossDihedral(utils::Property &pp) const;
    /**
     * Calculates the Bond interaction energy
     * @param val The bondlength
     * @param par A vector containing the optimum bond length @f$(b_{0_n})@f$
     *            and force constant @f$(K_{b_n})@f$, respectively
     * @return The potential energy for this bond
     */ 
    double calcBond(double val, const std::vector<double> & par) const;
    /**
     * Calculates the Angle interaction energy
     * @param val The value of the angle
     * @param par A vector containing the optimum angle @f$(\theta_{0_n})@f$ 
     * and force constant @f$(K_{\theta_n})@f$, respectively
     * @return The potential energy for this angle
     */ 
    double calcAngle(double val, const std::vector<double> & par) const;
    /**
     * Calculates the Dihedral interaction energy
     * @param val The value of the dihedral
     * @param par A vector containing the Phase @f$(\cos(\delta_n))@f$, 
     *            multiplicity @f$(m_n)@f$ and force constant 
     * @f$(K_{\phi_n})@f$ for this dihedral, respectively
     * @return The potential energy for this dihedral
     */
    double calcDihedral(double val, const std::vector<double> & par) const;
    /**
     * Calculates the CrossDihedral interaction energy
     * @param val The value of the cross dihedral
     * @param par A vector containing the Phase @f$(\cos(\delta_n))@f$,
     *            multiplicity @f$(m_n)@f$ and force constant
     * @f$(K_{\phi_n})@f$ for this dihedral, respectively
     * @return The potential energy for this dihedral
     */
    double calcCrossDihedral(double val, const std::vector<double> & par) const;
    /**
     * Calculates the Improper interaction energy
     * @param val The value of the improper
     * @param par A vector containing the optimum value @f$(\xi_{0_n})@f$ and 
     *            force constant @f$(K_{\xi_n})@f$ for this improper, 
     *            respectively
     * @return The potential energy for this improper
     */
    double calcImproper(double val, const std::vector<double> & par) const; 
};
//inline functions and methods

inline double Energy::tot() const
{
  return this->nb() + this->cov();
}
inline double Energy::nb() const
{
  return this->vdw() + this->el();
}
inline void Energy::setCutOff(double cut)
{
  d_cut=cut;
}
inline void Energy::setSoft(utils::AtomSpecifier &soft, double lam, double alj, double ac)
{
  d_soft=&soft;
  d_lam=lam;
  d_alj=alj;
  d_ac=ac;
}
inline void Energy::setRF(double eps, double kap)
{
  d_eps=eps;
  d_kap=kap;
}
}
#endif
