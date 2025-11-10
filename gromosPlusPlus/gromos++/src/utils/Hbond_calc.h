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

/**
 * @file
 * @author siggj, Update: @ref ms
 *
 * Created on October 24, 2012, 3:31 PM
 */

#ifndef HBOND_H
#define	HBOND_H

#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <map>
#include <vector>

#include "AtomSpecifier.h"
#include "CubeSystem.hcc"
#include "groTime.h"

namespace gcore {
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace gmath {
  class Vec;
}

namespace args {
  class Arguments;
}

namespace bound {
  class Boundary;
}

namespace utils {

/**
* Struct Key2c
* This struct is used as Key for a two-centered H-bond map. It consist of two atoms (donor and acceptor) and has comparison operators.
* @author @ref ms
* @ingroup utils
* @struct Key2c
*/
  struct Key2c{ //key for 2c hb map
    private:
    int don;
    int acc;

    public:
    /**
    * Constructor with 2 input atoms.
    */
    Key2c(int a, int b):don(a),acc(b){}
    /**
    * Constructor with an input Key.
    */
    Key2c(const Key2c& right){
        right.get_atom_num(don,acc);
    }
    /**
    * Default Constructor.
    */
    Key2c(){}
    /**
    * Method, which returns donor and acceptor atom by setting the two input arguments.
    */
    void get_atom_num(int& a, int& b) const{
        a=don;
        b=acc;
    }
    /**
    * Method to set donor and acceptor atom.
    */
    void set_index(int a, int b){
        don=a;
        acc=b;
    }
    /**
    * Member operator which returns if both Keys are equal.
    */
    bool operator== (const Key2c& right) const{
        return don == right.don && acc == right.acc;
    }
    /**
    * Member operator which returns if both Keys are not equal.
    */
    bool operator!= (const Key2c& right) const{
          return !(*this == right);
    }
    /**
    * Member operator for map construction.
    */
    bool operator< (const Key2c& right) const{ //for std::map & std::sort
      if (don < right.don)  return true;
      if (don > right.don)  return false;

      if (acc < right.acc)  return true;
      if (acc > right.acc)  return false;

      return false;
    }
  };

/**
* Struct Key3c
* This struct is used as Key for a three-centered or solvent-bridge H-bond map. It consist of four atoms (donor1, acceptor1, donor2, acceptor2) and has comparison operators.
* @author @ref ms
* @ingroup utils
* @struct Key3c
*/
  struct Key3c{  //key for 3c & bridges map
    private:
    int d_1, a_1, d_2, a_2; //donor & acceptor of first 2c hbond, and donor & acc of 2nd 2c hbond

    public:
    /**
    * Constructor with four atoms.
    */
    Key3c(int w, int x, int y, int z):d_1(w),a_1(x),d_2(y),a_2(z){}
    /**
    * Constructor with two two-centered Keys.
    */
    Key3c(const Key2c& left, const Key2c& right){
        left.get_atom_num(d_1,a_1);
        right.get_atom_num(d_2,a_2);
    }
    /**
    * Default Constructor.
    */
    Key3c(){}
    /**
    * Method, which returns all four atoms by setting the four input arguments.
    */
    void get_atom_num(int& w, int& x, int& y, int& z) const{
        w = d_1;
        x = a_1;
        y = d_2;
        z = a_2;
    }
    /**
    * Method, which sets all four atoms by specifying 4 atoms.
    */
    void set_index(int w, int x, int y, int z){
        d_1 = w;
        a_1 = x;
        d_2 = y;
        a_2 = z;
    }
    /**
    * Method, which sets all four atoms by specifying two two-centered Keys.
    */
    void set_index(const Key2c& left, const Key2c& right){
        left.get_atom_num(d_1,a_1);
        right.get_atom_num(d_2,a_2);
    }
    /**
    * Member operator which returns if both Keys are equal.
    */
    bool operator== (const Key3c& r) const{
        return d_1 == r.d_1 && a_1 == r.a_1 && d_2 == r.d_2 && a_2 == r.a_2;
    }
    /**
    * Member operator which returns if both Keys are not equal.
    */
    bool operator!= (const Key3c& r){
        return !(*this == r);
    }
    /**
    * Member operator for map construction.
    */
    bool operator< (const Key3c& r) const{
        if (d_1 < r.d_1)  return true;
        if (d_1 > r.d_1)  return false;
        //l.i==r.i
        if (a_1 < r.a_1)  return true;
        if (a_1 > r.a_1)  return false;
        //l.j==r.j
        if (d_2 < r.d_2)  return true;
        if (d_2 > r.d_2)  return false;

        if (a_2 < r.a_2)  return true;
        if (a_2 > r.a_2)  return false;
        //l.k==r.k
        return false;
    }
  };

  /**
   * Class Timeseries
   * Stores all relevant informations for generating the timeseries ouput: the current timestep, number of H-bonds at the current timestep,
   * and a vector of all H-bonds occurring at this timestep.
   * @author @ref ms
   * @ingroup utils
   * @class Timeseries
   */
  template <typename T1>
  class Timeseries{
    double tyme; //current timestep
    int number; //number of hbonds in total at this timestep

    std::vector<T1> key_list; //vector of hbond keys

    public:
    /**
    * Default constructor.
    */
    Timeseries() : tyme(0), number(0)
    {}
    /**
    * Constructor with time argument.
    */
    Timeseries(double t) : tyme(t), number(0)
    {}
    /**
    * Method to set the current time.
    */
    void set_time(double t){
        tyme = t;
    }
    /**
    * Method to add an H-bond Key ONCE to the vector storing the Keys.
    */
    void add_once(const T1& key){
        typename std::vector<T1>::iterator it = std::lower_bound(key_list.begin(), key_list.end(), key); //binary search for the iterator put

        if(it == key_list.end() || *it != key) //if it is not found in the vector, insert it
            key_list.insert(it, key);
    }
    /**
    * Method to add an H-bond Key to the vector storing the Keys.
    */
    void add_all(const T1& key){
        key_list.push_back(key);
    }
    /**
    * Method to sort the H-bond key_list in ascending order of the keys.
    */
    void sort_keys(){
        std::sort(key_list.begin(), key_list.end());
    }
    /**
    * Method to set the number of H-bonds.
    */
    void set_num(int n){
        number = n;
    }
    /**
    * Method returning the timestep.
    */
    double time() const{
        return tyme;
    }
    /**
    * Method returning the number of H-bonds at the current timestep.
    */
    double num() const{
        return number;
    }
    /**
    * Method returning a reference to the vector containing all H-bond Keys at the current timestep.
    */
    const std::vector<T1>& keys() const {
        return key_list;
    }
    /**
    * Method returning the estimated size of the Timeseries object.
    */
    double size() const {
        return (sizeof(double)+sizeof(int)+key_list.capacity()*sizeof(T1))/1000.0;
    }
  };

  /**
  * Function used for sorting a Container of Timeseries objects in ascending order of the timestep.
  */
  template <typename T1>
  inline bool CompTime(const Timeseries<T1>& left, const Timeseries<T1>& right){
      return left.time() < right.time();
  }
  /**
  * Function used for sorting iterators pointing to H-bond maps in descending order of occurence of the H-bond.
  */
  template <typename T1, typename T2>
  inline bool sort_rev_by_occ(typename std::map<T1, T2>::iterator left, typename std::map<T1, T2>::iterator right){
    return left->second.num() > right->second.num(); //compare the number of hbonds between map hbond entries
  }

  /**
   * Class HB_calc
   * purpose: calculate (native) intra-, intermolecular, solute-solvent and solvent-solvent
   * hydrogen bonds over the trajectory.
   *
   * Description:
   * The HB_calc class calculates (native) intra-, intermolecular, solute-solvent and solvent-solvent
   * hydrogen bonds over the trajectory. It prints out two timeseries (total hydrogen bonds per frame
   * and the occurrence of a specific hydrogen bond at a given time) and statistics. It makes use of
   * the AtomSpecifier if specific donor or acceptor atoms want to be specified or it can read in a
   * file containing the masses of hydrogens, donor, and acceptor atoms and thereby make the selection.
   * @author J. Sigg, M.Setz
   * @ingroup utils
   * @class HB_calc
   */
  class HB_calc {
  protected:
    double max_distance2, min_angle, time, time_start, time_dt;
    gcore::System *sys;
    args::Arguments *args;
    bound::Boundary *pbc;
    utils::AtomSpecifier donors, bound, acceptors;
    std::vector<int> donAB,donA,donB, accAB,accA,accB;
    std::vector<double> mass_hydrogens, mass_acceptors, mass_donors;
    int dummy;
    unsigned int frames, num_A_donors, num_A_acceptors, numHb;
    std::ofstream timeseriesHB, timeseriesHBtot;
    bool reduce, read_time, excl_dummy;
    std::vector<int> solv_donor, solv_acc;

    /**
     * Method that stores the system and all the arguments for further use.
     */
    void setval(gcore::System& sys, args::Arguments& args, int dummyIAC);
    /**
     * Method that reads the massfile, in which the hydrogen, donors,
     * and acceptors masses are stored.
     */
    void readinmasses(std::string filename);

    /**
     * Method that populates the vectors donAB, donA, donB, accAB, accA, accB:
     * <ul>
     * <li>donAB contains all atoms that occur in DonorAtomsA AND DonorAtomsB.</li>
     * <li>donA contains all atoms that occur ONLY in DonorAtomsA.</li>
     * <li>donB contains all atoms that occur ONLY in DonorAtomsB.</li>
     * <li>accAB contains all atoms that occur in AcceptorAtomsA AND AcceptorAtomsB.</li>
     * <li>accA contains all atoms that occur ONLY in AcceptorAtomsA.</li>
     * <li>accB contains all atoms that occur ONLY in AcceptorAtomsB.</li>
     * </ul>
     * The following donor-acceptor combinations <b>MUST BE USED</b> in hbond calculation to ensure that all donor-acceptor pairs are included:
     * <ol>
     * <li>donAB and accA</li>
     * <li>donB and (accAB + accA)</li>
     * <li>(donAB + donA) and (accAB + accB)</li>
     * </ol>
     * These combinations guarantee that no H-bonds will be counted twice (which was the case in the old hbond program), even if the user specifies it.
     * So a user input of
 @verbatim
DonorAtomsA=1:a
AcceptorAtomsA=1:a
DonorAtomsB=1:a
AcceptorAtomsB=1:a @endverbatim
     * will result in that only donAB and accAB are filled (all other vectors are empty), so only point 3. of the above list will take effect
     */
    void set_subatomspecs();
    /**
     * Method that opens the timeseries files.
     * timeseriesHB contains the number of H-bonds at a given timestep.
     * timeseriesHBtot contains all H-bonds at a given timestep.
     */
    void opents(string fi1, string fi2) {
      timeseriesHB.open(fi1.c_str());
      timeseriesHBtot.open(fi2.c_str());
    }//end opents()
    /**
     * Method to read in all atoms chosen by the input file. Calls add_atoms internally.
     */
    void determineAtoms();
    /**
    * Method that adds all atoms in the specified Argument names to the given AtomSpecifiers my_donor & my_acceptor. It also adds donor-bound atoms to the AtomSpecifier my_bound.
    * If \@excludedummy is requested by the user, it also removes all Dummy atoms from my_donor, my_acceptor and my_bound.
    */
    void add_atoms(const string donor_arg_name, const string acc_arg_name, AtomSpecifier& my_donor, AtomSpecifier& my_acc, AtomSpecifier& my_bound);
    /**
     * Method that cross-reference all atoms from the input file with the massfile.
     */
    void determineAtomsbymass();
    /**
     * Method to read in a frame.
     */
    void readframe();
    /**
     * Method that returns true if the atom i and j are NOT neighbours.
     */
    bool neighbour(int i, int j) {
      return (bound.atom(i) != acceptors.atom(j) ||
              bound.mol(i) != acceptors.mol(j));
    }//end HB_calc::neighbour()
    /**
     * Method that gives back true if the atom i and j have a distance smaller
     * than the maximal distance given.
     */
    bool distances(double &dist, gmath::Vec &tmpA) {
      dist = tmpA.abs2();
      return (max_distance2 >= dist);
    }//end HB_calc::distances()
    /**
     * Method that gives back true if the atom i and j have a angle bigger
     * than the minimal angle given.
     */
    bool angle(int i, double &angles, gmath::Vec &bound_i, gmath::Vec &tmpA) {
      gmath::Vec tmpB;
      tmpB = bound_i - donors.pos(i);
      angles = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;
      return (min_angle <= angles);
    }//end HB_calc::angle()
    /**
     * Method that populates the vectors solv_donor and solv_acc, which store information about the solvent, IF \@reducesolvent was requested.
     * solv_donor stores all donor atoms of the first solvent molecule and the position of the first donor solvent atom.
     * solv_acc stores all acceptor atoms of the first solvent molecule and the position of the first acceptor solvent atom.
     * This information is later used to construct the same Key for any solvent atom.
     */
    void set_reduce();
    /**
    * Pure virtual method to initialize H-bond calculation.
    */
    virtual void init_calc() =0;

  public:

    /**
     * Constructor that stores the maximal distance, the minimal angle, and if solvent should be reduced.
     */
    HB_calc(bool red, double max_distance2 = 0, double min_angle = 0) :
    max_distance2(max_distance2), min_angle(min_angle), time(0.0), frames(0), numHb(0), reduce(red) {
    }
    /**
    * Virtual destructor.
    */
    virtual ~HB_calc(){
        //delete pbc;
    }

    /**
     * Method setting the current time.
     */
    void settime(double times) {
      time = times;
    }//end HB_calc::settime()
    /**
     * Method returning a reference to the donors AtomSpecifier, so it only needs to be created once by HB2c_calc and can be re-used by HB3c_calc and HB_bridges_calc.
     */
    const AtomSpecifier& get_donors() const{
       return donors;
    }
    /**
     * Method returning a reference to the acceptors AtomSpecifier, so it only needs to be created once by HB2c_calc and can be re-used by HB3c_calc and HB_bridges_calc.
     */
    const AtomSpecifier& get_acceptors() const{
        return acceptors;
    }
    /**
     * Method returning a reference to the donors AtomSpecifier, so it only needs to be created once by HB2c_calc and can be re-used by HB3c_calc and HB_bridges_calc.
     */
    const AtomSpecifier& get_bound() const{
        return bound;
    }
    /**
     * Method which copies accAB, accA, accB, so they only need to be created once by HB2c_calc and can be re-used by HB3c_calc and HB_bridges_calc.
     */
    void get_all_acc(std::vector<int>& A, std::vector<int>& B, std::vector<int>& C) const{
        A = accAB;
		B = accA;
		C = accB;
    }
    /**
     * Method which copies donAB, donA, donB, so they only need to be created once by HB2c_calc and can be re-used by HB3c_calc and HB_bridges_calc.
     */
	void get_all_don(std::vector<int>& A, std::vector<int>& B, std::vector<int>& C) const{
        A = donAB;
        B = donA;
        C = donB;
    }
    /**
     * Method which returns the number of atoms in DonorsA and AcceptorsA.
     */
    void get_num_A(unsigned int& a_don, unsigned int& a_acc) const{
        a_don = num_A_donors;
        a_acc = num_A_acceptors;
    }
  }; //end class HB_calc
}

#endif	/* HBOND_H */

