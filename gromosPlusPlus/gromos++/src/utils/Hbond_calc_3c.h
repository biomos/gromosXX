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
 * @file:   Hbond_calc_3c.h
 * Author: siggj
 *
 * Created on October 10, 2012, 2:17 PM
 */

#ifndef INCLUDED_3C_H
#define	INCLUDED_3C_H

#include <vector>
#include <map>
#include <iterator>

#include "Hbond_calc.h"
#include "Hbond_calc_2c.h"

namespace utils {

  /**
   * Struct HBPara3c
   * purpose: to store all the paramaters (maximal distance, minimal angle, minimal
   * angle sum and maximal dihedral angle) given by the input file.
   * @struct HBPara3c
   * @author J.Sigg
   * @ingroup utils
   */
  struct HBPara3c {
    double maxdist, maxdist2, minangle;
    double minanglesum, maxdihedral;
  }; //end struct HBPara3c


  /**
   * Class HB3c
   * purpose: to store all 3-centred H-bonds, which have been found, with all
   * the informations, (distances, angles, sum of angels, dihedral angles, ID, H-bond count) needed.
   * @author J.Sigg, M.Setz
   * @ingroup utils
   * @class HB3c
   */
  class HB3c {
    int number, _id;
    vector<double> distance, angle;
    double dihed, angletot;
  public:

    /**
     * Constructor
     */
    HB3c() :number(0), _id(0) ,dihed(0), angletot(0) {
      distance.resize(2,0);
      angle.resize(2,0);
    }
    /**
     * Method, which increases the H-bond count.
     */
    void add(){
        ++number;
    }
    /**
     * Method, which adds the given values to all parameters, and increases the H-bond count.
     */
    void add(double dist1, double dist2, double a1, double a2, double as, double dih){
        ++number;
        distance[0] += sqrt(dist1);
        distance[1] += sqrt(dist2);
        angle[0] += a1;
        angle[1] += a2;
        angletot += as;
        dihed += dih;
    }
    /**
     * Method, which sets the H-bond ID.
     */
    void set_id(int i){
        _id = i;
    }
     /**
     * Method, which returns the H-bond count.
     */
    int num() const{
      return number;
    }//end HB3c::getnum()
    /**
     * Method, which returns the mean distance of the first or second H-bond contributing to this three-centered H-bond.
     */
    double meandist(unsigned int i) const{
      if(i>1)
        throw gromos::Exception("Hbond_calc_3c","getmeandist(i) does not exist. Must be 0<=i<=1");
      if(number <= 0)
        throw gromos::Exception("Hbond_calc_3c","number of 3c Hbonds <= 0");
      return distance[i] / number;
    }//end HB3c::getmeandist()
    /**
     * Method, which returns the mean angle of the first or second H-bond contributing to this three-centered H-bond.
     */
    double meanangle(unsigned int i) const{
      if(i>1)
        throw gromos::Exception("Hbond_calc_3c","getmeanangle(i) does not exist. Must be 0<=i<=1");
      if(number <= 0)
        throw gromos::Exception("Hbond_calc_3c","number of 3c Hbonds <= 0");
      return angle[i]/number;
    }//end HB3c::getmeanangle()
    /**
     * Method, which returns the mean of the sum of angles.
     */
    double meanangle_sum() const{
      if(number <= 0)
        throw gromos::Exception("Hbond_calc_3c","number of 3c Hbonds <= 0");
      return angletot/number;
    }//end HB3c::getmeanangle_sum()
    /**
     * Method, which returns the mean dihedral angle.
     */
    double meandihedral() const{
      if(number <= 0)
        throw gromos::Exception("Hbond_calc_3c","number of 3c Hbonds <= 0");
      return dihed/number;
    }//end HB3c::getmeandihedral()
    /**
     * Method, which returns the H-bond ID.
     */
    int id() const {
        return _id;
    }
    /**
     * Method, which adds two HB3c objects.
     */
    void add(const HB3c& rightop){ //HB2c
        distance[0] += rightop.distance[0]; //distances in dist are already sqrt, so we can simply add them
        distance[1] += rightop.distance[1];
        angle[0] += rightop.angle[0];
        angle[1] += rightop.angle[1];
        angletot += rightop.angletot;
        dihed += rightop.dihed;
        number += rightop.number;
    }
  }; //end class HB3c

  /**
   * Class HB3c_calc
   * purpose: Class, which inherit from HB_calc, to calculate the 3-centred H-bonds.
   * All found bonds are stored in a map with Key3c as key and HB3c as value.
   * @author J. Sigg, M.Setz
   * @ingroup utils
   * @class HB2c_calc
   */
  class HB3c_calc : public HB_calc {

    typedef std::map<Key3c, HB3c> HB3cContainer;
    typedef std::vector< HB3cContainer::iterator > HB3cContainerIteratorList;
    typedef std::vector<Timeseries<Key3c> > TimeseriesContainer;

    double min_angle_sum, max_dihedral;
    std::vector<Key3c> native_key_storage;
    TimeseriesContainer ts;
    HB3cContainer hb3cc, hb3cc_tmp;
    
    typedef std::map<int, TimeseriesContainer > TrajMap;
    TrajMap traj_map;

    /**
     * Method that gives back true if the sum of the angles between the atoms
     * i, j and k is as large or larger than the minimal sum of angles given.
     */
    bool anglesum(int i, double &angle_sum, gmath::Vec &acceptor_j, gmath::Vec &acceptor_k) {
      double angles3;
      gmath::Vec tmpA, tmpB;
      tmpA = acceptor_j - donors.pos(i);
      tmpB = acceptor_k - donors.pos(i);
      angles3 = acos((tmpA.dot(tmpB)) / (tmpA.abs() * tmpB.abs()))*180 / M_PI;
      angle_sum += angles3;
      return (angle_sum >= min_angle_sum);
    }

    /**
     * Method that gives back true if the dihedral angle between the atoms i,
     * the one bound to i, j and k is as small or smaller than the maximal dihedral angle given.
     */
    bool dihedrals(int i, double &dihedral, gmath::Vec &bound_i, gmath::Vec &acceptor_j, gmath::Vec &acceptor_k) {
      gmath::Vec tmpA, tmpB, tmpC, p1, p2, p3;
      tmpA = bound_i - acceptor_j;
      tmpB = donors.pos(i) - acceptor_k;
      tmpC = acceptor_k - acceptor_j;

      p1 = tmpA.cross(tmpC);
      p2 = tmpB.cross(tmpC);

      dihedral = acos((p1.dot(p2)) / (p1.abs() * p2.abs()))*180 / M_PI;
      p3 = p1.cross(p2);
      if (p3.dot(tmpC) < 0)
        dihedral = -dihedral;
      return (abs(dihedral) <= max_dihedral);
    }
    /**
     * Method to copy or reference all AtomSpecifiers from HB2c_calc.
     */
    void initialise(const HB2c_calc& to_copy){
        donors = to_copy.get_donors();
        bound = to_copy.get_bound();
        acceptors = to_copy.get_acceptors();
        to_copy.get_num_A(num_A_donors, num_A_acceptors);
        to_copy.get_all_acc(accAB,accA,accB);
        to_copy.get_all_don(donAB,donA,donB);
    }
    /**
     * Method to initialize the H-bond calculation of each frame.
     */
    void init_calc(){
        ++frames;
        numHb=0;
        hb3cc_tmp.clear();
        ts.push_back(Timeseries<Key3c>(time));
    }
    /**
    * Method which prints a generic header.
    */
    void print_header() const;
    /**
    * Method to print one H-bond (identified by its key) to standard output.
    */
    inline void print(const Key3c&);

    /**
     * Method that loops over a vector of k atoms and calculates if there is a three-centered H-bond between
     * atoms i, j and k.
     */
    inline void calc(int, int, const std::vector<int>&);
    /**
    * Method that loops over all cubes and calls calc to calculate the H-bonds.
    */
    inline void go_through_cubes(CubeSystem<int>&, CubeSystem<int>&);

  public:

    /**
     * Constructor, which stores all parameters given from the input file.
     */
    HB3c_calc(HBPara3c para, bool red) :  HB_calc(red, para.maxdist2, para.minangle),
                                min_angle_sum(para.minanglesum),
                                max_dihedral(para.maxdihedral){
    }

    /**
     * Destructor, which closes the timeseries files.
     */
    ~HB3c_calc() {
      timeseriesHB.close();
      timeseriesHBtot.close();
      //delete pbc;
    }

    /**
     * Method to store the system file and the argument file for further use, and
     * opens the timeseries file
     */
    void setval(const HB2c_calc&, gcore::System& , args::Arguments&);

    /**
     * Method to clear the parameters calculated during the native H-bond calculation.
     */
    void clear();

    /**
     * Method, which only calculates H-bonds, that were present in the first frame of the reference file.
     */
    void calc_native();

    /**
     * Method that prints all 3-centered H-bonds and the timeseries files.
     */
    void printstatistics(bool,double);
    /**
     * Method which stores the native H-bonds.
     */
    void store_index();
    /**
     * Method that merges two HB3c_calc objects.
     */
    void merge(HB3c_calc&, int traj_num);
    /**
     * Method to calculate all H-bonds in a frame. A CubeSystem for donors and acceptors provides a grid-based pairlist.
     */
    void calc_hb(CubeSystem<int>& , CubeSystem<int>&);
    /**
     * Method to calculate all H-bonds of a frame in vacuum.
     */
    void calc_vac();

  }; //end class HB3c_calc
}
#endif	/* NEW_3C_H */
