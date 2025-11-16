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
 * @file:   Hbond_calc_2c.h
 * Author: siggj
 *
 * Created on October 10, 2012, 2:16 PM
 */

#ifndef INCLUDED_2C_H
#define	INCLUDED_2C_H

#include <map>
#include <vector>
#include <iterator>

#include "CubeSystem.hcc"
#include "Hbond_calc.h"

namespace utils {

  /**
   * Struct HBPara2c
   * purpose: to store all the paramaters (maximal distance, minimal angle) given
   * by the input file.
   * @struct HBPara2c
   * @author J.Sigg
   * @ingroup utils
   */
  struct HBPara2c {
    double maxdist, maxdist2, minangle;
  }; //end struct HBPara2c


  /**
   * Class HB2c
   * Stores information of a single H-bond: how often it has occurred, its ID, the distance, and angle.
   * @author J.Sigg, M.Setz
   * @ingroup utils
   * @class HB2c
   */
  class HB2c {
    int _num, _id;
    double _dist, _angle;
  public:
    /**
     * Constructor
     */
    HB2c() : _num(0), _id(0), _dist(0), _angle(0)
    { }
    /**
     * Method, which adds the given values to all parameters, and increases the H-bond count.
     */
    void add(double distance, double ang){
        ++_num;
        _dist += sqrt(distance);
        _angle += ang;
    }
    /**
     * Method, which increases the H-bond count.
     */
    void add(){
        ++_num;
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
      return _num;
    }//end HB2c::getnum()
    /**
     * Method, which returns the H-bond ID.
     */
    int id() const {
        return _id;
    }
    /**
     * Method, which returns the mean distance.
     */
    double meandist() const{
      if(_num <= 0)
        throw gromos::Exception("Hbond_calc_2c","number of 3c Hbonds <= 0");
      return _dist/_num;
    }

    /**
     * Method, which returns the mean angle.
     */
    double meanangle() const{
      if(_num <= 0)
        throw gromos::Exception("Hbond_calc_2c","number of 3c Hbonds <= 0");
      return _angle/double(_num);
    }//end HB2c::getmeanangle()
    /**
     * Method, which adds two HB2c objects.
     */
    void add(const HB2c& rightop){ //HB2c
        _dist += rightop._dist; //distances in dist are already sqrt, so we can simply add them
        _angle += rightop._angle;
        _num += rightop._num;
    }
  }; //end class HB2c

  typedef std::map<Key2c, HB2c> HB2cContainer;

  /**
   * Class HB2c_calc
   * purpose: Class, which inherit from HB_calc, to calculate the 2-centred H-bonds.
   * If a bond is found it is stored in a map with Key3c as key and HB2c as value.
   * @class HB2c_calc
   * @author J. Sigg, M.Setz
   * @ingroup utils
   */
  class HB2c_calc : public HB_calc {

    typedef std::vector< HB2cContainer::iterator > HB2cContainerIteratorList;
    typedef std::vector<Timeseries<Key2c> > TimeseriesContainer;
    typedef std::map<int, TimeseriesContainer > TrajMap;

    HB2cContainer hb2cc, hb2cc_tmp;
    std::vector<Key2c> native_key_storage;
    TimeseriesContainer ts;
    TrajMap traj_map;

    /**
     * Method to initialise calculation.
     */
    void init_calc(){
        ++frames;
        numHb=0;
        hb2cc_tmp.clear();
        ts.push_back(Timeseries<Key2c>(time));
    }
    /**
     * Method to calculate if there is a H-bond between atom i and j.
     */
    inline void calc(int i, int j);//, const std::vector<int>& = std::vector<int>()); //default value for last argument= default constructor
    /**
     * Method which prints a generic header.
     */
    void print_header() const;
    /**
     * Method to print one H-bond (identified by its key) to standard output.
     */
    inline void print(const Key2c&);

    /**
    * Method that loops over all cubes and calls calc to calculate the H-bonds.
    */
    inline void go_through_cubes(CubeSystem<int>&, CubeSystem<int>&);
  public:
    /**
     * Constructor, calling parent constructor, which stores all parameters given from the input file.
     */
    HB2c_calc(HBPara2c para, bool red) : HB_calc(red, para.maxdist2, para.minangle)
    { }
    /**
     * Destructor, which closes the timeseries files.
     */
    ~HB2c_calc() {
      timeseriesHB.close();
      timeseriesHBtot.close();
      //delete pbc;
    }
    /**
     * Method to store the system file and the argument file for further use, and
     * opens the timeseries file
     */
    void setval(gcore::System &sys, args::Arguments &args, int);

    /**
     * Method to clear the parameters calculated during the native H-bond calculation.
     */
    void clear();

    /**
     * Method, which only calculates H-bonds, that were present in the first frame of the reference file.
     */
    void calc_native();

    /**
     * Method that prints all 2-centered H-bonds and the timeseries files.
     */
    void printstatistics(bool,double);

    /**
     * Method that merges two HB2c_calc objects. 
     */
    void merge(HB2c_calc&, int traj_num);
    /**
     * Method which stores the native H-bonds.
     */
    void store_index();

    /**
     * Method to calculate all H-bonds in a frame. A CubeSystem for donors and acceptors provides a grid-based pairlist.
     */
    void calc_hb(CubeSystem<int>&, CubeSystem<int>&);
    /**
     * Method to calculate all H-bonds of a frame in vacuum.
     */
    void calc_vac();
    /**
     * Method which returns a reference to the map that contains all H-bonds of all frames.
     */
    const HB2cContainer& get_whole_map() const {
        return hb2cc;
    }
    /**
     * Method which return a reference to a map that contains all H-bonds OF A SINGLE FRAME.
     */
    const HB2cContainer& get_tmp_map() const {
        return hb2cc_tmp;
    }
  }; //end class HB2c_calc

}
#endif	/* NEW_2C_H */
