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


#ifndef HB_CALC_BRIDGES
#define	HB_CALC_BRIDGES

#include <map>
#include <vector>
#include <iterator>

#include "Hbond_calc.h"
#include "Hbond_calc_2c.h"
#include "CubeSystem.hcc"

namespace utils {
  /**
   * Class Bridge
   * Stores information of a single solute-solvent-solute H-bond-bridge: how often it has occurred and its ID.
   * @author M.Setz
   * @ingroup utils
   * @class Bridge
   */
  class Bridge{
    int number, _id;
    public:
    /**
    * Default Constructor
    */
    Bridge() : number(0), _id(0)
    {}
    /**
     * Method, which increases the H-bond count.
     */
    void add(){
        ++number;
    }
    /**
     * Method, which sets the H-bond ID.
     */
    void set_id(int i){
        _id = i;
    }
    /**
     * Method, which returns the H-bond ID.
     */
    int id() const {
        return _id;
    }
    /**
     * Method, which adds two Bridge objects.
     */
    void add(const Bridge& rightop){
        number += rightop.number;
    }
    /**
     * Method, which returns the H-bond count.
     */
    int num() const{
        return number;
    }
  };

  /**
   * Class HB_bridges
   * Class, which inherits from HB_calc, to calculate the solute-solvent-solute H-bond-bridges.
   * If a bond is found it is stored in a map with Key3c as key and Bridge as value. It uses the two-centered H-bonds as input.
   * If two two-centered H-bonds share a common atom and this common atom is from a solvent molecule,
   * they are considered to form a solute-solvent-solute H-bond-bridge.
   * @author @ref ms
   * @ingroup utils
   * @class HB_bridges
   */
  class HB_bridges : public HB_calc {

    typedef std::map<Key3c, Bridge> BridgeContainer;
    typedef std::vector< BridgeContainer::iterator > BridgeContainerIteratorList;
    typedef std::vector<Timeseries<Key3c> > TimeseriesContainer;

    BridgeContainer bridges;
    std::vector<Key3c> native_key_storage;
    TimeseriesContainer ts;
    
    typedef std::map<int, TimeseriesContainer > TrajMap;
    TrajMap traj_map;

    /**
     * Method to copy or reference all AtomSpecifiers from HB2c_calc.
     */
    void initialise(const HB2c_calc& to_copy){
        donors = to_copy.get_donors();
        bound = to_copy.get_bound();
        acceptors = to_copy.get_acceptors();
        to_copy.get_num_A(num_A_donors, num_A_acceptors);
    }
    /**
     * Method to calculate if there is a solute-solvent-solute H-bond-bridge between two two-centered H-bonds.
     */
    inline void calc(const Key2c&, const Key2c&);
    /**
     * Method to print one solute-solvent-solute H-bond-bridge (identified by its key) to standard output.
     */
    inline void print(const Key3c&);
    /**
     * Method which prints a generic header.
     */
    void print_header() const;
    /**
     * Method to initialise calculation.
     */
    void init_calc(){
        numHb = 0;
        ++frames;
        ts.push_back(Timeseries<Key3c>(time));
    }


  public:
    /**
     * Constructor, which stores all parameters given from the input file.
      */
    HB_bridges(bool red): HB_calc(red)
    { }

    /**
     * Destructor, which closes the timeseries files.
     */
    ~HB_bridges() {
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
     * Method that prints all H-bond-bridges and the timeseries files.
     */
    void printstatistics(bool,double);
    /**
     * Method to clear the parameters calculated during the native H-bond calculation.
     */
    void clear();
    /**
     * Method that merges two HB_bridges objects .
     */
    void merge(const HB_bridges&, int traj_num);
    /**
     * Method to calculate all solute-solvent-solute H-bond-bridges in a frame from two-centered H-bonds.
     * A CubeSystem that stores the two-centered H-bonds provides a grid-based pairlist.
     */
    void calc_hb(const HB2c_calc&, CubeSystem<Key2c>&);
    /**
     * Method to calculate all solute-solvent-solute H-bond-bridges of a frame in vacuum.
     */
    void calc_vac(const HB2c_calc&);
    /**
     * Method, which only calculates solute-solvent-solute H-bond-bridges, that were present in the first frame of the reference file (="native" H-bond-bridges).
     */
    void calc_native();
    /**
     * Method which stores the native solute-solvent-solute H-bond-bridges.
     */
    void store_index();
  }; //end class HB_bridges
}
#endif
