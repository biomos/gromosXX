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



#ifndef INCLUDED_UTILS_DSSP
#define INCLUDED_UTILS_DSSP

#include <vector>
#include <string>
#include <fstream>

#include "AtomSpecifier.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../gcore/System.h"

namespace gcore{
  class System;
  class Molecule;
  class MoleculeTopology;
}
namespace gmath
{
  class Vec;
}
namespace args
{
  class Arguments;
}
namespace bound
{
  class Boundary;
}
namespace utils
{

  /**
   * Class Dssp
   * purpose: (calculate (intramolecular) backbone-backbone 
   * hydrogen bonds, and) define secondary structure of protein (DSSP)
   * over the trajectory.
   *
   * Description:
   * Dssp within gromos++ defines secondary structure for one single (protein) solute 
   * molecule, according to the DSSP rules defined by W. Kabsch and C. Sander 
   * (Biopolymers, 22, pp2577-2637 (1983)). Within these rules it may occur that one 
   * residue is defined as two different secondary-structure elements. In order to 
   * avoid duplicates in the output, the following (ad hoc) priority rules are applied 
   * here: Beta Sheet/Bridge > 5-helix > 4-helix > 3-helix > H-bonded turn > Bend. 
   * As a consequence, there may be, for instance, helices that are shorter than 
   * their minimal length. Thes rules are easily changed in Dssp::filter_SecStruct().
   * 
   * @class Dssp
   * @author U B
   * @ingroup utils
   * @sa utils::AtomSpecifier
   */
  class Dssp{
    std::vector<int> d_mol, d_atom, num, tsnum;
    std::vector<double> dist, ang, tstime;
    std::vector<int> acc_res, don_res;
    std::vector<int> helix3, helix4, helix5, Helix, helix;
    std::vector<int> bridge, extended;
    std::vector<int> Turn, turn, Bend, Beta;
    std::vector<std::vector< int > > summary;
    int d_numFrames;
    int numres;
    std::vector<int> d_resnum;
    std::vector<int> d_resOffSets;
    args::Arguments *d_args;
    gcore::System *d_sys;
    utils::AtomSpecifier d_H, d_N, d_O, d_C, d_CA;
    bound::Boundary *d_pbc; 
    std::ofstream timeseriesTurn, timeseries3Helix, timeseries4Helix, timeseries5Helix;
    std::ofstream timeseriesBBridge, timeseriesBStrand, timeseriesBend;
    int d_nummol;
    bool d_omit_self_species;
   
   
  public: 
    // Constructor
    /**
     * Dssp Constructor
     * @param sys The Dssp needs to know about the system. It 
     *            does not know about any atoms yet.
     * @param args all arguments are passed into Dssp. 
     */
    Dssp(gcore::System &sys, args::Arguments &args);

    /**
     * Dssp Deconstructor
     */
    virtual ~Dssp(){
      timeseriesTurn.close();
      timeseries3Helix.close();
      timeseries4Helix.close();
      timeseries5Helix.close();
      timeseriesBBridge.close();
      timeseriesBStrand.close();
      timeseriesBend.close();
	}

     /**
     * Method to determine the atoms specified using the AtomSpecifier.
     */
     void determineAtoms(utils::AtomSpecifier &protein);
     /**
     * Method to initialize the calculation for intramolecular hydrogen bonds.
     */
     void calcHintra_init(utils::AtomSpecifier &protein);
     /**
     * Method to calculate intramolecular hydrogen bonds over one frame.
     */
     virtual void calcHb_Kabsch_Sander();
     /**
      * Method to define 3-, 4-, and 5-helices over one frame.
      */
     void calc_Helices();
     /**
      * Method to define beta structures over one frame.
      */
     void calc_Betas();
     /**
      * Method to define bends over one frame.
      */
     void calc_Bends();
     /**
      * Method to filter secondary-structure output.
      */
     void filter_SecStruct();
     /**
      * Method to write secondary-structure output to files.
      */
    void writeToFiles(double time);     
    /**
     * Method to keep the statistics for later output
     */
    void keepStatistics();
    /**
     * Method to print the statistics
     */
    void writeSummary(std::ostream & of);
    /**
     * Method to calculate de number of residues
     */
    void calcnumres(utils::AtomSpecifier &protein, const gcore::System & sys);
    
    typedef void (Dssp::*MemPtr)();
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */ 
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says Dssp, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what): gromos::Exception("Dssp", what){}
    };
  protected:
    /**
     * Method that opens the timeseries files.
     * 
     */
    void opents(std::string fi1, std::string fi2, std::string fi3, 
		std::string fi4, std::string fi5, std::string fi6, 
		std::string fi7);
    /**
     * Method that reads a frame from either a reference coordinate file
     * or the first frame of the first trajectory file.
     */    
    void readframe();
  }; //end class Dssp
}

#endif
