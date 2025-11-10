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

#ifndef INCLUDED_UTILS_DISICL
#define INCLUDED_UTILS_DISICL

#include <map>
#include <vector>
#include <string>

#include "../gromos/Exception.h"
#include "../bound/Boundary.h"
#include "../gio/Ginstream.h"
#include "AtomSpecifier.h"
#include "CheckTopo.h"

using namespace std;
namespace utils
{
  /**
   * @class Dscl
   * @author M. Pechlaner
   * @ingroup utils
   *
   * @anchor DisiclLibrary
   * @section Disicl Disicl
   * Disicl classifies protein and nucleic acid secondary structure elements
   * based on dihedral angles. 
   *
   * - Nagy & Oostenbrink 2014, JCIM 54(1),278 
   * - Nagy & Oostenbrink 2014, JCIM 54(1),266
   *
   * @subsection LibraryFormat Library format
   * In addition to a TITLE block the DISICL library has to contain the 
   * following three blocks:
   * 
   * - DSCLANG (dihedral angle definitions)
   * An angle is defined by a unique name, followed by four atom names and 
   * relative residue numbers, i.e. in the example below angle PHI is defined by
   * the C of residue i-1, and by N, CA and C of residue i.
   *
   * @verbatim
   DSCLANG
   PHI      C    N    CA   C     -1   0   0   0 
   PSI      N    CA   C    N      0   0   0   1 
   END   @endverbatim
   * The number of dihedral angles is not limited to two.
   *
   * The selection of different atom names for different 
   * residues is possible in the following format: D;RES1,RES2:B;RES3:C , 
   * where RES1-3 are different residue names for which atoms B, B and C should
   * be used and D is the default atom name that should be used for all 
   * remaining residues.
   * @verbatim 
   CHI      O4*  C1*  N1;GUA,ADE:N9 C2;GUA,ADE:C4   0   0   0   0
   @endverbatim
   *
   *
   * - DSCLREG (region definitions)
   * A DISICL region is defined by a name and the lower and upper limits for 
   * each angle defined in block DSCLANG.
   *
   * @verbatim
   DSCLREG
   alfa1         -95      -40      -70      -32
   alfa2        -107      -40      -32      -12
   beta1        -135      -87       95      150
   beta2        -175     -135       95      136
   beta2        -180     -135      136      180
   END   @endverbatim
   * As many regions and subregions can be added as necessary.
   * 
   * - DSCLCLASS (class definitions)
   * A DISICL class is defined by a name, two DISICL region names and a short 
   * name.
   *
   @verbatim
   DSCLCLASS
   3/10-helix           alfa2    alfa2    3H    
   Alpha-helix          alfa1    alfa2    ALH   
   Normalbeta-strand    beta1    beta2    NBS  
   Mormalbeta-strand    beta2    beta1    NBS      
   Helix-cap            beta2    alfa1    HC     
   END   @endverbatim
   * As many different class definitions can be added as necessary.
   */
   
  class Dscl {
    unsigned int numAng, numReg, numClass, numResTot, numFrags;
    int minResShift, maxResShift;
    string const unclassString;
    
    vector<int> periodic;
    
    std::ofstream dihts, stats;
    std::string tsFile, statFile;
    map<string,ofstream* > classts;
    
    vector<string> libAngNames;
    vector<vector<string> > libAngAtoms;
    vector<vector<int> > libAngResShifts;
    vector<string> angAtomsUniq;
    
    vector<string> libRegionNames;
    vector<vector<double> >libRegions;
    map<string, vector<string> > libClassNames;
    map<string, string> classNameMap;
    map<string, int> classNumMap;
    vector<string> classShortnUniq;
    
    map<std::string, std::vector< int > > summary;
    
    vector<map<string,vector<string> > > vecAtomExists;
    vector<int> maxResSel; // maximal residue number for each molecule in the selection
    vector<int> maxAtom; // maximal atom number of each molecule
    vector<int> molSel; // molecules that occur in the atom selection
    vector<int> molStartRes; // starting residue for each molecule
    
    // a fragment contains the information for each set of consecutive residues 
    // from a single molecule in the atom selection
    struct fragment {
      unsigned int numRes;
      vector<PropertyContainer> propStore;
      vector<vector<string> > propNames;
      vector<int> resIds;
      vector<int> resNums;
      vector<int> molNums;
      vector<utils::AtomSpecifier> resAtoms;
      vector<vector<double> > torsions;
      vector<string> regions;
      vector<string> classes;
      vector<double> bfactors;
    };
    vector<fragment*> fragments;
    
    struct multAtom {
      string defaultAtom;
      string id;
      vector<string> atoms;
      vector<string> residues;
    };
    vector<multAtom> multAtoms;
    
    /**
     * modulo function transforming values to the region between min and max
     */
    double modulo(double value, int min, int max);
    vector<string> &split(const string &s, char delim, vector<string> &elems);
    vector<string> split(const string &s, char delim);
    
    void readLibAng(std::vector<std::string> &buffer);
    void readLibReg(std::vector<std::string> &buffer);
    void readLibClass(std::vector<std::string> &buffer);
    
    /**
     * collect all atoms of each selected residue
     */ 
    void getAtoms(gcore::System &sys);
  
  public:
    /**
     * default constructor
     */
    Dscl() : numAng(0), numReg(0), numClass(0), minResShift(9999), maxResShift(-9999), 
             unclassString("UC"), tsFile("ts_disicl.dat"), statFile("stat_disicl.out")
             {} 
       
    /**
     * Destructor.
     */
    ~Dscl() { 
      stats.close();
      for (unsigned int i=0; i<fragments.size(); i++) {
        delete fragments[i];
      }
    }
    
     
    // accessors
    int get_numAng() { return numAng; };
    int get_numClass() { return numClass; };
    int get_numReg() { return numReg; };
    vector<int> getPeriodic() { return periodic; };
    
    void print_angles();
    void print_regions();
    void print_classes();
    
    
    /**
     * read angles, regions and classes from library
     */
    void readLibrary(gio::Ginstream &lib);
    
    /**
     * set the periodic range for angle values
     */ 
    void setPeriodic(int min,int max);
    
    /**
     * assign DISICL regions based on the calculated dihedrals
     */ 
    void classifyRegions();
    
    /**
     * assign DISICL classes based on the regions 
     * and write class timeseries
     */ 
    void classifyClasses(double const &time);
    
    
    /**
     * write header for dihedral timeseries file
     */ 
    void writeHeader();
    
    /**
     *  initialize class classification statistics
     */ 
    void initSummary();
    
    /**
     * find those atoms from the given atom selection which are 
     * potentially needed for dihedral calculation
     */
    void determineAtoms(utils::AtomSpecifier &protein, gcore::System &sys);
    
    /**
     * create one property container for each angle specified in the library
     * and fill it with the dihedral angle specifications for each selected residue
     */ 
    void getPropSpec(gcore::System &sys, bound::Boundary * bound);
    
    /**
     * set B-factors for the output pdb according to class classification
     */ 
    void getBfactorValues(gcore::System &sys);
    
    /**
     * calculate dihedral angles
     */ 
    void calcDih();
    
    /**
     * write timeseries of dihedrals  
     * region and class classification are also written out as comments
     */ 
    void writeDihTs(double const &time);
    
    /**
     * write output pdb header
     */ 
    string pdbTitle();
    
    /**
     * build class classification statistics
     */ 
    void keepStatistics();
    
    /**
     * write class classification statistics to stat file
     * and dihedral time series statistics to the end of the dihedral
     * timeseries file
     */ 
    void writeStatistics(unsigned int  frameNum, bool do_tser);
    
    /**
     * initialize files for the class timeseries, one file per class
     */ 
    void initTimeseries();
    
    /**
     * close class timeseries files
     */ 
    void closeTimeseries();
    
    /**
     * write a pdb that can serve as a legend for the 
     * b-factor color code in the output pdbs
     */ 
    void writePdbColorLegend();
  };
}
#endif
