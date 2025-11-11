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
 * @file cgLJpot.cc
 * used to develop a coarse-grained potential from an all-atom simulation
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor cgLJpot
 * @section cgLJpot calculates coarse-grained LJ potentials form fine grained simulations
 * @author @ref ae
 * @date 07.11.2011
 *
 * Program cgLJpot calculates the Lennard-Jones potential energy function for a coarse-grained system based on a fine-grained (all-atom) simulation trajectory
 * (V_fg2cg). The resulting potential energy is not a 12/6-Lennard-Jones potential. However, the program also calculates the coarse-grained 12/6-Lennard-Jones
 * potential energy function as an approximation to V_fg2cg, both having a minimum at the same energy and r value, i.e. the sigma- and epsilon-value for the
 * approximation to V_fg2cg is read from the minimum of V-fg2cg. Therefor, in the case the beads do not include atoms of different molecules, the calculated
 * 12/6-Lennard-Jones potential is expected to reproduce the density and heat of vaporization of the fine-grained simulation in a coarse-grained simulation.
 * However, practice showed that this is normally not the case and parameterization is still necessary, but the Lennar-Jones parameters are reasonable to start
 * the parameterization.
 * 
 * NOTE: the current program version does only work for solute
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@method</td><td>&lt;method to goarse grain: atomic or molecular&gt; </td></tr>
 * <tr><td> [\@dist</td><td>&lt;min max ngrid&gt;] </td></tr>
 * <tr><td> \@beads</td><td>&lt;number of atoms per bead (atomic)&gt; or </td></tr>
 * <tr><td>        </td><td>&lt;sequence of bead size within one molecule (molecular)&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> \@trc</td><td>&lt;simulation trajectory or coordinate file&gt;</td></tr>
 * </table>
 *
 * 
 * Example:
 * @verbatim
   rdf
     @topo   ex.top
 * @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <ctime>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomPair.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Distribution.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/System.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"

namespace cgLJpot {

  /**
   * A class to store a pair of integers, e.g. the two IAC numbers of two particles interacting via a LJ potential
   */
  class IJ {
  private:
    /**
     * the first integer number 
     */
    int I;
    /**
     * the second integer number
     */
    int J;

  public:
    /**
     * constructor
     * @param i first integer number; standard: -1
     * @param j second integer number; standard: -1
     */
    IJ(int i = -1, int j = -1);
    /**
     * copy constructor 
     * @param ij the IJ class to be copied
     */
    IJ(const IJ &ij);
    void setValues(int i, int j);
    int i() const;
    int j() const;
  };
  
  /**
   * A class to store a quartet of integers, e.g. the two IAC numbers of two particles interacting via a LJ potential
   */
  class IJK {
  private:
    /**
     * the first integer number 
     */
    int I;
    /**
     * the second integer number
     */
    int J;
    /**
     * the third integer number
     */
    int K;


  public:
    /**
     * constructor
     * @param i first integer number; standard: -1
     * @param j second integer number; standard: -1
     * @param k third integer number; standard: -1
     */
    IJK(int i = -1, int j = -1, int k = -1);
    /**
     * copy constructor 
     * @param ijk the IJK class to be copied
     */
    IJK(const IJK &ijk);
    void setValues(int i, int j, int k);
    int i() const;
    int j() const;
    int k() const;
  };
  
  /**
   * A class to store a quartet of integers, e.g. the two IAC numbers of two particles interacting via a LJ potential
   */
  class IJKL {
  private:
    /**
     * the first integer number 
     */
    int I;
    /**
     * the second integer number
     */
    int J;
    /**
     * the third integer number
     */
    int K;
    /**
     * the fourth integer number
     */
    int L;

  public:
    /**
     * constructor
     * @param i first integer number; standard: -1
     * @param j second integer number; standard: -1
     * @param k third integer number; standard: -1
     * @param l fourth integer number; standard: -1
     */
    IJKL(int i = -1, int j = -1, int k = -1, int l = -1);
    /**
     * copy constructor 
     * @param ijkl the IJKL class to be copied
     */
    IJKL(const IJKL &ijkl);
    void setValues(int i, int j, int k, int l);
    int i() const;
    int j() const;
    int k() const;
    int l() const;
  };

  /**
   * just a < operator to make the class more complete... 
   * this is needed to iterate over a map using @ref IJs as the key value
   * @param ij1 first @ref IJ to be compared with...
   * @param ij2 ... the second @ref IJ
   */
  bool operator<(const IJ &ij1, const IJ &ij2);

  /** @class LJpot
   * stores the Lennard-Jones-potential of multiple pairs of particles on a grid */
  class LJpot {
  private:
    /**
     * the summed up LJ potential energies at each grid
     */
    vector<double> lj;
    /**
     * number of summed up potential energies at each grid
     */
    vector<int> count;
    /**
     * grid size of the LJ potential energy
     */
    double dgrid;
    /**
     * the minimum distance between two particles, pairs of particles with a smaller distance are not considered
     */
    double min;
    /**
     * the maximum distance between two particles, pairs of particles with a longer distance are not considered
     */
    double max;
    /**
     * the Lennard-Jones C12 parameter
     */
    double C12;
    /**
     * the Lennard-Jones C6 parameter
     */
    double C6;
    /**
     * indicates if there was any value added tot the potential
     */
    bool used;

  public:
    /**
     * constructor
     * @param min_ minimum interatomic distance
     * @param max_ maximum interatomic distance
     * @param grid number of grid points for interatomic distance
     * @param c12 Lennard-Jones C12 parameter
     * @param c6 Lennard-Jones C6 parameter
     */
    LJpot(double min_, double max_, int grid, double c12, double c6);
    /**
     * constructor
     * @param min_ minimum interatomic distance
     * @param max_ maximum interatomic distance
     * @param grid number of grid points for interatomic distance
     */
    LJpot(double min_ = 0.0, double max_ = 2.0, int grid = 200);
    /**
     * copy constructor
     * @param ljp a Lennard-Jones potential class of type @ref LJpot
     */
    LJpot(const LJpot &ljp);
    /**
     * initializer
     * @param min_ minimum interatomic distance
     * @param max_ maximum interatomic distance
     * @param grid number of grid points for interatomic distance
     */
    void init(double min_ = 0.1, double max_ = 2.0, int grid = 200);
    /**
     * adds the Lennard-Jones energy of two particles to the stored potential energy
     * @param pos inter-particle distance
     * @param val Lennard-Jones energy
     */
    void add(double pos, double val);
    /**
     * unifies two @ref LJpot classes into one
     * @param ljp the LJpot class which will be unified with the current one
     */
    void unify(const LJpot &ljp);
    /**
     * accessor: returns the inter-atomic distance r corresponding to grid point i
     * @param i grid point number
     * @return interatomic distance r
     */
    double r(unsigned int i);
    /**
     * accessor: returns the inter-atomic distance r corresponding to grid point i
     * @param i grid point number
     * @return interatomic distance r
     */
    double r(unsigned int i) const;
    /**
     * accessor: returns the (averaged) Lennard-Jones potential corresponding to grid point i, i.e. the summed up Lennard-Jones potential energy divided by
     * the number if summed up terms at grid point i
     * @param i grid point number
     * @return averaged Lennard-Jones potential energy
     */
    double pot(unsigned int i);
    /**
     * accessor: returns the minimum inter-particle distance to be considered
     * @return minimum inter-particle distance to be considered
     */
    double get_min();
    /**
     * accessor: returns the maximum inter-particle distance to be considered
     * @return maximum inter-particle distance to be considered
     */
    double get_max();
    /**
     * accessor: returns the number of grid points
     * @return number of grid points
     */
    int get_grid();
    /**
     * accessor: searches the inter-particle distance r_max for which the Lennard-Jones potential is maximal and returns the inter-particle distance r > r_max
     * for which the Lennard-Jones potential is minimal
     * @return inter-particle distance r for which the Lennard-Jones energy is minimal
     */
    double rmin();
    /** 
     * accessor: searches the inter-particle distance r_max for which the Lennard-Jones potential is maximal and returns the minimum Lennard-Jones potential
     * energu V_LJ(r) with r > r_max
     * @return inter-particle distance r for which the Lennard-Jones energy is minimal
     */
    double potmin();
    /*
     * returns true if the potential was used, false otherwise
     */
    bool wasUsed() const;
  };

  /**
   * stores multiple atoms to act as a bead as well as the different interaction Lennard-Jones potential energies to other beads
   */
  class bead {
  private:
    /**
     * a reference to the gromos force field
     */
    gcore::GromosForceField *gff;
    /**
     * the atoms which ar building this bead
     */
    utils::AtomSpecifier atoms;
    /**
     * the center of geometry or center of mass of all the atoms of the bead -> center of the bead
     */
    gmath::Vec centre;
    /**
     * in case all atoms of this bead belong to the same molecule this number indicates the molecule number
     */
    int memberOfMol;
    /**
     * a sequential number to tag the beads, either within a molecule or within the whole system
     */
    int beadnum;
    /**
     * stores the total Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> totLJ_ee;
    map<IJ, LJpot> totLJ_em;
    map<IJ, LJpot> totLJ_mm;
    /**
     * stores the total inter-molecular Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J,
     * respectively
     */
    map<IJ, LJpot> totinterLJ_ee;
    map<IJ, LJpot> totinterLJ_em;
    map<IJ, LJpot> totinterLJ_mm;
    /**
     * stores the total intra-molecular Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> totintraLJ_ee;
    map<IJ, LJpot> totintraLJ_em;
    map<IJ, LJpot> totintraLJ_mm;
    /**
     * stores the inter-molecular 12-neighboring Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> intra12LJ_ee;
    map<IJ, LJpot> intra12LJ_em;
    map<IJ, LJpot> intra12LJ_mm;
    /**
     * stores the inter-molecular 13-neighboring Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> intra13LJ_ee;
    map<IJ, LJpot> intra13LJ_em;
    map<IJ, LJpot> intra13LJ_mm;
    /**
     * stores the inter-molecular 14-neighboring Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> intra14LJ_ee;
    map<IJ, LJpot> intra14LJ_em;
    map<IJ, LJpot> intra14LJ_mm;
    /**
     * a bool to know if this bead is a head- or tail bead
     */
    bool istail;

  public:
    /**
     * constructor
     * @param sys the system
     * @param groff the GROMOS force field
     * @param mom the molecule the bead is a member of
     * @param bnum the sequential bead number
     * @param ij the possibly different IJ combinations (to other beads)
     * @param min the minimum distance to another bead for which the Lennard-Jones potential energy is calculated and stored
     * @param max the maximum distance to another bead for which the Lennard-Jones potential energy is calculated and stored 
     * @param grid the number of grid points of the Lennard-Jones potentials
     */
    bead(gcore::System &sys, gcore::GromosForceField &groff, int mom, int bnum, set<IJ> &ij, double min = 0.0, double max = 2.0, double grid = 200);
    /**
     * copy constructor
     * @param b another bead to be copied
     */
    bead(bead const &b);
    /**
     * destructor
     */
    ~bead() {
    };
    /**
     * adds an atom to the current bead
     * @param m molecule number of the atom to be added
     * @param a atom number of the atom to be added
     * @return total number of atoms within the bead
     */
    int addAtom(int m, int a);
    /**
     * returns the number of atoms within that bead
     */
    int size();
    /**
     * calculates, defines and returns the center of geometry of the bead
     * @param pbc periodic boundary, see @ref bound::Boundary
     * @param sys the system
     */
    gmath::Vec cog(bound::Boundary *pbc, gcore::System &sys);
    /**
     * calculates, defines and returns the center of mass
     * @param pbc periodic boundary, see @ref bound::Boundary
     * @param sys the system
     */
    gmath::Vec com(bound::Boundary *pbc, gcore::System &sys);
    /**
     * accessor to the (center) position of the bead (com or cog)
     * */
    gmath::Vec pos();
    /**
     * add a value to the total LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJtot_ee(const IJ &ij, double r, const double &lj);
    void addLJtot_em(const IJ &ij, double r, const double &lj);
    void addLJtot_mm(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJtotinter_ee(const IJ &ij, double r, const double &lj);
    void addLJtotinter_em(const IJ &ij, double r, const double &lj);
    void addLJtotinter_mm(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total intramolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJtotintra_ee(const IJ &ij, double r, const double &lj);
    void addLJtotintra_em(const IJ &ij, double r, const double &lj);
    void addLJtotintra_mm(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total 12-intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJintra12_ee(const IJ &ij, double r, const double &lj);
    void addLJintra12_em(const IJ &ij, double r, const double &lj);
    void addLJintra12_mm(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total 13-intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJintra13_ee(const IJ &ij, double r, const double &lj);
    void addLJintra13_em(const IJ &ij, double r, const double &lj);
    void addLJintra13_mm(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total 14-intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJintra14_ee(const IJ &ij, double r, const double &lj);
    void addLJintra14_em(const IJ &ij, double r, const double &lj);
    void addLJintra14_mm(const IJ &ij, double r, const double &lj);
    /** calculate the LJ interaction to another bead
     * @param b a @ref bead
     * @param pbc
     * @param sys
     * @return 
     */
    double calcLJ(bead &b, bound::Boundary *pbc, gcore::System &sys);
    /**
     * accessor: returns the total LJpot energy
     * @return total LJpot energy
     */
    map<IJ, LJpot> get_totLJ_ee();
    map<IJ, LJpot> get_totLJ_em();
    map<IJ, LJpot> get_totLJ_mm();
    /**
     * accessor: returns the total intermolecular LJpot energy
     * @return total intermolecular LJpot energy
     */
    map<IJ, LJpot> get_totinterLJ_ee();
    map<IJ, LJpot> get_totinterLJ_em();
    map<IJ, LJpot> get_totinterLJ_mm();
    /**
     * accessor: returns the total intramolecular LJpot energy
     * @return total intramolecular LJpot energy
     */
    map<IJ, LJpot> get_totintraLJ_ee();
    map<IJ, LJpot> get_totintraLJ_em();
    map<IJ, LJpot> get_totintraLJ_mm();
    /**
     * accessor: returns the 12-intermolecular LJpot energy
     * @return 12-intermolecular LJpot energy
     */
    map<IJ, LJpot> get_intra12LJ_ee();
    map<IJ, LJpot> get_intra12LJ_em();
    map<IJ, LJpot> get_intra12LJ_mm();
    /**
     * accessor: returns the 13-intermolecular LJpot energy
     * @return 13-intermolecular LJpot energy
     */
    map<IJ, LJpot> get_intra13LJ_ee();
    map<IJ, LJpot> get_intra13LJ_em();
    map<IJ, LJpot> get_intra13LJ_mm();
    /**
     * accessor: returns the 14-intermolecular LJpot energy
     * @return 14-intermolecular LJpot energy
     */
    map<IJ, LJpot> get_intra14LJ_ee();
    map<IJ, LJpot> get_intra14LJ_em();
    map<IJ, LJpot> get_intra14LJ_mm();
    /**
     * accessor: returns the molecule number the bead/atoms of the bead belong to
     */
    int mol();
    /**
     * accessor: returns the bead number within the molecule (sequential bead number)
     */
    int beadNum();
    /*
     * set this bead as head or tail bead
     */
    void setAsTail();
    /*
     * accessor to know if this bead is an endbead or not
     */
    bool isTail();
  };

  /**
   * prints the different potentials in columns where the first column is the inter-particle distance r
   * @param os the output stream
   * @param totLJpot the total Lennard-Jones potential
   * @param totLJinter the total inter-particle Lennard-Jones potential
   * @param totLJintra the total intra-particle Lennard-Jones potential
   * @param totLJ12 the total 12-Lennard-Jones potential
   * @param totLJ13 the total 13-Lennard-Jones potential
   * @param totLJ14 the total 14-Lennard-Jones potential
   */
  void printPot(ostream &os, std::map<cgLJpot::IJ, LJpot> &totLJpot,
          std::map<IJ, LJpot> &totLJinter,
          std::map<IJ, LJpot> &totLJintra,
          std::map<IJ, LJpot> &totLJ12,
          std::map<IJ, LJpot> &totLJ13,
          std::map<IJ, LJpot> &totLJ14);
  
  void printBeadBeadDist(std::string fname, std::map<IJ, gmath::Distribution> &beadbeadDist, std::set<IJ> IJs, double rmin, double rmax, int grid);

  /**
   * prints some header information to remember what the program was analyzing
   */
  void printTitleBlock(std::vector<int> &beadsizes, utils::AtomSpecifier &allatoms);
  /**
   * calculates the sigma and epsilon values of a Lennard-Jones potential @ref LJpot
   */
  void calcEpsSigma(std::map<IJ, double> &epsilons, std::map<IJ, double> &sigmas, std::map<IJ, LJpot> &LJ);
  /**
   * calculates the Lennard-Jones C12 and C6 parameters based on the corresponding epsilon- and sigma values
   */
  void calcC12C6(std::map<IJ, double> &C12s, std::map<IJ, double> &C6s, std::map<IJ, double> &epsilons, std::map<IJ, double> &sigmas);
  /**
   * prints the various Lennard-Jones parameters
   */
  void printLennardJonesParamters(string title, map<IJ, double> &epsilons_tot,
          map<IJ, double> &epsilons_totinter,
          map<IJ, double> &epsilons_totintra,
          map<IJ, double> &epsilons_intra12,
          map<IJ, double> &epsilons_intra13,
          map<IJ, double> &epsilons_intra14,
          map<IJ, double> &sigmas_tot,
          map<IJ, double> &sigmas_totinter,
          map<IJ, double> &sigmas_totintra,
          map<IJ, double> &sigmas_intra12,
          map<IJ, double> &sigmas_intra13,
          map<IJ, double> &sigmas_intra14,
          map<IJ, double> &C12_tot,
          map<IJ, double> &C12_totinter,
          map<IJ, double> &C12_totintra,
          map<IJ, double> &C12_intra12,
          map<IJ, double> &C12_intra13,
          map<IJ, double> &C12_intra14,
          map<IJ, double> &C6_tot,
          map<IJ, double> &C6_totinter,
          map<IJ, double> &C6_totintra,
          map<IJ, double> &C6_intra12,
          map<IJ, double> &C6_intra13,
          map<IJ, double> &C6_intra14);
  
}

using namespace args;
using namespace bound;
using namespace cgLJpot;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "method" << "beads" << "pbc" << "trc" << "dist";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@method        <method to goarse grain: atomic or molecular>\n";
  usage += "\t@beads         <number of atoms per bead (atomic)> or\n";
  usage += "\t               <sequence of bead size within one molecule (molecular)>\n";
  usage += "\t[@pbc          <boundary type (read from GENBOX block if not specified)> [<gather method>]]\n";
  usage += "\t@trc           <simulation trajectory or coordinate file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);
	  // this program is under development
	  args.underDevelopment();

    // read the method:
    // - atomic: the program does not care about the moleculs and where they start and end
    // - the program cares about the molecules (no beads over more than one molecule
    args.check("method", 1);
    string method = args["method"];
    if (method != "atomic" && method != "molecular") {
      stringstream msg;
      msg << "method \"" << method << "\" (@method) not implemented, chose \"atomic\" or \"molecular\"";
      throw gromos::Exception(argv[0], msg.str());
    }
    
    string fname_LJpot_FG2CG_ee = "LJpot_FG2CG_ee.dat";
    string fname_LJpot_FG2CG_em = "LJpot_FG2CG_em.dat";
    string fname_LJpot_FG2CG_mm = "LJpot_FG2CG_mm.dat";
    string fname_LJpot_CG_ee = "LJpot_CG_ee.dat";
    string fname_LJpot_CG_em = "LJpot_CG_em.dat";
    string fname_LJpot_CG_mm = "LJpot_CG_mm.dat";
    string fname_beadbead_dist = "bead-bead_dist.dat";
    string fname_beadbead_dist_ee = "bead-bead_dist_ee.dat";
    string fname_beadbead_dist_em = "bead-bead_dist_em.dat";
    string fname_beadbead_dist_mm = "bead-bead_dist_mm.dat";

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff = it.forceField();
    
    // read the bead size
    args.check("beads", 1);
    vector<int> beadsizes;
    {
      Arguments::const_iterator start = args.lower_bound("beads");
      Arguments::const_iterator stop = args.upper_bound("beads");
      Arguments::const_iterator it;
      for (it = start; it != stop; ++it) {
        stringstream ss;
        ss << it->second;
        int b;
        ss >> b;
        if (ss.fail() || ss.bad() || !ss.eof()) {
          stringstream msg;
          msg << "cannot use " << it->second << " as a bead size";
          throw gromos::Exception(argv[0], msg.str());
        }
        beadsizes.push_back(b);
      }
    }
    if (method == "atomic" && beadsizes.size() != 1) {
      throw gromos::Exception(argv[0], "method \"atomic\" (@atom) does not allow for more than one bead size");
    }

    // read the distribution parameters
    double distmin = 0.0;
    double distmax = 2.0;
    int distgrid = 200;
    if (args.count("dist") >= 0) {
      if (args.count("dist") == 3) {
        stringstream ss;
        Arguments::const_iterator it = args.lower_bound("dist");
        ss << it->second << endl;
        ++it;
        ss << it->second << endl;
        ++it;
        ss << it->second;
        ss >> distmin >> distmax >> distgrid;
        if (ss.fail() || ss.bad() || !ss.eof()) {
          stringstream msg;
          msg << "cannot convert the arguments of @dist to min, max and number of grid points";
          throw gromos::Exception(argv[0], msg.str());
        }
      } else {
        stringstream msg;
        msg << "cannot convert the arguments of @dist to min, max and number of grid points";
        throw gromos::Exception(argv[0], msg.str());
      }
      if (distmin >= distmax) {
        throw gromos::Exception(argv[0], "distmin >= distmax but should be distmin < distmax");
      }
    }

    // do some checks depending on the coarse-grain method
    int numMol = sys.numMolecules();
    int numAtMol = sys.mol(0).numAtoms();
    int checkMol = -1; // -1 if all molecules have the same number of atoms, otherwise
    // it gets is assigned to the molecule number of the first molecule
    // having a different atom number
    AtomSpecifier allAtoms(sys);
    for (int m = 0; m < numMol; ++m) {
      int numAt = sys.mol(m).numAtoms();
      if (numAt != numAtMol && checkMol == -1) {
        checkMol = m;
      }
      for (int a = 0; a < numAt; ++a) {
        allAtoms.addAtom(m, a);
      }
    }
    // checks for method = atomic:
    // is the total number of atoms "dividable" (without rest) by the bead size?
    int numBeads = 0;
    if (method == "atomic") {
      if (allAtoms.size() % beadsizes[0] != 0) {
        stringstream msg;
        msg << "the total number of atoms (" << allAtoms.size() << ") cannot be divided"
                " by the bead size (" << beadsizes[0] << "): " << allAtoms.size() % beadsizes[0]
                << " atoms left";
        throw gromos::Exception(argv[0], msg.str());
      }
      numBeads = allAtoms.size() / beadsizes[0];
    }
    // checks for method = molecular
    if (method == "molecular") {
      if (checkMol != -1) {
        stringstream msg;
        msg << "molecule " << checkMol + 1 << " has a different number of atoms"
                " than the previous molecules: method \"molecular\" (@method) "
                "does not work therefore";
        throw gromos::Exception(argv[0], msg.str());
      }
      int numAtMolBeads = 0;
      for (unsigned int b = 0; b < beadsizes.size(); ++b) {
        numAtMolBeads += beadsizes[b];
      }
      if (numAtMol != numAtMolBeads) {
        stringstream msg;
        msg << "the number of atoms per molecule is " << numAtMol << ", but the"
                " different bead sizes sum up to " << numAtMolBeads;
        throw gromos::Exception(argv[0], msg.str());
      }
      numBeads = allAtoms.size() / numAtMolBeads * beadsizes.size();
    }

    // get all possible combinations of beads (with respect to its size
    set<IJ> IJs;
    set<IJK> IJKs;
    set<IJKL> IJKLs;
    for (unsigned int bs1 = 0; bs1 < beadsizes.size(); ++bs1) {
      for (unsigned int bs2 = bs1; bs2 < beadsizes.size(); ++bs2) {
        IJ ij(beadsizes[bs1], beadsizes[bs2]);
        if (IJs.find(ij) == IJs.end()) {
          IJs.insert(ij);
        }
	  }
	  for (unsigned int bs2 = 0; bs2 < beadsizes.size(); ++bs2) {
        for (unsigned int bs3 = 0; bs3 < beadsizes.size(); ++bs3) {
          IJK ijk(beadsizes[bs1], beadsizes[bs2], beadsizes[bs3]);
          if (IJKs.find(ijk) == IJKs.end()) {
            IJKs.insert(ijk);
            for (unsigned int bs4 = 0; bs4 < beadsizes.size(); ++bs4) {
              IJKL ijkl(beadsizes[bs1], beadsizes[bs2], beadsizes[bs3], beadsizes[bs4]);
              if (IJKLs.find(ijkl) == IJKLs.end()) {
                IJKLs.insert(ijkl);
              }
            }
          }
        }
      }
    }

    // construct the beads and count the number of beads
    int endBeadNum = 0;
    int middleBeadNum = 0;
    vector<bead> beads;
    {
      int a = 0; // the current position in allAtoms
      while (a < allAtoms.size()) {
        int bnum = 0;
        for (int bs = 0; bs < (int) beadsizes.size(); ++bs, ++bnum) {
          bead B(sys, gff, allAtoms.mol(a), bnum, IJs, distmin, distmax, distgrid);
          for (int i = 0; i < beadsizes[bs]; ++i) {
            int mol = allAtoms.mol(a);
            int at = allAtoms.atom(a);
            B.addAtom(mol, at);
            ++a;
          }
          if(bs == 0 || bs == (int)(beadsizes.size() - 1)) {
            B.setAsTail();
            endBeadNum++;
          } else {
            middleBeadNum++;
          }
          beads.push_back(B);
        }
      }
    }
    
    printTitleBlock(beadsizes, allAtoms);

    // a map of distributions (for each IJ one) to remember the intramolecular
    // neighboring bead-bead distance (min and max automatically set based on the
    // first configuration of the trajectories
    map<IJ, Distribution> beadbeadDist, beadbeadDist_ee, beadbeadDist_em, beadbeadDist_mm;
    map<IJK, Distribution> angleDist;
    map<IJKL, Distribution> dihedralDist;
    double bondlength_min = sys.box().K().abs() + sys.box().L().abs() + sys.box().M().abs();
    double bondlength_max = 0.0;
    {
      Arguments::const_iterator trcfirs = args.lower_bound("trc");
      InG96 ic;
      ic.open(trcfirs->second.c_str());
      ic.select("SOLUTE");
      ic >> sys;
      ic.close();
      bound::Boundary *pbc;
      if (args.count("pbc") > 0) { // read from arguments
        pbc = args::BoundaryParser::boundary(sys, args);
      } else { // read from GENBOX block
        if (args.count("pbc") == 0) {
          cerr << "WARNING: @pbc given with no argument(s), reading boundary type "
                  "from GENBOX block of the trajectory/coordinate file\n";
        }
        pbc = args::BoundaryParser::boundary(sys);
      }
      // loop over the beads and get the min/max distance
      for (int b1 = 0; b1 < (int) beads.size() - 1; ++b1) {
        beads[b1].com(pbc, sys);
        for (int b2 = b1 + 1; b2 < (int) beads.size(); ++b2) {
          if (beads[b1].mol() == beads[b2].mol() &&
                  abs(beads[b1].beadNum() - beads[b2].beadNum()) == 1) {
            beads[b2].com(pbc, sys);
            double r = (beads[b1].pos() - pbc->nearestImage(beads[b1].pos(), beads[b2].pos(), sys.box())).abs();
            if (r < bondlength_min) {
              bondlength_min = r;
            }
            if (r > bondlength_max) {
              bondlength_max = r;
            }
          } else {
            continue;
          }
        }
      }
      double min = bondlength_min - (1.2 * bondlength_max - bondlength_max);
      bondlength_max *= 1.2;
      bondlength_min = min > 0 ? min * 0.8 : 0.0;
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        beadbeadDist.insert(pair<IJ, Distribution > (*it, Distribution(bondlength_min, bondlength_max, 1000)));
        beadbeadDist_ee.insert(pair<IJ, Distribution > (*it, Distribution(0.0, 2.0, 1000)));
        beadbeadDist_em.insert(pair<IJ, Distribution > (*it, Distribution(bondlength_min, bondlength_max, 1000)));
        beadbeadDist_mm.insert(pair<IJ, Distribution > (*it, Distribution(bondlength_min, bondlength_max, 1000)));
      }
      for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
        angleDist.insert(pair<IJK, Distribution> (*it, Distribution(0, 180, 180)));
      }
      for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
        dihedralDist.insert(pair<IJKL, Distribution> (*it, Distribution(0, 180, 180)));
      }
    }

    // a distribution for the RDFs
    double rdf_cut = 3;
    int rdf_grid = 300;
    vector<double> rdf_ee(rdf_grid);
    vector<double> rdf_em(rdf_grid);
    vector<double> rdf_me(rdf_grid);
    vector<double> rdf_mm(rdf_grid);
    double correct=4*acos(-1.0)*rdf_cut/double(rdf_grid);
    // a counter which counts all centre atoms in all frames, i.e. centre_count = #framse * #centre_atoms
    int count_frames = 0;
    
    // loop over the different trajectory files
    if (args.count("trc") < 1) {
      throw gromos::Exception(argv[0], "no coordinate or trajectory file specified (@trc)");
    }
    Arguments::const_iterator trcfirs = args.lower_bound("trc");
    Arguments::const_iterator trclast = args.upper_bound("trc");
    for (args::Arguments::const_iterator trc = trcfirs;
            trc != trclast; ++trc) {

      // the input coordinates
      InG96 ic;

      // the boundary
      bound::Boundary *pbc;

      // read boundary type, either from @pbc or GENBOX block
      if (args.count("pbc") > 0) { // read from arguments
        pbc = args::BoundaryParser::boundary(sys, args);
      } else { // read from GENBOX block
        if (args.count("pbc") == 0) {
          cerr << "WARNING: @pbc given with no argument(s), reading boundary type "
                  "from GENBOX block of the trajectory/coordinate file\n";
        }
        ic.open(trc->second.c_str());
        ic.select("SOLUTE");
        ic >> sys;
        pbc = args::BoundaryParser::boundary(sys);
        ic.close();
      }
      
      // loop over the configurations of the trajectory file
      ic.open(trc->second.c_str());
      ic.select("SOLUTE");
      while (!ic.eof()) {
        ic >> sys;
        count_frames++;

        // check if the box length is as least as big as twice the cut-off radius
        double L = sys.box().K().abs2();
        if (sys.box().L().abs2() < L) {
          L = sys.box().L().abs2();
        }
        if (sys.box().M().abs2() < L) {
          L = sys.box().M().abs2();
        }
        if (distmax > (sqrt(L) / 2.0)) {
          throw gromos::Exception(argv[0], "maximal @dist value bigger than "
                  "1/2 of the box length");
        }

        // calculate all centres of the beads (com)
        for (int b = 0; b < (int) beads.size(); ++b) {
          beads[b].com(pbc, sys);
        }

        // volume and volume correction needed for the rdf calculation
        double vol_corr = 1;
        if(pbc->type()=='t') vol_corr=0.5;
        double vol = sys.box().K_L_M() * vol_corr;
        
        // double loop over the beads
        int b2;
        double lj;
        double r2;
        double r;
        IJ ij;
#ifdef OMP
#pragma omp parallel for private(b2, ij, lj, r2, r) schedule(dynamic)
#endif
        for (int b1 = 0; b1 < (int) beads.size() - 1; ++b1) {
          
          // the distribution for the RDF calculation
          gmath::Distribution dist_ee(0, rdf_cut, rdf_grid);
          gmath::Distribution dist_em(0, rdf_cut, rdf_grid);
          gmath::Distribution dist_me(0, rdf_cut, rdf_grid);
          gmath::Distribution dist_mm(0, rdf_cut, rdf_grid);
          
          for (b2 = b1 + 1; b2 < (int) beads.size(); ++b2) {
            r2 = (beads[b1].pos() - pbc->nearestImage(beads[b1].pos(),
                    beads[b2].pos(), sys.box())).abs2();
            
            // add the distance to the rdf distribution
            if (beads[b1].mol() != beads[b2].mol()) {
#ifdef OMP
#pragma omp critical
#endif
              {
                if (beads[b1].isTail() && beads[b2].isTail()) {
                  dist_ee.add(sqrt(r2));
                  //cerr << "ee" << endl;
                } else if (beads[b1].isTail() && !beads[b2].isTail()) {
                  dist_em.add(sqrt(r2));
                  //cerr << "em" << endl;
                } else if (!beads[b1].isTail() && beads[b2].isTail()) {
                  dist_me.add(sqrt(r2));
                  //cerr << "me" << endl;
                } else {
                  dist_mm.add(sqrt(r2));
                  //cerr << "mm" << endl;
                }
                //cerr << endl;
              }
            } else {
#ifdef OMP
#pragma omp critical
#endif
              {
                // tail-tail contributions
                if (beads[b1].isTail() && beads[b2].isTail()) {
                  beadbeadDist_ee[ij].add(sqrt(r2));
                }
              }
            }
            
            // if the two beads are within the range of the distribution range,
            // calculate the LJ potential energy
            if (distmin * distmin <= r2 && distmax * distmax > r2) {
              if (method == "molecular") {
                ij.setValues(beads[b1].size(), beads[b2].size());
                // calculate the bead-bead LJ potential energy
                lj = beads[b1].calcLJ(beads[b2], pbc, sys);
                r = sqrt(r2);
                // and add the result to the two beads (corresponding potentials)
#ifdef OMP
#pragma omp critical
#endif               
                {
                  if(beads[b1].isTail() && beads[b2].isTail()) {
                    beads[b1].addLJtot_ee(ij, r, lj);
                    beads[b2].addLJtot_ee(ij, r, lj);
                  } else if (beads[b1].isTail() || beads[b2].isTail()) {
                    beads[b1].addLJtot_em(ij, r, lj);
                    beads[b2].addLJtot_em(ij, r, lj);
                  } else {
                    beads[b1].addLJtot_mm(ij, r, lj);
                    beads[b2].addLJtot_mm(ij, r, lj);
                  }
                  if (beads[b1].mol() != beads[b2].mol()) { // intermolecular LJ potential
                    if (beads[b1].isTail() && beads[b2].isTail()) {
                      beads[b1].addLJtotinter_ee(ij, r, lj);
                      beads[b2].addLJtotinter_ee(ij, r, lj);
                    }
                    else if (beads[b1].isTail() || beads[b2].isTail()) {
                      beads[b1].addLJtotinter_em(ij, r, lj);
                      beads[b2].addLJtotinter_em(ij, r, lj);
                    } else {
                      beads[b1].addLJtotinter_mm(ij, r, lj);
                      beads[b2].addLJtotinter_mm(ij, r, lj);
                    }
                  } else { // intramolecular LJ potential energy
                    if (beads[b1].isTail() && beads[b2].isTail()) {
                      beads[b1].addLJtotintra_ee(ij, r, lj);
                      beads[b2].addLJtotintra_ee(ij, r, lj);
                      //beadbeadDist_ee[ij].add(r);
                    } else if (beads[b1].isTail() || beads[b2].isTail()) {
                      beads[b1].addLJtotintra_em(ij, r, lj);
                      beads[b2].addLJtotintra_em(ij, r, lj);
                    } else {
                      beads[b1].addLJtotintra_mm(ij, r, lj);
                      beads[b2].addLJtotintra_mm(ij, r, lj);
                    }
                    // 12-LJ pot (neighboring beads)
                    if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 1) {
                      if (beads[b1].isTail() && beads[b2].isTail()) {
                        beads[b1].addLJintra12_ee(ij, r, lj);
                        beads[b2].addLJintra12_ee(ij, r, lj);
                        //beadbeadDist_ee[ij].add(r);
                      } else if (beads[b1].isTail() || beads[b2].isTail()) {
                        beads[b1].addLJintra12_em(ij, r, lj);
                        beads[b2].addLJintra12_em(ij, r, lj);
                        beadbeadDist_em[ij].add(r);
                      } else {
                        beads[b1].addLJintra12_mm(ij, r, lj);
                        beads[b2].addLJintra12_mm(ij, r, lj);
                        beadbeadDist_mm[ij].add(r);
                      }
                      beadbeadDist[ij].add(r);
                    } else if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 2) {
                      if (beads[b1].isTail() && beads[b2].isTail()) {
                        beads[b1].addLJintra13_ee(ij, r, lj);
                        beads[b2].addLJintra13_ee(ij, r, lj);
                      } else if (beads[b1].isTail() || beads[b2].isTail()) {
                        beads[b1].addLJintra13_em(ij, r, lj);
                        beads[b2].addLJintra13_em(ij, r, lj);
                      } else {
                        beads[b1].addLJintra13_mm(ij, r, lj);
                        beads[b2].addLJintra13_mm(ij, r, lj);
                      }
                    } else if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 3) {
                      if (beads[b1].isTail() && beads[b2].isTail()) {
                        beads[b1].addLJintra14_ee(ij, r, lj);
                        beads[b2].addLJintra14_ee(ij, r, lj);
                      } else if (beads[b1].isTail() || beads[b2].isTail()) {
                        beads[b1].addLJintra14_em(ij, r, lj);
                        beads[b2].addLJintra14_em(ij, r, lj);
                      } else {
                        beads[b1].addLJintra14_mm(ij, r, lj);
                        beads[b2].addLJintra14_mm(ij, r, lj);
                      }
                    }
                  }
                }
              } else {
                stringstream msg;
                msg << "method " << method << " not implemented";
                throw gromos::Exception(argv[0], msg.str());
              }
            }
          }
            
          double dens_ee = double(endBeadNum - 1) / vol;
          double dens_em = double(middleBeadNum) / vol;
          double dens_me = double(endBeadNum) / vol;
          double dens_mm = double(middleBeadNum -1) / vol;
          for (int k = 0; k < rdf_grid; k++) {
            // the factor 2 comes from the for loop over the with beads, which starts from the actual centre bead + 1 only
            // ==> we just calculate ij, but never ji bead combinations ==> factror 2
            double r = dist_ee.value(k);
            const double rdf_val_ee = 2 * double(dist_ee[k]) / (dens_ee * correct * r * r);
            r = dist_em.value(k);
            const double rdf_val_em = 2 * double(dist_em[k]) / (dens_em * correct * r * r);
            r = dist_me.value(k);
            const double rdf_val_me = 2 * double(dist_me[k]) / (dens_me * correct * r * r);
            r = dist_mm.value(k);
            const double rdf_val_mm = 2 * double(dist_mm[k]) / (dens_mm * correct * r * r);
#ifdef OMP
#pragma omp critical
#endif
            {
              rdf_ee[k] += rdf_val_ee;
              rdf_em[k] += rdf_val_em;
              rdf_me[k] += rdf_val_me;
              rdf_mm[k] += rdf_val_mm;
              //cerr << "r = " << r << endl << "rdf_val = " << rdf_val << endl << "rdf[k] = " << rdf[k] << endl << "dens = " << dens << endl << "correct = " << correct << endl << endl;
            }
          }

        }
        
        // bond, angle and dihedral distribution (actually, bonds are done above
        // and are not done here any more...)
        if (method == "molecular") {
          // do the angles
          for (int b = 0; b < (int) beads.size() - 2; ++b) {
            if((beads[b].beadNum() < beads[b + 1].beadNum()) &&
                    (beads[b].beadNum() < beads[b + 2].beadNum())) {
              Vec p2 = beads[b + 1].pos();
              Vec p1 = pbc->nearestImage(p2, beads[b].pos(), sys.box());
              Vec p3 = pbc->nearestImage(p2, beads[b + 2].pos(), sys.box());
              Vec v1 = p1 - p2;
              Vec v2 = p3 - p2;
              double a = acos(v1.dot(v2) / (v1.abs() * v2.abs())) / physConst.get_pi() * 180;
              IJK ijk(beads[b].size(), beads[b+1].size(), beads[b+2].size());
              angleDist[ijk].add(a);
            }
            // do the dihedrals
            if(b < int(beads.size() - 3) && (beads[b].beadNum() < beads[b + 1].beadNum()) && 
                    (beads[b].beadNum() < beads[b + 2].beadNum()) && 
                    (beads[b].beadNum() < beads[b + 3].beadNum())) {
              Vec tmpA = beads[b].pos() - pbc->nearestImage(beads[b].pos(), beads[b+1].pos(), sys.box());
              Vec tmpB = beads[b+3].pos() - pbc->nearestImage(beads[b+3].pos(), beads[b+2].pos(), sys.box());
              Vec tmpC = beads[b+2].pos() - pbc->nearestImage(beads[b+2].pos(), beads[b+1].pos(), sys.box());
              Vec p1 = tmpA.cross(tmpC);
              Vec p2 = tmpB.cross(tmpC);
              double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));
              if (cosphi > 1.0) cosphi = 1.0;
              if (cosphi <-1.0) cosphi = -1.0;
              double d = acos(cosphi) * 180 / physConst.get_pi();
              IJKL ijkl(beads[b].size(), beads[b+1].size(), beads[b+2].size(), beads[b+3].size());
              dihedralDist[ijkl].add(d);
            }
          }
        } else {
          stringstream msg;
          msg << "method " << method << " not implemented";
          throw gromos::Exception(argv[0], msg.str());
        }
        
      } // end of loop over the configurations of the trajectory file
      ic.close();

    } // end of loop over the different trajectory files
    
    // correct the rdf distribution for the number of frames and the number of centre atoms
    for (int i = 0; i < rdf_grid; i++) {  
      rdf_ee[i] /= double(endBeadNum * count_frames);
      rdf_em[i] /= double(endBeadNum * count_frames);
      rdf_me[i] /= double(middleBeadNum * count_frames);
      rdf_mm[i] /= double(middleBeadNum * count_frames); 
    }
    
    // now write the rdf
    {
      ofstream fout("rdf_ee.dat");
      for (int i = 0; i < rdf_grid; ++i) {
        double dr = rdf_cut / (rdf_grid);
        double r = 0 + i * dr;
        fout << r + dr/2 << " " << rdf_ee[i] << endl;
      }
      fout.close();
    }
    {
      ofstream fout("rdf_em.dat");
      for (int i = 0; i < rdf_grid; ++i) {
        double dr = rdf_cut / (rdf_grid);
        double r = 0 + i * dr;
        fout << r + dr/2 << " " << rdf_em[i] << endl;
      }
      fout.close();
    }
    {
      ofstream fout("rdf_me.dat");
      for (int i = 0; i < rdf_grid; ++i) {
        double dr = rdf_cut / (rdf_grid);
        double r = 0 + i * dr;
        fout << r + dr/2 << " " << rdf_me[i] << endl;
      }
      fout.close();
    }
    {
      ofstream fout("rdf_mm.dat");
      for (int i = 0; i < rdf_grid; ++i) {
        double dr = rdf_cut / (rdf_grid);
        double r = 0 + i * dr;
        fout << r + dr/2 << " " << rdf_mm[i] << endl;
      }
      fout.close();
    }

    // unify the calculated potential of all beads
    map<IJ, LJpot> totLJ_ee, totLJ_em, totLJ_mm;
    map<IJ, LJpot> totinterLJ_ee, totinterLJ_em, totinterLJ_mm;
    map<IJ, LJpot> totintraLJ_ee, totintraLJ_em, totintraLJ_mm;
    map<IJ, LJpot> intra12LJ_ee, intra12LJ_em, intra12LJ_mm;
    map<IJ, LJpot> intra13LJ_ee, intra13LJ_em, intra13LJ_mm;
    map<IJ, LJpot> intra14LJ_ee, intra14LJ_em, intra14LJ_mm;
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      totLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totinterLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totinterLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totinterLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totintraLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totintraLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totintraLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra12LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra12LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra12LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra13LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra13LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra13LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra14LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra14LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra14LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
    }
    for (unsigned int b = 0; b < beads.size(); ++b) {
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        totLJ_ee[*it].unify(beads[b].get_totLJ_ee()[*it]);
        totLJ_em[*it].unify(beads[b].get_totLJ_em()[*it]);
        totLJ_mm[*it].unify(beads[b].get_totLJ_mm()[*it]);
        totinterLJ_ee[*it].unify(beads[b].get_totinterLJ_ee()[*it]);
        totinterLJ_em[*it].unify(beads[b].get_totinterLJ_em()[*it]);
        totinterLJ_mm[*it].unify(beads[b].get_totinterLJ_mm()[*it]);
        totintraLJ_ee[*it].unify(beads[b].get_totintraLJ_ee()[*it]);
        totintraLJ_em[*it].unify(beads[b].get_totintraLJ_em()[*it]);
        totintraLJ_mm[*it].unify(beads[b].get_totintraLJ_mm()[*it]);
        intra12LJ_ee[*it].unify(beads[b].get_intra12LJ_ee()[*it]);
        intra12LJ_em[*it].unify(beads[b].get_intra12LJ_em()[*it]);
        intra12LJ_mm[*it].unify(beads[b].get_intra12LJ_mm()[*it]);
        intra13LJ_ee[*it].unify(beads[b].get_intra13LJ_ee()[*it]);
        intra13LJ_em[*it].unify(beads[b].get_intra13LJ_em()[*it]);
        intra13LJ_mm[*it].unify(beads[b].get_intra13LJ_mm()[*it]);
        intra14LJ_ee[*it].unify(beads[b].get_intra14LJ_ee()[*it]);
        intra14LJ_em[*it].unify(beads[b].get_intra14LJ_em()[*it]);
        intra14LJ_mm[*it].unify(beads[b].get_intra14LJ_mm()[*it]);
      }
    }

    // print the different potentials
    {
      {
        ofstream fout(fname_LJpot_FG2CG_ee.c_str());
        printPot(fout, totLJ_ee, totinterLJ_ee, totintraLJ_ee, intra12LJ_ee, intra13LJ_ee, intra14LJ_ee);
        fout.close();
      }
      {
        ofstream fout(fname_LJpot_FG2CG_em.c_str());
        printPot(fout, totLJ_em, totinterLJ_em, totintraLJ_em, intra12LJ_em, intra13LJ_em, intra14LJ_em);
        fout.close();
      }
      {
        ofstream fout(fname_LJpot_FG2CG_mm.c_str());
        printPot(fout, totLJ_mm, totinterLJ_mm, totintraLJ_mm, intra12LJ_mm, intra13LJ_mm, intra14LJ_mm);
        fout.close();
      }
    }

    
    // calculate and print the resulting LJ pot for the CG system
    map<IJ, double> sigmas_tot_ee, sigmas_tot_em, sigmas_tot_mm;
    map<IJ, double> sigmas_totinter_ee, sigmas_totinter_em, sigmas_totinter_mm;
    map<IJ, double> sigmas_totintra_ee, sigmas_totintra_em, sigmas_totintra_mm;
    map<IJ, double> sigmas_intra12_ee, sigmas_intra12_em, sigmas_intra12_mm;
    map<IJ, double> sigmas_intra13_ee, sigmas_intra13_em, sigmas_intra13_mm;
    map<IJ, double> sigmas_intra14_ee, sigmas_intra14_em, sigmas_intra14_mm;
    map<IJ, double> epsilons_tot_ee, epsilons_tot_em, epsilons_tot_mm;
    map<IJ, double> epsilons_totinter_ee, epsilons_totinter_em, epsilons_totinter_mm;
    map<IJ, double> epsilons_totintra_ee, epsilons_totintra_em, epsilons_totintra_mm;
    map<IJ, double> epsilons_intra12_ee, epsilons_intra12_em, epsilons_intra12_mm;
    map<IJ, double> epsilons_intra13_ee, epsilons_intra13_em, epsilons_intra13_mm;
    map<IJ, double> epsilons_intra14_ee, epsilons_intra14_em, epsilons_intra14_mm;
    map<IJ, double> C12_tot_ee, C12_tot_em, C12_tot_mm;
    map<IJ, double> C12_totinter_ee, C12_totinter_em, C12_totinter_mm;
    map<IJ, double> C12_totintra_ee, C12_totintra_em, C12_totintra_mm;
    map<IJ, double> C12_intra12_ee, C12_intra12_em, C12_intra12_mm;
    map<IJ, double> C12_intra13_ee, C12_intra13_em, C12_intra13_mm;
    map<IJ, double> C12_intra14_ee, C12_intra14_em, C12_intra14_mm;
    map<IJ, double> C6_tot_ee, C6_tot_em, C6_tot_mm;
    map<IJ, double> C6_totinter_ee, C6_totinter_em, C6_totinter_mm;
    map<IJ, double> C6_totintra_ee, C6_totintra_em, C6_totintra_mm;
    map<IJ, double> C6_intra12_ee, C6_intra12_em, C6_intra12_mm;
    map<IJ, double> C6_intra13_ee, C6_intra13_em, C6_intra13_mm;
    map<IJ, double> C6_intra14_ee, C6_intra14_em, C6_intra14_mm;
    calcEpsSigma(epsilons_tot_ee, sigmas_tot_ee, totLJ_ee);
    calcEpsSigma(epsilons_tot_em, sigmas_tot_em, totLJ_em);
    calcEpsSigma(epsilons_tot_mm, sigmas_tot_mm, totLJ_mm);
    calcEpsSigma(epsilons_totinter_ee, sigmas_totinter_ee, totinterLJ_ee);
    calcEpsSigma(epsilons_totinter_em, sigmas_totinter_em, totinterLJ_em);
    calcEpsSigma(epsilons_totinter_mm, sigmas_totinter_mm, totinterLJ_mm);
    calcEpsSigma(epsilons_totintra_ee, sigmas_totintra_ee, totintraLJ_ee);
    calcEpsSigma(epsilons_totintra_em, sigmas_totintra_em, totintraLJ_em);
    calcEpsSigma(epsilons_totintra_mm, sigmas_totintra_mm, totintraLJ_mm);
    calcEpsSigma(epsilons_intra12_ee, sigmas_intra12_ee, intra12LJ_ee);
    calcEpsSigma(epsilons_intra12_em, sigmas_intra12_em, intra12LJ_em);
    calcEpsSigma(epsilons_intra12_mm, sigmas_intra12_mm, intra12LJ_mm);
    calcEpsSigma(epsilons_intra13_ee, sigmas_intra13_ee, intra13LJ_ee);
    calcEpsSigma(epsilons_intra13_em, sigmas_intra13_em, intra13LJ_em);
    calcEpsSigma(epsilons_intra13_mm, sigmas_intra13_mm, intra13LJ_mm);
    calcEpsSigma(epsilons_intra14_ee, sigmas_intra14_ee, intra14LJ_ee);
    calcEpsSigma(epsilons_intra14_em, sigmas_intra14_em, intra14LJ_em);
    calcEpsSigma(epsilons_intra14_mm, sigmas_intra14_mm, intra14LJ_mm);
    calcC12C6(C12_tot_ee, C6_tot_ee, epsilons_tot_ee, sigmas_tot_ee);
    calcC12C6(C12_tot_em, C6_tot_em, epsilons_tot_em, sigmas_tot_em);
    calcC12C6(C12_tot_mm, C6_tot_mm, epsilons_tot_mm, sigmas_tot_mm);
    calcC12C6(C12_totinter_ee, C6_totinter_ee, epsilons_totinter_ee, sigmas_totinter_ee);
    calcC12C6(C12_totinter_em, C6_totinter_em, epsilons_totinter_em, sigmas_totinter_em);
    calcC12C6(C12_totinter_mm, C6_totinter_mm, epsilons_totinter_mm, sigmas_totinter_mm);
    calcC12C6(C12_totintra_ee, C6_totintra_ee, epsilons_totintra_ee, sigmas_totintra_ee);
    calcC12C6(C12_totintra_em, C6_totintra_em, epsilons_totintra_em, sigmas_totintra_em);
    calcC12C6(C12_totintra_mm, C6_totintra_mm, epsilons_totintra_mm, sigmas_totintra_mm);
    calcC12C6(C12_intra12_ee, C6_intra12_ee, epsilons_intra12_ee, sigmas_intra12_ee);
    calcC12C6(C12_intra12_em, C6_intra12_em, epsilons_intra12_em, sigmas_intra12_em);
    calcC12C6(C12_intra12_mm, C6_intra12_mm, epsilons_intra12_mm, sigmas_intra12_mm);
    calcC12C6(C12_intra13_ee, C6_intra13_ee, epsilons_intra13_ee, sigmas_intra13_ee);
    calcC12C6(C12_intra13_em, C6_intra13_em, epsilons_intra13_em, sigmas_intra13_em);
    calcC12C6(C12_intra13_mm, C6_intra13_mm, epsilons_intra13_mm, sigmas_intra13_mm);
    calcC12C6(C12_intra14_ee, C6_intra14_ee, epsilons_intra14_ee, sigmas_intra14_ee);
    calcC12C6(C12_intra14_em, C6_intra14_em, epsilons_intra14_em, sigmas_intra14_em);
    calcC12C6(C12_intra14_mm, C6_intra14_mm, epsilons_intra14_mm, sigmas_intra14_mm);
    printLennardJonesParamters("END-END", epsilons_tot_ee, epsilons_totinter_ee, epsilons_totintra_ee, epsilons_intra12_ee, epsilons_intra13_ee, epsilons_intra14_ee,
            sigmas_tot_ee, sigmas_totinter_ee, sigmas_totintra_ee, sigmas_intra12_ee, sigmas_intra13_ee, sigmas_intra14_ee,
            C12_tot_ee, C12_totinter_ee, C12_totintra_ee, C12_intra12_ee, C12_intra13_ee, C12_intra14_ee, C6_tot_ee,
            C6_totinter_ee, C6_totintra_ee, C6_intra12_ee, C6_intra13_ee, C6_intra14_ee);
    printLennardJonesParamters("END-MIDDLE", epsilons_tot_em, epsilons_totinter_em, epsilons_totintra_em, epsilons_intra12_em, epsilons_intra13_em, epsilons_intra14_em,
            sigmas_tot_em, sigmas_totinter_em, sigmas_totintra_em, sigmas_intra12_em, sigmas_intra13_em, sigmas_intra14_em,
            C12_tot_em, C12_totinter_em, C12_totintra_em, C12_intra12_em, C12_intra13_em, C12_intra14_em, C6_tot_em,
            C6_totinter_em, C6_totintra_em, C6_intra12_em, C6_intra13_em, C6_intra14_em);
    printLennardJonesParamters("MIDDLE-MIDDLE", epsilons_tot_mm, epsilons_totinter_mm, epsilons_totintra_mm, epsilons_intra12_mm, epsilons_intra13_mm, epsilons_intra14_mm,
            sigmas_tot_mm, sigmas_totinter_mm, sigmas_totintra_mm, sigmas_intra12_mm, sigmas_intra13_mm, sigmas_intra14_mm,
            C12_tot_mm, C12_totinter_mm, C12_totintra_mm, C12_intra12_mm, C12_intra13_mm, C12_intra14_mm, C6_tot_mm,
            C6_totinter_mm, C6_totintra_mm, C6_intra12_mm, C6_intra13_mm, C6_intra14_mm);

    // calculate the LJ potentials using the C12 and C6 values above...
    map<IJ, LJpot> ftotLJ_ee, ftotLJ_em, ftotLJ_mm;
    map<IJ, LJpot> ftotinterLJ_ee, ftotinterLJ_em, ftotinterLJ_mm;
    map<IJ, LJpot> ftotintraLJ_ee, ftotintraLJ_em, ftotintraLJ_mm;
    map<IJ, LJpot> fintra12LJ_ee, fintra12LJ_em, fintra12LJ_mm;
    map<IJ, LJpot> fintra13LJ_ee, fintra13LJ_em, fintra13LJ_mm;
    map<IJ, LJpot> fintra14LJ_ee, fintra14LJ_em, fintra14LJ_mm;
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      double c6 = 4 * epsilons_tot_ee[*it] * pow(sigmas_tot_ee[*it], 6);
      double c12 = 4 * epsilons_tot_ee[*it] * pow(sigmas_tot_ee[*it], 12);
      ftotLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_tot_em[*it] * pow(sigmas_tot_em[*it], 6);
      c12 = 4 * epsilons_tot_em[*it] * pow(sigmas_tot_em[*it], 12);
      ftotLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_tot_mm[*it] * pow(sigmas_tot_mm[*it], 6);
      c12 = 4 * epsilons_tot_mm[*it] * pow(sigmas_tot_mm[*it], 12);
      ftotLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totinter_ee[*it] * pow(sigmas_totinter_ee[*it], 6);
      c12 = 4 * epsilons_totinter_ee[*it] * pow(sigmas_totinter_ee[*it], 12);
      ftotinterLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totinter_em[*it] * pow(sigmas_totinter_em[*it], 6);
      c12 = 4 * epsilons_totinter_em[*it] * pow(sigmas_totinter_em[*it], 12);
      ftotinterLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totinter_mm[*it] * pow(sigmas_totinter_mm[*it], 6);
      c12 = 4 * epsilons_totinter_mm[*it] * pow(sigmas_totinter_mm[*it], 12);
      ftotinterLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totintra_ee[*it] * pow(sigmas_totintra_ee[*it], 6);
      c12 = 4 * epsilons_totintra_ee[*it] * pow(sigmas_totintra_ee[*it], 12);
      ftotintraLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totintra_em[*it] * pow(sigmas_totintra_em[*it], 6);
      c12 = 4 * epsilons_totintra_em[*it] * pow(sigmas_totintra_em[*it], 12);
      ftotintraLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totintra_mm[*it] * pow(sigmas_totintra_mm[*it], 6);
      c12 = 4 * epsilons_totintra_mm[*it] * pow(sigmas_totintra_mm[*it], 12);
      ftotintraLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra12_ee[*it] * pow(sigmas_intra12_ee[*it], 6);
      c12 = 4 * epsilons_intra12_ee[*it] * pow(sigmas_intra12_ee[*it], 12);
      fintra12LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra12_em[*it] * pow(sigmas_intra12_em[*it], 6);
      c12 = 4 * epsilons_intra12_em[*it] * pow(sigmas_intra12_em[*it], 12);
      fintra12LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra12_mm[*it] * pow(sigmas_intra12_mm[*it], 6);
      c12 = 4 * epsilons_intra12_mm[*it] * pow(sigmas_intra12_mm[*it], 12);
      fintra12LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra13_ee[*it] * pow(sigmas_intra13_ee[*it], 6);
      c12 = 4 * epsilons_intra13_ee[*it] * pow(sigmas_intra13_ee[*it], 12);
      fintra13LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra13_em[*it] * pow(sigmas_intra13_em[*it], 6);
      c12 = 4 * epsilons_intra13_em[*it] * pow(sigmas_intra13_em[*it], 12);
      fintra13LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra13_mm[*it] * pow(sigmas_intra13_mm[*it], 6);
      c12 = 4 * epsilons_intra13_mm[*it] * pow(sigmas_intra13_mm[*it], 12);
      fintra13LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra14_ee[*it] * pow(sigmas_intra14_ee[*it], 6);
      c12 = 4 * epsilons_intra14_ee[*it] * pow(sigmas_intra14_ee[*it], 12);
      fintra14LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra14_em[*it] * pow(sigmas_intra14_em[*it], 6);
      c12 = 4 * epsilons_intra14_em[*it] * pow(sigmas_intra14_em[*it], 12);
      fintra14LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra14_mm[*it] * pow(sigmas_intra14_mm[*it], 6);
      c12 = 4 * epsilons_intra14_mm[*it] * pow(sigmas_intra14_mm[*it], 12);
      fintra14LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
    }
    {
      ofstream fout(fname_LJpot_CG_ee.c_str());
      printPot(fout, ftotLJ_ee, ftotinterLJ_ee, ftotintraLJ_ee, fintra12LJ_ee, fintra13LJ_ee, fintra14LJ_ee);
      fout.close();
    }
    {
      ofstream fout(fname_LJpot_CG_em.c_str());
      printPot(fout, ftotLJ_em, ftotinterLJ_em, ftotintraLJ_em, fintra12LJ_em, fintra13LJ_em, fintra14LJ_em);
      fout.close();
    }
    {
      ofstream fout(fname_LJpot_CG_mm.c_str());
      printPot(fout, ftotLJ_mm, ftotinterLJ_mm, ftotintraLJ_mm, fintra12LJ_mm, fintra13LJ_mm, fintra14LJ_mm);
      fout.close();
    }
    
    // normalize and print the distribution
    printBeadBeadDist(fname_beadbead_dist, beadbeadDist, IJs, bondlength_min, bondlength_max, 1000);
    printBeadBeadDist(fname_beadbead_dist_ee, beadbeadDist_ee, IJs, 0.0, 2.0, 1000);
    printBeadBeadDist(fname_beadbead_dist_em, beadbeadDist_em, IJs, bondlength_min, bondlength_max, 1000);
    printBeadBeadDist(fname_beadbead_dist_mm, beadbeadDist_mm, IJs, bondlength_min, bondlength_max, 1000);

    // print the angle distribution
    {
      ofstream fout("angle.dist");
      double dgrid = 1.0;
      fout.precision(9);
      fout << "#" << setw(19) << "angle / degree";
      for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
        stringstream ss;
        ss << it->i() << "-" << it->j() << "-" << it->k();
        fout << scientific << setw(20) << ss.str();
      }
      fout << endl;
      for (int i = 0; i < 180; ++i) {
        double a = (i + 0.5) * dgrid;
        fout << scientific << setw(20) << a;
        for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
          if (angleDist[*it].nVal() > 0) {
            fout << scientific << setw(20) << (double) angleDist[*it][i] / (angleDist[*it].nVal() * dgrid);
          } else {
            fout << scientific << setw(20) << (double) angleDist[*it][i];
          }
        }
        fout << endl;
	  }
	  fout << "#" << setw(19) << "average:";
      for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
        if (angleDist[*it].nVal() > 9) {
		  fout << scientific << setw(20) << angleDist[*it].ave();
		} else {
		  fout << setw(20) << "-";
		}
      }
      fout << endl << "#" << setw(19) << "maximum value at:";
      for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
		if (angleDist[*it].nVal() > 0) {
		  fout << scientific << setw(20) << angleDist[*it].maxValAt();
		} else {
		  fout << setw(20) << "-";
		}
      }
	  fout << endl;
      fout.close();
    }
    
    // print the dihedral distribution
    {
      ofstream fout("dihedral.dist");
      double dgrid = 1.0;
      fout.precision(9);
      fout << "#" << setw(19) << "dihedal / degree";
      for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
        stringstream ss;
        ss << it->i() << "-" << it->j() << "-" << it->k() << "-" << it->l();
        fout << scientific << setw(20) << ss.str();
      }
      fout << endl;
      for (int i = 0; i < 180; ++i) {
        double a = (i + 0.5) * dgrid;
        fout << scientific << setw(20) << a;
        for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
          if (dihedralDist[*it].nVal() > 0) {
            fout << scientific << setw(20) <<  (double) dihedralDist[*it][i] / (dihedralDist[*it].nVal() * dgrid);
          } else {
            fout << scientific << setw(20) <<  (double) dihedralDist[*it][i];
          }
        }
		fout << endl;
	  }
      fout << endl << "#" << setw(19) << "average:";
      for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
        if (dihedralDist[*it].nVal() > 0) {
		  fout << scientific << setw(20) << dihedralDist[*it].ave();
		} else {
		  fout << setw(20) << "-";
		}
      }   
      fout << endl << "#" << setw(19) << "maximum value at:";
      for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
        if (dihedralDist[*it].nVal() > 0) {
		  fout << scientific << setw(20) << dihedralDist[*it].maxValAt();
		} else {
		  fout << setw(20) << "-";
		}
      }   
	  fout << endl;
      fout.close();

    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}

namespace cgLJpot {
  
  IJ::IJ(int i, int j) {
    if (j < i) {
      int t = i;
      i = j;
      j = t;
    }
    I = i;
    J = j;
  }

  IJ::IJ(const IJ &ij) {
    I = ij.I;
    J = ij.J;
  }

  int IJ::i() const {
    return I;
  }

  int IJ::j() const {
    return J;
  }

  void IJ::setValues(int i, int j) {
    if (j < i) {
      int t = i;
      i = j;
      j = t;
    }
    I = i;
    J = j;
  }
  
  IJK::IJK(int i, int j, int k) {
    if (k < i) {
      int t = i;
      i = k;
      k = t;
    }
    I = i;
    J = j;
    K = k;
  }

  IJK::IJK(const IJK &ijk) {
    I = ijk.I;
    J = ijk.J;
    K = ijk.K;
  }

  int IJK::i() const {
    return I;
  }

  int IJK::j() const {
    return J;
  }
  
  int IJK::k() const {
    return K;
  }

  void IJK::setValues(int i, int j, int k) {
    if (k < i) {
      int t = i;
      i = k;
      k = t;
    }
    I = i;
    J = j;
    K = k;
  }
  
  IJKL::IJKL(int i, int j, int k, int l) {
    if (k < j) {
      int t = j;
      j = k;
      k = t;
      t = i;
      i = l;
      l = t;
    }
    I = i;
    J = j;
    K = k;
    L = l;
  }

  IJKL::IJKL(const IJKL &ijkl) {
    I = ijkl.I;
    J = ijkl.J;
    K = ijkl.K;
    L = ijkl.L;
  }

  int IJKL::i() const {
    return I;
  }

  int IJKL::j() const {
    return J;
  }
  
  int IJKL::k() const {
    return K;
  }
  
  int IJKL::l() const {
    return L;
  }

  void IJKL::setValues(int i, int j, int k, int l) {
    if (k < j) {
      int t = j;
      j = k;
      k = t;
      t = i;
      i = l;
      l = t;
    }
    I = i;
    J = j;
    K = k;
    L = l;
  }

  LJpot::LJpot(double min_, double max_, int grid) {
    lj.resize(grid);
    count.resize(grid);
    min = min_;
    max = max_;
    dgrid = (max - min) / grid;
    C12 = -1.0;
    C6 = -1.0;
    for (int i = 0; i < grid; ++i) {
      lj[i] = 0.0;
      count[i] = 0;
    }
    used = false;
  }

  LJpot::LJpot(double min_, double max_, int grid, double c12, double c6) {
    lj.resize(grid);
    count.resize(grid);
    min = min_;
    max = max_;
    dgrid = (max - min) / grid;
    C12 = c12;
    C6 = c6;
    for (int i = 0; i < grid; ++i) {
      double R = r(i);
      double R3 = R * R * R;
      double R6 = R3 * R3;
      lj[i] = (C12 / R6 - C6) / R6;
      count[i] = 1;
    }
    used = false;
  }

  LJpot::LJpot(const LJpot &ljp) {
    lj = ljp.lj;
    count = ljp.count;
    min = ljp.min;
    max = ljp.max;
    dgrid = ljp.dgrid;
    C12 = ljp.C12;
    C6 = ljp.C6;
    used = ljp.used;
  }
  
  void LJpot::init(double min_, double max_, int grid) {
    lj.resize(grid);
    count.resize(grid);
    min = min_;
    max = max_;
    dgrid = (max - min) / grid;
    C12 = -1.0;
    C6 = -1.0;
    for (int i = 0; i < grid; ++i) {
      lj[i] = 0.0;
      count[i] = 0;
    }
    used = false;
  }

  void LJpot::add(double pos, double val) {
    if (min <= pos && pos < max) {
      int r = floor((pos - min) / dgrid);
      count[r] += 1;
      lj[r] += val;
    }
    if(!used) {
      used = true;
    }
  }

  void LJpot::unify(const LJpot & ljp) {
    assert(lj.size() == ljp.lj.size());
    for (unsigned int i = 0; i < lj.size(); ++i) {
      lj[i] += ljp.lj[i];
      count[i] += ljp.count[i];
    }
    if(ljp.wasUsed()) {
      used = true;
    }
  }

  double LJpot::r(unsigned int i) {
    assert(i < lj.size());
    return (i + 0.5) * dgrid + min;
  }

  double LJpot::r(unsigned int i) const {
    assert(i < lj.size());
    return (i + 0.5) * dgrid + min;
  }

  double LJpot::pot(unsigned int i) {
    assert(i < lj.size());
    return count[i] == 0 ? 0.0 : lj[i] / count[i];
  }

  double LJpot::get_min() {
    return min;
  }

  double LJpot::get_max() {
    return max;
  }

  int LJpot::get_grid() {
    return lj.size();
  }

  double LJpot::rmin() {
    // find the grid with maximal pot first
    double i_max = 0;
    double pot_max = pot(0);
    for (int i = 1; i < get_grid(); ++i) {
      if (pot(i) > pot_max) {
        i_max = i;
        pot_max = pot(i);
      }
    }
    // now find the minimum
    double r_min = r(i_max);
    double pot_min = pot(i_max);
    for (int i = i_max + 1; i < get_grid(); ++i) {
      if (pot(i) < pot_min) {
        r_min = r(i);
        pot_min = pot(i);
      }
    }
    return r_min;
  }

  double LJpot::potmin() {
    // find the grid with maximal pot first
    double i_max = 0;
    double pot_max = pot(0);
    for (int i = 1; i < get_grid(); ++i) {
      if (pot(i) > pot_max) {
        i_max = i;
        pot_max = pot(i);
      }
    }
    // now find the minimum
    double pot_min = pot(i_max);
    for (int i = i_max + 1; i < get_grid(); ++i) {
      if (pot(i) < pot_min) {
        pot_min = pot(i);
      }
    }
    return pot_min;
  }
  
  bool LJpot::wasUsed() const {
    return used;
  }

  bead::bead(gcore::System& sys, GromosForceField &groff, int mom, int bnum, set<IJ> &ij, double min, double max, double grid) {
    gff = &groff;
    atoms.setSystem(sys);
    memberOfMol = mom;
    beadnum = bnum;
    set<IJ>::const_iterator it = ij.begin();
    for (; it != ij.end(); ++it) {
      totLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totinterLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totinterLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totinterLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totintraLJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totintraLJ_em.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totintraLJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra12LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra12LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra12LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra13LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra13LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra13LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra14LJ_ee.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra14LJ_em.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra14LJ_mm.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
    }
    istail = false;
  }

  bead::bead(bead const &b) {
    gff = b.gff;
    atoms = b.atoms;
    centre = b.centre;
    memberOfMol = b.memberOfMol;
    beadnum = b.beadnum;
    totLJ_ee = b.totLJ_ee;
    totLJ_em = b.totLJ_em;
    totLJ_mm = b.totLJ_mm;
    totinterLJ_ee = b.totinterLJ_ee;
    totinterLJ_em = b.totinterLJ_em;
    totinterLJ_mm = b.totinterLJ_mm;
    totintraLJ_ee = b.totintraLJ_ee;
    totintraLJ_em = b.totintraLJ_em;
    totintraLJ_mm = b.totintraLJ_mm;
    intra12LJ_ee = b.intra12LJ_ee;
    intra12LJ_em = b.intra12LJ_em;
    intra12LJ_mm = b.intra12LJ_mm;
    intra13LJ_ee = b.intra13LJ_ee;
    intra13LJ_em = b.intra13LJ_em;
    intra13LJ_mm = b.intra13LJ_mm;
    intra14LJ_ee = b.intra14LJ_ee;
    intra14LJ_em = b.intra14LJ_em;
    intra14LJ_mm = b.intra14LJ_mm;
    istail = b.istail;
  }

  int bead::addAtom(int m, int a) {
    return atoms.addAtom(m, a);
  }

  int bead::size() {
    return atoms.size();
  }

  Vec bead::cog(Boundary *pbc, System &sys) {
    assert(atoms.size() > 0);
    Vec cog(atoms.pos(0));
    for (int i = 1; i < atoms.size(); ++i) {
      cog += pbc->nearestImage(atoms.pos(0), atoms.pos(i), sys.box());
    }
    cog /= atoms.size();
    centre = cog;
    return cog;
  }

  Vec bead::com(Boundary *pbc, System &sys) {
    assert(atoms.size() > 0);
    double sumMass = atoms.mass(0);
    Vec com(atoms.pos(0) * sumMass);
    for (int i = 1; i < atoms.size(); ++i) {
      double mass = atoms.mass(i);
      sumMass += mass;
      com += pbc->nearestImage(atoms.pos(0), atoms.pos(i), sys.box()) * mass;
    }
    com = com / (sumMass);
    centre = com;
    return com;
  }

  Vec bead::pos() {
    return centre;
  }

  map<IJ, LJpot> bead::get_totLJ_ee() {
    return totLJ_ee;
  }
  
  map<IJ, LJpot> bead::get_totLJ_em() {
    return totLJ_em;
  }
  
  map<IJ, LJpot> bead::get_totLJ_mm() {
    return totLJ_mm;
  }

  map<IJ, LJpot> bead::get_totinterLJ_ee() {
    return totinterLJ_ee;
  }
  
  map<IJ, LJpot> bead::get_totinterLJ_em() {
    return totinterLJ_em;
  }
  map<IJ, LJpot> bead::get_totinterLJ_mm() {
    return totinterLJ_mm;
  }

  map<IJ, LJpot> bead::get_totintraLJ_ee() {
    return totintraLJ_ee;
  }
  
  map<IJ, LJpot> bead::get_totintraLJ_em() {
    return totintraLJ_em;
  }
  
  map<IJ, LJpot> bead::get_totintraLJ_mm() {
    return totintraLJ_mm;
  }

  map<IJ, LJpot> bead::get_intra12LJ_ee() {
    return intra12LJ_ee;
  }
  
  map<IJ, LJpot> bead::get_intra12LJ_em() {
    return intra12LJ_em;
  }
  
  map<IJ, LJpot> bead::get_intra12LJ_mm() {
    return intra12LJ_mm;
  }

  map<IJ, LJpot> bead::get_intra13LJ_ee() {
    return intra13LJ_ee;
  }
  
  map<IJ, LJpot> bead::get_intra13LJ_em() {
    return intra13LJ_em;
  }
  
  map<IJ, LJpot> bead::get_intra13LJ_mm() {
    return intra13LJ_mm;
  }

  map<IJ, LJpot> bead::get_intra14LJ_ee() {
    return intra14LJ_ee;
  }

  map<IJ, LJpot> bead::get_intra14LJ_em() {
    return intra14LJ_em;
  }
  
  map<IJ, LJpot> bead::get_intra14LJ_mm() {
    return intra14LJ_mm;
  }
  
  int bead::mol() {
    return memberOfMol;
  }

  int bead::beadNum() {
    return beadnum;
  }
  
  void bead::setAsTail() {
    istail = true;
  }

  bool bead::isTail() {
    return istail;
  }
  
  double bead::calcLJ(bead &b, Boundary *pbc, System &sys) {
    double LJsum = 0.0;
    for (int i1 = 0; i1 < size(); ++i1) {
      for (int i2 = 0; i2 < b.size(); ++i2) {
        // the stuff for atom 1 is also defined in here since
        // the atoms at1 and a2 might be interchanged sometimes which has to be
        // undone before continuing with the nest at2
        int m1 = atoms.mol(i1);
        int a1 = atoms.atom(i1);
        int gnum1 = atoms.gromosAtom(i1);
        AtomTopology at1 = sys.mol(m1).topology().atom(a1);
        int m2 = b.atoms.mol(i2);
        int a2 = b.atoms.atom(i2);
        int gnum2 = b.atoms.gromosAtom(i2);
        AtomTopology at2 = sys.mol(m2).topology().atom(a2);
        // to memorize if the inter-bead atoms are excluded from each other
        bool excluded = false;
        bool excluded14 = false;
        // check if and which (normal, special 1,4; only happens if the atoms are
        // located in the same molecule) LJ interaction is calculated
        if (m1 == m2) {
          // make sure at1 is the atom first listed in the topology
          if (gnum2 < gnum1) {
            int atmp = a1;
            a1 = a2;
            a2 = atmp;
            int gnumtmp = gnum1;
            gnum1 = gnum2;
            gnum2 = gnumtmp;
            AtomTopology attmp(at1);
            at1 = at2;
            at2 = attmp;
          }
          // check if and which (normal, special 1,4) LJ interaction is calculated
          for (int e = 0; e < at1.exclusion().size(); ++e) {
            if (at1.exclusion().atom(e) == a2) {
              excluded = true;
              break;
            }
          }
          for (int e = 0; e < at1.exclusion14().size(); ++e) {
            if (at1.exclusion14().atom(e) == a2) {
              excluded14 = true;
              break;
            }
          }
        }
        // do the appropriate LJ potential energy calculation
        double c12;
        double c6;
        int iac1 = at1.iac();
        int iac2 = at2.iac();
        if (excluded) {
          continue;
        } else if (excluded14) {
          c12 = gff->ljType(AtomPair(iac1, iac2)).cs12();
          c6 = gff->ljType(AtomPair(iac1, iac2)).cs6();
        } else {
          c12 = gff->ljType(AtomPair(iac1, iac2)).c12();
          c6 = gff->ljType(AtomPair(iac1, iac2)).c6();
        }
        double r2a = (atoms.pos(i1) -
                (pbc->nearestImage(atoms.pos(i1), b.atoms.pos(i2), sys.box()))).abs2();
        double r6a = r2a * r2a * r2a;
        LJsum += (c12 / r6a - c6) / r6a;
      }
    }
    return LJsum;
  }

  void bead::addLJtot_ee(const IJ &ij, double r, const double &lj) {
    totLJ_ee[ij].add(r, lj);
  }
  
  void bead::addLJtot_em(const IJ &ij, double r, const double &lj) {
    totLJ_em[ij].add(r, lj);
  }
  
  void bead::addLJtot_mm(const IJ &ij, double r, const double &lj) {
    totLJ_mm[ij].add(r, lj);
  }

  void bead::addLJtotinter_ee(const IJ &ij, double r, const double &lj) {
    totinterLJ_ee[ij].add(r, lj);
  }
  
  void bead::addLJtotinter_em(const IJ &ij, double r, const double &lj) {
    totinterLJ_em[ij].add(r, lj);
  }
  
  void bead::addLJtotinter_mm(const IJ &ij, double r, const double &lj) {
    totinterLJ_mm[ij].add(r, lj);
  }

  void bead::addLJtotintra_ee(const IJ &ij, double r, const double &lj) {
    totintraLJ_ee[ij].add(r, lj);
  }
  
  void bead::addLJtotintra_em(const IJ &ij, double r, const double &lj) {
    totintraLJ_em[ij].add(r, lj);
  }
  
  void bead::addLJtotintra_mm(const IJ &ij, double r, const double &lj) {
    totintraLJ_mm[ij].add(r, lj);
  }

  void bead::addLJintra12_ee(const IJ &ij, double r, const double &lj) {
    intra12LJ_ee[ij].add(r, lj);
  }
  
  void bead::addLJintra12_em(const IJ &ij, double r, const double &lj) {
    intra12LJ_em[ij].add(r, lj);
  }
  
  void bead::addLJintra12_mm(const IJ &ij, double r, const double &lj) {
    intra12LJ_mm[ij].add(r, lj);
  }

  void bead::addLJintra13_ee(const IJ &ij, double r, const double &lj) {
    intra13LJ_ee[ij].add(r, lj);
  }
  
  void bead::addLJintra13_em(const IJ &ij, double r, const double &lj) {
    intra13LJ_em[ij].add(r, lj);
  }
  
  void bead::addLJintra13_mm(const IJ &ij, double r, const double &lj) {
    intra13LJ_mm[ij].add(r, lj);
  }

  void bead::addLJintra14_ee(const IJ &ij, double r, const double &lj) {
    intra14LJ_ee[ij].add(r, lj);
  }
  
  void bead::addLJintra14_em(const IJ &ij, double r, const double &lj) {
    intra14LJ_em[ij].add(r, lj);
  }
  
  void bead::addLJintra14_mm(const IJ &ij, double r, const double &lj) {
    intra14LJ_mm[ij].add(r, lj);
  }

  bool operator<(const IJ &ij1, const IJ &ij2) {
    if (ij1.i() < ij2.i() ||
		    (ij1.i() == ij2.i() && ij1.j() < ij2.j())) {
      return true;
    }
    return false;
  }
  
  bool operator<(const IJK &ijk1, const IJK &ijk2) {
    if (ijk1.i() < ijk2.i() ||
		    (ijk1.i() == ijk2.i() && ijk1.j() < ijk2.j()) ||
			(ijk1.i() == ijk2.i() && ijk1.j() == ijk2.j() && ijk1.k() < ijk2.k())) {
      return true;
    }
    return false;
  }
  
  bool operator<(const IJKL &ijkl1, const IJKL &ijkl2) {
    if (ijkl1.i() < ijkl2.i() || 
            (ijkl1.i() == ijkl2.i() && ijkl1.j() < ijkl2.j()) || 
            (ijkl1.i() == ijkl2.i() && ijkl1.j() == ijkl2.j() && ijkl1.k() < ijkl2.k()) ||
            (ijkl1.i() == ijkl2.i() && ijkl1.j() == ijkl2.j() && ijkl1.k() == ijkl2.k() && ijkl1.l() < ijkl2.l())) {
      return true;
    }
    return false;
  }

  void printPot(ostream &os, map<IJ, LJpot> &totLJpot, map<IJ, LJpot> &totLJinter,
          map<IJ, LJpot> &totLJintra, map<IJ, LJpot> &intraLJ12, map<IJ, LJpot> &intraLJ13,
          map<IJ, LJpot> &intraLJ14) {
    vector<string> header;
    vector<string> potentialnames;
    potentialnames.push_back("totLJpot");
    potentialnames.push_back("totLJinter");
    potentialnames.push_back("totLJintra");
    potentialnames.push_back("totLJ12");
    potentialnames.push_back("totLJ13");
    potentialnames.push_back("totLJ14");
    header.push_back("radius / nm");
    for (unsigned int i = 0; i < potentialnames.size(); ++i) {
      for (map<IJ, LJpot>::const_iterator it = totLJpot.begin(); it != totLJpot.end(); ++it) {
        stringstream ss;
        ss << potentialnames[i] << "_" << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    // print the header
    os << "#";
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        os << setw(19) << header[i];
      } else {
        os << setw(20) << header[i];
      }
    }
    os << endl;
    os.precision(9);
    map<IJ, LJpot>::const_iterator iter = totLJpot.begin();
    for (int i = 0; i < totLJpot.begin()->second.get_grid(); ++i) {
      os << scientific << setw(20) << iter->second.r(i);
      for (map<IJ, LJpot>::const_iterator it = totLJpot.begin(); it != totLJpot.end(); ++it) {
        os << scientific << setw(20) << totLJpot[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = totLJinter.begin(); it != totLJinter.end(); ++it) {
        os << scientific << setw(20) << totLJinter[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = totLJintra.begin(); it != totLJintra.end(); ++it) {
        os << scientific << setw(20) << totLJintra[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = intraLJ12.begin(); it != intraLJ12.end(); ++it) {
        os << scientific << setw(20) << intraLJ12[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = intraLJ13.begin(); it != intraLJ13.end(); ++it) {
        os << scientific << setw(20) << intraLJ13[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = intraLJ14.begin(); it != intraLJ14.end(); ++it) {
        os << scientific << setw(20) << intraLJ14[it->first].pot(i);
      }
      os << endl;
    }
  }
  
  void printBeadBeadDist(string fname, std::map<IJ, Distribution> &beadbeadDist, set<IJ> IJs, double rmin, double rmax, int grid) {
    ofstream fout(fname.c_str());
    double dgrid = (rmax - rmin) / grid;
    fout.precision(9);
    fout << "#" << setw(19) << "r / nm";
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      stringstream ss;
      ss << it->i() << "-" << it->j();
      fout << scientific << setw(20) << ss.str();
    }
    fout << endl;
    for (int i = 0; i < grid; ++i) {
      double r = (i + 0.5) * dgrid + rmin;
      fout << scientific << setw(20) << r;
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        if (beadbeadDist[*it].nVal() > 0) {
          fout << scientific << setw(20) << (double) beadbeadDist[*it][i] / (beadbeadDist[*it].nVal() * dgrid);
        } else {
          fout << scientific << setw(20) << beadbeadDist[*it][i];
        }
      }
      fout << endl;
    }
    fout << "#" << setw(19) << "average:";
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      if (beadbeadDist[*it].nVal() > 0) {
        fout << scientific << setw(20) << beadbeadDist[*it].ave();
      } else
        fout << setw(20) << "-";
    }
    fout << endl << "#" << setw(19) << "maximum value at:";
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      if (beadbeadDist[*it].nVal() > 0) {
        fout << scientific << setw(20) << beadbeadDist[*it].maxValAt();
      } else {
        fout << setw(20) << "-";
      }
    }
    fout << endl;
    fout.close();
  }

  void printTitleBlock(vector<int> &beadsizes, AtomSpecifier &allAtoms) {
    cout << "TITLE\n";
    cout << "   number of beads per molecule: " << beadsizes.size() << endl;
    cout << "   number of atoms per bead: ";
    for (unsigned int b = 0; b < beadsizes.size(); ++b) {
      cout << beadsizes[b] << " ";
    }
    cout << endl;
    cout << "   the molecule: |";
    int at = 0;
    for (unsigned int b = 0; b < beadsizes.size(); ++b) {
      for (int a = 0; a < beadsizes[b]; ++a) {
        cout << allAtoms.name(at);
        if (a < beadsizes[b] - 1) {
          cout << "--";
        } else if (b < beadsizes.size() - 1) {
          cout << "-|-";
        }
        at++;
      }
    }
    cout << "|\n";
    time_t rawtime;
    time(&rawtime);
    cout << "   Timestamp: " << ctime(&rawtime) << "END\n";
  }

  void calcEpsSigma(map<IJ, double> &epsilons, map<IJ, double> &sigmas, map<IJ, LJpot> &LJ) {
    if (sigmas.size() != 0) {
      sigmas.clear();
    }
    if (epsilons.size() != 0) {
      epsilons.clear();
    } 
    for (map<IJ, LJpot>::const_iterator it = LJ.begin(); it != LJ.end(); ++it) {
      double sigma = 0.0;
      double epsilon = 0.0;
      if (LJ[it->first].wasUsed()) {
        sigma = LJ[it->first].rmin() / pow(2.0, 1.0 / 6.0);
        epsilon = -LJ[it->first].potmin();
      }
      sigmas.insert(pair<IJ, double> (it->first, sigma));
      epsilons.insert(pair<IJ, double>(it->first, epsilon));
    }
  }
  
  void calcC12C6(map<IJ, double> &C12s, map<IJ, double> &C6s, map<IJ, double> &epsilons, map<IJ, double> &sigmas) {
    for(map<IJ, double>::const_iterator it =epsilons.begin(); it != epsilons.end(); ++it) {
      C6s.insert(pair<IJ, double>(it->first, 4 * epsilons[it->first] * pow(sigmas[it->first], 6)));
      C12s.insert(pair<IJ, double>(it->first, 4 * epsilons[it->first] * pow(sigmas[it->first], 12)));
    }
  }
  
  void printLennardJonesParamters(string title, map<IJ, double> &epsilons_tot,
          map<IJ, double> &epsilons_totinter,
          map<IJ, double> &epsilons_totintra,
          map<IJ, double> &epsilons_intra12,
          map<IJ, double> &epsilons_intra13,
          map<IJ, double> &epsilons_intra14,
          map<IJ, double> &sigmas_tot,
          map<IJ, double> &sigmas_totinter,
          map<IJ, double> &sigmas_totintra,
          map<IJ, double> &sigmas_intra12,
          map<IJ, double> &sigmas_intra13,
          map<IJ, double> &sigmas_intra14,
          map<IJ, double> &C12_tot,
          map<IJ, double> &C12_totinter,
          map<IJ, double> &C12_totintra,
          map<IJ, double> &C12_intra12,
          map<IJ, double> &C12_intra13,
          map<IJ, double> &C12_intra14,
          map<IJ, double> &C6_tot,
          map<IJ, double> &C6_totinter,
          map<IJ, double> &C6_totintra,
          map<IJ, double> &C6_intra12,
          map<IJ, double> &C6_intra13,
          map<IJ, double> &C6_intra14) {
    vector<string> header;
    vector<string> names;
    names.push_back("eps_tot_");
    names.push_back("eps_totinter_");
    names.push_back("eps_totintra_");
    names.push_back("eps_intra12_");
    names.push_back("eps_intra13_");
    names.push_back("eps_intra14_");
    stringstream Title;
    Title << title << "_LENNARD-JONES\n";
    cout << Title.str();
    cout << "# var_xy_i-j:\n";
    cout << "#   var ... eps:      epsilon\n";
    cout << "#       ... sig:      sigma\n";
    cout << "#   C12 ... C12 LJ-potential parameter\n";
    cout << "#   C6  ... C6 LJ-potential parameter\n";
    cout << "#   xy  ... tot:      total LJ-potential energy\n";
    cout << "#       ... totinter: total inter-molecular LJ potential\n";
    cout << "#       ... totintra: total intra-molecular LJ potential\n";
    cout << "#       ... intra12:  total 12-intra-molecular LJ potential\n";
    cout << "#       ... intra13:  total 13-intra-molecular LJ potential\n";
    cout << "#       ... intra14:  total 14-intra-molecular LJ potential\n";
    cout << "#   i   ... interaction from bead of size i to ...\n";
    cout << "#   j   ... bead of size j\n";
    cout << "#\n";
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::const_iterator it = epsilons_tot.begin(); it != epsilons_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    cout.precision(9);
    cout.setf(ios::scientific);
    for (map<IJ, double>::const_iterator it = epsilons_tot.begin(); it != epsilons_tot.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_totinter.begin(); it != epsilons_totinter.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_totintra.begin(); it != epsilons_totintra.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_intra12.begin(); it != epsilons_intra12.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_intra13.begin(); it != epsilons_intra13.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_intra14.begin(); it != epsilons_intra14.end(); ++it) {
      cout << setw(20) << it->second;
    }
    cout << "\n#\n";
    header.clear();
    names.clear();
    names.push_back("sig_tot_");
    names.push_back("sig_totinter_");
    names.push_back("sig_totintra_");
    names.push_back("sig_intra12_");
    names.push_back("sig_inter13_");
    names.push_back("sig_inter14_");
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::iterator it = sigmas_tot.begin(); it != sigmas_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    for (map<IJ, double>::const_iterator it = sigmas_tot.begin(); it != sigmas_tot.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_totinter.begin(); it != sigmas_totinter.end(); ++it) {
      cout << setw(20) << it->second;
    }
   for (map<IJ, double>::const_iterator it = sigmas_totintra.begin(); it != sigmas_totintra.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_intra12.begin(); it != sigmas_intra12.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_intra13.begin(); it != sigmas_intra13.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_intra14.begin(); it != sigmas_intra14.end(); ++it) {
      cout << setw(20) << it->second;
    }
    cout << "\n#\n";
    header.clear();
    names.clear();
    names.push_back("C12_tot_");
    names.push_back("C12_totinter_");
    names.push_back("C12_totintra_");
    names.push_back("C12_intra12_");
    names.push_back("C12_inter13_");
    names.push_back("C12_inter14_");
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::iterator it = C12_tot.begin(); it != C12_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    for (map<IJ, double>::const_iterator it = C12_tot.begin(); it != C12_tot.end(); ++it) {
      cout << setw(20) << C12_tot[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_totinter.begin(); it != C12_totinter.end(); ++it) {
      cout << setw(20) << C12_totinter[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_totintra.begin(); it != C12_totintra.end(); ++it) {
      cout << setw(20) << C12_totintra[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_intra12.begin(); it != C12_intra12.end(); ++it) {
      cout << setw(20) << C12_intra12[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_intra13.begin(); it != C12_intra13.end(); ++it) {
      cout << setw(20) << C12_intra13[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_intra14.begin(); it != C12_intra14.end(); ++it) {
      cout << setw(20) << C12_intra14[it->first];
    }
    cout << "\n#\n";
    header.clear();
    names.clear();
    names.push_back("C6_tot_");
    names.push_back("C6_totinter_");
    names.push_back("C6_totintra_");
    names.push_back("C6_intra12_");
    names.push_back("C6_inter13_");
    names.push_back("C6_inter14_");
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::iterator it = C12_tot.begin(); it != C12_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    for (map<IJ, double>::const_iterator it = C6_tot.begin(); it != C6_tot.end(); ++it) {
      cout << setw(20) << C6_tot[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_totinter.begin(); it != C6_totinter.end(); ++it) {
      cout << setw(20) << C6_totinter[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_totintra.begin(); it != C6_totintra.end(); ++it) {
      cout << setw(20) << C6_totintra[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_intra12.begin(); it != C6_intra12.end(); ++it) {
      cout << setw(20) << C6_intra12[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_intra13.begin(); it != C6_intra13.end(); ++it) {
      cout << setw(20) << C6_intra13[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_intra14.begin(); it != C6_intra14.end(); ++it) {
      cout << setw(20) << C6_intra14[it->first];
    }
    cout << "\nEND\n";
  }

}

