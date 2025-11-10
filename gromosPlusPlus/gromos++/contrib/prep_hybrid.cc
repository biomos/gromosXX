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
 * @file prep_hybrid.cc
 * Generates the files for a FG/CG hybrid simulation with a FG solvent layer
 */

/**
 * @page contrib Contrib Documentation
 *
 * @anchor prep_hybrid
 * @section prep_hybrid Generates the files for a FG/CG hybrid simulation with a FG solvent layer
 * @author @ref sr
 * @date 07-08-12
 *
 * Program prep_hybrid generates the files for a fine-grained (FG)/coarse-grained (CG) hybrid 
 * simulation with a FG solvent layer. Input coordinates are the gathered (!) solute in FG
 * solvent. The FG solvent layer can either be chosen (\@cuttype) as a sphere with \@radius around 
 * the center-of-mass (COM) of the FG solute, or as a layer with a thickness \@radius around the 
 * FG solute. For the layer, two values have to be given for \@radius. The first one is as 
 * mentioned the thickness of the FG solvent layer, the second is the radius of the distance 
 * restraining sphere around the COM of the FG solute.
 *  
 * Outside of the FG solvent layer, CG solvent is added (boxsize defined through \@minwall). 
 * Files generated: coordinate file ("hybrid_coordinates.cnf") and topology of the hybrid system
 * ("hybrid_topology.top"), and distance restraints file ("hybrid_disres.dat") of the FG solvent 
 * to COM of FG solute. The distance restraints are between the first atom of each FG solvent
 * molecule and the COM of the protein approximated by the 4 solute atoms given in \@com.
 * 
 * Note: In case a layer of FG solvent is chosen, the following tcl script allows to check in VMD
 * if all FG solvent molecules are within the specified sphere for distance restraining:
 * 
 * set selection [atomselect top "index x1 x2 x3 x4"]  // where x1,..,x4 are the 4 solute atoms given in \@com\n
 *                                                     // Note: index numbers start with 0 in VMD (not 1 as in GROMOS)!\n
 * set com [measure center $selection weight mass]\n
 * set mat Transparent\n
 * graphics top material $mat                          // sphere will be transparent (easier to check...)\n
 * graphics top sphere $com radius X resolution 80     // where X is the second number given in \@radius
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;solute topology with FG solvent&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions (r or t)&gt; </td></tr>
 * <tr><td> \@coord</td><td>&lt;coordinate file (cnf) of FG solute in FG solvent (gathered!)&gt; </td></tr>
 * <tr><td> \@fg_topo</td><td>&lt;topology of a single FG solvent molecule&gt; </td></tr>
 * <tr><td> \@cg_topo</td><td>&lt;topology of a single CG bead&gt; </td></tr>
 * <tr><td> \@cg_coord</td><td>&lt;coordinate file (cnf) of a box with CG beads&gt; </td></tr>
 * <tr><td> \@cuttype</td><td>&lt;type of FG solvent domain: layer or sphere (default: sphere)&gt; </td></tr>
 * <tr><td> \@radius</td><td>&lt;cutoff for FG solvent layer around FG solute: 1 value for sphere, 2 values for layer&gt; </td></tr>
 * <tr><td> [\@minwall</td><td>&lt;minimum FG region to wall distance [nm]&gt;] </td></tr>
 * <tr><td> [\@boxsize</td><td>&lt;(use boxsize specified in input file)&gt;] </td></tr>
 * <tr><td> [\@thresh</td><td>&lt;minimum solvent to solute distance (default 0.23 nm)&gt;] </td></tr>
 * <tr><td> \@com</td><td>&lt;@ref AtomSpecifier atomspecifier of four solute atoms approximating the COM of the solute (for distance restraints)&gt; </td></tr>
 * <tr><td> [\@disres</td><td>&lt;parameters for distance restraints: \<DISH DISC w0\> (default: 0.1 0.153 1.0)&gt;] </td></tr>
 * <tr><td> [\@solute</td><td>&lt;solute atoms&gt; (only, if cuttype == swd)] </td></tr>
 * <tr><td> [\@exponent</td><td>&lt;the exponent for calculating the swd&gt; (only, if cuttype == swd; default: 6)] </td></tr>
 * </table>
 *
 * Example 1: Sphere
 * @verbatim
  prep_hybrid
    @topo       protein.top
    @pbc        r
    @coord      protein.cnf
    @fg_topo    spc.top
    @cg_topo    cg_water.top
    @cg_coord   cg_waterbox.cnf
    @cuttype    sphere
    @radius     2.2
    @minwall    1.0
    @thresh     0.32
    @com        1:15,53,108,210
    @disres     0.1 0.153 1.0
 @endverbatim
 * 
 * Example 2: Layer
 * @verbatim
  prep_hybrid
    @topo       protein.top
    @pbc        r
    @coord      protein.cnf
    @fg_topo    spc.top
    @cg_topo    cg_water.top
    @cg_coord   cg_waterbox.cnf
    @cuttype    layer
    @radius     0.8 2.2
    @minwall    1.0
    @thresh     0.32
    @com        1:15,53,108,210
    @disres     0.1 0.153 1.0
 @endverbatim
 *
 * Example 3: Solute Weighted Distance (swd)
 * @verbatim
  prep_hybrid
    @topo       protein.top
    @pbc        r
    @coord      protein.cnf
    @fg_topo    spc.top
    @cg_topo    cg_water.top
    @cg_coord   cg_waterbox.cnf
    @cuttype    layer
    @radius     0.8 2.2
    @minwall    1.0
    @thresh     0.32
    @com        1:15,53,108,210
    @solute     1:C,N,CA,O       // These should be the same as for program swd
 @endverbatim
 * 
 * See also @ref AtomSpecifier and @ref swd
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace gmath;
using namespace bound;
using namespace utils;

double calc_max_size(System &sys, AtomSpecifier &as, int dim);
void write_disres(ofstream &out, AtomSpecifier &com, std::vector<double> disresinfo,
        int num, int firstatom, int num_atoms_per_solvent);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "coord" << "cg_topo" << "cg_coord" << "radius" << "minwall"
          << "boxsize" << "thresh" << "disres" << "fg_topo" << "com" << "cuttype"
          << "solute" << "fgsolv" << "cgsolv" << "exponent";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <solute topology with FG solvent>\n";
  usage += "\t@pbc       <periodic boundary conditions (r or t)>\n";
  usage += "\t@coord     <coordinate file (cnf) of FG solute in FG solvent (gathered!)>\n";
  usage += "\t@fg_topo   <topology of a single FG solvent molecule>\n";
  usage += "\t@cg_topo   <topology of a single CG bead>\n";
  usage += "\t@cg_coord  <coordinate file (cnf) of a box with CG beads>\n";
  usage += "\t@cuttype   <type of FG solvent domain: layer, sphere or swd [solute weighted distance] (default: sphere)>\n";
  usage += "\t@radius    <cutoff for FG solvent layer around FG solute: 1 value for sphere, swd, 2 values for layer>\n";
  usage += "\t[@minwall  <minimum FG region to wall distance [nm]>]\n";
  usage += "\t[@boxsize  (use boxsize specified in input file)]\n";
  usage += "\t[@thresh   <minimum solvent-solute distance (default 0.23 nm)>]\n";
  usage += "\t@com       <atomspecifier of four solute atoms approximating the COM of the solute>\n";
  usage += "\t[@disres   <parameters for distance restraints: <DISH DISC w0> (default: 0.1 0.153 1.0)>]\n";
  usage += "\t[@solute     <solute atoms> (only, if cuttype == swd)]\n";
  usage += "\t[@exponent   <the exponent for calculating the swd> (only, if cuttype == swd; default: 6)]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topologies
    InTopology it(args["topo"]);
    System sys(it.system());

    InTopology fgit(args["fg_topo"]);
    System fgsolv(fgit.system());

    InTopology cgit(args["cg_topo"]);
    System cgsys(cgit.system());

    // some variables for statistics
    int num_FG_solvent = 0;
    int num_CG_solvent = 0;

    // Read cuttype and radius
    string cuttype = "sphere";
    if (args.count("cuttype") > 0) {
      cuttype = args["cuttype"];
      if (!(cuttype == "sphere"
              || cuttype == "layer"
              || cuttype == "swd"))
        throw gromos::Exception("prep_hybrid",
              "@cuttype must be either sphere or layer");
    }
    std::vector<double> radius;
    if (cuttype == "sphere" || cuttype == "swd") {
      radius = args.getValues<double>("radius", 1, true);
      radius.push_back(radius[0]); // dummy value to make implementation easier
    } else { // layer: requires two values! 
      radius = args.getValues<double>("radius", 2, true);
    }
    if (radius[0] < 0.0 || radius[1] < 0.0)
      throw gromos::Exception("prep_hybrid",
            "@radius must be greater than 0!");

    // read the minimum solute to wall distances.
    // three possibilities:
    // 1. nothing specified: the user will specify the box size
    //    using the @boxsize flag
    // 2. one number specified: the box will be cubic (pbc: r) or a
    //    trunc. oct. (pbc: t) with the size defined by maximum solute
    //    atom-atom distance plus twice this number
    // 3. three numbers specified: the box will be rectangular.
    //    the specified numbers are added to the
    //    maximum distances in x, y and z (after possibly rotating the solute)

    vector<double> minwall;
    {
      Arguments::const_iterator iter = args.lower_bound("minwall"),
              to = args.upper_bound("minwall");
      while (iter != to) {
        minwall.push_back(atof(iter->second.c_str()));
        ++iter;
      }
    }

    // check for the boxsize flag
    // if it is given, the box from the solute coordinates is used
    bool boxsize = false;
    if (args.count("boxsize") != -1)
      boxsize = true;

    // read the minimum solvent-solute distance
    double minsol = args.getValue<double>("thresh", false, 0.23);
    double minsol2 = minsol * minsol;

    // read distance restraints info
    std::vector<double> disresinfo = args.getValues<double>("disres", 3, false,
            Arguments::Default<double>() << 0.1 << 0.153 << 1.0);
    disresinfo.push_back(radius[1]); // second radius is for distance restraining!

    // read FG coordinates into the system
    InG96 ic;
    ic.open(args["coord"]);
    // we also read in any solvent that is already in the file
    ic.select("ALL");
    ic >> sys;
    ic.close();

    // move the solute to the centre of geometry
    Vec shiftcog = fit::PositionUtils::shiftToCog(&sys);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // Construct a system to write out the coordinates
    System outsys;
    outsys.box() = sys.box();
    outsys.hasBox = true;

    int num_atoms_per_FGsolv = fgsolv.mol(0).numAtoms();
    int firstatom = 1; // info for distance restraints

    // add FG solute and coordinates to outsys  // same as before
    int actualNumSoluteMol = sys.numMolecules();
    for (int m = 0; m < sys.numMolecules(); m++) {
      outsys.addMolecule(sys.mol(m));
      firstatom += sys.mol(m).numAtoms();
      if (m == 0) {
        outsys.addPressureGroup(sys.mol(m).numAtoms());
        outsys.addTemperatureGroup(sys.mol(m).numAtoms());
      } else {
        outsys.addPressureGroup(outsys.pressureGroup(outsys.numPressureGroups() - 1) + sys.mol(m).numAtoms());
        outsys.addTemperatureGroup(outsys.temperatureGroup(outsys.numTemperatureGroups() - 1) + sys.mol(m).numAtoms());
      }
      outsys.mol(m).initPos();
      for (int a = 0; a < sys.mol(m).numAtoms(); a++) {
        outsys.mol(m).pos(a) = sys.mol(m).pos(a);
      }
    }

    if (cuttype == "sphere" || cuttype == "layer") {

      // Calculate center of mass of solute
      Vec com(0.0, 0.0, 0.0);
      double molmass = 0;
      for (int m = 0; m < sys.numMolecules(); m++) {
        for (int a = 0; a < sys.mol(m).numAtoms(); a++) {
          com += sys.mol(m).pos(a) * sys.mol(m).topology().atom(a).mass();
          molmass += sys.mol(m).topology().atom(a).mass();
        }
      }
      com /= molmass;

      // loop over FG solvent molecules and check if they are within radius
      // dependent on cuttype !!
      if (cuttype == "sphere") {
        for (int m = 0; m < sys.sol(0).numAtoms(); m += num_atoms_per_FGsolv) {
          Vec dist = pbc->nearestImage(sys.sol(0).pos(m), com, sys.box()) - sys.sol(0).pos(m);
          if (dist.abs() <= radius[0]) { // within sphere
            outsys.addMolecule(fgsolv.mol(0));
            num_FG_solvent++;
            int lastmol = outsys.numMolecules() - 1;
            if (outsys.numPressureGroups() == 0) {
              outsys.addPressureGroup(outsys.pressureGroup(num_atoms_per_FGsolv));
              outsys.addTemperatureGroup(outsys.temperatureGroup(num_atoms_per_FGsolv));
            } else {
              outsys.addPressureGroup(outsys.pressureGroup(lastmol - 1) + num_atoms_per_FGsolv);
              outsys.addTemperatureGroup(outsys.temperatureGroup(lastmol - 1) + num_atoms_per_FGsolv);
            }
            outsys.mol(lastmol).initPos();
            for (int a = 0; a < outsys.mol(lastmol).numAtoms(); a++) {
              outsys.mol(lastmol).pos(a) = sys.sol(0).pos(a + m);
            }
          }
        }
      } else { // layer
        double min_init = sys.box().K()[0] * sys.box().K()[0] + sys.box().L()[1] * sys.box().L()[1] +
                sys.box().M()[2] * sys.box().M()[2];
        Vec o(0.0, 0.0, 0.0);
        double radius2 = radius[0] * radius[0];
        for (int i = 0; i < sys.sol(0).numAtoms(); i += num_atoms_per_FGsolv) {
          double min2 = min_init;
          Vec check = pbc->nearestImage(o, sys.sol(0).pos(i), sys.box());
          for (int m = 0; m < sys.numMolecules(); m++) {
            for (int a = 0; a < sys.mol(m).numAtoms(); a++) {
              if (!sys.mol(m).topology().atom(a).isH() &&
                      (check - sys.mol(m).pos(a)).abs2() < min2)
                min2 = (check - sys.mol(m).pos(a)).abs2();
            }
          }

          if (min2 < radius2) {
            // yes! we keep this solvent
            outsys.addMolecule(fgsolv.mol(0));
            num_FG_solvent++;
            int lastmol = outsys.numMolecules() - 1;
            if (outsys.numPressureGroups() == 0) {
              outsys.addPressureGroup(outsys.pressureGroup(num_atoms_per_FGsolv));
              outsys.addTemperatureGroup(outsys.temperatureGroup(num_atoms_per_FGsolv));
            } else {
              outsys.addPressureGroup(outsys.pressureGroup(lastmol - 1) + num_atoms_per_FGsolv);
              outsys.addTemperatureGroup(outsys.temperatureGroup(lastmol - 1) + num_atoms_per_FGsolv);
            }
            outsys.mol(lastmol).initPos();
            for (int a = 0; a < outsys.mol(lastmol).numAtoms(); a++) {
              outsys.mol(lastmol).pos(a) = sys.sol(0).pos(a + i);
            }
          }
        }
      }
    } else { // swd
      if (args.count("solute") == -1) {
        throw (gromos::Exception("prep_hybrid",
                "For cuttype option 'swd', argument 'solute' is needed!"));
      }
      
      utils::AtomSpecifier solute(sys);
      for (args::Arguments::const_iterator it = args.lower_bound("solute"),
              to = args.upper_bound("solute");
           it != to; ++it) {
        solute.addSpecifier(it->second);
      }
      
      const double n = args.getValue("exponent", false, 6.0);
      for (int j = 0; j < sys.sol(0).numAtoms(); j += num_atoms_per_FGsolv){
        double sum = 0.0;
        gmath::Vec &solv_j = sys.sol(0).pos(j);
        for (int i = 0; i < solute.size(); i++){
          gmath::Vec nim_sol = pbc->nearestImage(solv_j, solute.pos(i), sys.box());
          gmath::Vec dist = solv_j - nim_sol;
          double r_ij = dist.abs();
          sum += pow(r_ij, -n);
        }
        double dist = pow(sum, -1.0 / n);
        if (dist < radius[0]) {
          outsys.addMolecule(fgsolv.mol(0));
          num_FG_solvent++;
          int lastmol = outsys.numMolecules() - 1;
          if (outsys.numPressureGroups() == 0) {
            outsys.addPressureGroup(outsys.pressureGroup(num_atoms_per_FGsolv));
            outsys.addTemperatureGroup(outsys.temperatureGroup(num_atoms_per_FGsolv));
          } else {
            outsys.addPressureGroup(outsys.pressureGroup(lastmol - 1) + num_atoms_per_FGsolv);
            outsys.addTemperatureGroup(outsys.temperatureGroup(lastmol - 1) + num_atoms_per_FGsolv);
          }
          outsys.mol(lastmol).initPos();
          for (int a = 0; a < outsys.mol(lastmol).numAtoms(); a++) {
            outsys.mol(lastmol).pos(a) = sys.sol(0).pos(a + j);
          }
        }
      }
    }

    // determine box shape, calculate relevant things and check the 
    // input for consistency
    vector<double> max_dist;
    AtomSpecifier as(outsys);

    if (minwall.size() == 3)
      throw gromos::Exception("prep_hybrid",
            "solute should be rotated in order to align largest extension with z axis");
    max_dist.push_back(calc_max_size(outsys, as, 3));

    enum boundary_enum {
      vacuum, rectangular, cubic, truncoct, triclinic
    } boundary = vacuum;

    double size_corr = 1.0;

    switch (pbc->type()) {
      case('t'):
        boundary = truncoct;
        size_corr = 2.0 * sqrt(3.0) / 3.0;
        if (minwall.size() > 1)
          throw (gromos::Exception("prep_hybrid",
                "For truncated octahedral boxes you can only specify one number for @minwall"));

        if (minwall.size() == 0) {
          if (outsys.box().K().abs() != outsys.box().L().abs() || outsys.box().K().abs() != outsys.box().M().abs())
            throw (gromos::Exception("prep_hybrid",
                  "For truncated octahedral boxes, the specified boxsize should be the same in all dimensions"));
        }
        break;
      case('r'):
        boundary = rectangular;
        if (minwall.size() == 1)
          boundary = cubic;
        else if (minwall.size() == 0 && (outsys.box().ntb() == gcore::Box::rectangular)
                && (outsys.box().K().abs() == outsys.box().L().abs()) && (outsys.box().K().abs() == outsys.box().M().abs()))
          boundary = cubic;
        break;
      case('c'):
        if (!boxsize)
          throw gromos::Exception("prep_hybrid", "boxsize has to be specified for triclinic");
        boundary = triclinic;
        break;

      case('v'):
        throw (gromos::Exception("prep_hybrid", "Why are you running this program if @pbc is vacuum?"));
        break;
    }

    // size of the solvent box
    vector<double> solvent_box(3, 0.0);

    if (boxsize) {
      if (minwall.size())
        throw (gromos::Exception("prep_hybrid", "Cannot specify both boxsize and minwall."));

      if (!outsys.hasBox)
        throw gromos::Exception("prep_hybrid",
              "If you specify boxsize, the box dimensions should be in "
              "the GENBOX block of the solute");

      if (pbc->type() == 't') {
        if (outsys.box().ntb() != gcore::Box::truncoct)
          throw gromos::Exception("prep_hybrid",
                "Specified truncated octahedral for @pbc, and try to read "
                "boxsize from file\nbut no truncated octahedral specified there");

        solvent_box[0] = outsys.box().K()[0];
        solvent_box[1] = outsys.box().L()[1];
        solvent_box[2] = outsys.box().M()[2];
      } else if (pbc->type() == 'r') {
        if (outsys.box().ntb() != gcore::Box::rectangular)
          throw gromos::Exception("prep_hybrid",
                "Specified rectangular box for @pbc, and try to read boxsize from "
                "file\nbut no rectangular box specified there");
        solvent_box[0] = outsys.box().K()[0];
        solvent_box[1] = outsys.box().L()[1];
        solvent_box[2] = outsys.box().M()[2];
      } else if (pbc->type() == 'c') {
        if (outsys.box().ntb() != gcore::Box::triclinic)
          throw gromos::Exception("prep_hybrid",
                "Specified triclinic box for @pbc, and try to read boxsize from "
                "file\nbut no triclinic box specified there");
        // calculate the dimension of a large enough box to encompass the triclinic box
        // first construct the four diagonals
        vector<Vec> diag(4);
        diag[0] = 0.5 * (outsys.box().K() + outsys.box().L() + outsys.box().M());
        diag[1] = 0.5 * (outsys.box().K() + outsys.box().L() - outsys.box().M());
        diag[2] = 0.5 * (outsys.box().K() - outsys.box().L() + outsys.box().M());
        diag[3] = 0.5 * (outsys.box().K() - outsys.box().L() - outsys.box().M());

        // find the maximum and minimum x,y,z, store these in boxsize
        for (int j = 0; j < 3; ++j) {
          double maxdim = diag[0][j], mindim = diag[0][j];
          for (int i = 1; i < 4; ++i) {
            if (diag[i][j] > maxdim) maxdim = diag[i][j];
            if (diag[i][j] < mindim) mindim = diag[i][j];
            if (-diag[i][j] > maxdim) maxdim = -diag[i][j];
            if (-diag[i][j] < mindim) mindim = -diag[i][j];
          }
          solvent_box[j] = maxdim - mindim;
        }
      }// triclinic
      else {
        throw gromos::Exception("prep_hybrid", "unknown boundary condition");
      }
      // check if all atoms are actually within the box.
      // If not, write a warning and put them in.
      bool warn_outside_box = false;
      for (int m = 0; m < outsys.numMolecules(); m++) {
        for (int a = 0; a < outsys.mol(m).numAtoms(); a++) {

          // are we inside the box
          Vec check = pbc->nearestImage(Vec(0.0, 0.0, 0.0), outsys.mol(m).pos(a), outsys.box());

          if (check[0] != outsys.mol(m).pos(a)[0] ||
                  check[1] != outsys.mol(m).pos(a)[1] ||
                  check[2] != outsys.mol(m).pos(a)[2]) {
            outsys.mol(m).pos(a) = check;
            warn_outside_box = true;
          }
        }
      }
      if (warn_outside_box)
        cerr << "WARNING: not all atoms were within the specified box !\n"
              << "         placed the atoms in the box before solvation\n"
              << "         according to the specified box.\n";

    }// boxsize
    else {
      if (minwall.size() == 0)
        throw gromos::Exception("prep_hybrid",
              "either use a specified boxsize "
              "or give a minimum distance from solute to the walls");
      if (boundary == truncoct) {
        outsys.box().setNtb(gcore::Box::truncoct);
        outsys.box().K()[0] = size_corr * (max_dist[0] + 2 * minwall[0]);
        solvent_box[0] = outsys.box().K()[0];
        outsys.box().L()[1] = size_corr * (max_dist[0] + 2 * minwall[0]);
        solvent_box[1] = outsys.box().L()[1];
        outsys.box().M()[2] = size_corr * (max_dist[0] + 2 * minwall[0]);
        solvent_box[2] = outsys.box().M()[2];
      } else {
        outsys.box().setNtb(gcore::Box::rectangular);
        if (minwall.size() == 1) {
          outsys.box().K()[0] = max_dist[0] + 2 * minwall[0];
          solvent_box[0] = outsys.box().K()[0];
          outsys.box().L()[1] = max_dist[0] + 2 * minwall[0];
          solvent_box[1] = outsys.box().L()[1];
          outsys.box().M()[2] = max_dist[0] + 2 * minwall[0];
          solvent_box[2] = outsys.box().M()[2];
        } else {
          // the solute has been rotated, 3 max_dist are known, 3 minwall distances are given
          outsys.box().K()[0] = max_dist[2] + 2 * minwall[0];
          solvent_box[0] = outsys.box().K()[0];
          outsys.box().L()[1] = max_dist[1] + 2 * minwall[1];
          solvent_box[1] = outsys.box().L()[1];
          outsys.box().M()[2] = max_dist[0] + 2 * minwall[2];
          solvent_box[2] = outsys.box().M()[2];
        }
      }
    }

    // Now add the CG solvent beads!!
    SolventTopology cgsolv;
    for (int a = 0; a < cgsys.mol(0).numAtoms(); a++) {
      cgsolv.addAtom(cgsys.mol(0).topology().atom(a));
    }
    cgsys.sol(0) = cgsolv;

    // read CG coordinates into the system
    InG96 cgic;
    cgic.open(args["cg_coord"]);
    cgic.select("SOLVENT");
    cgic >> cgsys;
    cgic.close();

    if (!cgsys.hasBox)
      throw gromos::Exception("prep_hybrid",
            "Could not read BOX block from CG solvent "
            "coordinates");

    // calculate the CG solvent COG
    Vec solv_cog(0.0, 0.0, 0.0);
    int num_solv_atoms_per_box = cgsys.sol(0).numPos();
    for (int i = 0; i < num_solv_atoms_per_box; i++)
      solv_cog += cgsys.sol(0).pos(i);
    solv_cog /= num_solv_atoms_per_box;
    // move to CG solvent cog
    for (int i = 0; i < num_solv_atoms_per_box; i++)
      cgsys.sol(0).pos(i) -= solv_cog;


    // set the solute box size, calculate how many solvent boxes we should
    // use; calculate where to move the initial box
    vector<int> needed_boxes;
    Vec move_solvent;
    needed_boxes.push_back(int(solvent_box[0] / cgsys.box().K()[0]) + 1);
    move_solvent[0] = -0.5 * (needed_boxes[0] - 1) * cgsys.box().K()[0];
    needed_boxes.push_back(int(solvent_box[1] / cgsys.box().L()[1]) + 1);
    move_solvent[1] = -0.5 * (needed_boxes[1] - 1) * cgsys.box().L()[1];
    needed_boxes.push_back(int(solvent_box[2] / cgsys.box().M()[2]) + 1);
    move_solvent[2] = -0.5 * (needed_boxes[2] - 1) * cgsys.box().M()[2];

    // move the initial box so that after multiplication the complete box
    // is centered around the origin.
    for (int i = 0; i < num_solv_atoms_per_box; i++)
      cgsys.sol(0).pos(i) += move_solvent;

    // do the multiplications
    for (int ix = 0; ix < needed_boxes[0]; ix++) {
      for (int iy = 0; iy < needed_boxes[1]; iy++) {
        for (int iz = 0; iz < needed_boxes[2]; iz++) {
          if (ix != 0 || iy != 0 || iz != 0) {
            Vec shift(ix * cgsys.box().K()[0], iy * cgsys.box().L()[1], iz * cgsys.box().M()[2]);

            for (int atom = 0; atom < num_solv_atoms_per_box; atom++) {
              cgsys.sol(0).addPos(cgsys.sol(0).pos(atom) + shift);
            }
          }
        } // iz
      } // iy
    } // ix

    int num_atoms_per_solvent = cgsys.mol(0).numAtoms();
    int num_solvent_molecules = cgsys.sol(0).numPos() / num_atoms_per_solvent;
    double minCGdist2 = radius[0] * radius[0];

    if (cuttype == "sphere" || cuttype == "layer") {
      // now we have to keep only those waters that are inside the box and 
      // far enough away from the solute
      // we look at the centre of geometry of the solvent molecule
      Vec o(0.0, 0.0, 0.0);
      double min_init = solvent_box[0] * solvent_box[0] + solvent_box[1] * solvent_box[1] +
              solvent_box[2] * solvent_box[2];

      for (int i = 0; i < num_solvent_molecules; i++) {
        // calculate the centre of geometry of this solvent
        Vec sol_i(0.0, 0.0, 0.0);
        for (int j = 0; j < num_atoms_per_solvent; j++)
          sol_i += cgsys.sol(0).pos(num_atoms_per_solvent * i + j);
        sol_i /= num_atoms_per_solvent;

        // are we inside the box
        Vec check = pbc->nearestImage(o, sol_i, outsys.box());

        if (check[0] == sol_i[0] && check[1] == sol_i[1] && check[2] == sol_i[2]) {
          // yes we are in the box
          // calculate the closest distance to any solute
          double min2 = min_init;
          for (int m = 0; m < actualNumSoluteMol; m++) {
            for (int a = 0; a < outsys.mol(m).numAtoms(); a++) {
              if (!outsys.mol(m).topology().atom(a).isH() &&
                      (check - outsys.mol(m).pos(a)).abs2() < min2)
                min2 = (check - outsys.mol(m).pos(a)).abs2();
            }
          }

          if (min2 > minsol2 && min2 > minCGdist2) {
            // yes! we keep this solvent
            outsys.addMolecule(cgsys.mol(0));
            num_CG_solvent++;
            int lastmol = outsys.numMolecules() - 1;
            outsys.addPressureGroup(outsys.pressureGroup(lastmol - 1) + num_atoms_per_solvent);
            outsys.addTemperatureGroup(outsys.temperatureGroup(lastmol - 1) + num_atoms_per_solvent);
            outsys.mol(lastmol).initPos();
            for (int a = 0; a < cgsys.mol(0).numAtoms(); a++) {
              outsys.mol(lastmol).pos(a) = cgsys.sol(0).pos(num_atoms_per_solvent * i + a);
            }
          }
        }
      }
    } else { // swd
      
      utils::AtomSpecifier solute(sys);
      for (args::Arguments::const_iterator it = args.lower_bound("solute"),
              to = args.upper_bound("solute");
           it != to; ++it) {
        solute.addSpecifier(it->second);
      }
      
      const double n = args.getValue("exponent", false, 6.0);
      Vec o(0.0, 0.0, 0.0);
      for (int j = 0; j < cgsys.sol(0).numAtoms(); j += num_atoms_per_solvent){
        double sum = 0.0;
        gmath::Vec &solv_j = cgsys.sol(0).pos(j);
        Vec check = pbc->nearestImage(o, solv_j, outsys.box());
        if (check == solv_j) {
          for (int i = 0; i < solute.size(); i++) {
            // check if cg is in box
            gmath::Vec nim_sol = pbc->nearestImage(solv_j, solute.pos(i), sys.box());
            gmath::Vec dist = solv_j - nim_sol;
            double r_ij = dist.abs();
            sum += pow(r_ij, -n);
          }
          double dist = pow(sum, -1.0 / n);
          if (dist > radius[0]) {
            outsys.addMolecule(cgsys.mol(0));
            num_CG_solvent++;
            int lastmol = outsys.numMolecules() - 1;
            outsys.addPressureGroup(outsys.pressureGroup(lastmol - 1) + num_atoms_per_solvent);
            outsys.addTemperatureGroup(outsys.temperatureGroup(lastmol - 1) + num_atoms_per_solvent);
            outsys.mol(lastmol).initPos();
            for (int a = 0; a < outsys.mol(lastmol).numAtoms(); a++) {
              outsys.mol(lastmol).pos(a) = cgsys.sol(0).pos(a + j);
            }
          }
        }
      }
      
    }

    // Write out coordinates
    OutG96S oc;
    ofstream outcoord("hybrid_coordinates.cnf");
    oc.open(outcoord);
    stringstream title;
    title << "Hybrid FG/CG system with FG solvent layer containing " << num_FG_solvent
            << " molecules surrounded by " << num_CG_solvent << " CG beads";
    oc.writeTitle(title.str());
    oc.select("ALL");
    oc << outsys;
    oc.close();

    // add solvent info
    outsys.addSolvent(sys.sol(0));

    // Write out topology
    ofstream outtopo("hybrid_topo.top");
    OutTopology ot(outtopo);
    ot.setTitle(title.str());
    ot.write(outsys, cgit.forceField()); // take CG force field in order to have the extra atom types!!!

    // Write distance restraints file

    // read 4 atoms to approximate COM of solute
    if (cuttype != "swd") {
      AtomSpecifier comatoms(sys);
      Arguments::const_iterator iter = args.lower_bound("com"),
              to = args.upper_bound("com");
      for (; iter != to; ++iter) {
        comatoms.addSpecifierStrict(iter->second);
      }
      if (comatoms.size() < 3)
        throw gromos::Exception("prep_hybrid",
              "not enough atoms specified in @com");
      if (comatoms.size() > 4)
        throw gromos::Exception("prep_hybrid",
              "too many atoms specified in @com");

      ofstream outdisres("hybrid_disres.dat");
      write_disres(outdisres, comatoms, disresinfo, num_FG_solvent, firstatom, num_atoms_per_FGsolv);
    }

    // Write summary
    int last_solute_atom = firstatom - 1;
    int last_fgsolv_atom = last_solute_atom + num_atoms_per_FGsolv*num_FG_solvent;
    int last_cgsolv_atom = last_fgsolv_atom + num_atoms_per_solvent*num_CG_solvent;
    cout << "# SUMMARY: " << endl
            << "# FG Solute: " << sys.numMolecules() << " (last atom: " << last_solute_atom << ")" << endl
            << "# FG Solvent Molecules: " << num_FG_solvent << " (last atom: " << last_fgsolv_atom << ")" << endl
            << "# CG Solvent Molecules: " << num_CG_solvent << " (last atom: " << last_cgsolv_atom << ")" << endl;
    if (cuttype == "swd") {
      cout << "# Files written: hybrid_coordinates.cnf and hybrid_topo.top" << endl;
    } else {
      cout << "# Files written: hybrid_coordinates.cnf, hybrid_topo.top and hybrid_disres.dat" << endl;
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);

  }

  return 0;
}

double calc_max_size(System &sys, AtomSpecifier &as, int dim) {
  // calculate the longest distance between solute atoms considering
  // the first dim dimensions.
  double max2 = 0.0;
  double d2 = 0;
  int max_m1 = 0, max_m2 = 0, max_a1 = 0, max_a2 = 0;

  for (int m1 = 0; m1 < sys.numMolecules(); m1++) {
    for (int a1 = 0; a1 < sys.mol(m1).numAtoms(); a1++) {
      Vec current = sys.mol(m1).pos(a1);
      for (int m2 = m1; m2 < sys.numMolecules(); m2++) {
        int start = 0;
        if (m1 == m2) start = a1;
        for (int a2 = start; a2 < sys.mol(m2).numAtoms(); a2++) {
          d2 = 0.0;
          for (int i = 0; i < dim; i++) {
            d2 += (current[i] - sys.mol(m2).pos(a2)[i]) *
                    (current[i] - sys.mol(m2).pos(a2)[i]);
          }
          if (d2 > max2) {
            max2 = d2;
            max_m1 = m1;
            max_m2 = m2;
            max_a1 = a1;
            max_a2 = a2;
          }
        }
      }
    }
  }
  as.addAtomStrict(max_m1, max_a1);
  as.addAtomStrict(max_m2, max_a2);

  return sqrt(max2);
}

void write_disres(ofstream &out, AtomSpecifier &com, std::vector<double> disresinfo,
        int num, int firstatom, int num_atoms_per_solvent) {
  // Write title
  out << "TITLE" << endl
          << "   distance restraints file for " << num << " FG solvent molecules" << endl
          << "   to the COM of the solute approximated by 4 solute atoms." << endl
          << "   FG solvent sphere around solute: " << disresinfo[3] << endl
          << "   To be used in hybrid FG/CG simulations." << endl
          << "END" << endl;

  // Write DISTANCERESSPEC block
  out << "DISTANCERESSPEC" << endl
          << "# DISH  DISC" << endl
          << "  " << disresinfo[0] << "   " << disresinfo[1] << endl
          << "# i    j    k    l    type     i    j    k    l    type     r0    w0    rah" << endl;

  for (int m = 0; m < num; m++) {
    out << setw(6) << firstatom << " 0    0    0    0        " << com.atom(0) + 1
            << setw(5) << com.atom(1) + 1 << setw(5) << com.atom(2) + 1 << setw(5) << com.atom(3) + 1
            << "  -2       " << disresinfo[3] << "   " << disresinfo[2] << "     1" << endl;
    // Note: rah is always 1, because the distance restraints should only be
    // attractive in FG/CG hybrid simulations!!
    firstatom += num_atoms_per_solvent;
  }
  out << "END" << endl;
}
