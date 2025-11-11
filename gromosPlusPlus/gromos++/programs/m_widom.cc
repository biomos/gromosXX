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
 * @file m_widom.cc
 * calculate particle insertion free energies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor m_widom
 * @section m_widom calculate particle insertion free energies
 * @author @ref sr @ref co
 * @date 14-6-10   21-6-07
 *
 * Program m_widom can calculate the free energy of inserting a test particle
 * into configurations of a molecular system. For every configuration in the given trajectory file, the
 * program places the particle at a user specified number of random positions
 * and evaluates the nonbonded interaction energy, @f$V^{nbd}@f$. The solvation free energy
 * is calculated as
 *
 * @f[ \Delta G_S = -k_BT \ln \frac{<V e^{-V^{nbd}/k_BT}>}{<V>} @f]
 *
 * with @f$k_B@f$ the Boltzmann constant and @f$T@f$ and @f$V@f$ the temperature and volume
 * of the system. The programme will also calculate the solute-solvent energies
 * according to 
 *
 * @f[ \Delta U_{uv} = \frac{<V^{nbd} V e^{-V^{nbd}/k_BT}>}{V e^{-V^{nbd}/k_BT}} @f]
 *
 * which equals the solute-solvent enthalpy, @f$H_{uv}@f$, as no volume change 
 * upon solvation is taking place. The solute-solvent entropy is subsequently
 * calculated from 
 *
 * @f[ T \Delta S_{uv} = \Delta G_S - \Delta H_{uv} @f]
 *
 * For a more complete description of these free energies, see e.g.
 * [J.Phys.Chem.B 108, 1056 - 1064 (2004)].
 * 
 * In addition to the energetics of the system, the program can also calculate
 * radial distribution functions for all inserted species, with respect to 
 * user-specified atoms in the original system. Each group of atoms to include 
 * in the rdf-calculations is preceded by the keyword "new" in the input 
 * string. The radial distribution function is calculated as in the program
 * @ref rdf (Vol. 5, Section 4.14), where all averages are weighted with the
 * Boltzmann probability of every insertion attempt.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@intopo</td><td>&lt;topology of the inserted particle&gt; </td></tr>
 * <tr><td> \@inpos</td><td>&lt;coordinates of the inserted particle&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time dt"&gt; </td></tr>
 * <tr><td> [\@stride</td><td>&lt;take every n-th frame&gt;] </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> [\@eps</td><td>&lt;epsilon for reaction field (default: 1)&gt;] </td></tr>
 * <tr><td> [\@kap</td><td>&lt;kappa for reaction field (default: 0)&gt;] </td></tr>
 * <tr><td> [\@rdf</td><td>&lt;rdf with atom types&gt;] </td></tr>
 * <tr><td> [\@rdfparam</td><td>&lt;rdf-cutoff&gt; &lt;grid&gt;] </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> \@ntry</td><td>&lt;number of insertion tries per frame&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  m_widom
    @topo     ex.top
    @pbc      r
    @intopo   namta45a3.dat
    @inpos    nasx.dat
    @time     0 1  
    @stride   1
    @cut      1.4
    @eps      61
    @kap      0
    @rdf      s:OW
    @rdfparam 1.5 100
    @temp     298
    @ntry     10000
    @traj     ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InPtTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/utils/Energy.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Distribution.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "intopo" << "inpos"
          << "time" << "stride" << "cut" << "eps" << "kap"
          << "rdf" << "rdfparam" << "temp" << "ntry" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc       <boundary type>\n";
  usage += "\t@intopo    <topology of the inserted particle>\n";
  usage += "\t@inpos     <coordinates of the inserted particle>\n";
  usage += "\t@time     <time> <dt>\n";
  usage += "\t[@stride   <take every n-th frame>]\n";
  usage += "\t@cut       <cut-off distance>\n";
  usage += "\t[@eps      <epsilon for reaction field (default: 1)>]\n";
  usage += "\t[@kap      <kappa for reaction field (default: 0)>]\n";
  usage += "\t[@rdf      <rdf with atom types>]\n";
  usage += "\t[@rdfparam <rdf-cutoff> <grid>]\n";
  usage += "\t@temp      <temperature>\n";
  usage += "\t@ntry      <number of insertion tries per frame>\n";
  usage += "\t@traj      <trajectory files>\n";


  try {
    clock();

    Arguments args(argc, argv, knowns, usage);

    //   get simulation time
    bool usertime = args.count("time") >= 0;
    vector<double> timearg = args.getValues<double>("time", 2, usertime,
            Arguments::Default<double>() << 0.0 << 1.0);
    double t0 = timearg[0], dt = timearg[1];

    // get the stride
    int stride = args.getValue<int>("stride", false, 1);

    //  read topologies, create a system that will contain only the topology.
    InTopology it(args["topo"]);
    System systop(it.system());
    GromosForceField gff(it.forceField());

    // count the number of atoms in the system
    int start = 0;
    for (int m = 0; m < systop.numMolecules(); ++m) start += systop.mol(m).numAtoms();

    InTopology it2(args["intopo"]);
    System insys(it2.system());

    // define input coordinate
    InG96 ic;
    ic.open(args["inpos"]);
    ic >> insys;
    ic.close();

    // move the insertion particle to its first atom
    // usually it will be only one atom, but we keep our possibilities
    // open
    for (int m = 0; m < insys.numMolecules(); m++)
      for (int a = 0; a < insys.mol(m).numAtoms(); a++)
        insys.mol(m).pos(a) -= insys.mol(0).pos(0);

    // get non-bonded parameters
    double cut = atof(args["cut"].c_str());
    double eps = args.getValue<double>("eps", false, 1.0);
    double kap = args.getValue<double>("kap", false, 0.0);

    // which atom types
    vector<AtomSpecifier> rdfatoms;
    {
      Arguments::const_iterator iter = args.lower_bound("rdf");
      Arguments::const_iterator to = args.upper_bound("rdf");
      for (int i = -1; iter != to; ++iter) {
        if (iter->second == "new") {
          AtomSpecifier as(systop);
          rdfatoms.push_back(as);
          i++;
        } else {
          rdfatoms[i].addSpecifier(iter->second);
        }
      }
    }

    // prepare rdf arrays
    vector<double> rdfparam = args.getValues<double>("rdfparam", 2, false,
            Arguments::Default<double>() << cut << 100.0);
    double rdfcut = rdfparam[0];
    int rdfgrid = int(rdfparam[1]);

    vector <gmath::Distribution> rdf;
    vector <vector<double> > s_rdf;

    for (unsigned int i = 0; i < rdfatoms.size(); i++) {
      gmath::Distribution d(0, rdfcut, rdfgrid);
      rdf.push_back(d);
    }
    s_rdf.resize(rdfatoms.size());
    for (unsigned int i = 0; i < s_rdf.size(); i++) {
      s_rdf[i].resize(rdfgrid, 0.0);
    }

    double rdfdist = rdfcut / double(rdfgrid);
    double rdfcorr = 1 / (4.0 * acos(-1.0) * rdfdist);

    // read in the temperature
    double temp = args.getValue<double>("temp");
    double beta = 1.0 / (gmath::physConst.get_boltzmann() * temp);

    // read in number of tries per frame
    int ntry = args.getValue<int>("ntry");

    // define some values for averaging
    int numframes = 0;
    double fexp = 0.0; // exp(-E/kT)

    double vol, s_vol = 0.0; // volume
    double v_exp;
    double s_v_exp = 0.0; // v.exp(-E/kt)
    double s_v_Eexp = 0.0; // v.E.exp(-E/kt)

    // initialize random seed
    srand(int(clock()));

    // do we need to correct the volume for trunc-oct
    double vcorr = 1.0;

    // print a fancy title
    cout.precision(4);
    cout << "#     ";
    for (int i = 0; i < 29; i++) cout << "-";
    cout << endl;
    cout << "# time ";
    cout << setw(10) << "DG  "
              << setw(10) << "DH  "
              << setw(10) << "TDS" << endl;

    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj"); iter != to; ++iter) {

      // open file
      ic.open(iter->second.c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        numframes++;

        // create a new system to which we can add the molecules of the
        // insys;
        System sys(systop);

        System refSys(sys);

        // parse boundary conditions

        Boundary *pbc = BoundaryParser::boundary(sys, args);
        Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
        if (pbc->type() == 't') vcorr = 0.5;

        // read in the coordinates
        ic >> sys;

        if (!(numframes % stride)) {
          // we have to gather to get covalent interactions
          // and charge-groups connected

          (*pbc.*gathmethod)();

          // add the molecules of the insys
          for (int m = 0; m < insys.numMolecules(); ++m) {
            sys.addMolecule(insys.mol(m));
            sys.mol(systop.numMolecules() + m).initPos();
          }

          // create the energy class for the new system
          Energy en(sys, gff, *pbc);

          // set the atoms for which to calculate the energy
          AtomSpecifier energyatoms(sys);
          for (int m = 0; m < insys.numMolecules(); ++m) {
            for (int a = 0; a < insys.mol(m).numAtoms(); ++a) {
              energyatoms.addAtom(systop.numMolecules() + m, a);
            }
          }
          en.setAtoms(energyatoms);
          en.setCutOff(cut);
          en.setRF(eps, kap);

          // now reset the rdfatoms to point at this system
          for (unsigned int i = 0; i < rdfatoms.size(); ++i) {
            rdfatoms[i].setSystem(sys);
          }

          // we need the volume, correct for truncated octahedron!!
          sys.box().update_triclinic();
          vol = vcorr * sys.box().K_L_M();
          s_vol += vol;

          // now we can go into the loop over the trial positions
          for (int i = 0; i < ntry; i++) {
            //get three random numbers between the box dimensions
            Vec move;

            int r = rand();
            move[0] = double(sys.box().K()[0] * r) / double(RAND_MAX);
            r = rand();
            move[1] = double(sys.box().L()[1] * r) / double(RAND_MAX);
            r = rand();
            move[2] = double(sys.box().M()[2] * r) / double(RAND_MAX);
            //for(int d=0; d<3; d++){
            //  int r=rand();
            //  move[d]=double(sys.box()[d]*r)/double(RAND_MAX);
            //}

            // move the inatoms
            for (int m = 0; m < insys.numMolecules(); ++m) {
              for (int a = 0; a < insys.mol(m).numAtoms(); ++a) {
                sys.mol(systop.numMolecules() + m).pos(a) =
                        insys.mol(m).pos(a) + move;
              }
            }

            // do the rdf's
            for (unsigned int j = 0; j < rdf.size(); j++) {

              // set them to zero
              rdf[j].clear();

              // loop over the 'with' atoms
              for (unsigned int i = 0; i < rdfatoms[j].size(); i++) {
                Vec tmp;
                tmp = pbc->nearestImage(move, *rdfatoms[j].coord(i),
                        sys.box());
                rdf[j].add((tmp - move).abs());
              }
            }

            // calculate the interactions
            en.calcNb();

            // store and sum everything to the appropriate arrays
            fexp = exp(-beta * en.tot());
            v_exp = vol*fexp;

            s_v_exp += v_exp;
            s_v_Eexp += en.tot() * v_exp;

            for (unsigned int j = 0; j < rdfatoms.size(); j++) {
              for (int i = 0; i < rdfgrid; i++) {
                double r = (i + 0.5) * rdfdist;

                s_rdf[rdfatoms.size() + j][i] +=
                        v_exp * rdf[j][i] * vol * rdfcorr / (r * r * rdfatoms[j].size());
              }
            }
          } // loop over trials

          // add the time
          t0 += dt*stride;

          cout << t0 << "  ";
          double DG = -log(s_v_exp / (s_vol * ntry)) / beta;
          double DH = s_v_Eexp / s_v_exp;
          double TDS = (DH - DG);

          cout << setw(10) << DG << "  "
                  << setw(10) << DH << "  "
                  << setw(10) << TDS;
          cout << endl;
        } // if stride
      } //loop over frames
    }
    // print out averages
    cout << "#\n# Summary per species:\n";

    double divide = double(numframes / stride);
      cout << "#\n# ---------------------\n";

      cout.setf(ios::right, ios::adjustfield);
      const double DG = -log(s_v_exp / (s_vol * ntry)) / beta;
      const double DH = s_v_Eexp / s_v_exp;
      const double ds = (DH - DG) / temp;
      // const double DS=(DH-DG)/temp;

      cout << "# " << endl;
      cout << "# Number of frames          : " << numframes << endl;
      cout << "# Number of tries per frame : " << ntry << endl;
      cout << "# <V.exp(-E/kT)>   : " << s_v_exp / (divide*ntry) << endl;
      cout << "# <V.E.exp(-E/kT)> : " << s_v_Eexp / (divide*ntry) << endl;
      cout << "# <V>              : " << s_vol / divide << endl;
      cout << "# DG    : " << DG << endl;
      cout << "# DH_uv : " << DH << endl;
      cout << "# DS_uv : " << ds << endl;


      ostringstream os;
      os << "rdf_widom.out";

      ofstream fout(os.str().c_str());
      fout << "# rdf's for insertion of:" << endl;
      fout << "#    " << args["inpos"] << endl;

      fout << "# into:" << endl;
      if (args.count("traj") == 1)
        fout << "#    " << args["traj"] << endl;
      else
        fout << "#    " << args.count("traj") << " trajectory files" << endl;
      fout << "#" << endl;
      fout << "# taking rdf of test position with "
              << rdfatoms.size() << " specified groups of atoms\n";
      fout << "#" << endl;
      fout << "#" << setw(11) << "r";
      for (unsigned int j = 0; j < rdfatoms.size(); j++)
        fout << " " << setw(12) << j;
      fout << endl;
      for (int i = 0; i < rdfgrid; i++) {

        fout << setw(12) << (i + 0.5) * rdfdist;
        for (unsigned int j = 0; j < rdfatoms.size(); j++)
          fout << " " << setw(12) << s_rdf[rdfatoms.size() + j][i] / s_v_exp;
        fout << endl;

      }
      fout.close();

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}









