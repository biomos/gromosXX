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
 * @file m_ener.cc
 * Calculates (non-bonded) interaction energies for specific atoms using 
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor m_ener
 * @section m_ener Calculates (non-bonded) interaction energies for specific 
 * atoms 
 * @author 
 * @date 
 *
 * The program m_ener calculates the (non-bonded) interaction energies over molecular 
 * trajectory files (as the program ener). It reads in the molecular 
 * topology file and the MPERTATOM block in the pertubation topology file. 
 * With that the (non-bonded) interaction energies are calculated for multiple sets
 * of parameters. 
 *
 *  Nonbonded interactions are calculated for all selected atoms with all other
 * atoms in the system. Some atoms can be specified as being soft, indicating 
 * that interactions involving any of these atoms have a specified softness
 * parameter, for all other atoms in the system, the softness parameter
 * @f$\alpha = 0@f$. Vanderwaals interactions between particles i and j are
 * calculated as
 * 
 * @f[ V^{LJ}_{ij}=\left[\frac{C_{12}(i,j)}{ ( r_{ij}^6 + \alpha_{LJ} \lambda ^2 C_{126})}-C_6(i,j)\right] \frac{1}{(r_{ij}^6 + \alpha_{LJ} \lambda ^2 C_{126})} @f]
 *
 * with @f$C_{126} = C_{12}/C_6 @f$ for @f$C_{12}@f$ and @f$C_6@f$ unequal 0,
 * @f$C_{126} = 0@f$ otherwise. @f$C_{12}@f$ and @f$C_6@f$ are the interaction
 * parameters taken from the topology, @f$\lambda@f$ and @f$\alpha_{LJ}@f$ are
 * specified by the user. Similarly, the electrostatic interaction, including
 * reaction field contribution for a homogeneous medium outside the cutoff
 * sphere is calculated as 
 *
 * @f[ V^{CRF}_{ij}=\frac{q_iq_j}{4\pi\epsilon_0}\left[\frac{1}{(r^2_{ij}+\alpha_{CRF}\lambda^2)^{1/2}} - \frac{\frac{1}{2}C_{rf}r_{ij}^2}{(R_{rf}^2+\alpha_{CRF}\lambda^2)^{3/2}} - \frac{(1-\frac{1}{2}C_{rf})}{R_{rf}}\right] @f]
 *
 * where @f$\epsilon_0@f$ is the dielectric permittivity of vacuum and 
 * @f$q_i@f$ and @f$q_j@f$ are the atomic partial charges. @f$R_{rf}@f$ is the
 * reaction field cutoff distance, here assumed to be the same as the
 * interaction cutoff. @f$\alpha_{CRF}@f$ and @f$\lambda@f$ are again user 
 * specified. @f$C_{rf}@f$ is calculated from the reaction field dielectric
 * constant @f$\epsilon_{rf}@f$ and @f$\kappa@f$ (user specified) as
 *
 * @f[ C_{rf} = \frac{ (2 - 2 \epsilon_{rf}) (1 + \kappa R_{rf}) - \epsilon_{rf} (\kappa R_{rf})^2 }{ (1 + 2 \epsilon_{rf}) (1 + \kappa R_{rf}) + \epsilon_{rf} (\kappa R_{rf})^2 } @f]
 * 
 * The program m_ener writes out the total Van-der-Waals interaction energies and the 
 * electrostatic interaction energies seperately as well as the overall non-bonded 
 * interaction energy. 
 *  
 *
 ** <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" for nonbonded interaction&gt; </td></tr>
 * <tr><td> \@props</td><td>&lt;@ref PropertySpecifier "properties" to be calculated&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field contribution&gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa for reaction field contribution&gt; </td></tr>
 * <tr><td> \@soft</td><td>&lt;soft @ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@softpar</td><td>&lt;lam&gt; &lt;a_lj&gt; &lt;a_crf&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@firstatom</td><td>&lt;first pertubed atom&gt; </td></tr>
 *<tr><td> \@pttopo</td><td>&lt;pertubation topology&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ener
    @topo      ex.top
    @pbc       r
    @atoms     1:3-13
    @props     d%1:1,2 a%1:1,2,3 t%1:1,2,4,6 t%1:4,2,5,6
    @time      0 0.2
    @cut       1.4
    @eps       61
    @kap       0.0
    @soft      1:4
    @softpar   0.5 1.51 0.5
    @traj      ex.tr
    @firstatom 1:4 
    @pttopo    ex.mpt
 @endverbatim
 *
 * <hr>
 */






#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/PtTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InPtTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/utils/Energy.h"
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
  knowns << "topo" << "pbc" << "atoms" << "props" << "time" << "cut" << "eps"
          << "kap" << "soft" << "softpar" << "traj" << "firstatom" << "pttopo";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@atoms <atoms for nonbonded interaction>\n";
  usage += "\t@props <properties to be calculated>\n";
  usage += "\t@time <time amd dt>\n";
  usage += "\t@cut <cut-off distance>\n";
  usage += "\t@eps <epsilon for reaction field correction>\n";
  usage += "\t@kap <kappa for reaction field correction>\n";
  usage += "\t@soft <soft atoms>\n";
  usage += "\t@softpar <lam> <a_lj> <a_c>\n";
  usage += "\t@traj  <trajectory files>\n";
  usage += "\t@firstatom <first perturbed atom>\n";
  usage += "\t@pttopo <perturbation topology>\n";



  try {
    Arguments args(argc, argv, knowns, usage);

    //   get simulation time
    Time time(args);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());


    //with which atom do we start?
    int start = 0;
    {
      Arguments::const_iterator iter = args.lower_bound("firstatom");
      if (iter != args.upper_bound("firstatom")) {
        utils::AtomSpecifier at(sys, iter->second.c_str());
        start = at.gromosAtom(0);
      }
    }

    // read in multiple-perturbation-topology-file
    gio::InPtTopology ipt(args["pttopo"]);

    // create perturbation class to contain all perturbation data
    gcore::PtTopology pt = ipt.ptTopo();

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // declare the energy class
    Energy en(sys, gff, *pbc);

    //  set atoms
    AtomSpecifier atoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atoms");
      Arguments::const_iterator to = args.upper_bound("atoms");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atoms.addSpecifier(spec);
      }
    }
    en.setAtoms(atoms);

    // set properties
    PropertyContainer props(sys, pbc);
    {
      Arguments::const_iterator iter = args.lower_bound("props");
      Arguments::const_iterator to = args.upper_bound("props");
      for (; iter != to; iter++) {
        string p = iter->second.c_str();
        props.addSpecifier(p);
      }
    }
    en.setProperties(props);

    // set non-bonded parameters
    //   get cut-off distance
    en.setCutOff(args.getValue<double>("cut", false, 1.4));
    en.setRF(args.getValue<double>("eps", false, 1.0),
            args.getValue<double>("kap", false, 0.0));
    // get soft atom list
    AtomSpecifier soft(sys);
    {
      bool lsoft = false;
      Arguments::const_iterator iter = args.lower_bound("soft");
      Arguments::const_iterator to = args.upper_bound("soft");
      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        soft.addSpecifier(spec);
        lsoft = true;
      }
      std::vector<double> softpar = args.getValues<double>("softpar", 3, lsoft,
              Arguments::Default<double>() << 0.0 << 0.0 << 0.0);
      if (lsoft)
        en.setSoft(soft, softpar[0], softpar[1], softpar[2]);

      bool warn = false;
      for(int i = 0; i < pt.numAtoms(); ++i) {
        if (softpar[1] != pt.alphaLJ(i) || softpar[2] != pt.alphaCRF(i)) {
          warn = true; break;
        }
      }
      if (warn) {
        std::cerr << "Warning: softness parameters from perturbation topology do "
                  << " not match the input parameters. The latter are taken for "
                  << " the energy evalulation." << endl;
      }
    }

    // define input coordinate
    InG96 ic;

    // declare some variables for averaging
    int num_frames = 0;
    vector<double> nb(pt.numPt(), 0.0);

    // open files for output
    vector<ofstream*> fout(pt.numPt());
    cout << "opened " << pt.numPt() << " files" << endl;

    for (int i = 0; i < pt.numPt(); ++i) {
      string name = "en_" + pt.pertName(i) + ".dat";
      cout << "  " << name << endl;

      fout[i]->open(name.c_str());

      *fout[i] << "# Time"
              << "           vanderwaals"
              << "         electrostatic"
              << "            non-bonded"
              << endl;
      fout[i]->precision(10);
      fout[i]->setf(ios::right, ios::adjustfield);
    }

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        // we have to gather with any method to get covalent interactions 
        // and charge-groups connected
        pbc->gathergr();

        // make the pairlists
        en.calcPairlist();

        // loop over the different perturbations
        for (int p = 0; p < pt.numPt(); ++p) {
          pt.apply(sys, p, start);
          // calculate the interactions
          en.calcNb_interactions();

          nb[p] += en.nb();

          *fout[p] << time << setw(22) << en.vdw()
                  << setw(22) << en.el() << setw(22) << en.nb() << endl;
        }
        cout << time << endl;

        num_frames++;
      }
    }
    // print out averages
    if (num_frames > 1) {
      for (int p = 0; p < pt.numPt(); ++p) {
        *fout[p] << "# ave " << setw(36) << nb[p] / num_frames << endl;
        fout[p]->close();
      }
    }
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
