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
 * @file swd.cc
 * Calculates the solute averaged distance
 */

/**
 * @page programs Program Documentation
 *
 * @anchor swd
 * @section swd Calculates time serie of solute averaged distances
 * @author @ref nb
 * @date 22. 11. 2004
 *
 * Program swd calculates the solute averaged distance for fg and cg solvent.
 * If the force constant and cutoff radius is supplied, the energy is
 * calculated as well.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@solute</td><td>&lt;@ref AtomSpecifier "solute atoms"&gt; </td></tr>
 * <tr><td> \@fgsolv</td><td> &lt;@ref AtomSpecifier "fg solvent atoms"&gt;</td></tr>
 * <tr><td> \@cgsolv</td><td> &lt;@ref AtomSpecifier "cg solvent atoms"&gt;</td></tr>
 * <tr><td> [\@energy</td><td>&lt;fg force constant&gt; &lt;fg force cut off&gt;  
 *                            &lt;cg force constant&gt; &lt;cg force cut off&gt;
 *                            &lt;cg energy file&gt;  ]</td></tr>
 * <tr><td> [\@exponent</td><td>&lt;the exponent for calculating the swd&gt;] </td></tr>
 * <tr><td> [\@measure</td><td>&lt;output file for the measure of the splitting&gt;] </td></tr>
 * <tr><td> [\@weights</td><td>&lt;output file for all the weights&gt;] </td></tr>
 * <tr><td> [\@mindist</td><td>&lt;output file for all the minimal distances&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 * 
 * <hr>
 */
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/groTime.h"
#include "../src/utils/SoluteWeightedDistance.h"
#include "../src/gromos/Exception.h"

int main(int argc, char **argv) {

  args::Argument_List knowns;
  knowns << "topo" << "pbc" << "time"  << "traj" << "energy"
          << "solute" << "fgsolv" << "cgsolv" << "exponent" 
          << "measure" << "weights" << "mindist";

  std::string usage = "# " + std::string(argv[0]) + "\n";
  usage += "\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gathermethod>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@solute     <solute atoms>\n";
  usage += "\t@fgsolv     <fg solvent atoms>\n";
  usage += "\t@cgsolv     <cg solvent atoms>\n";
  usage += "\t[@exponent   <the exponent for calculating the swd>]\n";
  usage += "\t[@energy     <fg force constant> <fg cutoff> <cg force constant> <cg cutoff> <energy file>]\n";
  usage += "\t[@measure    <output file for the measure of the splitting>]\n";
  usage += "\t[@weights    <output file for all the weights>]\n";
  usage += "\t[@mindist    <output file for all the minmal distances>]\n";
  usage += "\t@traj       <trajectory files>\n";

  std::ofstream eneout;
  int returnCode = 0;

  try {
    args::Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args);

    //  read topology
    args.check("topo", 1);
    gio::InTopology it(args["topo"]);

    gcore::System sys(it.system());

    utils::SoluteWeightedDistance swd(sys, args); // :-)

    // define input coordinate
    gio::InG96 ic;

    std::cout << swd.title();
    if (swd.withEnergy()) {
      std::string eneOutFile = swd.energyFile();
      eneout.open(eneOutFile.c_str());
      eneout << swd.title();
      std::cerr << "# Energy parameter given. Will print energies to "
              << eneOutFile << std::endl;
    } else {
      std::cerr << "# No energy parameters given, but don't panic!\n";
    }

    // loop over all trajectories
    for (args::Arguments::const_iterator iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        if (ic.stride_eof())
          break;

        swd.calculate(time.time());

        std::cout << time << " \t\t";
        std::cout << swd;
        std::cout << std::endl;
        
        if (swd.withEnergy()){
          eneout << time << " \t\t";
          swd.energies(eneout);
          eneout << std::endl;
        }
      }


      ic.close();
    }

  } catch (const gromos::Exception &e) {
    std::cerr << e.what() << std::endl;
    returnCode = 1;
  }
  if (eneout.is_open()) {
    eneout.close();
  }
  return returnCode;
}
