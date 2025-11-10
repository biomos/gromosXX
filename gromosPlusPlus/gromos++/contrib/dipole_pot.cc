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
 * @file dipole_pot.cc
 * calculate the potential of dipoles on a ion
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor dipole_pot
 * @section dipole_pot calculate the potential of the water dipoles generted
 *                     at the 'centre'
 * @author @ref nb
 * @date 5.9.2013
 *
 * Similar to the @ref iwdcf program
 * 
 * WARNING: This program works only for water, nothing else!
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@centre</td><td>&lt;@ref AtomSpecifier "atoms" to take as centre&gt; </td></tr>
 * <tr><td> \@nsm</td><td>&lt;number of solvent molecules; </td></tr>
 * <tr><td> \@cut</td><td>&lt;maximum distance&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;inital time&gt; &lt;time step&gt;] </td></tr>
 * <tr><td> [\@weighted</td><td>should the numbers be normed?] </td></tr>
 * <tr><td> [\@factor</td><td>&lt;multiply by factor, e.g. 1/(4 pi eps)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 * 
 * The following is calculated
 * 
 * @f[ \frac{f}{N} \cdot \sum_i \frac{\vec{\mu}_i \cdot \vec{r}_i}{\mu_i r_i} \cdot \frac{1}{r_i^2} @f]
 *
 * Here @f$N@f$ is the normalization (either 1 or @f$\sum_i 1/r_i^2 @f$).
 * @f$\vec{\mu}_i@f$ is the dipole moment of water @f$i@f$ and @f$\vec{r}_i@f$
 * the distance vector to the center.
 * @f$f@f$ is the arbitrarly chosen factor (default 1).
 *
 * Example:
 * @verbatim
  rdf
    @topo   ex.top
    @pbc    r
    @centre atom 1:1
    @nsm    7000
    @cut    2.0
    @factor 0.3
    @weighted
    @traj   ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cstdlib>
#include <string>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/bound/Boundary.h"
#include "../src/gromos/Exception.h"

int main(int argc, char** argv) {
  
  args::Argument_List knowns;
  knowns << "topo" << "pbc" << "centre" << "nsm" << "cut" << "traj" 
         << "weighted" << "factor" << "time";
  
  
  std::string usage = "# " + std::string(argv[0]);
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@centre <type> or\n";
  usage += "\t        <atoms> or\n";
  usage += "\t        <cog> or\n";
  usage += "\t        all\n";
  usage += "\t@nsm    <number of solvent molecules>\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t[@time   <inital time> <time step>]\n";
  usage += "\t[@factor <Multiply by factor, e.g. q / (4 pi eps)>]\n";
  usage += "\t[@weighted   normalize to +/- 1]\n";
  usage += "\t@traj   <trajectory files>\n";
  
  
  try {
    args::Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    gio::InTopology it(args["topo"]);
    gcore::System sys(it.system());
    
    // parse boundary conditions
    bound::Boundary *pbc = args::BoundaryParser::boundary(sys, args);

    // read in number of solvent molecules
    int nsm = args.getValue<int>("nsm");
    
    // the cutoff
    double cutoff = args.getValue<double>("cut");
    
    double factor = args.getValue("factor", false, 1.0);
    bool weighted = args.count("weighted") != -1;
    
    // set centre atoms
    int sol_c = 0;

    utils::AtomSpecifier centre(sys);

    {
      args::Arguments::const_iterator iter = args.lower_bound("centre");
      args::Arguments::const_iterator to = args.upper_bound("centre");
      int error = 1;

      if (iter != to) {
        std::string s = iter->second.c_str();
        iter++;
        if (s == "type") {
          error = 0;
          for (; iter != to; iter++) {
            std::string name = iter->second.c_str();
            for (int i = 0; i < sys.numMolecules(); i++)
              for (int j = 0; j < sys.mol(i).topology().numAtoms(); j++)
                if (name == sys.mol(i).topology().atom(j).name())
                  centre.addAtom(i, j);
            for (int i = 0; i < sys.sol(0).topology().numAtoms(); i++)
              if (name == sys.sol(0).topology().atom(i).name())
                for (int j = 0; j < nsm; j++) {
                  int off = j * sys.sol(0).topology().numAtoms();
                  centre.addAtom(-1, i + off);
                  sol_c++;
                }

          }

        }
        if (s == "atom") {
          error = 0;
          for (; iter != to; iter++) {
            std::string spec = iter->second.c_str();
            centre.addSpecifier(spec);
          }

        }
        if (s == "cog") {
          error = 0;
          utils::AtomSpecifier temp(sys);
          centre.addAtom(-2, 0);


          for (; iter != to; iter++) {
            std::string spec = iter->second.c_str();
            centre.addSpecifier(spec);
          }
        }
        if (s == "all") {
          error = 0;
          for (int i = 0; i < sys.numMolecules(); i++)
            for (int j = 0; j < sys.mol(i).numAtoms(); j++)
              centre.addAtom(i, j);
        }
        if (error == 1 || centre.size() == 0)
          throw gromos::Exception("iwdcf @centre", s +
                " unknown or no atoms specified. Give 'type', 'atom', 'cog' or 'all'");
      }
    }
    
    if (sys.sol(0).topology().atom(0).name() != "OW" ||
        sys.sol(0).topology().atom(1).name() != "HW1" ||
        sys.sol(0).topology().atom(2).name() != "HW2") {
      std::cerr << "Solvent is " << sys.sol(0).topology().solvName() << std::endl;
      for (int i = 0; i < 3; i++){
        std::cerr << "Atom " << i << " is " << sys.sol(0).topology().atom(0).name()
                << "; should be OW, HW1 or HW2\n";
      }
      throw gromos::Exception("dipole_pot",
              " works only for water (H2O or H2OE) as solvent...Abort!");
    }
    
    // define input coordinate
    gio::InG96 ic;
    for (args::Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      utils::Time time(args);

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;

        if (nsm > sys.sol(0).numPos() / sys.sol(0).topology().numAtoms())
          throw gromos::Exception("iwdcf",
                " nsm specified is more than in coordinate file");
        else
          sys.sol(0).setNumPos(nsm * sys.sol(0).topology().numAtoms());

        Vec cog(0.0, 0.0, 0.0);
        for (int i = 1; i < centre.size(); i++){
          Vec ni = pbc->nearestImage(centre.pos(0), centre.pos(i), sys.box());
          cog += ni - centre.pos(0);
        }
        cog /= (centre.size());
        cog += centre.pos(0);
        
        double potential = 0.0;
        double normalizeBy = weighted ? 0.0 : 1.0;

        //loop over the atoms to consider
        for (int j = 0; j < nsm; j++) {
          //calculate distance only if this atom is not the current centre
          Vec tmp;
          
          int solv_j = j * sys.sol(0).topology().numAtoms();

          tmp = pbc->nearestImage(cog,
                  sys.sol(0).pos(solv_j),
                  sys.box());
          Vec Ionwater(tmp - cog);

          double dis = Ionwater.abs();
          
          if (dis > cutoff)
            continue;

          Vec H1;
          Vec H2;

          H1 = pbc->nearestImage(tmp,
                  sys.sol(0).pos(solv_j + 1),
                  sys.box());

          H2 = pbc->nearestImage(tmp,
                  sys.sol(0).pos(solv_j + 2),
                  sys.box());

          Vec t1(H1 - tmp);
          Vec t2(H2 - tmp);

          Vec e(t1 + t2);
          e /= e.abs();
          
          double dis2 = dis * dis;
          if (weighted)
            normalizeBy += 1.0 / dis2;
          
          potential += Ionwater.dot(e) / (dis * dis2);

        }
        
        std:: cout << time << " \t" << factor * potential / normalizeBy << std::endl;
      }
      ic.close();
    } //end loop over trajectory
  }  catch (const gromos::Exception &e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  return 0;
}
