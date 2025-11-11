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
 * @file stacking.cc
 * Monitors the occurrence of stacking residues
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor stacking
 * @section stacking monitors the occurrence of stacking residues
 * @author @ref ns
 * @date 11-12-2007
 *
 * Program stacking monitors the occurrence of stacking residues over a molecular 
 * trajectory file through geometric criteria.
 *
 * Ring systems were considered to stack if the distance between
 * the centres of geometry of the rings is less than a given distance (typically
 * 0.5 nm) and the angle between the planes through the two rings is
 * less than a user specified angle (typically 30 degree).
 *
 * The user can specify two groups of atoms (donors and acceptors) between which
 * the stacking interactions are to be monitored. Ring systems are identified by
 * a list given in the library file. This library file has the following format:
 *  @verbatim
TITLE
my stacking library
END
STACKINGRINGS
# residuename atomname1 atomname2 atomname3 ...
PHE CG CE1 CE2 CD1 CD2 CZ
END@endverbatim
 * 
 * In every line of the STACKINGRINGS block, a residue name followed by at least
 * three atom names has to be given. The first three atoms of the list define 
 * the plane through the ring system. The further atoms are only used for the
 * centre of geometry calculation: see @ref utils::StackingProperty "here" for
 * details. For residues containing two or more ring systems one has to specify
 * a line for each ring.
 *
 * The program calculates occurances and prints out a time series file
 * (stackingts.dat) of the observed residue pairs.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@donor</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@acceptor</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> [\@paras</td><td>&lt;distance [nm] and angle [deg]; default: 0.5, 30&gt;] </td></tr>
 * <tr><td> \@library</td><td>&lt;stacking library file&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  stacking
    @topo             ex.top
    @pbc              r
    [@time            0 1]
    @donor            1:a
    @acceptor         2:a
    @paras            0.5 30
    @library          ../data/stacking.lib
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <map>
#include <map>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iterator>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;

/**
 * returns the atom specifiers for a set of atoms and a library
 * @param atoms the atoms which will be sorted and are scaned for rings
 * @param lib the library that maps the residue name
 */
vector<string> get_atom_specifiers(AtomSpecifier & atoms,
        const multimap<string, vector<string> > & lib);

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "donor" << "acceptor" << "paras"
          << "library" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@donor          <atoms>\n";
  usage += "\t@acceptor       <atoms>\n";
  usage += "\t[@paras          <distance [nm] and angle; default: 0.5, 135>]\n";
  usage += "\t@library        <stacking library file\n";
  usage += "\t@traj           <trajectory files>\n";


  try {
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // get the @time argument
    utils::Time time(args);

    // the library defines rings for given residue names
    multimap<string, vector<string> > library;
    Ginstream nf(args["library"]);
    vector<string> buffer;
    nf.getblock(buffer);
    if (buffer[0] != "STACKINGRINGS")
      throw gromos::Exception("stacking",
            "stacking library file does not contain a STACKINGRINGS block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("stacking", "library file file " + nf.name() +
            " is corrupted. No END in STACKINGRINGS"
            " block. Got\n"
            + buffer[buffer.size() - 1]);
    // read in the lib
    for (size_t i = 1; i < buffer.size() - 1; i++) {
      string residue;
      vector<string> atoms;
      istringstream ss(buffer[i]);
      ss >> residue;
      if (ss.fail()) {
        ostringstream msg;
        msg << "bad line in library (residue name): " << buffer[i];
        throw gromos::Exception(argv[0], msg.str());
      }
      for (string atom; ss >> atom;)
        atoms.push_back(atom);

      if (atoms.size() < 3) {
        ostringstream msg;
        msg << "bad line in library: " << buffer[i];
        throw gromos::Exception(argv[0], msg.str());
      }

      library.insert(pair<string, vector<string> >(residue, atoms));
    }

    /* // to debug the library parsing
    multimap<string, vector<string> >::const_iterator iter = library.begin(), to = library.end();
    for(;iter != to; iter++) {
      cerr << iter->first << ": ";
      copy(iter->second.begin(), iter->second.end(), ostream_iterator<string>(cerr));
      cerr << endl;
    } */

    AtomSpecifier donor_atoms = AtomSpecifier(sys);
    {
      Arguments::const_iterator to = args.upper_bound("donor");
      for (Arguments::const_iterator iter = args.lower_bound("donor"); iter != to; iter++)
        donor_atoms.addSpecifier(iter->second);
    }
    vector<string> donor = get_atom_specifiers(donor_atoms, library);
    if (donor.size() == 0)
      throw gromos::Exception("stacking", "No donor-atoms specified!");

    AtomSpecifier acceptor_atoms = AtomSpecifier(sys);
    {
      Arguments::const_iterator to = args.upper_bound("acceptor");
      for (Arguments::const_iterator iter = args.lower_bound("acceptor"); iter != to; iter++)
        acceptor_atoms.addSpecifier(iter->second);
    }
    vector<string> acceptor = get_atom_specifiers(acceptor_atoms, library);
    if (acceptor.size() == 0)
      throw gromos::Exception("stacking", "No acceptor-atoms specified!");

    ofstream timeseries("stackingts.dat");
    if (!timeseries.is_open()) {
      throw gromos::Exception("stacking", "Can't open time series file for writing!");
    }

    // get the parameters
    vector<double> paras = args.getValues<double>("paras", 2, false,
          Arguments::Default<double>() << 0.5 << 30.0);
    double dist_upper = paras[0];
    double angle_upper = paras[1];

    // now we have all the parameters. Let's create the properties
    PropertyContainer props(sys, pbc);
    // loop donor against acceptor
    for (vector<string>::const_iterator d_it = donor.begin(), d_to = donor.end();
            d_it != d_to; ++d_it) {
      for (vector<string>::const_iterator a_it = acceptor.begin(), a_to = acceptor.end();
              a_it != a_to; ++a_it) {
        // check if the residues are not equal
        if (*d_it != *a_it) {
          // create the property specifier
          ostringstream propspec;
          propspec << "st%" << *d_it << "%" << *a_it << "%" << dist_upper
                  << "%" << angle_upper;
          props.addSpecifier(propspec.str());

        } // not eqal
      } // acceptor
    } // donor

    // a vector that contains all times
    std::vector<double> times;

    // loop over all trajectories
    InG96 ic;
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());

      // loop over single trajectory
      while (!ic.eof()) {
        ic >> sys >> time;
        times.push_back(time.time());
        // gather system
        (*pbc.*gathmethod)();

        // calculate the stackings
        props.calc();
      } // while frames
      ic.close();
    } // for traj

    // now we have calculated all properties. Lets sort them and kick out those
    // that are not interesting (i.e. never stack)
    //
    // the sorting we do by a multimap which is sorted by key by default

    multimap<double, StackingProperty *> sorted_props;
    for (PropertyContainer::size_type i = 0; i < props.size(); ++i) {
      StackingProperty * prop = dynamic_cast<StackingProperty *> (props[i]);
      const double occurence = 100.0 * prop->getScalarStat().ave(); // in percent
      sorted_props.insert(pair<double, StackingProperty *>(occurence, prop));
    }

    // write title to the output
    cout << "#" << setw(4) << right << "ID" << " "
            << setw(5) << left << "MOL1" << setw(10) << left << "RES1"
            << setw(5) << left << "MOL2" << setw(10) << left << "RES2"
            << setw(6) << right << "OCCUR%" << endl;
    cout.precision(2);

    // loop over the sorted properties
    multimap<double, StackingProperty *>::const_reverse_iterator
    prop_it = sorted_props.rbegin(), prop_to = sorted_props.rend();
    for (unsigned int id = 1; prop_it != prop_to; ++prop_it, ++id) {
      StackingProperty * prop = prop_it->second;

      // don't print if occurance is zero.
      if (prop_it->first == 0.0)
        continue;

      // loop over the data and write the time series
      vector<double>::const_iterator
      data_it = prop->getScalarStat().data().begin(),
              data_to = prop->getScalarStat().data().end(),
              time_it = times.begin();

      for (; data_it != data_to; ++data_it, ++time_it) {
        // only write if they stack
        if (*data_it != 0.0)
          timeseries << setw(10) << *time_it << setw(5) << id << endl;
      }

      // write statistics to the output
      cout.setf(ios::fixed);
      cout << setw(5) << right << id << " "
              << setw(5) << left << prop->atoms1().mol(0) + 1
              << setw(5) << left << prop->atoms1().resnum(0) + 1
              << setw(5) << left << prop->atoms1().resname(0)
              << setw(5) << left << prop->atoms2().mol(0) + 1
              << setw(5) << left << prop->atoms2().resnum(0) + 1
              << setw(5) << left << prop->atoms2().resname(0)
              << setw(6) << right << prop_it->first << endl;

      // if there are multiple rings for this resname possible let's be a bit
      // more specific and name the atoms of the ring
      if (library.count(prop->atoms1().resname(0)) > 1) {
        cout << "#     ring 1:";
        for (int i = 0; i < prop->atoms1().size(); ++i)
          cout << " " << prop->atoms1().name(i);
        cout << endl;
      }
      if (library.count(prop->atoms2().resname(0)) > 1) {
        cout << "#     ring 2:";
        for (int i = 0; i < prop->atoms2().size(); ++i)
          cout << " " << prop->atoms2().name(i);
        cout << endl;
      }

      // finish and clean
      timeseries.close();
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

vector<string> get_atom_specifiers(AtomSpecifier & atoms,
        const multimap<string, vector<string> > & lib) {
  vector<string> result;

  int current_mol = -1;
  int current_res = -1;

  // make sure the atoms are sorted for the looping over molecules and residues
  atoms.sort();

  for (int i = 0; i < atoms.size(); ++i) {
    const int mol = atoms.mol(i);
    const int res = atoms.resnum(i);

    if (current_mol != mol || current_res != res) {
      // next residue
      current_mol = mol;
      current_res = res;

      const string resname = atoms.resname(i);
      // search the library for the residue name
      multimap<string, vector<string> >::const_iterator
      it = lib.lower_bound(resname),
              to = lib.upper_bound(resname);
      for (; it != to; ++it) {
        ostringstream spec;
        // create a specifier like 1:res(1:CA,CB,CG)
        spec << (mol + 1) << ":res(" << (res + 1) << ":";

        // add the atoms to the spec
        const unsigned int num_atoms = it->second.size();
        for (unsigned int i = 0; i < num_atoms;) {
          spec << it->second[i];

          if (++i != num_atoms) // add comma - but not at last atom
            spec << ",";
        }
        // close the residue bracket
        spec << ")";
        // found a ring - add it
        result.push_back(spec.str());
      } // for library
    } // if new residue
  } // for atoms

  return result;
}
