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
 * @file prep_noe.cc
 * Converts X-plor NOE data to GROMOS format
 */

#include <cstdlib>
#include <fstream>


/**
 * @page programs Program Documentation
 *
 * @anchor prep_noe
 * @section prep_noe Converts X-plor NOE data to GROMOS format
 * @author @ref mk @ref co
 * @date 11-8-2006, update 25.11.2017 @ref ms
 *
 * Program prep_noe converts NOE data from an X-plor like format to GROMOS
 * format, determining the proper choice of pseudo- or virtual atoms based on
 * the topology and a library file. The output can be used to apply distance
 * restraints during a simulation using programs promd or md, or to analyse a
 * molecular trajectory using program @ref noe "noe". For a definition of the
 * different types of pseudo- and virtual atoms see volume 2, section 9.4. In cases
 * where the library-file specifies a stereospecific CH2 atom (type 4), but
 * does not indicate which of the two protons is specified, NOE upper bounds
 * are created for both protons. Program @ref post_noe can process the output
 * of an NOE analysis to determine the best assignment.
 *
 * <b>Parse types</b><br>
 * The experimentally determined upper bounds are generally listed in a three
 * column format, with distances in Angstrom. Program prep_noe has three types of
 * parsing these three columns. Specify one of the following numbers for the <b>\@parsetype</b> flag:
 * - 1: take the first value as the upper bound
 * - 2: take the sum of the first and third values as the upper bound (default); or
 * - 3: take the difference between the first and second values (commonly the
 * lower bound).
 *
 *<b>Corrections</b><br>
 * The experimentally determined upper bounds can be corrected for pseudo-atom
 * distances (addition of a geometric constant) or multiplicity factors
 * (typically multiplication with @f$N^{1/p}@f$, where N is the multiplicity of
 * indistinguishable protons involved and p is the averaging power). Such
 * corrections can either be applied to the distances or can be taken out of a
 * set of distances. Specify the corrections file via <b>\@correction</b> flag.
 *
 * <b>NOE specification file example</b><br>
 * This file is using the described format for ambiguous and unambiguous NOEs. Use as input for <b>\@noe</b> flag.
 * @verbatim
TITLE
NOE specification file of 1d3z
END
NOESPEC
1  1 HA   1 HG1 3.567 3.567 0.892 1 2 2           # ambiguous NOE: 2 distances possible
2  1 HA   1 HG2 3.567 3.567 0.892 2 2 1
...
21 1 HA   2 HN  2.133 0.533 0.533                 # unambiguous NOE: Only a single distance
22 1 HB1  2 HN  3.032 3.032 0.758 22 4 23 24 25   # ambiguous NOE: 4 distances possible
23 2 HN  63 HB2 3.032 3.032 0.758 23 4 22 24 25
24 1 HB2  2 HN  3.032 3.032 0.758 24 4 22 23 25
25 2 HN  63 HB2 3.032 3.032 0.758 25 4 22 23 24
END @endverbatim
 *
 * <b>Description of "columns"</b><br>
 * A column refers not to actual columns, but words separated by whitespace.
 * - Column 1: The NOE number. Does not need to be consecutive but must be ascending (i.e. 3, 5, 6, 10 is fine)
 * - Column 2: Resid 1
 * - Column 3: Atom name 1
 * - Column 4: Resid 2
 * - Column 5: Atom name 2
 * - Column 6-8: Distances. Use as specified by <b>\@parsetype</b> flag
 * 
 * If <b>ambiguous NOEs</b> are present (i.e. linked by "OR" in XPLOR or CIF file), they can be marked by <i>additional</i> columns:
 * - Column 9: NOE number, must match number in first column
 * - Column 10: Number (count) of ambiguous NOEs linked to this NOE (N)
 * - Column 11 - 11+N-1: Specify the other NOE numbers that belong to the ambiguous NOE, but exclude the current NOE number.
 * 
 * The current implementation can handle any combination of type 4/0 (non-stereoassigned CH2 group, which creates 2 NOEs to monitor) and
 * other ambiguous NOEs, as specified via column 9.
 *
 * This can be exploited with @ref post_noe to remove all but one of the ambiguous NOE
 * distances by using the <b>\@minmax</b> flag and setting it to <b>min</b>.
 * prep_noe writes a <b>filter file</b>, which can be used to re-evaluate
 * a given analysis over a specific trajectory, without recalculating all
 * distances, through program @ref post_noe.
 *
 * <b>Output files</b><br>
 * - Stdout: NOESPEC file to use as input for @ref noe program via \@noe flag
 * - noe.filter: filter file. Can be used via \@filter flag for @ref post_noe program
 * - noe.dsr: distance restraints file. Can be used as \@disres input for MD program
 *
 *
 * <b>Arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@title</td><td>&lt;NOE title for output&gt; </td></tr>
 * <tr><td> \@noe</td><td>&lt;X-plor like NOE specification file&gt; </td></tr>
 * <tr><td> \@lib</td><td>&lt;NOE specification library&gt; </td></tr>
 * <tr><td> [\@dish</td><td>&lt;carbon-hydrogen distance; default: 0.1 nm&gt;] </td></tr>
 * <tr><td> [\@disc</td><td>&lt;carbon-carbon distance; default: 0.153 nm&gt;] </td></tr>
 * <tr><td> [\@parsetype</td><td>&lt;Upper bound parse type: 1, 2 or 3&gt; ] </td></tr>
 * <tr><td> [\@correction</td><td>&lt;correction file&gt; [&lt;correction type&gt;] ] </td></tr>
 * <tr><td> [\@action</td><td>&lt;add&gt; or &lt;sub&gt; correction from upper bound; default: add ] </td></tr>
 * <tr><td> [\@filter</td><td>&lt;discard NOE's above a certain distance [nm]; default 10000 nm&gt;] </td></tr>
 * <tr><td> [\@factor</td><td>&lt;conversion factor Ang to nm; default is 10&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  prep_noe
    @topo          ex.top
    @title         octa-alanine
    @noe           noe.prep
    @lib           ../data/noelib.45a3
    @dish          0.1
    @disc          0.153
    @parsetype     2
    @correction    ../data/noecor.gromos96
    @action        sub
    @filter        0.8
    @factor        10
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/StringTokenizer.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/VirtualAtoms.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/VirtualAtom.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Noe.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;
using namespace utils;


class Noeprep {
public:
  int residA, vacogsubtypeA, type_A, subtype_A;
  string atomA;
  int residB, vacogsubtypeB, type_B, subtype_B, count;
  string atomB;
  double dis;
  vector<int> links;
  string atname;
  vector<VirtualAtom*> vatomA, vatomB;

/*  Noeprep(){}

  Noeprep(int resA, string A, int resB, string B, double d) {
    residA = resA;
    atomA = A;
    residB = resB;
    atomB = B;
    dis = d;
  }

  ~Noeprep() {
  }
  */
  void set(int resA, string A, int resB, string B, double d){
  	residA = resA;
    atomA = A;
    residB = resB;
    atomB = B;
    dis = d;
  }

  void add_link(int link_number){
  	links.push_back(link_number);
  }

  const vector<int>& get_links() const{
  	return links;
  }
  const vector<VirtualAtom*>& get_vatomA() const{
  	return vatomA;
  }
  const vector<VirtualAtom*>& get_vatomB() const{
  	return vatomB;
  }
};

int main(int argc, char *argv[]) {
  // Usage string

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@title        <NOE title for output>\n";
  usage += "\t@noe          <NOE specification file>\n";
  usage += "\t@lib          <NOE specification library>\n";
  usage += "\t[@dish        <carbon-hydrogen distance; default: 0.1 nm>]\n";
  usage += "\t[@disc        <carbon-carbon distance; default: 0.153 nm>]\n";
  usage += "\t[@parsetype   <Upper bound parse type: 1, 2 or 3> ]\n";
  usage += "\t        Choices are:\n";
  usage += "\t        1: Upper bound == first number\n";
  usage += "\t        2: Upper bound == first + third number (most common, default)\n";
  usage += "\t        3: Upper bound == first - second number (commonly the lower bound)\n";
  usage += "\t[@correction  <correction file> [<correction type>] ]\n";
  usage += "\t[@action      <add> or <sub> correction from upper bound; default: add ]\n";
  usage += "\t[@filter      <discard NOE's above a certain distance [nm]; default 10000 nm>]\n";
  usage += "\t[@factor      <conversion factor Ang to nm; default is 10>]\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "title" << "filter" << "factor" << "noe" << "lib"
          << "parsetype" << "correction" << "dish" << "disc" << "action";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);

  try {

    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    //read in title
    string tit;
    {
      Arguments::const_iterator iter = args.lower_bound("title");
      for (; iter != args.upper_bound("title"); ++iter) {
        tit += iter->second + " ";
      }
    }

    //read in filter threshold
    double filt = args.getValue<double>("filter", false, 10000.0);

    //read in conversion factor
    double conv = args.getValue<double>("factor", false, 10.0);


    //try for disc and dish
    double dish = args.getValue<double>("dish", false, 0.1);
    double disc = args.getValue<double>("disc", false, 0.153);


    // Read in and create the NOE list
    Ginstream nf(args["noe"]);
    vector<string> buffer;
    nf.getblock(buffer);
    if (buffer[0] != "NOESPEC")
      throw gromos::Exception("main",
            "NOESPEC file does not contain an NOESPEC block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("prep_noe", "NOE file " + nf.name() +
            " is corrupted. No END in NOESPEC"
            " block. Got\n"
            + buffer[buffer.size() - 1]);

    // in noe all noes will be stored
    map<int, Noeprep > noevec;
    // get parsetype
    int ptype = args.getValue<int>("parsetype", false, 2);


    // Read in and create the NOE library
    Ginstream nff(args["lib"]);
    vector<string> lib_buffer;
    nff.getblock(lib_buffer);

    vector<Noelib> noelib;
    parse_noelib(lib_buffer, noelib);
    nff.close();

    for (unsigned int j = 1; j < buffer.size() - 1; j++) {
      StringTokenizer tok(buffer[j], " \t");
      vector<string> tokens = tok.tokenize();
      // check if the input file format has at least 8 columns
      if (tokens.size() < 8) {
        throw gromos::Exception("prep_noe", "Too few columns in \""
                + args["noe"] + "\". Did you set the sequential numbers in column 1 "
                "(see manual for further information)?\n");
      }

      // tokens[0] would be the sequential NOE number not used here but nice
      // to have it in the file for comparision with the output of post_noe
      int a = atoi(tokens[1].c_str());
      int b = atoi(tokens[3].c_str());
      unsigned int noe_number = atoi(tokens[0].c_str());
      double d = atof(tokens[5].c_str());
      double e = atof(tokens[6].c_str());
      double f = atof(tokens[7].c_str());

      // throw an exception if the noe has been used already
      if (noevec.count(noe_number)){
	    stringstream out;
        out << "NOE number " << noe_number << " is used twice. Please make sure every NOE number is used only once.";
      	throw gromos::Exception("prep_noe", out.str());
      }
      // the noe numbers must be consecutive due to downstream programs
      if (!noevec.empty()) {
          unsigned int last_noe_num = noevec.rbegin()->first;
          if (noe_number <= last_noe_num){ // check if the last element of this ordered map is bigger than the current noe number
            stringstream msg;
            msg <<"NOE number " << noe_number << " is smaller than its predecessor " << last_noe_num << ". This is not allowed";
            throw gromos::Exception("prep_noe", msg.str());
          }
      }

      if (tokens.size() > 8) {
        unsigned int g = atoi(tokens[8].c_str());
        if (g != noe_number)
          throw gromos::Exception("prep_noe",
                "Numbering in NOESPEC file (9th column) is not correct. Must match NOE number in 1st column.");
        int h = atoi(tokens[9].c_str());
        for (int ii = 0; ii < h - 1; ii++)
          noevec[noe_number].add_link(atoi(tokens[10 + ii].c_str()));
      }
      // apply parse type
      switch (ptype) {
        case 1: d = d;
          break;
        case 2: d = d + f;
          break;
        case 3: d = d - e;
          break;
        default:
          throw gromos::Exception("prep_noe", args["parsetype"] +
                  " unknown. Known types are 1, 2 and 3");
      }

      noevec[noe_number].set(a, tokens[2], b, tokens[4], d);

      // get the virtual atom type
      string resnameA, resnameB;
	  int molA = 0, molB = 0, resnumA = 0, resnumB = 0, mol = 0, atA = 0, atB = 0, atNumA = 0, atNumB = 0, count = 0;

      //Noeprep NOE = noevec[noe_number];


      //first find the residue-name corresponding to
      //your input-number in the topology
      mol = 0;
      atNumA = 0;
      while (noevec[noe_number].residA > (atNumA += sys.mol(mol).topology().numRes())) {
        ++mol;
        if (mol > sys.numMolecules())
          throw gromos::Exception("prep_noe",
                "Residue number too high in input line:\n");
      }
      atA = noevec[noe_number].residA;
      atA -= atNumA - sys.mol(mol).topology().numRes();
      resnumA = (atA - 1);
      molA = mol;
      resnameA = (sys.mol(mol).topology().resName(atA - 1));

      mol = 0;
      atNumB = 0;
      while (noevec[noe_number].residB > (atNumB += sys.mol(mol).topology().numRes())) {
        ++mol;
        if (mol > sys.numMolecules())
          throw gromos::Exception("prep_noe", +"Residue number too high in input line:\n");
      }
      atB = noevec[noe_number].residB;
      atB -= atNumB - sys.mol(mol).topology().numRes();
      resnumB = (atB - 1);
      molB = mol;
      resnameB = (sys.mol(mol).topology().resName(atB - 1));

      //then map the correct gromos-topology atomname
      //to your input-atomname based on the residue-name
      //and the input-atomname
      bool foundA = false;

      int p = 0;
      for (int k = 0; k< int (noelib.size()); ++k) {
        Noelib NOELIB = noelib[k];
        if (NOELIB.resname == resnameA && NOELIB.orgatomname == noevec[noe_number].atomA) {
          //back to topology to get the atom number
          for (int f = 0; f < sys.mol(molA).numAtoms() && foundA == false; ++f) {
            if (sys.mol(molA).topology().atom(f).name() == NOELIB.gratomname &&
                    sys.mol(molA).topology().resNum(f) == resnumA) {
              int addA = 0;
              foundA = true;
              p = k;

              for (int i = 0; i < molA; ++i) addA += sys.mol(i).numAtoms();

              if (noevec[noe_number].dis / conv > filt) cout << "#";
              noevec[noe_number].vatomA = getvirtual(f + addA, NOELIB.NOETYPE, NOELIB.NOESUBTYPE, sys,
                      dish, disc);
              // store the type and subtype. needed for filter file
              noevec[noe_number].type_A = NOELIB.NOETYPE;
              noevec[noe_number].subtype_A = NOELIB.NOESUBTYPE;
              // remember the subtype if found VA type is COM (-1)
              if (NOELIB.NOETYPE == -1) {
                noevec[noe_number].vacogsubtypeA  = NOELIB.NOESUBTYPE;
              } else {
                noevec[noe_number].vacogsubtypeA = 0;
              }
            }
          }
        }
      }

      if (!foundA) {
        std::stringstream ss;
        string a;
        ss << noevec[noe_number].residA;
        ss >> a;
        string b = noevec[noe_number].atomA;
        string c = " ";
        string d = a + c + b;
        throw gromos::Exception("prep_noe ", d +
                " Noe specification not found in library!");
      }
      Noelib NA = noelib[p];
      ostringstream atname;

      bool foundB = false;
      for (int z = 0; z< int (noelib.size()); ++z) {
        Noelib NOELIBB = noelib[z];
        if (NOELIBB.resname == resnameB && NOELIBB.orgatomname == noevec[noe_number].atomB) {
          //back to topology to get the atom number
          for (int g = 0; g < sys.mol(molB).numAtoms() && foundB == false; ++g) {
            if (sys.mol(molB).topology().atom(g).name() == NOELIBB.gratomname &&
                    sys.mol(molB).topology().resNum(g) == resnumB) {
              int addB = 0;
              foundB = true;
              count += 1;
              for (int i = 0; i < molB; ++i) addB += sys.mol(i).numAtoms();
              noevec[noe_number].vatomB = getvirtual(g + addB, NOELIBB.NOETYPE, NOELIBB.NOESUBTYPE, sys,
                      dish, disc);
              noevec[noe_number].type_B = NOELIBB.NOETYPE;
              noevec[noe_number].subtype_B = NOELIBB.NOESUBTYPE;
              // remember the subtype if found VA type is COM (-1)
              if (NOELIBB.NOETYPE == -1) {
                noevec[noe_number].vacogsubtypeB = NOELIBB.NOESUBTYPE;
              } else {
                noevec[noe_number].vacogsubtypeB = 0;
                }

              atname << setw(3) << molA + 1 << " "
                      << setw(5) << resnumA + 1 << " "
                      << setw(4) << NA.resname << " "
                      << setw(4) << NA.gratomname << " "
                      << setw(4) << NA.orgatomname << " "
                      << setw(3) << molB + 1 << " "
                      << setw(5) << (resnumB + 1) << " "
                      << setw(4) << NOELIBB.resname << " "
                      << setw(4) << NOELIBB.gratomname << " "
                      << setw(4) << NOELIBB.orgatomname << " ";

            }
          }
        }
      }
      noevec[noe_number].atname = atname.str();
      noevec[noe_number].count = count;

      if (!foundB) {
        std::stringstream s;
        string aa;
        s << noevec[noe_number].residB;
        s >> aa;
        string bb = noevec[noe_number].atomB;
        string cc = " ";
        string dd = aa + cc + bb;
        throw gromos::Exception("prep_noe ", dd +
                " Noe specification not found in library!");
      }

    }
    nf.close();

    //check whether to add or subtract correction
    bool add = true;
    bool sub = false;
    if (args.count("action") > 0) {
      if (args["action"] == "add" || args["action"] == "ADD") {
        add = true;
        sub = false;
      } else if (args["action"] == "sub" || args["action"] == "SUB") {
        add = false;
        sub = true;
      } else
        throw gromos::Exception("prep_noe",
              "action type " + args["action"] + " not known!");
    }

    //read in the correction file if it exists
    map<int, map<int, double> > pseudocorrectiondata; // first key: type
    map<int, map<int, double> > multiplicitydata; // second key: subtype
    bool pseudocorrection = false;
    bool multiplicitycorrection = false;
    try {
      args.check("correction");
      Ginstream corf(args["correction"]);
      //get NOECORGROMOS block
      buffer.clear();
      corf.getblock(buffer);

      if (buffer[0] != "NOECORGROMOS")
        throw gromos::Exception("main",
              "NOE correction file does not contain the "
              "NOECORGROMOS block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("prep_noe", "Correction file " + corf.name() +
              " is corrupted. No END in NOECORGROMOS"
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 1; j < buffer.size() - 1; j++) {
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
        pseudocorrectiondata[atoi(tokens[0].c_str())][atoi(tokens[1].c_str())]
                = atof(tokens[2].c_str());
      }


      //get MULTIPLICITY block
      buffer.clear();
      corf.getblock(buffer);

      if (buffer[0] != "MULTIPLICITY")
        throw gromos::Exception("main",
              "NOE correction file does not contain the"
              " MULTIPLICITY block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("prep_noe", "Correction file " + corf.name() +
              " is corrupted. No END in MULTIPLICITY"
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 1; j < buffer.size() - 1; j++) {
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
        multiplicitydata[atoi(tokens[0].c_str())][atoi(tokens[1].c_str())] = atof(tokens[2].c_str());
      }
      corf.close();

      //determine the correction type
      Arguments::const_iterator it = args.lower_bound("correction");
      ++it;

      string ctype;

      if (it == args.upper_bound("correction")) {
        pseudocorrection = true;
        multiplicitycorrection = true;
      } else ctype = it->second;

      if (ctype == "pseudo") pseudocorrection = true;
      else if (ctype == "multiplicity") multiplicitycorrection = true;
      else if (ctype != "")
        throw gromos::Exception("prep_noe",
              "Correction type " + ctype + " not known!" +
              "Use 'pseudo' or 'multiplicity' to apply " +
              "only one type of correction");

    }    catch (Arguments::Exception e) {
      cout << "# No correction file used!" << endl;
    }


    //cout the title and that kind of stuff
    cout << "TITLE" << endl;
    cout << "NOE specification file for: " << tit << endl;
    cout << "END" << endl;
    cout << "NOECALCSPEC" << endl;
    cout << "# DISH: carbon-hydrogen distance" << endl;
    cout << "# DISC: carbon-carbon distance" << endl;
    cout << "#" << setw(9) << "DISH" << setw(10) << "DISC" << endl;
    cout.precision(5);
    cout << setw(10) << dish << setw(10) << disc << endl;
    cout << "# IDR1, JDR1, KDR1, LDR1:\n";
    cout << "#       atom sequence numbers of the real atoms defining the geometric position\n";
    cout << "#       of the first atom of a distance restraint pair\n";
    cout << "# IDR2, JDR2, KDR2, LDR2:\n";
    cout << "#       atom sequence numbers of the real atoms defining the geometric position\n";
    cout << "#       of the second atom of a distance restraint pair\n";
    cout << "# ICDR: geometric code defining the positions of the atoms of a distance\n";
    cout << "#       restraint pair\n";
    cout << "# VACS: subtypes of virual atoms of type COG (-1). Possible subtypes are:\n";
    cout << "#         0: no subtype defined\n";
    cout << "#         1: aromatic flipping ring\n";
    cout << "#         2: non-stereospecific NH2 group\n";
    cout << "# R0:   upper bound\n";
    cout << "#" << setw(4) << "IDR1" << setw(5) << "JDR1"
            << setw(5) << "KDR1" << setw(5) << "LDR1"
            << setw(5) << "ICDR" << setw(5) << "VACS"
            << setw(5) << "IDR2" << setw(5) << "JDR2"
            << setw(5) << "KDR2" << setw(5) << "LDR2"
            << setw(5) << "ICDR" << setw(5) << "VACS"
            << setw(10) << "R0" << endl << "#" << endl;

    //open the filter file...
    ofstream filterfile;
    filterfile.open("noe.filter");
    filterfile << "TITLE" << endl;
    filterfile << "NOE filter file for: " << tit << endl;
    filterfile << "END" << endl;
    filterfile << "NOEFILTER" << endl;
    filterfile << "# noe"
            << setw(4) << "mol"
            << setw(11) << "residue"
            << setw(5) << "atom"
            << setw(5) << "atom"
            << setw(4) << "mol"
            << setw(11) << "residue"
            << setw(5) << "atom"
            << setw(5) << "atom"
            << setw(6) << "r0"
            << " filter noe" << endl;

    ofstream disresfile;
    disresfile.open("noe.dsr");
    disresfile << "TITLE" << endl;
    disresfile << "NOE distance restraints file for: " << tit << endl;
    disresfile << "END" << endl;
    disresfile << "DISTANCERESSPEC" << endl;
    disresfile << "#" << setw(9) << "DISH" << setw(10) << "DISC" << endl;
    disresfile.precision(5);
    disresfile << setw(10) << dish << setw(10) << disc << endl;
    disresfile << "#" << setw(4) << "IDR1" << setw(5) << "JDR1"
            << setw(5) << "KDR1" << setw(5) << "LDR1"
            << setw(5) << "ICDR" << setw(5) << "IDR2"
            << setw(5) << "JDR2" << setw(5) << "KDR2" << setw(5) << "LDR2"
            << setw(5) << "ICDR" << setw(10) << "R0" << setw(10) << "W0"
            << setw(10) << "NRAH" << endl << "#" << endl;

    //here goes the crappy code
    double filterbound = 0;
    int totalnoecount = 1;

	for(map<int, Noeprep>::iterator noe_it = noevec.begin(); noe_it != noevec.end(); ++noe_it) {
      int noe_num = noe_it->first;  // get the key (=noe number)
      Noeprep NOE = noe_it->second;  // get the map value (=Noeprep)

      vector<int> links = NOE.get_links();
      string atname = NOE.atname;
      vector<VirtualAtom*> vatomA = NOE.get_vatomA();
	  vector<VirtualAtom*> vatomB = NOE.get_vatomB();

      //spit out disresblock...
      int atomsA[4], atomsB[4];
      int creatednoe = 0;

      // print out "readable" constraints as a coment
      disresfile << "# " << NOE.count << atname << endl;

      for (int va = 0; va < (int) vatomA.size(); ++va) {
        int offsetA = 1;
        VirtualAtom VA(*vatomA[va]);

        int mol = VA.conf().mol(0);
        for (int l = 0; l < mol; ++l) offsetA += sys.mol(l).numAtoms();

        for (unsigned int aa = 0; aa < 4; ++aa) {
          int att;
          if (VA.conf().size() > aa)
            att = VA.conf().atom(aa);
          else
            att = -1;

          atomsA[aa] = att;
        }

        for (int vb = 0; vb < (int) vatomB.size(); ++vb) {
          int offsetB = 1;
          VirtualAtom VB(*vatomB[vb]);
          int mol = VB.conf().mol(0);
          for (int l = 0; l < mol; ++l) offsetB += sys.mol(l).numAtoms();

          for (unsigned int bb = 0; bb < 4; ++bb) {
            int att;
            if (VB.conf().size() > bb)
              att = VB.conf().atom(bb);
            else
              att = -1;


            atomsB[bb] = att;
          }

          ostringstream ss, noeline, disresline;
          ss.setf(ios::right, ios::adjustfield);
          ss.setf(ios::fixed, ios::floatfield);
          ss.precision(5);

          for (int kk = 0; kk < 4; ++kk) {
            if (atomsA[kk] == -1) ss << setw(5) << 0;
            else ss << setw(5) << atomsA[kk] + offsetA;
          }
          ss << setw(5) << VA.type();
          disresline << ss.str();
          ss << setw(5) << NOE.vacogsubtypeA;
          noeline << ss.str();
          ss.str(""); // clear it

          for (int kk = 0; kk < 4; ++kk) {
            if (atomsB[kk] == -1) ss << setw(5) << 0;
            else ss << setw(5) << atomsB[kk] + offsetB;
          }
          ss << setw(5) << VB.type();
          disresline << ss.str();
          ss << setw(5) << NOE.vacogsubtypeB;
          noeline << ss.str();
          ss.str(""); // clear it

          double bound = NOE.dis / conv;
          //set up corrections
          vector<int> type;
          vector<int> stype;
          type.push_back(VA.type());
          type.push_back(VB.type());
          stype.push_back(NOE.vacogsubtypeA);
          stype.push_back(NOE.vacogsubtypeB);

          // first do the multiplicity correction and then the
          // pseudo atom correction

          double mult = 1;
          double cor = 0;

          //check for multiplicity-correction
          if (multiplicitycorrection) {
            for (int i = 0; i < (int) type.size(); ++i) {
              if (multiplicitydata.find(type[i])
                      != multiplicitydata.end()) {
                if (multiplicitydata.find(type[i])->second.find(stype[i])
                        != multiplicitydata.find(type[i])->second.end()) {
                  mult *= multiplicitydata.find(type[i])->second.find(stype[i])->second;
                }
              }
            }
          } //end if (multiplicitycorrection)


          //check for gromos-correction
          if (pseudocorrection) {
            for (int i = 0; i < (int) type.size(); ++i) {
              if (pseudocorrectiondata.find(type[i])
                      != pseudocorrectiondata.end()) {
                if (pseudocorrectiondata.find(type[i])->second.find(stype[i])
                        != pseudocorrectiondata.find(type[i])->second.end()) {
                  cor += pseudocorrectiondata.find(type[i])->second.find(stype[i])->second;
                }
              }
            }
          } //end if (pseudocorrection)
          if (add) bound = bound * mult + cor;
          else if (sub) bound = (bound - cor) / mult;

          // in the filter file I also want the corrected bound
          filterbound = bound;

          disresline << ss.str();
          noeline << ss.str();

          cout << noeline.str() << setw(10) << bound << endl;
          disresfile << disresline.str() << setw(10) << bound
                  << setw(10) << 1.0 << setw(10) << 1 << endl; // half-harmonic attractive
          if (va > 0 || vb > 0) cout << "#automatically generated NOE distance to monitor!" << endl;

          ++creatednoe;
        }
      }

      //smack out the filterfile...
      for (int noe_index = 0; noe_index < creatednoe; ++noe_index) {
        filterfile << setw(5) << totalnoecount << " ";
        filterfile << atname;
        filterfile << setw(5) << filterbound << " ";

        map<int, int> type_4_counts;
		// for type 4 (unassigned methylene group) we create 2 NOEs, which are linked,
		// without expliticlty specifying the link in the @noe input file
		// check if other NOEs exist in this group that also have interally created NOEs:
		for (unsigned int links_index = 0; links_index < links.size(); ++links_index){
		    // check for type 4 0 for both atoms
		    int other_noe = links[links_index];
		    if (!noevec.count(other_noe)){
                stringstream msg;
                msg << "NOE link " << other_noe << " specified for NOE number " << noe_num << " does not exist";
                throw gromos::Exception("prep_noe", msg.str());
		    }
            if (noevec[other_noe].type_A == VirtualAtom::CH2  && noevec[other_noe].subtype_A == 0){
                if (type_4_counts.count(other_noe)) type_4_counts[other_noe] += 1;
                else type_4_counts[other_noe] = 1;
             }
            if (noevec[other_noe].type_B == VirtualAtom::CH2  && noevec[other_noe].subtype_B == 0){
                if (type_4_counts.count(other_noe)) type_4_counts[other_noe] += 1;
                else type_4_counts[other_noe] = 1;
            }
		}
		int type_4_count = 0; // only catches type 4s that are not of the current noe
		for(map<int,int>::iterator map_it=type_4_counts.begin(); map_it != type_4_counts.end(); ++map_it)
            type_4_count += (map_it->second)*2-1;

		int number_of_links = links.size() + 1 + type_4_count + creatednoe-1;  // +creatednoe-1 for the current noe
		filterfile << " " << number_of_links << " ";

        // get all links
        // first check if any links came before:
        vector<int> pre_noes;
        vector<int> post_noes;
        for(vector<int>::iterator lnk=links.begin(); lnk!=links.end(); ++lnk){
            if(*lnk < noe_num) pre_noes.push_back(*lnk);
            else if(*lnk > noe_num) post_noes.push_back(*lnk);
        }

        // create all link numbers based on previous and future links
        int neg_offset = 0;
        for(vector<int>::iterator pre_noe=pre_noes.begin(); pre_noe!=pre_noes.end(); ++pre_noe){
                if(type_4_counts.count(*pre_noe)) neg_offset += type_4_counts[*pre_noe]*2;
                else ++neg_offset;
        }
        // spit out past and current links
        for(int of=totalnoecount-neg_offset-noe_index; of<totalnoecount+creatednoe-noe_index; ++of)
            if (of != totalnoecount)
                filterfile << " " << of;
        // spit out future links
        int link_end = 0;
        for(vector<int>::iterator post_noe=post_noes.begin(); post_noe!=post_noes.end(); ++post_noe){
                if(type_4_counts.count(*post_noe)) link_end += type_4_counts[*post_noe]*2;
                else ++link_end;
        }
        for(int le=totalnoecount+creatednoe-noe_index; le<totalnoecount+creatednoe+link_end-noe_index; ++le)
                filterfile << " " << le;

        filterfile << endl;
        ++totalnoecount;
      }
    } //end looping over all lines in @noe input

    cout << "END" << endl;
    filterfile << "END" << endl;
    disresfile << "END" << endl;
    filterfile.close();
  } catch (gromos::Exception e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


