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
 * @file make_rdc_spec.cc
 * converts RDCs listed per residue into GROMOS input format
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor make_rdc_spec
 * @section make_rdc_spec converts RDCs listed per residue into GROMOS input format
 * @author gp, jra, lnw
 * @date 16.12.2014, improved mid-2015
 *
 * how to use
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@mol</td><td>&lt;molecule that the rdcs are for (default: 1)&gt;] </td></tr>
 * <tr><td> \@rdc</td><td>&lt;RDC data file. Be aware that numbers are with respect to the first residue of the current molecule. This is not always what you get from the NMR people.&gt; </td></tr>
 * <tr><td> [\@weights</td><td>&lt;rdc-specific weight factors&gt;] </td></tr>
 * <tr><td> [\@nmf</td><td>&lt;number of magnetic field vectors (default: 3)&gt;] </td></tr>
 * <tr><td> [\@mass</td><td>&lt;mass of magnetic field representation particles (default: 10)&gt;] </td></tr>
 * <tr><td> \@type</td><td>&lt;type of rdc (as in library file)&gt; </td></tr>
 * <tr><td> \@lib</td><td>&lt;rdc library file&gt; </td></tr>
 * <tr><td> \@group</td><td>&lt;rdc groups file, or 'ONE'&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   atominfo
     @topo         topo.top
     @mol          2
     @rdc          NH.rdc
     @type         1
     @lib          rdc.lib
     @group        ONE
   @endverbatim

 * <hr>
 *
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gromos/Exception.h"

#if (__cplusplus > 199711L) // we have c++11 or newer
#include <chrono>
#include <ctime>
#endif

using namespace std;

struct rdc_struct {
  int type;
  string name_i;
  string name_j;
  string name_k;
  string name_l;
  double gyr_i;
  double gyr_j;
  double rij;
  double rik;
};

struct rdc_data_struct {
  int residue;
  int num_i;
  int num_j;
  int num_k;
  int num_l;
  double rdc;
  double weight;
};

int main(int argc, char **argv) {
  // print the time and command for future reference
#if (__cplusplus > 199711L) // we have c++11 or newer
  std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::cout << "time: " << std::ctime(&now_time) ;
#endif
  cout <<  "command: ";
  for (int i=0; i<argc; i++) cout << argv[i] << " ";
  cout << endl;

  args::Argument_List knowns;
  knowns << "topo"
         << "mol"
         << "rdc"
         << "weights"
         << "nmf"
         << "mass"
         << "type"
         << "lib"
         << "group";

  string usage = "# " + string(argv[0]) + "\n";
  usage += "\t@topo       <molecular topology file>\n";
  usage += "\t[@mol       <molecule that the rdcs are for (default: 1)>]\n";
  usage += "\t@rdc        <RDC data file. Be aware that numbers are with respect to the first residue of the current molecule. This is not always what you get from the NMR people.>\n";
  usage += "\t[@weights   <RDC-specific weight factors>]\n";
  usage += "\t[@nmf       <number of magnetic field vectors (default: 3)>]\n";
  usage += "\t[@mass      <mass of magnetic field representation particles (default: 10.0)>]\n";
  usage += "\t@type       <type of RDC>\n";
  usage += "\t@lib        <RDC library file>\n";
  usage += "\t@group      <file with RDC groups, or 'ONE' for one group with all RDCs>\n";

  try {
    args::Arguments args(argc, argv, knowns, usage);

    // read topology
    gio::InTopology it(args["topo"]);
    gcore::System sys(it.system());

    // read molecule number
    unsigned molnum = 0;
    if (args.count("mol") == 1) {
      istringstream is(args["mol"]);
      int tmp_molnum = 0;
      is >> tmp_molnum;
      if (tmp_molnum < 1) {
        throw gromos::Exception("make_rdc_spec", "Molecule indices should be positive, starting at 1.");
      }
      // convert to internal GROMOS numbering
      molnum = tmp_molnum - 1;
    }
    // check number of molecules
    if (molnum + 1 > unsigned(sys.numMolecules())) {
#if (__cplusplus > 199711L)
      throw gromos::Exception("make_rdc_spec", "You asked for molecule " + to_string(molnum + 1) + ", but the topology only contains " + to_string(sys.numMolecules()) + " molecules.");
#else
      throw gromos::Exception("make_rdc_spec", "The topology does not contain the specified molecule.");
#endif
    }

    // read the type of RDC
    int rdctype;
    if (args.count("type") != 1) {
      throw gromos::Exception("make_rdc_spec", "No or too many RDC type(s) given");
    } else {
      istringstream is(args["type"]);
      is >> rdctype;
    }

    // get information for this rdc type from the library
    if (args.count("lib") != 1)
      throw gromos::Exception("make_rdc_spec", "No RDC library file");
    gio::Ginstream lib_file(args["lib"]);
    vector<string> lib_buf;
    lib_file.getblock(lib_buf);
    if (lib_buf[0] != "RDCSPEC")
      throw gromos::Exception("make_rdc_spec", "The RDC library file does not contain an RDCSPEC block.");
    if (lib_buf[lib_buf.size() - 1].find("END") != 0)
      throw gromos::Exception("make_rdc_spec", "RDC library file " + lib_file.name() + " is corrupted. No END in RDCSPEC block. Got\n" +
                                                   lib_buf[lib_buf.size() - 1]);

    rdc_struct rdc_spec;
    bool found_rdc_type = false;
    vector<string>::const_iterator il = lib_buf.begin() + 1, lo = lib_buf.end() - 1;
    for (; il != lo; ++il) {
      istringstream line(*il);
      line >> rdc_spec.type;
      if (rdc_spec.type == rdctype) {
        found_rdc_type = true;
        line >> rdc_spec.name_i >> rdc_spec.name_j >> rdc_spec.name_k >> rdc_spec.name_l >> rdc_spec.gyr_i >> rdc_spec.gyr_j >>
            rdc_spec.rij >> rdc_spec.rik;
        if (line.fail())
          throw gromos::Exception("make_rdc_spec", "bad line in RDCSPEC block of library.");
        break;
      }
    } // RDCSPEC block
    if (!found_rdc_type) {
      throw gromos::Exception("make_rdc_spec", "Unknown RDC type given in library");
    }

    // read in the rdc data
    vector<rdc_data_struct> rdc_data;
    string rdc_filename;
    if (args.count("rdc") != 1) {
      throw gromos::Exception("make_rdc_spec", "No RDC data file");
    } else {
      rdc_filename = args["rdc"].c_str();
      ifstream rdc_file(args["rdc"].c_str());
      string line;
      std::string::size_type iterator;
      istringstream is;
      while (true) {
        std::getline(rdc_file, line);
        if (rdc_file.eof())
          break;
        if ((iterator = line.find('#')) != std::string::npos) {
          if (iterator == 0)
            continue;
          line = line.substr(0, iterator);
        }
        rdc_data_struct data;
        data.num_i = -1, data.num_j = -1, data.num_k = -1, data.num_l = -1;
        data.weight = 1.0;
        is.clear();
        is.str(line);
        if (!(is >> data.residue >> data.rdc))
          throw gromos::Exception("make_rdc_spec", "Bad line or non-existent RDC data file.");
        // set weights to one
        rdc_data.push_back(data);
      }
    }

    // read the rdc-specific weights, if present
    if (args.count("weights") == 1) {
      ifstream weights_file(args["weights"].c_str());
      string line;
      std::string::size_type iterator;
      istringstream is;
      while (true) {
        std::getline(weights_file, line);
        if (weights_file.eof())
          break;
        if ((iterator = line.find('#')) != std::string::npos) {
          if (iterator == 0)
            continue;
          line = line.substr(0, iterator);
        }
        is.clear();
        is.str(line);
        int resnum;
        double weight;
        if (!(is >> resnum >> weight)) {
          throw gromos::Exception("make_rdc_spec", "bad line in RDC weight file.");
        } else {
          // find right place in data struct
          bool found = false;
          for (unsigned i = 0; i < rdc_data.size(); i++) {
            int residue = rdc_data[i].residue;
            if (residue == resnum) {
              rdc_data[i].weight = weight;
              found = true;
            }
          }
          if (!found)
            throw gromos::Exception("make_rdc_spec", "residues with weights do not match residues with RDCs");
        }
      }
    }

    // get the number of magnetic field vectors
    int nmf = 3;
    if (args.count("nmf") == 1) {
      istringstream is(args["nmf"]);
      is >> nmf;
    }

    // get the number of magnetic field representations
    double mass = 10.0;
    if (args.count("nmf") == 1) {
      istringstream is(args["mass"]);
      is >> mass;
    }

    // initialise the output
    ostringstream out;
    // title block
    out << "TITLE\n"
        << "RDC specifications created from " << rdc_filename << " using make_rdc_spec\n"
        << "END\n";
    // conversion factors block
    out << "CONVERSION\n"
        << "# factors to convert the frequency from [RDC]=(s-1) to (ps-1)\n"
        << "# and to convert the gyromagnetic ratios from [gamma]=10^6*(rad/T s) to (e/u)\n"
        << "0.000000000001 0.010364272\n"
        << "END\n";
    // magnetic field vectors block (currently all vectors are initialised identically by default)
    out << "INPUTMODE\n"
        << "# input mode: choose '0' for basic and '1' for advanced.\n"
        << "0\n"
        << "END\n";
    out << "MAGFIELDC\n"
        << "# number and mass of magnetic field vectors\n" << nmf << " " << mass << endl;
    out << "END\n";
    // alignment tensor block (masses initialised to 10 by default)
    out << "ALIGNT\n"
        << "# mass of all components\n" << mass << "\n"
        << "END\n";
    // spherical harmonics block (masses initialised to 10 by default)
    out << "SPHERICALHARMONICS\n"
        << "# mass of all components\n" << mass << "\n"
        << "END\n";
    // RDC specification block
    out << "RDCRESSPEC\n"
        << "# For each RDC restraint the following should be specified:\n"
        << "# IPRDCR, JPRDCR, KPRDCR, LPRDCR, atom numbers (defining the vector that forms the angle with the magnetic field)\n"
        << "# WRDCR                           weight factor (for weighting some RDCs higher than others)\n"
        << "# PRDCR0                          RDC restraint value (i.e. the experimental data)\n"
        << "# RDCGI, RDCGJ                    gyromagnetic ratios for atoms i and j\n"
        << "# RDCRIJ, RDCRIK                  distance between atoms i and j or i and k (RIJ = RCH for CA:HA)\n"
        << "# TYPE                            code to define type of RDC\n"
        << "# IPRDCR JPRDCR KPRDCR LPRDCR   WRDCR      PRDCR0      RDCGI      RDCGJ     RDCRIJ     RDCRIK    RDCTYPE\n";

    // find index of first atom in the molecule of interest
    unsigned atom_offset = 0;
    for (unsigned i = 0; i < molnum; i++) {
      atom_offset += sys.mol(i).topology().numAtoms();
    }

    // loop over all atoms in the selected molecule
    for (unsigned i = 0; i < sys.mol(molnum).topology().numAtoms(); ++i) {
      // for C:N RDCs, the N is from the next residue
      int thisRes = sys.mol(molnum).topology().resNum(i) + 1;
      int prevRes = 1;
      if (thisRes >= 1)
        prevRes = thisRes - 1;

      // loop over RDC data
      for (unsigned j = 0; j < rdc_data.size(); ++j) {
        switch (rdctype) {
        // N:H
        case 1:
        // CA:C
        case 2:
          if (thisRes == rdc_data[j].residue) {
            // check atom name
            string atom_name = sys.mol(molnum).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i + 1 + atom_offset;
            }
          }
          rdc_data[j].num_k = 0;
          rdc_data[j].num_l = 0;
          break;
        // C:N
        case 3:
          // find the C
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(molnum).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i + 1 + atom_offset;
            }
          } else if (prevRes == rdc_data[j].residue) {
            // find the N (from residue i+1)
            string atom_name = sys.mol(molnum).topology().atom(i).name();
            if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i + 1 + atom_offset;
            }
          }
          rdc_data[j].num_k = 0;
          rdc_data[j].num_l = 0;
          break;
        // CA:HA - need to get CA, N, C and CB (all from same residue)
        case 4:
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(molnum).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_k) {
              rdc_data[j].num_k = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_l) {
              rdc_data[j].num_l = i + 1 + atom_offset;
            }
          }
          break;
        // side-chain N:H (two possible hydrogens)
        case 5:
        case 6:
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(molnum).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_k) {
              rdc_data[j].num_k = i + 1 + atom_offset;
            }
          }
          rdc_data[j].num_l = 0;
          break;
        // side-chain H:H
        case 7:
        case 8:
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(molnum).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i + 1 + atom_offset;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i + 1 + atom_offset;
            }
          }
          rdc_data[j].num_k = 0;
          rdc_data[j].num_l = 0;
          break;
        } // switch
      }   // RDC data
    }     // atoms

    // loop over RDC data again to write output
    for (unsigned i = 0; i < rdc_data.size(); ++i) {
      if (rdc_data[i].num_i == -1 || rdc_data[i].num_j == -1 || rdc_data[i].num_k == -1 || rdc_data[i].num_l == -1) {
        throw gromos::Exception("make_rdc_spec", "You're inputting RDCs of non-existing atoms. Check if RDCs are numbered correctly and if the sequence is the same for RDCs and topology.");
      }
      out << setw(8) << rdc_data[i].num_i << setw(7) << rdc_data[i].num_j << setw(7) << rdc_data[i].num_k << setw(7) << rdc_data[i].num_l
          << setw(8) << rdc_data[i].weight << setprecision(6) << setw(12) << rdc_data[i].rdc << setw(11) << rdc_spec.gyr_i << setw(11)
          << rdc_spec.gyr_j << setw(11) << rdc_spec.rij << setw(11) << rdc_spec.rik << setw(11) << rdc_spec.type << endl;
    }
    out << "END\n";

    // read rdc group file
    if (args.count("group") != 1)
      throw gromos::Exception("make_rdc_spec", "no RDC groups chosen");
    vector<vector<int> > rdcgroups;
    if (args["group"] == "ONE") {
      rdcgroups.push_back(vector<int>());
      for (unsigned i = 0; i < rdc_data.size(); i++) {
        rdcgroups[0].push_back(i + 1);
      }
    } else {
      gio::Ginstream group_file(args["group"]);
      vector<string> group_buf;
      group_file.getblock(group_buf);
      if (group_buf[0] != "RDCGROUPS")
        throw gromos::Exception("make_rdc_spec", "RDC group file does not contain an RDCGROUPS block.");
      if (group_buf[group_buf.size() - 1].find("END") != 0)
        throw gromos::Exception("make_rdc_spec", "RDC group file " + group_file.name() + " is corrupted. No END in RDCGROUPS block. Got\n" +
                                                     group_buf[group_buf.size() - 1]);

      // we are reusing il from above
      il = group_buf.begin() + 1, lo = group_buf.end() - 1;
      for (; il != lo; ++il) {
        istringstream line(*il);
        int tmp_int;
        vector<int> tmp_vec;
        while (line >> tmp_int) {
          tmp_vec.push_back(tmp_int);
        }
        rdcgroups.push_back(tmp_vec);
      } // RDCGROUPS block
    }

    out << "RDCGROUPS\n"
        << "# one line per group\n";
    for (unsigned k = 0; k < rdcgroups.size(); k++) {
      for (unsigned l = 0; l < rdcgroups[k].size(); l++) {
        out << rdcgroups[k][l] << " ";
      }
      out << "\n";
    }
    out << "END\n";

    cout << out.str();

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
