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
 * @file sasa_hasel.cc
 * compute sasa using Hasel formula (as done in the MD++ SASA implicit solvent model)
 */

/**
 * @page programs Program Documentation
 *
 * @anchor sasa_hasel
 * @section sasa_hasel compute sasa using hasel formula
 * @author @ref ja
 * @date 23. 4. 2009
 *
 * Program sasa_hasel computes the solvent-accessible surface area (sasa)
 * of all atoms in the solute part of the molecular system according to the
 * method of Hasel et al. [Tetra. Comput. Method., 1, 103-116, (1988)]. This is
 * the same method implemented in the SASA/VOL implicit solvent model. If a
 * single conformation is given, either the atomic sasa values or the total sasa,
 * along with the hydrophilic, hydrophobic and "other" contributions (defined by the sign
 * of the sigma values given in the sasaspec file) may be printed. If multiple
 * conformations are given, the averaged totals and, if requested,
 * the time-series of the total sasa values may be printed.
 * 
 * Note that the sasa_hasel program uses the Neighbours algorithm in prep_bb.cc 
 * to generate lists of bonded atoms, and then uses these to generate lists of first, 
 * second, third and higher neighbours that are then used in the sasa calculation algorithm.
 * The atoms which belong to different molecules will ALWAYS be higher neighbours.

 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt" (optional and only if time-series)&gt;] </td></tr>
 * <tr><td> [\@timeseries</td><td>&lt;if you want the time-series as well as the average (if not atomic)&gt;] </td></tr>
 * <tr><td> [\@timespec</td><td>&lt;timepoints at which to compute the sasa: ALL (default), EVERY or SPEC (if time-series)&gt;] </td></tr>
 * <tr><td> [\@timepts</td><td>&lt;timepoints at which to compute the sasa (if time-series and timespec EVERY or SPEC)&gt;] </td></tr>
 * <tr><td> [\@atomic</td><td>&lt;print atomic sasa (if not time-series - only for single frame)&gt;] </td></tr>
 * <tr><td> \@sasaspec</td><td>&lt;sasa specification library file&gt; </td></tr>
 * <tr><td> \@probe</td><td>&lt;IAC of central atom of solvent molecule and radius of solvent molecule (e.g. 5 0.14 nm for 53A6 SPC water)&gt; </td></tr>
 * <tr><td> [\@noH</td><td>&lt;do not include hydrogen atoms in the sasa calculation (default: include)&gt;] </td></tr>
 * <tr><td> [\@p_12</td><td>&lt;overlap parameter for bonded atoms (default: 0.8875)&gt;] </td></tr>
 * <tr><td> [\@p_13</td><td>&lt;overlap parameter for atoms separated by two bonds (default: 0.3516)&gt;] </td></tr>
 * <tr><td> [\@p_1x</td><td>&lt;overlap parameter for atoms separated by more than one bond (default: 0.3516)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory file(s)&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   atominfo
     @topo       ex.top
     @pbc        v
     @timeseries
     @timespec   EVERY
     @timepts    100
     @sasaspec   sasaspec53b6.lib
     @probe      5 0.14
     @traj       ex.trj
   @endverbatim

 * <hr>
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace utils;

bool compute_sasa(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done);

void calculate_sasa(gcore::System & sys, bool higher, vector<double> surfaces,
        vector<double> & sasa_areas, unsigned int ii, unsigned int mi, unsigned int i,
        unsigned int jj, unsigned int mj, unsigned int j,
        double Ri_Rsolv, const double Rj_Rsolv, const double sum_of_radii,
        const double p_i, const double p_j, const double pij, const double pi);

struct sasa_parameter {
  double radius;
  double probability;
  double sigma;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "timeseries" << "timespec"
          << "timepts" << "atomic" << "sasa_spec" << "probe" << "noH"
          << "p_12" << "p_1x" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gather method>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t[@timeseries <if you want the time-series as well as the average>]\n";
  usage += "\t[@timespec   <timepoints at which to compute the SASA: ALL (default), EVERY or SPEC>])\n";
  usage += "\t[@timepts    <timepoints at which to compute the SASA (if timespec EVERY or SPEC)>]\n";
  usage += "\t[@atomic     <print atomic sasa (only if not time-series)>]\n";
  usage += "\t@sasa_spec   <sasa specification file>\n";
  usage += "\t@probe       <IAC and radius of solvent molecule>\n";
  usage += "\t[@noH        <do not include hydrogen atoms in the sasa calculation (default: include)]\n";
  usage += "\t[@p_12       <overlap parameter for bonded atoms> (default: 0.8875)]\n";
  usage += "\t[@p_13       <overlap parameter for atoms separated by two bonds> (default: 0.3516)]\n";
  usage += "\t[@p_1x       <overlap parameter for atoms separated by more than one bond> (default: 0.3516)]\n";
  usage += "\t@traj        <trajectory file(s)>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    gcore::System sys(it.system());

    System refSys(it.system());

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // get time
    Time time(args);

    // parse timespec
    string timespec = "ALL";
    vector<int> timepts;
    if (args.count("timespec") > 0) {
      timespec = args["timespec"];
      if (timespec != "ALL" && timespec != "EVERY" && timespec != "SPEC")
        throw gromos::Exception("sasa_hasel",
              "timespec format " + timespec + " unknown.\n");
      if (timespec == "EVERY" || timespec == "SPEC") {
        for (Arguments::const_iterator it = args.lower_bound("timepts");
                it != args.upper_bound("timepts"); ++it) {
          int bla = atoi(it->second.c_str());
          timepts.push_back(bla);
        }
        if (timepts.size() == 0) {
          throw gromos::Exception("sasa_hasel",
                  "if you give EVERY or SPEC you have to use "
                  "@timepts as well");
        }
        if (timepts.size() != 1 && timespec == "EVERY") {
          throw gromos::Exception("sasa_hasel",
                  "if you give EVERY you have to give exactly"
                  " one number with @timepts");
        }
      }
    } else if (args.count("timepts") > 0) {
      throw gromos::Exception("sasa_hasel",
              "timepts is only read if you give @timespec too");
    }
    // check if we want to print the sasa of every atom
    bool sasa_at = false;
    if (args.count("atomic") != -1) {
      sasa_at = true;
      // can't check ALL here because that is the default
      if (timespec == "EVERY" ||
              (timespec == "SPEC" && timepts.size() > 1)) {
        throw gromos::Exception("sasa_hasel",
                "printing the atomic SASA for multiple frames is not allowed");
      }
    }
    // check if we want to print the time-series
    bool sasa_ts = false;
    if (args.count("timeseries") != -1) {
      sasa_ts = true;
      if (sasa_at) {
        throw gromos::Exception("sasa_hasel",
                "you cannot print a time-series of the atomic SASA");
      }
      if (timespec == "SPEC" && timepts.size() == 1) {
        throw gromos::Exception("sasa_hasel",
                "you cannot have a time-series for one frame");
      }
    }

    // store sasa specifications according to IAC
    if (args.count("sasa_spec") != 1)
      throw gromos::Exception("sasa_hasel", "No sasa specification file");

    map<int, sasa_parameter> sasa_spec;
    {
      Ginstream spec_file(args["sasa_spec"]);
      vector<string> buffer;
      spec_file.getblock(buffer);
      if (buffer[0] != "SASASPEC")
        throw gromos::Exception("sasa_hasel",
              "sasa specification file does not contain a SASASPEC block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("sasa_hasel", "sasa specification file " + spec_file.name() +
              " is corrupted. No END in SASASPEC"
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      vector<string>::const_iterator it = buffer.begin() + 1, to = buffer.end() - 1;
      for (; it != to; ++it) {
        sasa_parameter param;
        istringstream line(*it);
        unsigned int num;
        line >> param.radius >> param.probability >> param.sigma >> num;
        if (line.fail())
          throw gromos::Exception("sasa_spec",
                "bad line in SASASPEC block!");

        for (unsigned int i = 0; i < num; ++i) {
          int iac;
          line >> iac;
          if (line.fail()) {
            ostringstream msg;
            msg << "bad line in SASASPEC block: could not read " << num
                    << " IACs from line.";
            throw gromos::Exception("sasa_spec", msg.str());
          }
          sasa_spec[iac - 1] = param;
        } // for iacs
      } // SASASPEC block
    }

    // get the solvent IAC and radius
    vector<double> probe = args.getValues<double>("probe", 2, true);
    int probe_iac = int(probe[0]);
    double R_solv = probe[1];

    utils::compute_atomic_radii_vdw(probe_iac, R_solv, sys, it.forceField());

    // check whether we include hydrogens or not
    bool noH = false;
    if (args.count("noH") != -1)
      noH = true;

    // check if we want to change the p_12, p_13 and p_1x values
    double p_12 = args.getValue<double>("p_12", false, 0.8875);
    double p_13 = args.getValue<double>("p_13", false, 0.3516);
    double p_1x = args.getValue<double>("p_1x", false, 0.3516);

    // define input coordinate
    InG96 ic;

    // check we have trajectories
    if (args.count("traj") == -1) {
      throw gromos::Exception("sasa_hasel", "No coordinate file");
    }

    // INITIALISE NEIGHBOUR LISTS
    // read in first frame
    ic.open(args.lower_bound("traj")->second);
    ic.select("SOLUTE");
    ic >> sys;
    ic.close();
    // gather
    (*pbc.*gathmethod)();

    // find number of "sasa" atoms and total number of atoms
    // and store the true atom number of each sasa atom in a list
    // of size numSasaAtoms
    vector<unsigned int> sasa_atoms;
    vector<unsigned int> sasa_mols;
    unsigned int numSasaAtoms = 0;
    unsigned int numAtoms = 0;
    for (int m = 0; m < sys.numMolecules(); ++m) {
      numAtoms += sys.mol(m).numAtoms();
      for (int i = 0; i < sys.mol(m).numAtoms(); ++i) {
        if (!noH || !sys.mol(m).topology().atom(i).isH()) {
          numSasaAtoms += 1;
          sasa_atoms.push_back(i);
          sasa_mols.push_back(m);
        }
      }
    }

    // define a bondnumber storage list
    vector<unsigned int> pathlength(numSasaAtoms * numSasaAtoms);
    // and fill with large numbers ("infinity")
    const unsigned int infinity = numeric_limits<unsigned int>::max();
    for (unsigned int i = 0; i < pathlength.size(); ++i) {
      pathlength[i] = infinity;
    }

    // define neighbour lists
    vector<vector<unsigned int> > first_neighbours(numSasaAtoms);
    vector<vector<unsigned int> > second_neighbours(numSasaAtoms);
    vector<vector<unsigned int> > third_neighbours(numSasaAtoms);
    vector<vector<unsigned int> > higher_neighbours(numSasaAtoms);
    vector<double> surfaces(numSasaAtoms);
    vector<double> sasa_areas(numSasaAtoms);

    // define physical constants
    const double pi = gmath::physConst.get_pi();

    // loop through sasa atoms. have to loop over all so that surface is initialised.
    for (unsigned int ii = 0; ii < numSasaAtoms; ++ii) {

      // temporary first neighbour storage
      vector<unsigned int> fnb;

      // first compute and store total surface area of atom ii
      // and get radius, probability and sigma for atom ii
      unsigned int mol_i = sasa_mols[ii];
      unsigned int atom_i = sasa_atoms[ii];
      map<int, sasa_parameter>::const_iterator para_ii =
              sasa_spec.find(sys.mol(mol_i).topology().atom(atom_i).iac());
      const sasa_parameter & s = para_ii->second;
      const double Rsum = s.radius + R_solv;
      double S = 4.0 * pi * (Rsum) * (Rsum);
      surfaces[ii] = S;

      // get the bonded neighbours of atom i in molecule m
      Neighbours neighbours(sys, mol_i, atom_i);
      // loop over neighbours
      Neighbours::const_iterator itn = neighbours.begin(), ton = neighbours.end();
      for (; itn != ton; ++itn) {
        
        // now loop over all other  sasa atoms jj!=ii
        for (unsigned int jj = ii + 1; jj < numSasaAtoms; ++jj) {

          unsigned int atom_j = sasa_atoms[jj];
          unsigned int mol_j = sasa_mols[jj];
          
          // if the atoms are not from the same molecule, they are not considered here
          if (mol_j == mol_i) {        
            if (*itn == int(atom_j)) { // j is bonded to i
              // assign ii and jj a pathlength of 1
              pathlength[ii * numSasaAtoms + jj] = 1;
              pathlength[jj * numSasaAtoms + ii] = 1;
              // and store in tmp first neighbours list
              // note we store indexes not actual atom numbers
              fnb.push_back(jj);
            } // end *itn==j
          }
        } // end jj
      } // end neighbours
      // store first neighbours for atom ii
      first_neighbours[ii] = fnb;
    } // end ii

    // do Floyds algorithm to find minimum path lengths
    for (unsigned int k = 0; k < numSasaAtoms; ++k) {
      for (unsigned int i = 0; i < numSasaAtoms - 1; ++i) {
        // stored path between i and k
        const unsigned int pathlength_ik = pathlength[i * numSasaAtoms + k];
        if (pathlength_ik < infinity) {
          // stored path between k and j
          for (unsigned int j = i + 1; j < numSasaAtoms; ++j) {
            const unsigned int pathlength_kj = pathlength[k * numSasaAtoms + j];
            if (pathlength_kj < infinity) {
              // path between i and j via k
              const unsigned int complength = pathlength_ik + pathlength_kj;
              // compare to stored path between i and j
              if (complength < pathlength[i * numSasaAtoms + j]) {
                pathlength[i * numSasaAtoms + j] = complength;
              }
            }
          }
        }
      }
    } // end floyd's

    // now make the third and higher lists
    for (unsigned int ii = 0; ii < numSasaAtoms - 1; ++ii) {

      // temporary neighbour lists
      vector<unsigned int> snb;
      vector<unsigned int> tnb;
      vector<unsigned int> hnb;
      for (unsigned int jj = ii + 1; jj < numSasaAtoms; ++jj) {
          
        unsigned int mol_i = sasa_mols[ii];  
        unsigned int mol_j = sasa_mols[jj];
          
        // if the atoms are not from the same molecule, they are not considered only as higher neighbour atoms 
        if (mol_j != mol_i) {
          hnb.push_back(jj);   
        }
        
        const unsigned int pathlength_ij = pathlength[ii * numSasaAtoms + jj];
        if (pathlength_ij < infinity) {
          if (pathlength_ij == 2) {
            // store in 2nd neighbours
            snb.push_back(jj);
          }
          if (pathlength_ij == 3) {
            // store in 3rd neighbours
            tnb.push_back(jj);
          } else if (pathlength_ij > 3) {
            // note we are storing the indexes not the actual atom numbers
            hnb.push_back(jj);
          }
        }
      }
      // now put tmp neighbour lists into overall neighbour lists
      second_neighbours[ii] = snb;
      third_neighbours[ii] = tnb;
      higher_neighbours[ii] = hnb; // note last space will not be filled...
    } // end ii

    // declare some variables for averaging
    double ave_phobic_sasa = 0.0;
    double ave_philic_sasa = 0.0;
    double ave_other_sasa = 0.0;
    double ave_tot_sasa = 0.0;

    // start at -1 to get times right
    int num_frames = -1;
    // number of time-points for which SASAs have been written
    unsigned int times_written = 0;
    // for SPEC: so that we stop trying when all requested timepoints are written
    bool done = false;

    // print titles for time-series (if wanted)
    if (sasa_ts) {
      cout << "#    Time          hydrophobic       hydrophilic             "
              "other             total" << endl;
    }

    // loop over all trajectories
    for (Arguments::const_iterator
      iter = args.lower_bound("traj"),
            to = args.upper_bound("traj");
            iter != to; ++iter) {

      // open file
      ic.open((iter->second).c_str());
      ic.select("SOLUTE");

      // loop over single trajectory
      while (!ic.eof()) {

        // declare variables for this frame
        double phobic_sasa = 0.0;
        double philic_sasa = 0.0;
        double other_sasa = 0.0;
        double tot_sasa = 0.0;

        ic >> sys >> time;
        // gather
        (*pbc.*gathmethod)();

        // check whether to skip or not
        num_frames++;

        if (compute_sasa(num_frames, timespec, timepts, times_written, done)) {

          // write headers for atomic sasa
          if (sasa_at) {
            cout.precision(10);
            cout.setf(ios::right, ios::adjustfield);
            cout << setw(6) << "# Atom" << setw(18) << "SASA" << endl;
          }

          // initialise sasa array
          sasa_areas = surfaces;

          // loop through sasa atoms
          for (unsigned int ii = 0; ii < numSasaAtoms; ++ii) {

            unsigned int atom_i = sasa_atoms[ii];
            unsigned int mol_i = sasa_mols[ii];

            // get radius, probability and sigma for atom ii
            map<int, sasa_parameter>::const_iterator para_ii =
                    sasa_spec.find(sys.mol(mol_i).topology().atom(atom_i).iac());

            const sasa_parameter & s = para_ii->second;
            const double R_i = s.radius;
            const double p_i = s.probability;
            const double Ri_Rsolv = R_i + R_solv;

            // not higher neighbours yet
            bool higher = false;

            // first neighbours
            for (unsigned int j = 0; j < first_neighbours[ii].size(); ++j) {

              // get atom number
              unsigned int jj = first_neighbours[ii][j];
              if (jj < ii) {
                throw gromos::Exception("sasa_hasel",
                        "neighbours must have higher sequence numbers");
              } else {
                unsigned int atom_j = sasa_atoms[jj];
                unsigned int mol_j = sasa_mols[jj];

                // get radius (and probability and sigma) for atom jj
                map<int, sasa_parameter>::const_iterator para_jj =
                        sasa_spec.find(sys.mol(mol_j).topology().atom(atom_j).iac());

                const sasa_parameter & t = para_jj->second;
                const double R_j = t.radius;
                double p_j = t.probability;

                // compute sum of radii
                const double Rj_Rsolv = R_j + R_solv;
                const double sum_of_radii = Ri_Rsolv + Rj_Rsolv;
                
                calculate_sasa(sys, higher, surfaces, sasa_areas, ii, mol_i, atom_i,
                        jj, mol_j, atom_j, Ri_Rsolv, Rj_Rsolv, sum_of_radii, p_i, p_j, p_12, pi);
              }

            } // end j

            // second neighbours
            for (unsigned int j = 0; j < second_neighbours[ii].size(); ++j) {

              // get atom number
              unsigned int jj = second_neighbours[ii][j];
              if (jj < ii) {
                throw gromos::Exception("sasa_hasel",
                        "neighbours must have higher sequence numbers");
              } else {
                unsigned int atom_j = sasa_atoms[jj];
                unsigned int mol_j = sasa_mols[jj];

                // get radius (and probability and sigma) for atom jj
                map<int, sasa_parameter>::const_iterator para_jj =
                        sasa_spec.find(sys.mol(mol_j).topology().atom(atom_j).iac());

                const sasa_parameter & t = para_jj->second;
                const double R_j = t.radius;
                double p_j = t.probability;

                // compute sum of radii
                const double Rj_Rsolv = R_j + R_solv;
                const double sum_of_radii = Ri_Rsolv + Rj_Rsolv;

                calculate_sasa(sys, higher, surfaces, sasa_areas, ii, mol_i, atom_i,
                        jj, mol_j, atom_j, Ri_Rsolv, Rj_Rsolv, sum_of_radii, p_i, p_j, p_13, pi);
              }
            }

            higher = true;
            // third neighbours
            for (unsigned int j = 0; j < third_neighbours[ii].size(); ++j) {

              // get atom number
              unsigned int jj = third_neighbours[ii][j];
              if (jj < ii) {
                throw gromos::Exception("sasa_hasel",
                        "neighbours must have higher sequence numbers");
              } else {
                unsigned int atom_j = sasa_atoms[jj];
                unsigned int mol_j = sasa_mols[jj];

                // get radius (and probability and sigma) for atom jj
                map<int, sasa_parameter>::const_iterator para_jj =
                        sasa_spec.find(sys.mol(mol_j).topology().atom(atom_j).iac());

                const sasa_parameter & t = para_jj->second;
                const double R_j = t.radius;
                double p_j = t.probability;

                // compute sum of radii
                const double Rj_Rsolv = R_j + R_solv;
                const double sum_of_radii = Ri_Rsolv + Rj_Rsolv;

                calculate_sasa(sys, higher, surfaces, sasa_areas, ii, mol_i, atom_i,
                        jj, mol_j, atom_j, Ri_Rsolv, Rj_Rsolv, sum_of_radii, p_i, p_j, p_1x, pi);
              }
            }

            // higher neighbours
            for (unsigned int j = 0; j < higher_neighbours[ii].size(); ++j) {

              // get atom number
              unsigned int jj = higher_neighbours[ii][j];
              if (jj < ii) {
                throw gromos::Exception("sasa_hasel",
                        "neighbours must have higher sequence numbers");
              } else {
                unsigned int atom_j = sasa_atoms[jj];
                unsigned int mol_j = sasa_mols[jj];

                // get radius (and probability and sigma) for atom jj
                map<int, sasa_parameter>::const_iterator para_jj =
                        sasa_spec.find(sys.mol(mol_j).topology().atom(atom_j).iac());

                const sasa_parameter & t = para_jj->second;
                const double R_j = t.radius;
                double p_j = t.probability;

                // compute sum of radii
                const double Rj_Rsolv = R_j + R_solv;
                const double sum_of_radii = Ri_Rsolv + Rj_Rsolv;

                calculate_sasa(sys, higher, surfaces, sasa_areas, ii, mol_i, atom_i,
                        jj, mol_j, atom_j, Ri_Rsolv, Rj_Rsolv, sum_of_radii, p_i, p_j, p_1x, pi);
              }
            }

          } // ii

          // new loop over sasa atoms now that all back/forward calculations are done
          for (unsigned int ii = 0; ii < numSasaAtoms; ++ii) {

            unsigned int atom_i = sasa_atoms[ii];
            unsigned int mol_i = sasa_mols[ii];
            // get radius (and probability and sigma) for atom ii
            map<int, sasa_parameter>::const_iterator para_ii =
                    sasa_spec.find(sys.mol(mol_i).topology().atom(atom_i).iac());

            const sasa_parameter & s = para_ii->second;
            const double sigma_i = s.sigma;

            // add to total
            tot_sasa += sasa_areas[ii];

            // assign sasa to various types according to sigma
            if (sigma_i > 0)
              phobic_sasa += sasa_areas[ii];
            else if (sigma_i < 0)
              philic_sasa += sasa_areas[ii];
            else if (sigma_i == 0)
              other_sasa += sasa_areas[ii];

            // print atomic SASA
            if (sasa_at) {
              cout.precision(10);
              cout.setf(ios::right, ios::adjustfield);
              fprintf(stdout, "%6d %17.5f\n", atom_i + 1, sasa_areas[ii]);
            }

          } // ii

          // if time-series, print areas for this frame
          if (sasa_ts) {
            cout.precision(10);
            cout.setf(ios::right, ios::adjustfield);
            double timeout = time.time();
            fprintf(stdout, "%12.5f %17.5f %17.5f %17.5f %17.5f\n",
                    timeout, phobic_sasa, philic_sasa, other_sasa, tot_sasa);
          }

          //store values for averaging averages
          ave_phobic_sasa += phobic_sasa;
          ave_philic_sasa += philic_sasa;
          ave_other_sasa += other_sasa;
          ave_tot_sasa += tot_sasa;

        } // if compute sasa
        if (done)
          break;

      } // while not eof
    } // loop over traj

    // print out averages (remember num_frames starts from -1)
    if (num_frames > 0) {
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      if (sasa_ts) cout << endl;
      cout << "#   Averages:      hydrophobic       hydrophilic             "
              "other             total" << endl;
      fprintf(stdout, "#            %17.5f %17.5f %17.5f %17.5f\n",
      ave_phobic_sasa / times_written, ave_philic_sasa / times_written,
              ave_other_sasa / times_written, ave_tot_sasa / times_written);
      // otherwise print totals
    } else {
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      if (sasa_at) cout << endl;
      cout << "#   Totals:        hydrophobic       hydrophilic             "
              "other             total" << endl;
      fprintf(stdout, "#            %17.5f %17.5f %17.5f %17.5f\n",
              ave_phobic_sasa, ave_philic_sasa, ave_other_sasa, ave_tot_sasa);
    }


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

bool compute_sasa(int i, std::string const & timespec, vector<int> const & timepts,
        unsigned int & timesWritten, bool & done) {
  if (timespec == "ALL") {
    ++timesWritten;
    return true;
  } else if (timespec == "EVERY" && i % timepts[0] == 0) {
    ++timesWritten;
    return true;
  } else if (timespec == "SPEC") {
    for (unsigned int j = 0; j < timepts.size(); ++j) {
      if (timepts[j] == i) {
        ++timesWritten;
        if (timesWritten == timepts.size())
          done = true;
        return true;
      } // compute
    } // times
  }
  return false;
}

void calculate_sasa(gcore::System & sys, bool higher, vector<double> surfaces,
        vector<double> & sasa_areas, unsigned int ii, unsigned int mi, unsigned int i,
        unsigned int jj, unsigned int mj, unsigned int j,
        const double Ri_Rsolv, const double Rj_Rsolv, const double sum_of_radii,
        const double p_i, const double p_j, const double pij, const double pi) {

  // compute distance between atoms i and j
  const double rdist = (sys.mol(mi).pos(i) - sys.mol(mj).pos(j)).abs();

  if (!higher) {
    // compute components of area function
    const double c1 = (sum_of_radii - rdist) * pi;
    const double c2 = (Rj_Rsolv - Ri_Rsolv) / rdist;
    // note that "bij", "bji" are actually pi*pij*bij/Si
    const double bij = (c1 * Ri_Rsolv * (1.0 + c2) * pij * p_i) /
            surfaces[ii];
    const double bji = (c1 * Rj_Rsolv * (1.0 - c2) * pij * p_j) /
            surfaces[jj];
    // modify areas
    sasa_areas[ii] *= (1.0 - bij);
    sasa_areas[jj] *= (1.0 - bji);
  } else {
    if (rdist <= sum_of_radii) {
      // compute components of area function
      const double c1 = (sum_of_radii - rdist) * pi;
      const double c2 = (Rj_Rsolv - Ri_Rsolv) / rdist;
      // note that "bij", "bji" are actually pi*pij*bij/Si
      const double bij = (c1 * Ri_Rsolv * (1.0 + c2) * pij * p_i) /
              surfaces[ii];
      const double bji = (c1 * Rj_Rsolv * (1.0 - c2) * pij * p_j) /
              surfaces[jj];
      // modify areas
      sasa_areas[ii] *= (1.0 - bij);
      sasa_areas[jj] *= (1.0 - bji);
    }
  }
}
