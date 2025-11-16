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
 * @file cry_rms.cc
 * Investigate symmetry of molecules
 */

/**
 * @page programs Program Documentation
 *
 * @anchor cry_rms
 * @section cry_rms Investigate symmetry of unit cells
 * @author @ref ns
 * @date 08.10.2010
 *
 * Program cry_rms is used to compute atom positional RMSDs and RMSFs between the
 * asymmetric units within a unit cell of a crystalline system. The symmetry operations
 * are either specified using a special file (\@spec, \@factor) or by the space group
 * (\@spacegroup). In order to identify the individual asymmetric units (ASUs)
 * an @ref AtomSpecifier to the first atom of every ASU have to
 * be given (\@asuspec).
 * If an RMSD is requested (\@atomsrmsd), the atom positional RMSD between all
 * the asymmetric units is printed in separate columns.
 * The RMSF is calculated for the requested atoms (\@atomsrmsf) while taking
 * also the fluctuations of the symmetry related copies into account.
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory coordinate files&gt; </td></tr>
 * <tr><td>[\@spec</td><td>&lt;specification file for the symmetry transformations]</td></tr>
 * <tr><td>[\@factor</td><td>&lt;conversion factor for distances&gt;]</td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup symbol&gt;]</td></tr>
 * <tr><td>\@asuspec</td><td>&lt;@ref AtomSpecifier to the first atom of every asymetric unit&gt;</td></tr>
 * <tr><td>[\@atomsrmsd</td><td>&lt;@ref AtomSpecifier used for RMSD calculation&gt;]</td></tr>
 * <tr><td>[\@atomsrmsf</td><td>&lt;@ref AtomSpecifier used for RMSF calculation&gt;]</td></tr>
 * </table>
 *
 * Example using a specification file:
 * @verbatim
  ccry_rms
    @topo      ex.top
    @pbc       r
    @traj      ex.trc.gz
    @spec      cry.spec
    @factor    0.1
    @asuspec   1:1 2:1 3:1 4:1
    @atomsrmsd 1:CA,N,C
 @endverbatim
 * Example using a spacegroup
 * @verbatim
  ccry_rms
    @topo       ex.top
    @pbc        r
    @traj       ex.trc.gz
    @spacegroup P 21 21 21
    @asuspec    1:1 2:1 3:1 4:1
    @atomsrmsf  1:CA
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <iostream>
#include <map>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/groTime.h"

#include "../config.h"
#include "../src/gromos/Exception.h"

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#endif

using namespace std;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace bound;
using namespace args;
using namespace utils;

void read_spec(std::string name,
        vector<Matrix> & rotation,
        vector<Vec> &translation,
        double factor);

void set_box_dimensions(const System &sys);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "traj" << "time" << "pbc" << "spec" << "factor" << "spacegroup" << "asuspec" << "atomsrmsd" << "atomsrmsf";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pbc         <boundary type>\n";
  usage += "\t@traj        <trajectory coordinate files>\n";
  usage += "\t[@time       <t0 dt>]\n";
  usage += "\t[@spec       <specification file for the symmetry transformations]\n";
  usage += "\t[@factor     <conversion factor for distances>]\n";
  usage += "\t[@spacegroup <spacegroup symbol, Hall or Hermann-Mauguin>]\n";
  usage += "\t@asuspec     <atom specifiers to first atom of every ASU>]\n";
  usage += "\t[@atomsrmsd  <atoms used for RMSD calculation>]\n";
  usage += "\t[@atomsrmsf  <atoms used for RMSF calculation>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    // get simulation time
    Time time(args);
    System sys(it.system());

    System refSys(it.system());

    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // prepare the rotation matrices and translation vectors
    vector<Matrix> rotation;
    vector<Vec> translation;

    // this program can either use a specification file OR a spacegroup
    if (args.count("spec") > 0 && args.count("spacegroup") > 0) {
      throw gromos::Exception(argv[0], "Either give a specification file or a spacegroup.");
    }

    unsigned int num_symop = 0;

    // read the specification file
    if (args.count("spec") > 0) {
      // read the conversion factor
      double factor = args.getValue<double>("factor", false, 1.0);
      // check input for consistency
      if (args.count("cell") > 0)
        throw gromos::Exception(argv[0], "No cell needed when using specification file.");
      read_spec(args["spec"], rotation, translation, factor);
      num_symop = rotation.size();
    }

    // create the matrices and vectors from the spacegroup
#ifdef HAVE_CLIPPER
    clipper::Spacegroup spacegroup;
    bool use_spacegroup = false;
#endif
    
    if (args.count("spacegroup") > 0) {
#ifdef HAVE_CLIPPER
      // check input for consistency
      if (args.count("factor") > 0)
        throw gromos::Exception(argv[0], "No conversion factor needed when using spacegroup.");

      // concatenate the spacegroup symbol from the arguments
      std::string spacegroup_symbol;
      {
        Arguments::const_iterator iter = args.lower_bound("spacegroup");
        Arguments::const_iterator to = args.upper_bound("spacegroup");
        spacegroup_symbol = iter->second;
        for (++iter; iter != to; ++iter) {
          spacegroup_symbol += " ";
          spacegroup_symbol += iter->second;
        }
      }
      // get the spacegroup

      try {
        const clipper::Spgr_descr spacegroup_descriptor(spacegroup_symbol);
        spacegroup.init(spacegroup_descriptor);
      } catch (clipper::Message_fatal & msg) {
        throw gromos::Exception(argv[0], msg.text());
      }

      num_symop = spacegroup.num_symops();
      use_spacegroup = true;
#else
      throw gromos::Exception(argv[0], "You have to compile GROMOS++ with CCP4"
              " and clipper libraries in order to use the spacegroup feature"
              " of cry.\nUse --with-ccp4 and --with-clipper for configuration.");
#endif

    } // if spacegroup

    AtomSpecifier asu_pointer(sys);
    //get structure_factor atoms
    {
      Arguments::const_iterator iter = args.lower_bound("asuspec");
      Arguments::const_iterator to = args.upper_bound("asuspec");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        asu_pointer.addSpecifier(spec);
      }
    }
    if (asu_pointer.size() != (unsigned int)(num_symop)) {
      ostringstream os;
      os << "Please give the first atom of every ASU for @asuspec. Got "
              << asu_pointer.size() << " atoms but need " << num_symop << ":" << endl;
      const vector<string> & str = asu_pointer.toString();
      for(vector<string>::const_iterator it = str.begin(), to = str.end();
              it != to; ++it)
              os << *it << " ";
      os << endl;
      throw gromos::Exception(argv[0], os.str());
    }

    int num_atoms = 0;
    for(int m = 0; m < sys.numMolecules(); ++m)
      num_atoms += sys.mol(m).numAtoms();

    // this vector contains the atomspecifiers for every asu.
    vector<AtomSpecifier> atomsrmsd(num_symop, AtomSpecifier(sys));
    {
      Arguments::const_iterator iter = args.lower_bound("atomsrmsd");
      Arguments::const_iterator to = args.upper_bound("atomsrmsd");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atomsrmsd[0].addSpecifier(spec);
      }
      for (unsigned int i = 1; i < num_symop; ++i) {
        for (unsigned int j = 0; j < atomsrmsd[0].size(); ++j) {
          int rel_pointer = atomsrmsd[0].gromosAtom(j) - asu_pointer.gromosAtom(0);
          int atom = asu_pointer.gromosAtom(i) + rel_pointer;
          if (rel_pointer < 0 || atom >= num_atoms) {
            throw gromos::Exception(argv[0], "Atoms out of range. There is something wrong with @atomsrmsd or @asuspec.");
          }
          atomsrmsd[i].addGromosAtom(atom);
        }
      }
    }

    // this vector contains the atomspecifiers for every asu.
    vector<AtomSpecifier> atomsrmsf(num_symop, AtomSpecifier(sys));
    {
      Arguments::const_iterator iter = args.lower_bound("atomsrmsf");
      Arguments::const_iterator to = args.upper_bound("atomsrmsf");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        atomsrmsf[0].addSpecifier(spec);
      }
      atomsrmsf[0].sort();
      for (unsigned int i = 1; i < num_symop; ++i) {
        for (unsigned int j = 0; j < atomsrmsf[0].size(); ++j) {
          int rel_pointer = atomsrmsf[0].gromosAtom(j) - asu_pointer.gromosAtom(0);
          int atom = asu_pointer.gromosAtom(i) + rel_pointer;
          if (rel_pointer < 0 || atom >= num_atoms) {
            throw gromos::Exception(argv[0], "Atoms out of range. There is something wrong with @atomsrmsf or @asuspec.");
          }
          atomsrmsf[i].addGromosAtom(atom);
        }
      }
    }
    
    if (atomsrmsd[0].empty() && atomsrmsf[0].empty()) {
      throw gromos::Exception(argv[0], "Give at least atoms for either RMSD or RMSF calculation.");
    }

    int num_rmsf_atoms = atomsrmsf[0].size();
    vector<Vec> rmsf_sum;
    vector<double> rmsf_sum2;
    if (num_rmsf_atoms != 0) {
      rmsf_sum.resize(num_rmsf_atoms, Vec(0.0, 0.0, 0.0));
      rmsf_sum2.resize(num_rmsf_atoms, 0.0);
    }

    //===========================
    // loop over all trajectories
    InG96 ic;

    unsigned int frame_number = 0;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);
      ic.select("ALL");

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys >> time;
        ++frame_number;
        if (!sys.hasPos)
          throw gromos::Exception(argv[0],
                "Unable to read POSITION(RED) block from "
                "trajectory file.");

        if (!sys.hasBox)
          throw gromos::Exception(argv[0],
                "Cannot proceed without a box.");

        (*pbc.*gathmethod)();

#ifdef HAVE_CLIPPER
        if (use_spacegroup) {
          // create the cell
          const double a = sys.box().K().abs();
          const double b = sys.box().L().abs();
          const double c = sys.box().M().abs();
          const double alpha = sys.box().alpha();
          const double beta = sys.box().beta();
          const double gamma = sys.box().gamma();

          if (!a || !b || !c || !alpha || !beta || !gamma)
            throw gromos::Exception(argv[0], "Box has zero volume!");

          clipper::Cell_descr cellinit(a, b, c,
                  alpha, beta, gamma);
          const clipper::Cell cell(cellinit);
          // loop over symmetry operations
          rotation.clear();
          translation.clear();
          for (int i = 0; i < spacegroup.num_symops(); ++i) {
            // get the rotation/translation operator from the cell and the sym op.
            const clipper::RTop_orth & rt_operator = spacegroup.symop(i).rtop_orth(cell);
            // get the matrix
            const clipper::Mat33<double> & clipper_matrix = rt_operator.rot();
            Matrix rot;
            // convert the matrix to the GROMOS format
            for (unsigned int row = 0; row < 3; ++row) {
              for (unsigned int col = 0; col < 3; ++col) {
                rot(row, col) = clipper_matrix(row, col);
              }
            }
            rotation.push_back(rot);

            // get the translation vector and convert to nm.
            const clipper::Vec3<double> & clipper_vec = rt_operator.trn();
            const Vec trans(clipper_vec[0], clipper_vec[1], clipper_vec[2]);
            translation.push_back(trans);
          } // loop over symmetry operations
        }
#endif

        // calculate the inverse rotation
        vector<Matrix> inv_rotation;
        for (unsigned int i = 0; i < num_symop; ++i) {
          inv_rotation.push_back(rotation[i].invert());
        }
        \
        if (atomsrmsd[0].size()) {
          if (frame_number == 1) {
            cout << "#" << setw(14) << "time";
            for (unsigned int i = 0; i < num_symop; ++i) {
              for (unsigned int j = i + 1; j < num_symop; ++j) {
                ostringstream os; os << (i+1) << ":" << (j+1);
                cout << setw(15) << os.str(); 
              }
            }
            cout << endl;
          }

          cout << time;
          for (unsigned int i = 0; i < num_symop; ++i) {
            for (unsigned int j = i + 1; j < num_symop; ++j) {
              // calculate RMSD between each pair.
              double sum = 0.0;
              for (unsigned int a = 0; a < atomsrmsd[i].size(); ++a) {
                Vec pos_i(inv_rotation[i] * (atomsrmsd[i].pos(a) - translation[i]));
                Vec pos_j(inv_rotation[j] * (atomsrmsd[j].pos(a) - translation[j]));
                Vec d(pos_i - pbc->nearestImage(pos_i, pos_j, sys.box()));
                sum += d.abs2();
              }
              double rmsd = sqrt(sum / atomsrmsd[i].size());
              cout.precision(8);
              cout << setw(15) << rmsd;
            }
          }
          cout << endl;
        }

        for (int a = 0; a < num_rmsf_atoms; ++a) {
          Vec pos_ref(inv_rotation[0] * (atomsrmsf[0].pos(a) - translation[0]));
          for (unsigned int i = 0; i < num_symop; ++i) {
            Vec pos_i(inv_rotation[i] * (atomsrmsf[i].pos(a) - translation[i]));
            Vec gath_pos(pbc->nearestImage(pos_ref, pos_i, sys.box()));
            rmsf_sum[a] += gath_pos;
            rmsf_sum2[a] += gath_pos.abs2();
          }
        }
      } // while frames
    } // for traj files

    if (num_rmsf_atoms) {
      cout << "#" << setw(4) << "i" << setw(15) << "rmsf" << setw(5) << "name" << endl;
    }
    for (int i=0; i < num_rmsf_atoms; ++i) {
      cout.precision(8);
      rmsf_sum[i] /= frame_number * num_symop;
      rmsf_sum2[i] /= frame_number * num_symop;

      const double rmsf = sqrt(rmsf_sum2[i] - rmsf_sum[i].abs2());
      cout << setw(5) << i+1
	   << setw(15) << rmsf
	   << setw(5) << atomsrmsf[0].name(i)
	   << endl;
    }


  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void read_spec(std::string name,
        vector<Matrix> & rotation,
        vector<Vec> &translation,
        double factor) {
  Ginstream file(name);
  vector<string> buffer;
  file.getblock(buffer);
  file.close();

  if (buffer[0] != "TRANSFORM")
    throw gromos::Exception("cry_cms", "Could not read TRANSFORM block in specification file");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("cry", "Specification file " + file.name() +
          " is corrupted. No END in " + buffer[0] +
          " block. Got\n"
          + buffer[buffer.size() - 1]);
  int num = 0;
  vector<string>::iterator iter = buffer.begin() + 1;

  istringstream is(*iter);
  ++iter;

  if (!(is >> num) || num <= 0)
    throw gromos::Exception("cry", "Need some transformations");
  if (buffer.size() - 3 != unsigned (num * 3))
    throw gromos::Exception("cry", "Line count wrong in " + file.name());

  Matrix rot(3, 3);
  Vec v;

  for (int i = 0; i < num; i++) {
    for (int j = 0; j < 3; j++, ++iter) {
      is.clear();
      is.str(*iter);
      for (int k = 0; k < 3; ++k) {
        if (!(is >> rot(j, k)))
          throw gromos::Exception("cry", "error reading file");
      }
      if (!(is >> v[j]))
        throw gromos::Exception("cry", "error reading file");
    }
    rotation.push_back(rot);
    translation.push_back(factor * v);
  }

}





