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
 * @file rmsdmat.cc
 * Calculates the rmsd-matrix for given structures
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rmsdmat
 * @section rmsdmat Calculates the rmsd-matrix for given structures
 * @author @ref mc @ref co
 * @date 22-8-06
 *
 * Program rmsdmat calculates the root-mean-square deviation
 * between all pairs of structures in a given trajectory file. This matrix of
 * RMSDs can subsequently be used by program @ref cluster to perform a
 * conformational clustering. The matrix can be written out in human readable
 * form, or -to save disk space- in binary format. For efficiency reasons, the
 * RMSD values are written in an integer format. The user can specify the
 * required precision of the RMSD values that are stored. In the case of the
 * binary format, the values are stored as unsigned short int if the precision is
 * less or equal to 4, or otherwise as unsigned int.
 *
 * For an atom-positional RMSD matrix different sets of atoms can be selected to perform a rotational
 * least-squares-fit and to calculate the RMS deviation from. The RMSD matrix
 * can also be calculated from differences in internal coordinates defined by
 * a set of properties (e.g. torsional angles or hydrogen bonds).
 * A selection of
 * structures in the trajectory file to consider can be made using the options
 * skip and stride. Structure pairs may occur for which the least-squares
 * rotational fit fails for numerical reasons. In these cases both structures
 * are fit to the reference structure. If no user specified reference structure
 * is available, the first structure in the trajectory is taken as such.
 * Specifying a reference structure allows the program @ref cluster to perform a
 * forced clustering as well, requiring that the first cluster contains the
 * reference structure, regardless of the cluster size.
 *
 * The program allows calculation of rmsds between related structures with 
 * different topologies if multiple topologies are given under \@topo and 
 * sets of trajectories belonging to each topology are separated by keyword
 * "newgroup". Whenever the keyword is found, the program will simply move on to
 * the next topology. This will only give useful results if atom selections 
 * are chosen such that in each system equivalent atoms are selected.  
 * Optionally an individual topology can also be specified for
 * the reference (\@reftopo), otherwise the first topology is used.
 * 
 * If argument \@dist is specified the distribution of rmsd values is printed to a file.
 *
 * <B>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@reftopo</td><td>&lt;molecular topology file for the reference; default: first topology in \@topo &gt;] </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary conditions&gt; &lt;gather type&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;]</td></tr>
 * <tr><td> [\@atomsrmsd</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsd&gt; </td></tr>
 * <tr><td> [\@prop</td><td>&lt;@ref PropertySpecifier "properties" to be used for rmsd computation.&gt;]</td></tr>
 * <tr><td> [\@skip</td><td>&lt;skip frames at beginning&gt;] </td></tr>
 * <tr><td> [\@stride</td><td>&lt;use only every step frame&gt;] </td></tr>
 * <tr><td> [\@human</td><td>(write the matrix in human readable form)] </td></tr>
 * <tr><td> [\@precision</td><td>&lt;number of digits in the matrix (default 4)&gt;] </td></tr>
 * <tr><td> [\@dist</td><td>&lt;binsize for distribution histogram (default:0.02)&gt;]</td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates&gt;] </td></tr>
 * <tr><td> [\@printatoms</td><td>print list of selected atoms] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rmsdmat
    @topo         sys1.top sys2.top
    @reftopo      ref.top
    @pbc          r
    @atomsfit     1:a
    @atomsrmsd    1:CA
    @skip         5
    @stride       10
    @human
    @precision    4
    @ref          exref.coo
    @dist         0.02 
    @traj         sys1_1.trc sys1_2.trc newgroup sys2_1.trc sys2_2.trc 
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/Value.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/FastRotationalFit.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace bound;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

double props_rmsd(const vector<vector<Value> > &props, int i, int j) {
  double rmsd = 0.0;
  for (unsigned int k=0; k< props[i].size(); k++) {
    // get nearestimagedistance for dihedrals+angles
    // assuming distances or Hbonds will never be >|180|, they should be ok
    double diff = props[j][k].scalar() - props[i][k].scalar(); 
    while (diff >= 180.0) diff -= 360.0;
    while (diff < -180.0) diff += 360.0;

    rmsd += Value(diff).abs2();
  }
  return std::sqrt(rmsd / props[i].size());
}

void update_bins(map<int, int > & bins, double bin_size, double x) {
  int i = std::floor(x / bin_size);
  if ( bins.find(i) == bins.end() ) {
    bins.insert(map<int,int >::value_type(i,1));
  } else {
    bins[i]++;
  }
}

void write_bins(map<int, int > bins, double bin_size) {
  //int maxbin=bins.rbegin()->first;
  //int minbin=bins.begin()->first;
  ofstream distfile;
  distfile.open("distr.out");
  distfile << "# rmsd distribution, bin size " <<bin_size <<" nm" << endl;
  
  for (std::map<int,int>::iterator it=bins.begin(); it!=bins.end(); ++it)
    distfile << it->first*bin_size << " " << it->second << '\n';
}

void getAtomSpecs(System &sys, Arguments &args, AtomSpecifier &atomspecs, AtomSpecifier &fitatoms, AtomSpecifier &rmsdatoms) {
    //clear atomspecs
    atomspecs.clear();
    fitatoms.clear();
    rmsdatoms.clear();
    
    atomspecs.setSystem(sys);
    // read the fit atoms
    fitatoms.setSystem(sys);
    if (args.count("atomsfit")>0) {
      Arguments::const_iterator iter = args.lower_bound("atomsfit"),
              to = args.upper_bound("atomsfit");
      for (; iter != to; ++iter) {
        fitatoms.addSpecifier(iter->second);
      }
      if (fitatoms.size() == 0)
          throw gromos::Exception("rmsdmat",
          "No atoms in atomspecifier (@atomsfit)\n");
    }

    // read the rmsd atoms
    rmsdatoms.setSystem(sys);
    if (args.count("atomsrmsd")>0) {
      Arguments::const_iterator iter = args.lower_bound("atomsrmsd"),
              to = args.upper_bound("atomsrmsd");
      for (; iter != to; ++iter) {
        rmsdatoms.addSpecifier(iter->second);
      }
      if (rmsdatoms.size() == 0) {
          throw gromos::Exception("rmsdmat",
          "No atoms in atomspecifier (@atomsrmsd)\n");
        }
      } else {
          std::cerr << "# WARNING: no @atomsrmsd given, using "
          << "atoms specified by @atomsfit for fit and rmsd.\n";
          rmsdatoms = fitatoms;
    }

    // for which atoms do we want to keep the coordinates
    atomspecs=fitatoms + rmsdatoms;
}

void renewPropSpecs(System &sys, Arguments &args, Boundary* pbc, 
                    PropertyContainer &propspecs) {
    // clear the property container and refill it with the specs so
    //  that they belong to the new system 
    propspecs.reinitialize(sys, pbc);  
    for (Arguments::const_iterator p_iter = args.lower_bound("prop"); 
           p_iter != args.upper_bound("prop"); p_iter++) {
      propspecs.addSpecifier(p_iter->second);
    }   
}

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "traj" << "pbc" << "ref" << "atomsrmsd" << "atomsfit"
            << "skip" << "stride" << "human" << "precision" << "prop" 
            << "reftopo" << "dist" << "printatoms";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file(s), one for each trajectory group>\n";
  usage += "\n\t[@reftopo     <molecular topology file for the reference, default: first topology from @topo>]\n";
  usage += "\t@pbc          <boundary conditions> <gather type>\n";
  usage += "\t[@prop        <PropertySpecifier>]\n";
  usage += "\t[@atomsfit    <atoms to consider for fit>]\n";
  usage += "\t[@atomsrmsd   <atoms to consider for rmsd>\n";
  usage += "\t              (only required if different from atomsfit)]\n";
  usage += "\t[@skip        <skip frames at beginning>]\n";
  usage += "\t[@stride      <use only every step frame>]\n";
  usage += "\t[@human       (write the matrix in human readable form)]\n";
  usage += "\t[@precision   <number of digits in the matrix (default 4)>]\n";
  usage += "\t[@ref         <reference coordinates>]\n";
  usage += "\t[@dist        < binsize(default:0.02) >]\n";
  usage += "\t[@printatoms  < print list of selected atoms >]\n";
  usage += "\t@traj         <groups of trajectory files, separated by keyword 'newgroup'>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try {
    Arguments args(argc, argv, knowns, usage);
    
    int skip = args.getValue<int>("skip", false, 0);
    int stride = args.getValue<int>("stride", false, 1);

    // read the precision
    int ii = args.getValue<int>("precision", false, 4);
    int precision = 1;
    for (int i = 0; i < ii; ++i) {
      precision *= 10;
    }
    
    //print selected atoms or not
    bool printatoms=false, do_atomrmsd=false, do_proprmsd=false;
    if (args.count("printatoms") >= 0) printatoms=true;

    if (args.count("atomsrmsd") > 0 || args.count("atomsfit") > 0) do_atomrmsd=true;
    if (args.count("prop") > 0) do_proprmsd=true;
    if (!do_proprmsd && !do_atomrmsd)
          throw gromos::Exception(argv[0],
          "You need to specify either atoms (@atomsrmsd or @fitatoms) or properties (@prop)!\n");

    if (do_proprmsd && do_atomrmsd)
        throw gromos::Exception(argv[0],
          "specify either atoms (@atomsrmsd or @fitatoms) or properties (@prop), not both!\n");
   
    // create the vector to store the trajectory or properties
    vector< vector < Vec > > traj;
    vector< vector < Value > > props;
    
    // read reference topology
    System refSys;
    if (args.count("reftopo") == 1 && args.count("ref") > 0) {
      cout << "# reference topology: " << args["reftopo"] << endl;
      InTopology refit(args["reftopo"]);
      refSys=refit.system();
    } else if (args.count("reftopo") >1) {
          throw gromos::Exception(argv[0],
          "@reftopo only takes one argument\n");
    } else if (args.count("reftopo") >0 && args.count("ref")<=0) {
          throw gromos::Exception(argv[0],
          "@reftopo can only be used in combination with @ref\n"); 
    } else {
      cout << "# reference topology: " << args.lower_bound("topo")->second << endl;
      InTopology refit(args.lower_bound("topo")->second);
      refSys=refit.system();
    }
    
    // read reference coordinates...
    InG96 ic;
    if (args.count("ref") > 0) {
      cout << "# reference " << args["ref"]  << endl;
      ic.open(args["ref"]);
    } else {
      cout << "# reference is first frame" << endl;
      ic.open(args.lower_bound("traj")->second);
    }
    ic >> refSys;
    ic.close();
     
    // make sure no gathering method with a reference is used
    // because this will mess things up when using different topologies 
    if (args.count("topo") > 1 || args.count("reftopo") > 0)
    {
      Arguments::const_iterator iter = args.lower_bound("pbc"),
              to = args.upper_bound("pbc");
      for (; iter != to; ++iter) {
        if (iter->second == "refg" || iter->second == "gfit" || iter->second == "gtime" || iter->second == "grtime" || iter->second == "gltime") {
          throw gromos::Exception(argv[0],
            "ERROR: When using multiple topologies you can not use a gathering method that needs a reference or the previous frame's coordinates!");
        }
      }
    }
    
    // Parse boundary conditions without reading args
    // -> box information from the incoords
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(refSys, refSys, args);

    AtomSpecifier fitatoms(refSys), rmsdatoms(refSys), atomspecs(refSys);

    // get atom specifications for refsys from args
    if (do_atomrmsd) 
        getAtomSpecs(refSys, args, atomspecs, fitatoms, rmsdatoms);
    unsigned int atomspecnum=atomspecs.size();
    
    if (printatoms) {
        for (unsigned int i=0; i < atomspecs.size(); i++) {
            cout << atomspecs.toString(i) << " "<< atomspecs.mol(i)+1 <<":"<< atomspecs.resname(i) << atomspecs.resnum(i)+1<<":" << atomspecs.name(i) <<  endl;
        }
    }

    // if fitatoms != rmsdatoms keep lists of what to do
    vector<bool> fit_spec, rmsd_spec;
    if (atomspecs.size() != fitatoms.size() || atomspecs.size() != rmsdatoms.size()) {
      fit_spec.resize(atomspecs.size(), true);
      rmsd_spec.resize(atomspecs.size(), true);
      for (unsigned int i = 0; i < atomspecs.size(); ++i) {
        if (fitatoms.findAtom(atomspecs.mol(i), atomspecs.atom(i)) < 0)
          fit_spec[i] = false;
        if (rmsdatoms.findAtom(atomspecs.mol(i), atomspecs.atom(i)) < 0)
          rmsd_spec[i] = false;
      }
    }

    FastRotationalFit frf(fit_spec, rmsd_spec);

    // read in properties
    PropertyContainer propspecs(refSys, pbc);
    if (do_proprmsd) {
      Arguments::const_iterator iter = args.lower_bound("prop");
      Arguments::const_iterator to = args.upper_bound("prop");
      for (int i = 0; iter != to; iter++, ++i) {
        propspecs.addSpecifier(iter->second);
      }
      if (propspecs.empty()) {
        throw gromos::Exception(argv[0],
            "Empty property specifier (@prop).");
      } else {
        propspecs.calc();
        vector < Value > propvalues;
        for (PropertyContainer::const_iterator it = propspecs.begin(), to = propspecs.end();
             it != to; ++it) {
        propvalues.push_back((*it)->getValue(0));
      }
      props.push_back(propvalues);
      }
    }

    (*pbc.*gathmethod)();

    // calculate the centre of geometry of the relevant atoms
    Vec cog;
    for (unsigned int i = 0; i < fitatoms.size(); ++i) {
      cog += *fitatoms.coord(i);
    }
    cog /= fitatoms.size();

    // put it in the trajectory
    vector< Vec > frame(atomspecs.size());
    for (unsigned int i = 0; i < atomspecs.size(); ++i) {
      frame[i] = *atomspecs.coord(i) - cog;
    }
    traj.push_back(frame);
    
    
    // parse parameters for printing a distribution
    bool do_dist = false;
    double bin_size=0.02;
    map<int, int > bins;
    if (args.count("dist") >= 0) {
      do_dist = true;
      Arguments::const_iterator iter=args.lower_bound("dist");
      if(iter!=args.upper_bound("dist")){
	    std::istringstream is(iter->second);
	    is >> bin_size;
	    ++iter;
	    if (bin_size <= 0){
	      throw Arguments::Exception("distribution: invalid bin size" +
				     iter->second);
	    }
      }
    }
    
    
    // initialize first new system from first topology
    Arguments::const_iterator topo_iter = args.lower_bound("topo"),
            topo_to = args.upper_bound("topo");
            
    InTopology it(topo_iter->second);
    cout << "# " << topo_iter->second << endl;
    System sys(it.system());
    System refs(it.system());
    pbc = BoundaryParser::boundary(sys, args);
    gathmethod = args::GatherParser::parse(sys, refs, args);
    
    // get atom specifications for sys from args
    if (do_atomrmsd)
      getAtomSpecs(sys, args, atomspecs, fitatoms, rmsdatoms);
    if (printatoms) {
        for (unsigned int i=0; i < atomspecs.size(); i++) {
            cout << atomspecs.toString(i) << " "<< atomspecs.mol(i)+1 <<":"
                 << atomspecs.resname(i) << atomspecs.resnum(i)+1<<":" 
                 << atomspecs.name(i) <<  endl;
        }
    }
    if (atomspecs.size() != atomspecnum) {
      throw gromos::Exception("rmsdmat",
            "Atom selection does not refer to the same number of atoms in all the topologies.");
    }
    
    if (do_proprmsd) renewPropSpecs(sys, args, pbc, propspecs);    
    topo_iter++; 


    int framenum = 0;
    int propcnt = 0;
    // loop over all trajectories
    // put rmsdatom coordinates or property values in one long list each
    Arguments::const_iterator iter = args.lower_bound("traj"),
                               to = args.upper_bound("traj");
    for (; iter != to; ++iter) {
         
      if (iter->second=="newgroup") { 
        if (topo_iter == topo_to) {
          throw gromos::Exception("rmsdmat",
            "Not enough topologies for the number of trajectory groups.");   
        }
        InTopology it(topo_iter->second);
        cout << "# " <<topo_iter->second<< endl;
        sys=it.system();
        refs=it.system();
        pbc = BoundaryParser::boundary(sys, args);
        gathmethod = args::GatherParser::parse(sys, refs, args);
        if (do_atomrmsd) getAtomSpecs(sys, args, atomspecs, fitatoms, rmsdatoms);
        
        if (printatoms) {
          for (unsigned int i=0; i < atomspecs.size(); i++) {
            cout << atomspecs.toString(i) << " "<< atomspecs.mol(i)+1 <<":"<< atomspecs.resname(i) << atomspecs.resnum(i)+1<<":" << atomspecs.name(i) <<  endl;
          }
        }
        if (atomspecs.size() != atomspecnum) {
        throw gromos::Exception("rmsdmat",
            "Atom selection does not refer to the same number of atoms in all the topologies.");
        }
    
        if (do_proprmsd) renewPropSpecs(sys, args, pbc, propspecs);
        topo_iter++;;
        iter++;
        propcnt=0;
      }
      ic.open(iter->second);

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys;

        if (!((framenum - skip) % stride)) {

          //pbc call
          (*pbc.*gathmethod)();
          
          if (propspecs.size()) {
            // calculate props for this frame, this adds them to the propertycontainer
            propspecs.calc();
            vector < Value > propvalues;
            for (PropertyContainer::const_iterator it = propspecs.begin(), to = propspecs.end();
                  it != to; ++it) {
              propvalues.push_back((*it)->getValue(propcnt));
            }
            props.push_back(propvalues);
            propcnt++;
          } else {
            Vec cog;
            for (unsigned int i = 0; i < fitatoms.size(); ++i) {
              cog += *fitatoms.coord(i);
            }
            cog /= fitatoms.size();
            for (unsigned int i = 0; i < atomspecs.size(); ++i) {
              frame[i] = *atomspecs.coord(i) - cog;
            }
            int err = frf.fit(traj[0], frame);

            if (err) {
              ostringstream os;
              os << "Error while fitting to the reference structure\n"
                      << "Error code " << err << " in frame number " << framenum + 1;

              throw gromos::Exception(argv[0], os.str());
            }

            // store coordinates from sys in traj
            traj.push_back(frame);
          }
        }
        framenum++;
      }
      ic.close();
    }
    if (topo_iter != topo_to) {
      cout << "# WARNING: number of topologies is larger than the number of " <<
               "trajectory groups. Some topologies were not used." << endl;
    }

    // everything is in the thing; create the thingy
    // open a file
    ofstream fout;
    bool human = false;
    if (args.count("human") >= 0) {
      fout.open("RMSDMAT.dat");
      human = true;
    } else {
      fout.open("RMSDMAT.bin", ios::out | ios::binary);
    }

    // make a double loop
    int num = traj.size();
    if (props.size())
//      num = props.front()->num();
      num=props.size();
    Matrix rot(3, 3, 0);
    Matrix unit(3, 3, 0);
    for (size_t i = 0; i < 3; i++) unit(i, i) = 1;

    cout << "Read " << num << " out of " << framenum
            << " structures from trajectory" << endl;

    if (human) { // human format
      fout << "TITLE\n"
              << "\trmsd-matrix for " << num - 1 << " + 1 (ref) = "
              << num << " structures\n"
              << "END\n"
              << "RMSDMAT\n"
              << "# number of frames   skip   stride\n"
              << num << "\t" << skip << "\t" << stride << "\n"
              << "# precision\n"
              << precision << "\n";

      for (int i = 0; i < num; ++i) {
        for (int j = i + 1; j < num; ++j) {
          double rmsd = 0.0;
          if (props.size()) {
            rmsd = props_rmsd(props, i, j);
          } else {
            if (frf.fit(rot, traj[i], traj[j])) {
              cout << "error rotational fit on frames " << i + 1 << " and "
                      << j + 1
                      << "\nfitting to reference structure instead" << endl;
              rot = unit;
            }

            rmsd = frf.rmsd(rot, traj[i], traj[j]);
          }
          if (do_dist)
            update_bins(bins, bin_size, rmsd);

          rmsd *= precision;
          if (rmsd > std::numeric_limits<unsigned int>::max()) {
            std::cout << "frame " << i << " - " << j << " rmsd " << rmsd << endl;
            throw gromos::Exception(argv[0], "RMSD value is too big for a 'int'. Adjust @precision.");
          }
          fout << setw(8) << i
                  << setw(8) << j
                  << setw(8) << unsigned(rmsd)
                  << endl;
        }
      }
      fout << "END\n";
    } else { // binary format
      fout.write((char*) &num, sizeof (int));
      fout.write((char*) &skip, sizeof (int));
      fout.write((char*) &stride, sizeof (int));
      fout.write((char*) &precision, sizeof (int));

      if (precision < 1e5) { // small precision -> short format
        std::cout << "using 'unsigned short' as format" << std::endl;

        typedef unsigned short ushort;
        ushort irmsd;
        for (int i = 0; i < num; ++i) {
          for (int j = i + 1; j < num; ++j) {
            double rmsd = 0.0;
            if (props.size()) {
              rmsd = props_rmsd(props, i, j);
            } else {
              if (frf.fit(rot, traj[i], traj[j])) {
                cout << "error rotational fit on frames " << i + 1 << " and "
                        << j + 1
                        << "\nfitting to reference structure instead" << endl;
                rot = unit;
              }

              rmsd = frf.rmsd(rot, traj[i], traj[j]);
            }
            if (do_dist)
              update_bins(bins, bin_size, rmsd);

            rmsd *= precision;
            if (rmsd > std::numeric_limits<unsigned short>::max()) {
              std::cout << "frame " << i << " - " << j << " rmsd " << rmsd << endl;
              throw gromos::Exception(argv[0], "RMSD value is too big for a 'short'. Adjust @precision.");
            }
            irmsd = ushort(rmsd);
            fout.write((char*) &irmsd, sizeof (ushort));
          }
        }
      } else { // higher precision -> int format
        std::cout << "using 'unsigned int' as format" << std::endl;
        unsigned irmsd;
        for (int i = 0; i < num; ++i) {
          for (int j = i + 1; j < num; ++j) {
            double rmsd = 0.0;
            if (props.size()) {
              rmsd = props_rmsd(props, i, j);
            } else {
              if (frf.fit(rot, traj[i], traj[j])) {
                cout << "error rotational fit on frames " << i + 1 << " and "
                        << j + 1
                        << "\nfitting to reference structure instead" << endl;
                rot = unit;
              }

              rmsd = frf.rmsd(rot, traj[i], traj[j]);
            }
            if (do_dist)
              update_bins(bins, bin_size, rmsd);

            rmsd *= precision;
            if (rmsd > std::numeric_limits<unsigned int>::max())
              throw gromos::Exception(argv[0], "RMSD value is too big for a 'int'. Adjust @precision.");

            irmsd = unsigned(rmsd);
            fout.write((char*) &irmsd, sizeof (unsigned int));
          }
        }
      } // if precision
    } // if human
    if (do_dist)
       write_bins(bins, bin_size);
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
