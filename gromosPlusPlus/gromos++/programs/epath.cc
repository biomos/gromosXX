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
 * @file epath.cc
 * Calculate the electron-tunneling pathway in a protein
 */

/**
 * @page programs Program Documentation
 *
 * @anchor epath
 * @section epath calculates electron-tunneling pathways in a protein
 * @author @ref as @ref co
 * @date 06-10-11
 *
 * Program epath finds electron-tunneling pathways in proteins.
 * It uses Dijkstra's graph search algorithm, 
 * <I> Dijkstra, E. W. (1959), "A note on two problems in connexion with graphs",
 * Numerische Mathematik 1: 269–271 </I>, to find the pathway with the
 * highest product of the decay factors, corresponding to the "shortest path".
 *
 * The decay factor @f$\epsilon_{ij}@f$ for the electron transfer
 * between atoms i and j is calculated according to:
 * @f[ \epsilon_{ij} = A e^{B(r_{ij} - R)}  @f]
 * where @f$r_{ij}@f$ is the distance between the atoms and different 
 * parameters A, B and R are specified for jumps through covalent bonds,
 * hydrogen bonds and space, as described in <I>Beratan D. N. (1990),
 * "Electron-Tunneling Pathways in Ruthenated Proteins",
 * J. Am. Chem. Soc. 112: 7915-7921 </I>.
 *
 * For every atom of the system, the neighbouring atoms within a user-specified
 * cutoff are determined.
 * For every neighbouring atom its connectivity is classified as covalent, 
 * H-bonded or through space.
 * The decay factor for the neighbouring atoms is calculated using the
 * appropriate parameters and stored if it is higher than the decay factor that
 * was already stored from a previous cycle.
 * Additionally the jump type and the atom from where this jump occured are
 * stored. When all atoms have been visisted, the "shortest path" is backtraced
 * from the acceptor to the donor.
 *
 * The program outputs per default a pdb file containing the coordinates of
 * all atoms that have been part of a path throughout the different frames and
 * how often they have been part of a path as percentage in the b-factor column.
 *
 * Furthermore the program can produce a detailed output consiting off a summary
 * per frame to stdout and a pdb file containing the coordinates of the system
 * as well as the coordinates of the path corresponding to that frame linked
 * together in order to make the path visualisable.
 *
 * Also the program can output to a file a timeseries of the product of the
 * decay factor of the diffrent frames as well as its log, and the avarages
 * of both at the end.
 *
 * The parameters A, B and R are configurable via a parameter file, as well as
 * the parameters for the Hbond detection.
 *
 * All filenames are configurable.
 *
 * The program can produce a verbose output to stderr giving runtime details.
 * 
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@donor</td><td>&lt;@ref AtomSpecifier electron donor&gt; </td></tr>
 * <tr><td> \@acceptor</td><td>&lt;@ref AtomSpecifier electron acceptor&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input trajectory file(s)&gt; </td></tr>
 * <tr><td> [\@cutoff]</td><td>&lt;cut-off distance (nm, default: 0.6)&gt; </td></tr>
 * <tr><td> [\@param]</td><td>&lt;parameter file&gt; </td></tr>
 * <tr><td> [\@outfile]</td><td>&lt;name of the output file&gt; </td></tr>
 * <tr><td> [\@details]</td><td>&lt;detailed output&gt; </td></tr>
 * <tr><td> [\@detailsfile]</td><td>&lt;name of the output file of the details&gt; </td></tr>
 * <tr><td> [\@timseries]</td><td>&lt;print timeseries&gt; </td></tr>
 * <tr><td> [\@timseriesfile]</td><td>&lt;filename of the output file of the timeseries&gt; </td></tr>
 * <tr><td> [\@verbose]</td><td></td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 epath
    @topo      ex.top
    @pbc       r
    @donor     1:1
    @acceptor  2:1
    [@param     param.dat]
    [@cutoff    0.6]
    @traj      ex.tr
    @verbose
 @endverbatim
 *
 * <hr>
 */
#include <cmath>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/OutformatParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Property.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/Ginstream.h"
#include "../src/gromos/Exception.h"



using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace std;
using namespace utils;

class atom_map {
public:
    double decayparam;
    int previous;
    int jumptype;
    double jumpvalue;

    atom_map() {
        decayparam = 0.0;
        previous = -1;
        jumptype = -1;
        jumpvalue = 0.0;
    }
};

int main(int argc, char **argv) {

    // Argument list
    Argument_List knowns;
    knowns << "topo" << "pbc" << "donor" << "acceptor" << "traj" << "cutoff" << "param" << "details" << "detailsfile" << "timeseries" << "timeseriesfile" << "outfile" << "verbose";

    // usage
    string usage = "# " + string(argv[0]);
    usage += "\n\t@topo               <molecular topology file that you want to compare to>\n";
    usage += "\t@pbc                <boundary type> [<gathermethod>]\n";
    usage += "\t@donor              <electron donor> \n";
    usage += "\t@acceptor           <electron acceptor>\n";
    usage += "\t@traj               <trajectories>\n";
    usage += "\t[@cutoff            <cutoff (default 0.6)>]\n";
    usage += "\t[@param             <parameter file>]\n";
    usage += "\t[@outfile           <filename of the outputfile>]\n";
    usage += "\t[@details           <print out a detailed summary and a detailed pdb of the structure and the paths>]\n";
    usage += "\t[@detailsfile       <filename of the pdb file containing the structure and the paths>]\n";
    usage += "\t[@timeseries        <print out a timeseries of the product of decay factors>]\n";
    usage += "\t[@timeseriesfile    <filename of the timeseries file>]\n";
    usage += "\t[@verbose           <verbose runtime output>]\n";

    // try to do something
    try {
        ////////////////
        // start code //
        ////////////////

        Arguments args(argc, argv, knowns, usage);

        ///// debug switch /////
        bool debug = false;
        if (args.count("verbose") >= 0) {
            debug = true;
        }
        ////////////////////////
        //// detail switch /////
        bool details = false;
        if (args.count("details") >= 0) {
            details = true;
        }
        ////////////////////////////
        //// timeseries switch /////
        bool times = false;
        if (args.count("timeseries") >= 0) {
            times = true;
        }

        ////////////////////////////////////////////
        ///// Sanity checks /////

        // do we have arguments?
        if (!(args.count("topo") >= 0)) {
            throw gromos::Exception(string(argv[0]),
                    "no topology file given");
        }

        if (!(args.count("pbc") >= 0)) {
            throw gromos::Exception(string(argv[0]),
                    "no boundary type given");
        }

        if (!(args.count("donor") >= 0)) {
            throw gromos::Exception(string(argv[0]),
                    "no electron donor given");
        }

        if (!(args.count("acceptor") >= 0)) {
            throw gromos::Exception(string(argv[0]),
                    "no electron acceptor given");
        }

        if (!(args.count("traj") >= 0)) {
            throw gromos::Exception(string(argv[0]),
                    "no trajectory file given");
        }

        // declare filenames
        string outfile;
        string detailsfile;
        string timesfile;

        if (args.count("outfile") >= 0) {
            outfile = args["outfile"];
        } else {
            outfile = "et_occurrence.pdb";
        }

        if (args.count("detailsfile") >= 0) {
            detailsfile = args["detailsfile"];
        } else {
            detailsfile = "et_paths.pdb";
        }

        if (args.count("timeseriesfile") >= 0) {
            timesfile = args["timeseriesfile"];
        } else {
            timesfile = "et_decay.out";
        }

        // declare parameter variables with default values //
        // cutoff
        double cutoff = 0.6;
        //covalent
        double Ac = 0.6;
        double Bc = 0;
        double Rc = 0;
        //hbond
        double Ah = 0.36;
        double Bh = -17;
        double Rh = 0.28;
        //space
        double As = 0.6;
        double Bs = -17;
        double Rs = 0.14;
        //masses
        vector<double> masses;
        masses.push_back(15.99940);
        masses.push_back(14.00670);
        masses.push_back(32.06000);
        // hbparams
        double hbmaxdist = 0.25;
        double hbminangle = 135.0;


        // actually filll them from a file//
        // cutoff
        if (args.count("cutoff") >= 0) {
            cutoff = args.getValue<double>("cutoff");
        }
        // parmeters for the hbond and the decay calculation
        if (args.count("param") >= 0) {

            if (debug) {
                cerr << "Found parameter file " << args["param"] << endl;
            }

            // read in block
            Ginstream param(args["param"]);
            // vector to hold a block
            vector<string> linebuffer;
            // vector to hold all blocks
            vector<vector<string> > blockbuffer;
            // block counter
            while (!param.stream().eof()) {
                param.getblock(linebuffer);
                blockbuffer.push_back(linebuffer);
            }
            for (unsigned int i = 0; i < blockbuffer.size()-1; i++) {
                if (blockbuffer[i][0] == "HBONDPARAM") {

                    if (debug) {
                        cerr << "     Found HBONDPARAM block" << endl;
                    }
                    if (blockbuffer[i].size() != 3) {
                        throw gromos::Exception(string(argv[0]), "bad format for HBONDPARAM block");
                    }
                    // vector to hold the params
                    vector<double> hbparam;
                    stringstream ss(blockbuffer[i][1]);
                    double dum;

                    if (!(ss >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read HB max distance");
                    }
                    hbparam.push_back(dum);

                    if (!(ss >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read HB min angle");
                    }
                    hbparam.push_back(dum);

                    // now set the values
                    hbmaxdist = hbparam[0];
                    hbminangle = hbparam[1];

                } else if (blockbuffer[i][0] == "ACCEPTORMASS") {

                    if (debug) {
                        cerr << "     Found ACCEPTORMASS block" << endl;
                    }
                    // reset vector masses
                    masses.clear();
                    // loop over it
                    for (unsigned int j = 1; j < blockbuffer[i].size() - 1; j++) {
                        stringstream ss(blockbuffer[i][j]);
                        double dum;
                        if (!(ss >> dum))
                            throw gromos::Exception(string(argv[0]), "could not read mass");
                        masses.push_back(dum);
                    }
                } else if (blockbuffer[i][0] == "DECAYPARAM") {

                    if (debug) {
                        cerr << "     Found DECAYPARAM block" << endl;
                    }
                    // 3 lines 3 values per line
                    double dum;
                    // covalent
                    stringstream ss1(blockbuffer[i][1]);
                    if (!(ss1 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Ac");
                    }
                    Ac = dum;
                    if (!(ss1 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Bc");
                    }
                    Bc = dum;
                    if (!(ss1 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Rc");
                    }
                    Rc = dum;
                    // hb
                    stringstream ss2(blockbuffer[i][2]);
                    if (!(ss2 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Ah");
                    }
                    Ah = dum;
                    if (!(ss2 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Bh");
                    }
                    Bh = dum;
                    if (!(ss2 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Rh");
                    }
                    Rh = dum;
                    // space
                    stringstream ss3(blockbuffer[i][3]);
                    if (!(ss3 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read As");
                    }
                    As = dum;
                    if (!(ss3 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Bs");
                    }
                    Bs = dum;
                    if (!(ss3 >> dum)) {
                        throw gromos::Exception(string(argv[0]), "could not read Rs");
                    }
                    Rs = dum;
                    // done
                } else {
                    throw gromos::Exception(string(argv[0]), "Unknown block: " + blockbuffer[i][0]);
                }
            }
        }
        /////                                  /////
        ////////////////////////////////////////////
        ///// Initialise /////

        // read topology
        InTopology it(args["topo"]);
        System sys(it.system());
        System refSys(it.system());

        ///// check whether start and stop are hydrogens and add them to a specifier /////
        AtomSpecifier start(sys);
        start.addSpecifier(args["donor"]);
        if (sys.mol(start.mol(0)).topology().atom(start.atom(0)).isH()) {
            throw gromos::Exception("ASTAR",
                    "electron donor cannot be a hydrogen");
        }

        AtomSpecifier stop(sys);
        stop.addSpecifier(args["acceptor"]);
        if (sys.mol(stop.mol(0)).topology().atom(stop.atom(0)).isH()) {
            throw gromos::Exception("ASTAR",
                    "electron acceptor cannot be a hydrogen");
        }
        //////////////////////////////////////////////////////

        if (debug) {
            cerr << "Read in topology from " << args["topo"] << endl;
        }

        // parse boundary conditions
        Boundary *pbc = BoundaryParser::boundary(sys, args);
        //parse gather method
        Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

        //////////////////////////////////
        ///// Create atom specifier /////

        // create specifierer "allatoms" which contains all solutes
        AtomSpecifier allatoms(sys);

        {
            AtomSpecifier withH(sys);
            // put in all molecules
            withH.addSpecifier("a:a");

            // remove H
            for (unsigned int i = 0; i < withH.size(); i++) {
                if (!sys.mol(withH.mol(i)).topology().atom(withH.atom(i)).isH()) {
                    allatoms.addAtom(withH.mol(i), withH.atom(i));
                }
            }

        }

        if (debug) {
            cerr << "The system has " << allatoms.size() << " heavy atoms" << endl;
            cerr << "##### Parameters used for the calculation" << endl;
            cerr << "# cutoff = " << cutoff << endl;
            cerr << "# HB max distance = " << hbmaxdist << endl;
            cerr << "# HB min angle = " << hbminangle << endl;
            cerr << "# Using masses: ";
            for (unsigned int i = 0; i < masses.size(); i++) {
                cerr << masses[i] << " ";
            }
            cerr << "as potential hydrogenbond acceptors" << endl;
            cerr << setw(13) << "# covalent_A = " << setw(6) << Ac << setw(16) << " ; covalent_B = " << setw(6) << Bc << setw(16) << " ;  covalent_R = " << setw(6) << Rc << endl;
            cerr << setw(13) << "# hbond_A    = " << setw(6) << Ah << setw(16) << " ; hbond_B    = " << setw(6) << Bh << setw(16) << " ;  hbond_R    = " << setw(6) << Rh << endl;
            cerr << setw(13) << "# space_A    = " << setw(6) << As << setw(16) << " ; space_B    = " << setw(6) << Bs << setw(16) << " ;  space_R    = " << setw(6) << Rs << endl;
            cerr << "#####" << endl;
        }
        /////////////////////////////////

        // prepare occurence vector
        vector<int> occ(allatoms.size(), 0);

        // prepair coordinatin input
        InG96 ic;
        // fancy output
        if (details) {
            cout << "Summary output of the electron transfer path calculation" << endl;
            cout << "Electron donor: " << start.mol(0) + 1 << ":" << start.atom(0) + 1 << endl;
            cout << "Electron acceptor: " << stop.mol(0) + 1 << ":" << stop.atom(0) + 1 << endl;
            cout << endl;
        }
        // fancy title only once
        bool title = true;
        // table header once
        bool header = true;

        // prepare fancy pdb output file, this makes sure we have an empty file
        if (details) {
            ofstream pdbfile;
            pdbfile.open(detailsfile.c_str());
            pdbfile << "";
            pdbfile.close();
        }
        if (times) {
            // prepare fancy timeseries output file, this makes sure we have an empty file
            ofstream timefile;
            timefile.open(timesfile.c_str());
            timefile << "# Timeseries of the product of the decay factors, k" << endl;
            timefile << "#" << setw(10) << "time" << setw(20) << "log_k" << setw(20) << "k" << endl;
            timefile.close();
        }

         // initialise some helper variables
        int globframecounter = 1;
        vector<double> k_decay;
        vector<double> log_k_decay;
        ////// Loop over trajectory //////
        for (Arguments::const_iterator
            iter = args.lower_bound("traj"),
                to = args.upper_bound("traj");
                iter != to; ++iter) {

            ic.open(iter->second);
            ic.select("SOLUTE");

            while (!ic.eof()) {
                ic >> sys;
                (*pbc.*gathmethod)();

                if (debug) {
                    cerr << "Read in coordinates from " << iter->second << endl;
                    cerr << "Frame number: " << globframecounter << endl;
                }

                // initialise map vector
                vector<atom_map> map(allatoms.size());

                // create specifier current atom
                AtomSpecifier current(sys);
                // fill it with the start atom
                current.addSpecifier(args["donor"]);

                // create set visited and insert the starting point
                set<int> visited;
                visited.insert(allatoms.findAtom(current.mol(0), current.atom(0)));

                // insert decay for starting point default 1
                map[allatoms.findAtom(current.mol(0), current.atom(0))].decayparam = 1;

                ////////////////////
                ///// BIG LOOP /////
                ////////////////////
                if (debug) {
                    cerr << "Calculating decay factors..." << endl;
                }

                bool domore = true;
                while (visited.size() != allatoms.size() && domore) {
                    ///////////////////

                    // create specifier pairlist
                    AtomSpecifier pairs(sys);
                    // create speicifier covalent
                    AtomSpecifier cov(sys);
                    // create speicifier hbond
                    AtomSpecifier hb(sys);


                    // create pairlist of atoms around the current atom
                    SimplePairlist around(sys, *pbc, cutoff);
                    // set the starting point
                    around.setAtom(current.mol(0), current.atom(0));
                    around.setType("ATOMIC");
                    around.calc();

                    Neighbours nb(sys, current.mol(0), current.atom(0));
                    // put the pairlist into the pair specifier
                    for (unsigned int i = 0; i < around.size(); ++i) {
                        // check if atom is visited
                        if (!sys.mol(around.mol(i)).topology().atom(around.atom(i)).isH()
                                && !visited.count(allatoms.findAtom(around.mol(i), around.atom(i)))) {
                            // then add it
                            pairs.addAtom(around.mol(i), around.atom(i));
                        }
                    }

                    // determine bond typ: covalent, hbond, space
                    // covalent //

                    // put the neigbours in the covalent specifier
                    for (unsigned int i = 0; i < nb.size(); i++) {
                        if (!sys.mol(current.mol(0)).topology().atom(nb[i]).isH()
                                && !visited.count(allatoms.findAtom(current.mol(0), nb[i]))) {
                            cov.addAtom(current.mol(0), nb[i]);
                            pairs.removeAtom(current.mol(0), nb[i]);
                        }
                    }

                    // first case: current is the donor
                    for (unsigned int i = 0; i < nb.size(); i++) {
                        if (sys.mol(current.mol(0)).topology().atom(nb[i]).isH()) {
                            // now search in the pairs list for a suitable acceptor
                            for (unsigned int j = 0; j < pairs.size(); j++) {
                                // loop over masses
                                bool couldbeh = false;
                                for (unsigned int h = 0; h < masses.size(); h++) {
                                    if (pairs.mass(j) == masses[h]) {
                                        couldbeh = true;
                                    }
                                }
                                if (couldbeh) {
                                    utils::HBProperty hbprop(sys, pbc);
                                    stringstream istr;
                                    istr << current.toString(0) << ";" << current.mol(0) + 1 << ":" << nb[i] + 1 << ";"
                                            << pairs.toString(j) << "%" << hbmaxdist << "%" << hbminangle;
                                    vector<string> pprop;
                                    pprop.push_back(istr.str());
                                    hbprop.parse(pprop, 1);

                                    if (hbprop.calc().scalar() != 0) {
                                        hb.addAtom(pairs.mol(j), pairs.atom(j));
                                    }
                                }
                            }
                        }
                    }
                    // second case: current is the acceptor
                    // loop over masses
                    bool couldbeh = false;
                    for (unsigned int h = 0; h < masses.size(); h++) {
                        if (current.mass(0) == masses[h]) {
                            couldbeh = true;
                        }
                    }
                    if (couldbeh) {
                        for (unsigned int i = 0; i < pairs.size(); i++) {
                            Neighbours pnb(sys, pairs.mol(i), pairs.atom(i));
                            for (unsigned int j = 0; j < pnb.size(); j++) {
                                if (sys.mol(pairs.mol(i)).topology().atom(pnb[j]).isH()) {
                                    utils::HBProperty hbprop(sys, pbc);
                                    stringstream istr;
                                    istr << pairs.toString(i) << ";" << pairs.mol(i) + 1 << ":" << pnb[j] + 1 << ";"
                                            << current.toString(0) << "%" << hbmaxdist << "%" << hbminangle;
                                    vector<string> pprop;
                                    pprop.push_back(istr.str());
                                    hbprop.parse(pprop, 1);

                                    if (hbprop.calc().scalar() != 0) {
                                        hb.addAtom(pairs.mol(i), pairs.atom(i));
                                    }
                                }
                            }
                        }
                    }

                    // remove the hb from the pairs
                    for (unsigned int i = 0; i < hb.size(); i++) {
                        pairs.removeAtom(hb.mol(i), hb.atom(i));
                    }

                    // space - everything that is neither covalent nor hbond
                    // this is what is left in the pairs atom specifier

                    // calculate the decay parameters ec, eh, es for every item in pair
                    // global values
                    double ec;
                    double eh;
                    double es;

                    // local values
                    double ec_l;
                    double eh_l;
                    double es_l;

                    double dist;
                    double currentdecay = map[allatoms.findAtom(current.mol(0), current.atom(0))].decayparam;

                    // jump
                    int jtype;

                    // covalent
                    for (unsigned int i = 0; i < cov.size(); i++) {
                        // current in the donor or acceptor range?
                        if (((start.findAtom(cov.mol(i), cov.atom(i)) != -1) && (start.findAtom(current.mol(0), current.atom(0)) != -1))|| ((stop.findAtom(cov.mol(i), cov.atom(i)) != -1) && (stop.findAtom(current.mol(0), current.atom(0)) != -1))) {
                            ec = currentdecay * 1.00;
                            ec_l = 1.00;
                            jtype = 3;
                        } else {
                            cov.pos(i) = pbc->nearestImage(current.pos(0), cov.pos(i), sys.box());
                            dist = (current.pos(0) - cov.pos(i)).abs();
                            // calculatec ec
                            ec = currentdecay * Ac * exp(Bc * (dist - Rc));
                            ec_l = Ac * exp(Bc * (dist - Rc));
                            jtype = 0;
                        }
                        // store it if it is bigger
                        int index = allatoms.findAtom(cov.mol(i), cov.atom(i));
                        if (index != -1 && map[index].decayparam < ec) {
                            map[index].decayparam = ec;
                            map[index].previous = allatoms.findAtom(current.mol(0), current.atom(0));
                            map[index].jumptype = jtype;
                            map[index].jumpvalue = ec_l;
                        }
                    }

                    // hbond
                    for (unsigned int i = 0; i < hb.size(); i++) {
                        // current in the donor or acceptor range?
                        if (((start.findAtom(hb.mol(i), hb.atom(i)) != -1) && (start.findAtom(current.mol(0), current.atom(0)) != -1))|| ((stop.findAtom(hb.mol(i), hb.atom(i)) != -1) && (stop.findAtom(current.mol(0), current.atom(0)) != -1))) {
                            eh = currentdecay * 1.00;
                            eh_l = 1.00;
                            jtype = 3;
                        } else {
                            hb.pos(i) = pbc->nearestImage(current.pos(0), hb.pos(i), sys.box());
                            dist = (current.pos(0) - hb.pos(i)).abs();
                            // calculatec eh
                            eh = currentdecay * Ah * exp(Bh * (dist - Rh));
                            eh_l = Ah * exp(Bh * (dist - Rh));
                            jtype = 1;
                        }
                        // store it if it is bigger
                        int index = allatoms.findAtom(hb.mol(i), hb.atom(i));
                        if (index != -1 && map[index].decayparam < eh) {
                            map[index].decayparam = eh;
                            map[index].previous = allatoms.findAtom(current.mol(0), current.atom(0));
                            map[index].jumptype = jtype;
                            map[index].jumpvalue = eh_l;
                        }
                    }

                    // space
                    for (unsigned int i = 0; i < pairs.size(); i++) {
                        // current in the donor or acceptor range?
                        if (((start.findAtom(pairs.mol(i), pairs.atom(i)) != -1) && (start.findAtom(current.mol(0), current.atom(0)) != -1))|| ((stop.findAtom(pairs.mol(i), pairs.atom(i)) != -1) && (stop.findAtom(current.mol(0), current.atom(0)) != -1))) {
                            es = currentdecay * 1.00;
                            es_l = 1.00;
                            jtype = 3;
                        } else {
                            pairs.pos(i) = pbc->nearestImage(current.pos(0), pairs.pos(i), sys.box());
                            dist = (current.pos(0) - pairs.pos(i)).abs();
                            // calculatec es
                            es = currentdecay * As * exp(Bs * (dist - Rs));
                            es_l = As * exp(Bs * (dist - Rs));
                            jtype = 2;
                        }
                        // store it if it is bigger
                        int index = allatoms.findAtom(pairs.mol(i), pairs.atom(i));
                        if (index != -1 && map[index].decayparam < es) {
                            map[index].decayparam = es;
                            map[index].previous = allatoms.findAtom(current.mol(0), current.atom(0));
                            map[index].jumptype = jtype;
                            map[index].jumpvalue = es_l;
                        }
                    }

                    double max = 0.0;
                    int next = -1;
                    for (unsigned int i = 0; i < map.size(); i++) {
                        if (!visited.count(i) && map[i].decayparam > max) {
                            max = map[i].decayparam;
                            next = i;
                        }
                    }
                    if (next == -1) {
                        domore = false;

                        if (debug) {
                            int outofcutoff = allatoms.size() - visited.size();
                            cerr << endl;
                            cerr << "   " << outofcutoff << " atoms are not within the cutoff distance of a possible path";
                        }

                        break;
                    }

                    // break if acceptor reached

                    if (stop.atom(0) == allatoms.atom(next) &&
                            stop.mol(0) == allatoms.mol(next)) {
                        domore = false;

                     if (debug) {
                        cerr << "   " << "breaking, acceptor reached: " << stop.atom(0) << endl;
                     }
                        break;
                    }

                    // Clean up and remove current from allatoms
                    current.clear();
                    current.addAtom(allatoms.mol(next), allatoms.atom(next));

                    // add next to visisted
                    visited.insert(allatoms.findAtom(current.mol(0), current.atom(0)));

                    // output progress
                    if (debug) {
                        cerr << "   " << visited.size() << "/" << allatoms.size() << " atoms" << "\r";
                    }


                    ///// END BIG LOOP /////
                }

                // regather
                (*pbc.*gathmethod)();

                ////////////////////////
                if (debug) {
                    cerr << endl;
                    cerr << "...done" << endl;
                }

                //////////////////////////////////////
                //          FANCY OUTPUT            //
                //////////////////////////////////////
                if (details) {
                    cout << "## Frame " << globframecounter << endl;
                    cout << "###################" << endl;
                    cout << "Path:" << endl;
                }

                // calculate path
                if (debug) {
                    cerr << "Calculating path...";
                }
                // create path AtomSpeficier
                AtomSpecifier path(sys);
                path.addAtom(stop.mol(0), stop.atom(0));
                // ending point
                int pindex = allatoms.findAtom(path.mol(0), path.atom(0));
                // decay of this path
                double current_decay = map[pindex].decayparam;

                // safe it for the average
                k_decay.push_back(current_decay);
                log_k_decay.push_back(log(current_decay));

                if (details) {
                    cout << "  Product of the decay factors of this path: " << current_decay << endl;
                }
                if (times) {
                    // timeseries out
                    ofstream timefile;
                    timefile.open(timesfile.c_str(), ofstream::app);
                    timefile << setw(11) << globframecounter << setw(20) << log(current_decay) << setw(20) << current_decay << endl;
                    timefile.close();
                }
                // starting point
                int pend = allatoms.findAtom(start.mol(0), start.atom(0));
                // the start atom occurs also in the path :)
                occ[pend]++;
                // number of covalent jumps
                int countcov = 0;
                // number of hbond jumps
                int counthb = 0;
                // number of jumps through space
                int countspace = 0;

                // backtrace of the path
                while (pindex != pend) {
                    if (map[pindex].previous == -1) {
                        //bool nopath = true;
                        break;
                    }

                    if (details) {
                        if (header) {
                            cout << "|" << setfill ('-') << setw(54) << "|" << setfill (' ') << endl;
                            cout << "|" << setw(10) << "   From   " << "|" << setw(10) << "    To    " << "|" << setw(15) << "      Type     " << "|" << setw(15) << "     Value     " << "|" << endl;
                            cout << "|" << setfill ('-') << setw(54) << "|" << setfill (' ') << endl;
                            header = false;
                        }
                        cout << "|" << setw(9) << allatoms.toString(pindex) << " |" << setw(9) << allatoms.toString(map[pindex].previous) << " |" ;
                        if (map[pindex].jumptype == 0) {
                            cout << setw(15) << " covalent bond ";
                        }
                        if (map[pindex].jumptype == 1) {
                            cout << setw(15) << " hydrogen bond ";
                        }
                        if (map[pindex].jumptype == 2) {
                            cout << setw(15) << "     space     ";
                        }
                        if (map[pindex].jumptype == 3) {
                            cout << setw(15) << "  inside group ";
                        }
                        cout << "|" << setw(14) << map[pindex].jumpvalue << " |";
                        cout << endl;
                    }

                    path.addAtomStrict(allatoms.mol(map[pindex].previous), allatoms.atom(map[pindex].previous));
                    if (map[pindex].jumptype == 0) countcov++;
                    if (map[pindex].jumptype == 1) counthb++;
                    if (map[pindex].jumptype == 2) countspace++;
                    // add an occurence
                    occ[pindex]++;
                    // next one
                    pindex = map[pindex].previous;
                }

                if (debug) {
                    cerr << "done" << endl;
                }
                // done with path backtracing
                // reset header
                header = true;

                // if there is no path...
                // path length is 1 atom algorithm can handle it

                if (details) {
                    // table footer
                    cout << "|" << setfill ('-') << setw(54) << "|" << setfill (' ') << endl;

                    // continue output
                    cout << endl;
                    cout << "This path goes over: " << endl;
                    cout << "     " << countcov << " covalent bonds" << endl;
                    cout << "     " << counthb << " hydrogen bonds" << endl;
                    cout << "     " << countspace << " jumps through space" << endl;
                    cout << "This path has a total of " << countcov + counthb + countspace + 1 << " atoms" << endl;
                }

                if (details) {
                    cout << "###################" << endl;
                }

                // PDB OUT
                if (details) {
                    // open output file
                    ofstream pdbfile;
                    pdbfile.open(detailsfile.c_str(), ofstream::app);
                    // output title once
                    if (title) {
                        pdbfile.setf(ios::fixed, ios::floatfield);
                        pdbfile.setf(ios::unitbuf);
                        pdbfile.precision(3);
                        pdbfile << "TITLE";
                        pdbfile.setf(ios::right, ios::adjustfield);
                        pdbfile << "  " << setw(7) << "  This file shows the molecule strucutre in chain A and the corresponding calculated electron transfer path in chain B" << endl;
                        title = false;
                    }

                    // write summary as remark 3
                   pdbfile.setf(ios::fixed, ios::floatfield);
                   pdbfile.setf(ios::unitbuf);
                   pdbfile.precision(3);
                   pdbfile << "REMARK   ";
                   pdbfile.setf(ios::right, ios::adjustfield);
                   pdbfile << "3";
                   pdbfile.setf(ios::left, ios::adjustfield);
                   pdbfile << "  " << setw(4) << "This corresponds to frame: " << globframecounter << endl;

                    // data
                    // start a frame
                    pdbfile << "MODEL" << endl;

                    // output allatoms as chain A
                    int resoffset = 0;
                    int atomoffset = 0;
                    int m = 0;
                    int atomsoff = 0;
                    for (unsigned int i = 0; i < allatoms.size(); i++) {

                        if (allatoms.mol(i) > m) {
                            resoffset += sys.mol(allatoms.mol(i) - 1).topology().numRes();
                            atomoffset += sys.mol(allatoms.mol(i) - 1).topology().numAtoms();
                            m++;
                        }
                        Vec pos = allatoms.pos(i);
                        // constrcut pdb line
                        pdbfile.setf(ios::fixed, ios::floatfield);
                        pdbfile.setf(ios::unitbuf);
                        pdbfile.precision(3);
                        pdbfile << "ATOM  ";
                        pdbfile.setf(ios::right, ios::adjustfield);
                        pdbfile << setw(5) << atomoffset + allatoms.atom(i) + 1;
                        pdbfile.setf(ios::left, ios::adjustfield);
                        pdbfile << "  " << setw(4) << allatoms.name(i).substr(0, 4);
                        pdbfile << setw(3) << allatoms.resname(i).substr(0, 3);
                        pdbfile << setw(2) << " A"; // chain number
                        pdbfile.setf(ios::right, ios::adjustfield);
                        pdbfile << setw(4) << resoffset + allatoms.resnum(i) + 1 << "    "
                                << setw(8) << pos[0]*10
                                << setw(8) << pos[1]*10
                                << setw(8) << pos[2]*10
                                << "  " << "1.00" << " ";
                        pdbfile.precision(2);
                        pdbfile << " " // some thing later
                                << endl;
                        // end pdb line
                        // clear vector
                    }
                    // end atoms
                    pdbfile << "TER" << endl;


                    // Chain B
                    for (unsigned int i = 0; i < path.size(); i++) {
                        atomsoff = 0;
                        // complicated way of doing the resoffset
                        resoffset = 0;
                        if (allatoms.mol(allatoms.findAtom(path.mol(i), path.atom(i))) > 0) {
                            int j = allatoms.mol(allatoms.findAtom(path.mol(i), path.atom(i))) - 1;
                            while (j >= 0) {
                                resoffset += sys.mol(j).topology().numRes();
                                atomsoff += sys.mol(j).topology().numAtoms();
                                j = j - 1;
                            }
                        }
                        // path atomofset is last atom of the system
                        atomoffset = allatoms.gromosAtom(allatoms.size()-1) + 1;
                        atomoffset += atomsoff;
                        // done with offset calc

                        Vec pos = path.pos(i);
                        // constrcut pdb line
                        pdbfile.setf(ios::fixed, ios::floatfield);
                        pdbfile.setf(ios::unitbuf);
                        pdbfile.precision(3);
                        pdbfile << "ATOM  ";
                        pdbfile.setf(ios::right, ios::adjustfield);
                        pdbfile << setw(5) << allatoms.atom(allatoms.findAtom(path.mol(i), path.atom(i))) + atomoffset + 1;
                        pdbfile.setf(ios::left, ios::adjustfield);
                        pdbfile << "  " << setw(4) << path.name(i).substr(0, 4);
                        pdbfile << setw(3) << path.resname(i).substr(0, 3);
                        pdbfile << setw(2) << " B"; // chain number
                        pdbfile.setf(ios::right, ios::adjustfield);
                        pdbfile << setw(4) << resoffset + allatoms.resnum(allatoms.findAtom(path.mol(i), path.atom(i))) + 1 << "    "
                                << setw(8) << pos[0]*10
                                << setw(8) << pos[1]*10
                                << setw(8) << pos[2]*10
                                << "  " << "1.00" << "  ";
                        pdbfile.precision(2);
                        pdbfile << " " // some thing later
                                << endl;
                        // end pdb line
                    }
                    // end atoms
                    pdbfile << "TER" << endl;

                    // connect from OutPdb.cc for the molecule
                    for (int i = 0, offatom = 1; i < sys.numMolecules(); ++i) {
                        BondIterator bit(sys.mol(i).topology());
                        for (int count = 0; bit; ++bit) {
                            pdbfile << setw(6) << "CONECT"
                                    << setw(5) << bit()[0] + offatom
                                    << setw(5) << bit()[1] + offatom
                                    << endl;


                            ++count;
                        }
                        offatom += sys.mol(i).numAtoms();
                    }
                    // end connect

                    // connect block forthe path
                    // two offsets needed as i and i+1 can be two diffrent molecules
                    int atomoffset_connect1 = 0;
                    int atomoffset_connect1plus = 0;
                    int atomoffset_plus1 = 0;

                    for (unsigned int i = 0; i < path.size() - 1; i++) {
                        // complicated way of doing the offset
                        // reset
                        atomoffset_connect1 = 0;
                        atomoffset_connect1plus = 0;

                        if (allatoms.mol(allatoms.findAtom(path.mol(i), path.atom(i))) > 0) {
                            int j = allatoms.mol(allatoms.findAtom(path.mol(i), path.atom(i))) - 1;
                            while (j >= 0) {
                                atomoffset_connect1 += sys.mol(j).topology().numAtoms();
                                j = j - 1;
                            }
                        }
                        if (allatoms.mol(allatoms.findAtom(path.mol(i+1), path.atom(i+1))) > 0) {
                            int j = allatoms.mol(allatoms.findAtom(path.mol(i+1), path.atom(i+1))) - 1;
                            while (j >= 0) {
                                atomoffset_connect1plus += sys.mol(j).topology().numAtoms();
                                j = j - 1;
                            }
                        }
                        // path atomofset is last atom of the system
                        atomoffset = allatoms.gromosAtom(allatoms.size()-1) + 1;
                        atomoffset_plus1 = allatoms.gromosAtom(allatoms.size()-1) + 1;
                        atomoffset += atomoffset_connect1;
                        atomoffset_plus1 += atomoffset_connect1plus;
                        // done with offset calc

                        pdbfile << setw(6) << "CONECT"
                                << setw(5) << allatoms.atom(allatoms.findAtom(path.mol(i), path.atom(i))) + atomoffset + 1
                                << setw(5) << allatoms.atom(allatoms.findAtom(path.mol(i + 1), path.atom(i + 1))) + atomoffset_plus1 + 1
                                << endl;
                    }
                    // end frame
                    pdbfile << "ENDMDL" << endl;
                    //close
                    pdbfile.close();
                }
                /////////////////////////////////////////////////
                //            END FANCY OUTPUT                 //
                /////////////////////////////////////////////////

                /////// End looping over trajectories ///////////
                globframecounter++;
            }
            ic.close();
        }
        ////////////////////////////////////////////////
        // simple output
        // open file stream
        ofstream simplefile;
        simplefile.open(outfile.c_str());
        // output title
        simplefile.setf(ios::fixed, ios::floatfield);
        simplefile.setf(ios::unitbuf);
        simplefile.precision(3);
        simplefile << "TITLE";
        simplefile.setf(ios::right, ios::adjustfield);
        simplefile << "  " << setw(7) << "  This file indicates the percantage of how often an atom is part of a path in the B-factor column" << endl;
        // write summary as remark 3
        simplefile.setf(ios::fixed, ios::floatfield);
        simplefile.setf(ios::unitbuf);
        simplefile.precision(3);
        simplefile << "REMARK   ";
        simplefile.setf(ios::right, ios::adjustfield);
        simplefile << "3";
        simplefile.setf(ios::left, ios::adjustfield);
        simplefile << "  " << setw(4) << "Electron donor: " << start.mol(0) + 1 << ":" << start.atom(0) + 1 << endl;

        simplefile.setf(ios::fixed, ios::floatfield);
        simplefile.setf(ios::unitbuf);
        simplefile.precision(3);
        simplefile << "REMARK   ";
        simplefile.setf(ios::right, ios::adjustfield);
        simplefile << "3";
        simplefile.setf(ios::left, ios::adjustfield);
        simplefile << "  " << setw(4) << "Electron acceptor: " << stop.mol(0) + 1 << ":" << stop.atom(0) + 1 << endl;

        // count atoms that are part of a path
        int partof = 0;
        for (unsigned int i = 0; i < occ.size(); i++) {
            if (occ[i] != 0) {
                partof++;
            }
        }
        simplefile.setf(ios::fixed, ios::floatfield);
        simplefile.setf(ios::unitbuf);
        simplefile.precision(3);
        simplefile << "REMARK   ";
        simplefile.setf(ios::right, ios::adjustfield);
        simplefile << "3";
        simplefile.setf(ios::left, ios::adjustfield);
        simplefile << "  " << setw(4) << partof << " atoms were part of a possible path" << endl;

        // loop over atoms
        int resoffset = 0;
        int m = 0;
        for (unsigned int i = 0; i < occ.size(); i++) {

            if (allatoms.mol(i) > m) {
                resoffset += sys.mol(allatoms.mol(i) - 1).topology().numRes();
                m++;
            }

            if (occ[i] != 0) {

                Vec pos = allatoms.pos(i);
                // constrcut pdb line
                simplefile.setf(ios::fixed, ios::floatfield);
                simplefile.setf(ios::unitbuf);
                simplefile.precision(3);
                simplefile << "ATOM  ";
                simplefile.setf(ios::right, ios::adjustfield);
                simplefile << setw(5) << allatoms.atom(i) + 1;
                simplefile.setf(ios::left, ios::adjustfield);
                simplefile << "  " << setw(4) << allatoms.name(i).substr(0, 4);
                simplefile << setw(3) << allatoms.resname(i).substr(0, 3);
                simplefile << setw(2) << " A"; // chain number
                simplefile.setf(ios::right, ios::adjustfield);
                simplefile << setw(4) << resoffset + allatoms.resnum(i) + 1 << "    "
                        << setw(8) << pos[0]*10
                        << setw(8) << pos[1]*10
                        << setw(8) << pos[2]*10
                        << "  " << "1.00" << " ";
                simplefile.precision(2);
                //helper double
                double bf_out = double(occ[i])*100.00 / (double(globframecounter) - 1.00);
                simplefile << bf_out // b-factor column as %
                        << endl;
                // end pdb line
                //  cout << "Atom " << allatoms.toString(i) << " occurs in " << occ[i] << " paths" << endl;
            }
        }
        // end pdb
        simplefile << "TER" << endl;

        // close file stream
        simplefile.close();

        // end simple output

        // average over timeseries
        if (times) {
            //average over log_k_decay and k_decay vector
            double log_k_decay_av = 0.00;
            double k_decay_av = 0.00;

            for (unsigned int i = 0; i < log_k_decay.size() ; i++) {
                log_k_decay_av += log_k_decay[i];
            }
            log_k_decay_av = log_k_decay_av/log_k_decay.size();

            for (unsigned int i = 0; i < k_decay.size() ; i++) {
                k_decay_av += k_decay[i];
            }
            k_decay_av = k_decay_av/k_decay.size();

            // timeseries out
            ofstream timefile;
            timefile.open(timesfile.c_str(), ofstream::app);
            timefile << "#######################" << endl;
            timefile << "# Averages:" << setw(20) << "log_k" << setw(20) << "k" << endl;
            timefile << "           " << setw(20) << log_k_decay_av << setw(20) << k_decay_av << endl;
            timefile.close();
          }

        //////////////
        // end code //
        //////////////
    }    // catch all
    catch (const gromos::Exception &e) {
        cerr << e.what() << endl;
        exit(1);
    }
    return 0;
}


