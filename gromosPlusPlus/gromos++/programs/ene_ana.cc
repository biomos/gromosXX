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
 * @file ene_ana.cc
 * extracts time series from (energy) trajectory files
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ene_ana
 * @section ene_ana analyse (energy) trajectories
 * @author @ref mc @ref co
 * @date 26. 7. 2006
 *
 * GROMOS can write energies, free-energy derivatives and block averages of 
 * these to separate trajectory files for later analysis. Program ene_ana 
 * extracts individual values from such files and can perform simple 
 * mathematical operations on them. The format for (free) energy trajectory 
 * files as written by promd and md are known to the program. In addition, the
 * user can define custom made formats of any trajectory file that comes in a
 * block-format through a library file. ene_ana is able to read and interpret 
 * series of two types of such files simultaneously, typically referred to as 
 * the "energy file" and the "free energy file".
 * 
 * Using the same library file one can define properties to be calculated from 
 * the values that are listed in them. For the selected properties, ene_ana 
 * will calculate the time series, averages, root-mean-square fluctuations and 
 * a statistical error estimate. The error estimate is calculated from block 
 * averages of different sizes, as described in Allen and Tildesley: "Computer 
 * Simulation of Liquids", 1987. The time for the time series is taken from the
 * trajectory files, unless a different time interval between blocks is 
 * specified through an input parameter. If a topology is supplied, the ene_ana 
 * uses this to define the total solute mass (MASS) and the total number of 
 * solute molecules (NUMMOL).
 * 
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@en_files</td><td>&lt;energy files&gt; (and/or) </td></tr>
 * <tr><td> \@fr_files</td><td>&lt;free energy files&gt; </td></tr>
 * <tr><td> \@prop</td><td>&lt;"properties" to monitor&gt; </td></tr>
 * <tr><td> \@library</td><td>&lt;library for property names&gt; [print] </td></tr>
 * <tr><td> [\@topo</td><td>&lt;molecular topology file&gt; (for MASS and NUMMOL)] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt; (overwrites TIME in the trajectory files)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ene_ana
    @topo       ex.top
    @en_files   ex.tre
    @prop       densit
    @library    ene_ana.lib

   @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cctype>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include "../src/utils/EnergyTraj.h"
#include "../src/utils/debug.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace args;
using namespace gio;
using namespace gcore;

bool file_exists(const std::string &filename) {
    std::ifstream file(filename);
    return file.good();
}

void print(gmath::Stat<double> &p, const string &s, const vector<double> &time);
void set_standards(utils::EnergyTraj &e);
void read_library(const string &name, utils::EnergyTraj &e);
void check_version(const string &library, const utils::EnergyTraj &etrj, 
                   const Ginstream &gin, bool& version_checked);

int main(int argc, char **argv) {
    Argument_List knowns; 
    knowns << "topo" << "time" << "en_files" << "fr_files" << "prop" << "library";

    string usage = "# " + string(argv[0]);
    usage += "\n\t@en_files    <energy files> (and/or)\n";
    usage += "\t@fr_files    <free energy files>\n";
    usage += "\t@prop        <properties to monitor>\n";
    usage += "\t@library     <library for property names> [print]\n";
    usage += "\t[@topo       <molecular topology file> (for MASS and NUMMOL)]\n";
    usage += "\t[@time       <t and dt> (overwrites TIME in the trajectory files)]\n";

    try {
        Arguments args(argc, argv, knowns, usage);

        // Get simulation time
        bool usertime = false;
        vector<double> time;
        double t0 = 0, dt = 1; 
        {
            Arguments::const_iterator iter = args.lower_bound("time");
            if (iter != args.upper_bound("time")) {
                t0 = atof(iter->second.c_str());
                ++iter;
            }
            if (iter != args.upper_bound("time")) {
                dt = atof(iter->second.c_str());
                usertime = true;
                t0 -= dt; // Adjust initial time
            }
        }

        // Check arguments
        if (args.count("en_files") <= 0 && args.count("fr_files") <= 0)
            throw gromos::Exception("ene_ana", "no data specified:\n" + usage);
        if (args.count("prop") <= 0)
            throw gromos::Exception("ene_ana", "no properties to follow:\n" + usage);
        if (args.count("library") <= 0)
            throw gromos::Exception("ene_ana", "no library file specified:\n" + usage);
        
        // Read library file
        string library = "";
        int print_library = 0;
        {
            Arguments::const_iterator iter = args.lower_bound("library"), 
                                       to = args.upper_bound("library");
            if (iter != to) {
                library = iter->second;
                ++iter;
            }
            if (iter != to && iter->second == "print") print_library = 1;
        }
        
        // Read properties
        vector<string> prop;
        int num_prop = 0;
        {
            Arguments::const_iterator iter = args.lower_bound("prop"), 
                                       to = args.upper_bound("prop");
            while (iter != to) {
                prop.push_back(iter->second);
                ++iter;
            }
            num_prop = prop.size();
        }

        // Define an energy trajectory
        utils::EnergyTraj etrj;

        // Read topology for mass and number of molecules
        double mass = 0;
        if (args.count("topo") > 0) {
            InTopology it(args["topo"]);
            System sys(it.system());
            etrj.addConstant("NUMMOL", sys.numMolecules());
            for (int m = 0; m < sys.numMolecules(); m++) {
                for (int a = 0; a < sys.mol(m).topology().numAtoms(); a++) {
                    mass += sys.mol(m).topology().atom(a).mass();
                }
            }
        }
        etrj.addConstant("MASS", mass);
        
        // Read library for variable names
        read_library(library, etrj);
        
        if (print_library) etrj.write_map();

        // Prepare for statistical information
        vector<gmath::Stat<double>> s(num_prop);

        // Define two input streams
        Ginstream gin_en;
        Ginstream gin_fr;
        bool do_energy_files = (args.count("en_files") > 0);
        bool do_free_energy_files = (args.count("fr_files") > 0);
        
        Arguments::const_iterator it_en = args.lower_bound("en_files"),
                                  to_en = args.upper_bound("en_files"),
                                  it_fr = args.lower_bound("fr_files"),
                                  to_fr = args.upper_bound("fr_files");
        bool version_checked = false;
        // Open energy files
        if (do_energy_files) {
            if (!file_exists(it_en->second)) {
                cerr << "ERROR: Energy file " << it_en->second << " does not exist.\n";
                return 1; // Exit with error
            }
            gin_en.open(it_en->second.c_str()); 
            check_version(library, etrj, gin_en, version_checked);
            ++it_en; 
            version_checked = false;
        }   
        
        // Open free energy files
        if (do_free_energy_files) {
            if (!file_exists(it_fr->second)) {
                cerr << "ERROR: Free energy file " << it_fr->second << " does not exist.\n";
                return 1; // Exit with error
            }
            gin_fr.open(it_fr->second.c_str());
            check_version(library, etrj, gin_fr, version_checked);
            ++it_fr;
            version_checked = false;
        }
        
        
        while (true) {
            // Version check
            if (!version_checked) {
                version_checked = true;
            }

            // Read the numbers into the energy trajectory
            if (do_energy_files) {
                int end_of_file = etrj.read_frame(gin_en, "ENERTRJ");
                if (end_of_file) {
                    if (it_en != to_en) {
                        gin_en.close();
                        if (file_exists(it_en->second)) {
                            gin_en.open(it_en->second.c_str());
                            check_version(library, etrj, gin_en, version_checked);
                            ++it_en;
                            version_checked = false;
                        } else {
                            cerr << "ERROR: Energy file " << it_en->second << " does not exist.\n";
                            break;
                        }
                    } else {
                        if (do_free_energy_files) {
                            int end_of_file = etrj.read_frame(gin_fr, "FRENERTRJ");
                            if (!end_of_file) {
                                std::cerr << "# WARNING: frames left over in free energy trajectory,\n"
                                          << "#   they will be disregarded.\n";
                            }
                        }
                        break;
                    }
                }
            }

            if (do_free_energy_files) {
                int end_of_file = etrj.read_frame(gin_fr, "FRENERTRJ");
                if (end_of_file) {
                    if (it_fr != to_fr) {
                        gin_fr.close();
                        if (file_exists(it_fr->second)) {
                            gin_fr.open(it_fr->second.c_str());
                            check_version(library, etrj, gin_fr, version_checked);
                            ++it_fr;
                            version_checked = false;
                        } else {
                            cerr << "ERROR: Free energy file " << it_fr->second << " does not exist.\n";
                            break;
                        }
                    } else {
                        if (do_energy_files) {
                            std::cerr << "# WARNING: frames left over in energy trajectory,\n"
                                      << "#   they will be disregarded.\n";
                        }
                        break;
                    }
                }
            }

            // Calculate and store the necessary number in the stat-classes
            for (int i = 0; i < num_prop; i++)
                s[i].addval(etrj[prop[i]]);
            if (usertime)
                t0 += dt;
            else      
                t0 = etrj["TIME[2]"];
            time.push_back(t0);
        }

        bool flag_error = false;
        for (int i = 0; i < num_prop; i++) {
            if (std::isnan(s[i].ee())) { 
                flag_error = true;
            }      
        }
        if (flag_error) {
            cerr << "# WARNING: One of the values is a NaN,\n"
                 << "#   the data provided are not enough to \n"
                 << "#   give a sensible error estimate" << endl;
        }

        // Print out the statistical information
        cout << setw(10) << "property"
             << " "
             << setw(15) << "average"
             << " "
             << setw(15) << "rmsd"
             << " "
             << setw(15) << "error est."
             << endl;
        for (int i = 0; i < num_prop; i++)
            print(s[i], prop[i], time);
    } 
    catch (const gromos::Exception &e) {
        cerr << e.what() << endl;
        exit(1);
    }
    return 0;
}

void print(gmath::Stat<double> &p, const string &s, const vector<double> &time) {
    if (p.n() != int(time.size())) 
        throw gromos::Exception("ene_ana", "number of time steps read is not equal"
                                 " to number of data points to print out");
    
    ostringstream os;
    os << s << ".dat";
    ofstream fout(os.str().c_str());
    fout.precision(9); // Set precision of numbers going to ofstream
    fout << "#"
         << setw(14) << "time"
         << " "
         << setw(15) << s
         << endl;
    for (int i = 0; i < p.n(); i++) {
        fout << setw(15) << time[i]
             << " "
             << setw(15) << p.val(i)
             << endl;
    }
    fout.close();

    // Print the averages etc. to cout
    cout.precision(9); // Set precision of number going to cout
    cout << setw(10) << s
         << " "
         << setw(15) << p.ave()
         << " "
         << setw(15) << p.rmsd()
         << " "
         << setw(15) << p.ee()
         << endl;
}

void set_standards(utils::EnergyTraj &e) {  
    e.addConstant("BOLTZ", gmath::physConst.get_boltzmann());
}

void read_library(const string &name, utils::EnergyTraj &e) {
    Ginstream gin;
    
    try {
        gin.open(name);
    }
    catch (gromos::Exception ex) {
        throw gromos::Exception("read_library", "failed to open library file " + name);
    }
    
    while (true) {
        vector<string> buffer;
        gin.getblock(buffer);
        if (gin.stream().eof()) break;
        if (buffer[buffer.size() - 1].find("END") != 0)
            throw gromos::Exception("ene_ana", "Library file " + gin.name() + " is corrupted. No END in " + buffer[0] + " block. Got\n" + buffer[buffer.size() - 1]);
        
        string sdum;
        if (buffer[0] == "ENERTRJ" || buffer[0] == "FRENERTRJ") {
            for (unsigned int i = 1; i < buffer.size() - 1; i++) {
                e.addBlock(buffer[i], buffer[0]);
            }
        }
        
        if (buffer[0] == "ENEVERSION") {
            std::cerr << "Found ENEVERSION block" << std::endl;
            if (gin.has_version() || e.has_version()) {
                throw gromos::Exception("ene_ana", "Library file " + gin.name() + " is corrupted. Two ENEVERSION blocks found.");
            }
            string ver;
            gio::concatenate(buffer.begin() + 1, buffer.end() - 1, ver);
            ver.erase(std::remove_if(ver.begin(), ver.end(), ::isspace), ver.end());
            std::cerr << "Version string before setting: " << ver << std::endl;
            e.set_version(ver);
            std::cerr << "Setting EnergyTraj version to: " << ver << std::endl;
        }
        vector<string> data;
        if (buffer[0] == "VARIABLES") {
            set_standards(e);
            
            string bufferstring;
            gio::concatenate(buffer.begin() + 1, buffer.end(), bufferstring);
            
            istringstream iss(bufferstring);
            while (sdum != "END") {
                iss >> sdum;
                data.push_back(sdum);
            }
            
            for (unsigned int i = 0; i < data.size(); i++) {
                if (data[i] == "=") {
                    unsigned int to = i + 1;
                    for (; to < data.size(); to++) if (data[to] == "=") break;
                    
                    ostringstream os;
                    for (unsigned int j = i + 1; j < to - 1; j++) os << " " << data[j]; 
                    e.addKnown(data[i - 1], os.str());
                }
            }
        }
    }
    if (gin.has_version()) {
        e.set_version(gin.version());
    }
}

void check_version(const string &library, const utils::EnergyTraj& etrj, 
                   const Ginstream &gin, bool& version_checked) {
    if (!version_checked){
      version_checked = true;
    }
    DEBUG(1, "Library version: " << etrj.get_version());
    DEBUG(1, "Energy Trajectory version: " << gin.version());

    if (etrj.has_version()) {
        if (gin.has_version()) {
            if (!etrj.version_match(gin.version())) {
                cerr << "WARNING: Version number mismatch!\n"
                     << "         Library " << library << " version: "
                     << etrj.get_version() << std::endl
                     << "         Energy Trajectory " << gin.name() << " version: "
                     << gin.version() << std::endl;
            }
        } else {
            cerr << "WARNING: Version number missing!\n"
                 << "         Library " << library << " version: "
                 << etrj.get_version() << std::endl
                 << "         Energy Trajectory " << gin.name() << " version: "
                 << "NONE" << std::endl;
        }
    } else {
        if (gin.has_version()) {
            cerr << "WARNING: Version number missing!\n"
                 << "         Library " << library << " version: "
                 << "NONE" << std::endl
                 << "         Energy Trajectory " << gin.name() << " version: "
                 << gin.version() << std::endl;
        } else {
            cerr << "WARNING: Version number missing!\n"
                 << "         Library " << library << " version: "
                 << "NONE" << std::endl
                 << "         Energy Trajectory " << gin.name() << " version: "
                 << "NONE" << std::endl;
        }
    }
}

