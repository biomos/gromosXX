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
 * @file mk_script.cc
 * generate scripts to run MD simulations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor mk_script
 * @section mk_script generate scripts to run MD simulations
 * @author @ref co
 * @date 10-6-07
 *
 * A molecular dynamics simulation is usually performed by executing a
 * small script that combines all the necessary files and redirects the
 * output to the appropriate places. When simulations are
 * performed on a queue, such scripts become indispensable. Additionally, in
 * many simulation projects the user prepares similar input files and
 * scripts over and over again.
 * Program mk_script can either create a series of similar scripts that
 * run a sequential list of simulations (keyword \@script) or it can
 * create scripts for a more complex set of simulations that perform a
 * specific task (start-up, perturbation; keyword \@joblist).
 *
 * GROMOS does not require specific filenames for specific types of files.
 * However, most users find it useful to retain some order in
 * their filenames. mk_script has a standard way of constructing
 * filenames that depends on the script number and the system name. The
 * user can specify another set of rules to create filenames through
 * the mk_script library file (keyword \@template). In this file,
 * machine-dependent modifications to the scripts that are to be written can also be
 * specified, such as the job submission command, the MPI command, the stopcommand
 * (to deleta all subsequent jobs from the queue in case the current job fails) and
 * which queue to use (keyword \@queue). A standard location of the mk_script
 * library file can be specified through the environment variable MK_SCRIPT_TEMPLATE.
 *
 * Program mk_script can write input files for promd or md++ (keyword \@version).
 * The md input file (keyword \@files->input) should also be of the correct
 * format: mk_script cannot convert program-specific md input blocks into the
 * analogous blocks for the other version of GROMOS. promd blocks in the input
 * file will be ignored when preparing a md run and vice versa (printing a warning).
 *
 * In addition to write out scripts and input files, mk_script performs a small number
 * of tests on the given input files to prevent the user from submitting
 * a simulation that will fail within the first few seconds. In the messages produced by
 * these tests, a distinction is made between warnings and errors. A
 * warning is given for inconsistencies in the inputs that may lead to an erroneous simulation,
 * but could also be intentional. An error is produced by inconsistencies that will definitely
 * result in the program crashing. Note that passing the tests carried out by mk_script
 * does not guarantee that a simulation will work, as these checks are not exhaustive.
 * All performed tests are listed below (Warnings, Errors).
 *
 * The mentioned tests are done for every script since some variables may change
 * due to a joblist file. If there are no errors, the input file and script will
 * be written to disc. If there are errors, the script and input file will not be
 * written, unless the user forces this (keyword \@force).
 *
 *
 * <b>Warnings:</b>
 * <ol>
 * <li>  the specified coordinate file is not found</li>
 * <li>  the GROMOS binary specified in the mk_script input file cannot be found
 * <li>  the highest LAST atom number in the MULTIBATH block in the md input file is not
 *       equal to the total number of atoms calculated from the topology file and SYSTEM block
 * <li>  DT in the STEP block is too large in combinations with the geometric constraints
 *       in the CONSTRAINT/GEOMCONSTRAINTS block (md/promd). Suggested step sizes are:
 *       0.0005 ps: no constraints on solvent or bonds involving hydrognes
 *       0.001 ps : no constraints on bons not involving hydrogens
 *       0.002 ps : all bonds constraint
 * <li>  the FORCE for bonds that are SHAKEn is calculated
 * <li>  the FORCE for bonds that are not SHAKEn is not calculated
 * <li>  smallest box dimension (length) of the periodic box is less than twice the long-range
 *       cut-off RCUTL in the PAIRLIST/NEIGHBOURLIST block of the md/promd input file
 * <li>  the reaction field cut-off distance RCRF in the NONBONDED block of the md/promd
 *       input file is not equal to the long-range cut-off RCUTL/RCUTI in the
 *       PAIRLIST/NEIGHBOURLIST block (md/promd)
 * <li>  a perturbation topology was specified in the mk_script input file but no
 *       perturbation was requested in the md input file
 * <li>  the combination of RLAM and DLAMT in the PERTURBATION block and the number of steps
 *       from the STEP block in the md/promd input file will lead to a lambda value larger than 1
 * </ol>
 *
 *
 * <b>Errors:</b>
 * <ol>
 * <li>  any of the specified input files (except the coordinate file) is not found</li>
 * <li>  one of the essential blocks is missing (md/promd):
 *       STEP, BOUNDCOND, , INITIALISE, FORCE, CONSTRAINT/GEOMCONSTRAINTS, PAIRLIST/NEIGHBOURLIST, NONBONDED
 * <li>  there is no VELOCITY block in the coordinate file, but NTIVEL in the
 *       INITIALISE block of the md/promd input file specifies that the velocities should be read from file
 * <li>  non-zero NTISHI in the INITIALISE block of the md input file specifies that the lattice shifts
 *       should be initialised, but zero NTB in the BOUNDCOND block specifies a vacuum simulation
 * <li>  there is no LATTICESHIFT block in the coordinate file, but NTISHI in the
 *       INITIALISE block of the md input file specifies that the lattce shifts should be read from file
 * <li>  there is no GENBOX block in the coordinate file, but non-zero NTB in the BOUNDCOND block
 *       specifies a non-vacuum simulation
 * <li>  the number of the last atom given in the FORCE block of the md/promd input file is
 *       not equal to the total number of atoms calculated from the topology and SYSTEM block
 * <li>  the number of atoms calculated from the topology and SYSTEM block of the md/promd input
 *       file is not equal to the number of atoms in the POSITION block of the coordinate file
 * <li>  in the PAIRLIST/NEIGHBOURLIST block, the short-range cutoff RCUTP/RCUTS is larger than the
 *       long-range cutoff RCUTL/RCUTI (md/promd)
 * <li>  no position restraints specification file is specified in the mk_script input file,
 *       but position restraining is switched on in the md input file
 * <li>  no perturbation topology file is specified in the mk_script input file,
 *       but perturbation is switched on in the md input file
 * </ol>
 *
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@sys</td><td>&lt;system name&gt; </td></tr>
 * <tr><td> \@bin</td><td>&lt;GROMOS binary to use&gt; </td></tr>
 * <tr><td> \@version</td><td>&lt;md++ / promd&gt; </td></tr>
 * <tr><td> \@dir</td><td>&lt;directory where mk_script is executed&gt; </td></tr>
 * <tr><td> \@files</td><td>(give all files relative to \@dir)</td></tr>
 * <tr><td></td><td>topo</td><td>&lt;molecular topology file&gt;</td></tr>
 * <tr><td></td><td>input</td><td>&lt;input file&gt;</td></tr>
 * <tr><td></td><td>coord</td><td>&lt;initial coordinates&gt;</td></tr>
 * <tr><td></td><td>[refpos</td><td>&lt;reference positions&gt;]</td></tr>
 * <tr><td></td><td>[posresspec</td><td>&lt;position restraints specifications&gt;]</td></tr>
 * <tr><td></td><td>[disres</td><td>&lt;distance restraints&gt;]</td></tr>
 * <tr><td></td><td>[dihres</td><td>&lt;dihedral restraints&gt;]</td></tr>
 * <tr><td></td><td>[angres</td><td>&lt;bond-angle restraints&gt;]</td></tr>
 * <tr><td></td><td>[colvar</td><td>&lt;collective variable restraints&gt;]</td></tr>
 * <tr><td></td><td>[jvalue</td><td>&lt;j-value restraints&gt;]</td></tr>
 * <tr><td></td><td>[order</td><td>&lt;order parameter restraints&gt;]</td></tr>
 * <tr><td></td><td>[sym</td><td>&lt;symmetry restraints&gt;]</td></tr>
 * <tr><td></td><td>[ledih</td><td>&lt;local elevation dihedrals&gt;]</td></tr>
 * <tr><td></td><td>[friction</td><td>&lt;friction coefficients&gt;]</td></tr>
 * <tr><td></td><td>[leumb</td><td>&lt;local elevation umbrellas&gt;]</td></tr>
 * <tr><td></td><td>[bsleus</td><td>&lt;B&S-LEUS topology file&gt;]</td></tr>
 * <tr><td></td><td>[jin</td><td>&lt;J formatted input file&gt;]</td></tr>
 * <tr><td></td><td>[jtrj</td><td>&lt;J formatted trajectory file&gt;]</td></tr>
 * <tr><td></td><td>[pttopo</td><td>&lt;perturbation topology&gt;]</td></tr>
 * <tr><td></td><td>[gamd</td><td>&lt;gamd input file&gt;]</td></tr>
 * <tr><td></td><td>[xray</td><td>&lt;xray restraints file&gt;]</td></tr>
 * <tr><td></td><td>[repout</td><td>&lt;replica exchange output file&gt;]</td></tr>
 * <tr><td></td><td>[repdat</td><td>&lt;replica exchange data file&gt;]</td></tr>
 * <tr><td> [\@script</td><td>&lt;first script&gt; &lt;number of scripts&gt;] </td></tr>
 * <tr><td> [\@joblist</td><td>&lt;joblist file&gt;] </td></tr>
 * <tr><td> [\@template</td><td>&lt;mk_script library file&gt;] </td></tr>
 * <tr><td> [\@queue</td><td>&lt;which queue to use in mk_script library file&gt;] </td></tr>
 * <tr><td> [\@cmd</td><td>&lt;overwrite last command in mk_script library file&gt;] </td></tr>
 * <tr><td> [\@force</td><td>(write script regardless of errors)] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  mk_script
    @sys         ex
    @bin         /usr/local/gromos/md++/bin/md
    @version     md++
    @dir         /home/user
    @joblist     joblist.startup
    @files
       topo      ex.top
       coord     exref.coo
       input     imd.dat
    @template    mk_script.lib
    @queue       igcpc
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cmath>
#include <unistd.h>

#include <sys/syscall.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mk_script.h"

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InG96.h"
#include "../src/gmath/Vec.h"
#include "../src/gromos/Exception.h"

#include "../config.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace gmath;


void printInput(string ofile, input gin);
void readLibrary(string file, vector<directive> &directives,
        vector<filename> &names, vector<filename> &misc,
        vector<string> &linknames, vector<int> &linkadditions,
        string system, string queue, double t,
        double dt, int ns);
void readJobinfo(string file, map<int, jobinfo> &ji);
void setParam(input &gin, jobinfo const &job);
bool file_exists(string file);

string word_expansion(const string & arg);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "sys" << "script" << "bin" << "dir" << "queue"
          << "files" << "template" << "version" << "cmd" << "joblist"
          << "force" << "mail" << "putdev" << "debug";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@sys           <system name>\n";
  usage += "\t@bin           <gromos binary to use>\n";
  usage += "\t@dir           <directory where mk_script is executed>\n";
  usage += "\t@version       <md++ / promd>\n";
  usage += "\t[@script       <first script> <number of scripts>]\n";
  usage += "\t[@joblist      <joblist file>]\n";
  usage += "\t@files         <specify all files relative to @dir>\n";
  usage += "\t\ttopo         <molecular topology file>\n";
  usage += "\t\tinput        <input file>\n";
  usage += "\t\tcoord        <initial coordinates>\n";
  usage += "\t\t[refpos      <reference positions>]\n";
  usage += "\t\t[posresspec  <position restraints specifications>]\n";
  usage += "\t\t[disres      <distance restraints>]\n";
  usage += "\t\t[dihres      <dihedral restraints>]\n";
  usage += "\t\t[angres      <angle restraints>]\n";
  usage += "\t\t[colvarres      <collective variable restraints>]\n";
  usage += "\t\t[jvalue      <j-value restraints>]\n";
  usage += "\t\t[order       <order parameter restraints>]\n";
  usage += "\t\t[sym         <symmetry restraints>]\n";
  usage += "\t\t[ledih       <local elevation dihedrals>]\n";
  usage += "\t\t[friction    <friction coefficients>]\n";
  usage += "\t\t[leumb       <local elevation umbrellas>]\n";
  usage += "\t\t[bsleus      <B&S-LEUS topology file>]\n";
  usage += "\t\t[qmmm        <QMMM specification file>]\n";
  usage += "\t\t[jin         <J formatted input file>]\n";
  usage += "\t\t[jtrj        <J formatted trajectory file>]\n";
  usage += "\t\t[pttopo      <perturbation topology>]\n";
  usage += "\t\t[gamd        <gamd input file>]\n";
  usage += "\t\t[xray        <xray restraints file>]\n";
  usage += "\t\t[repout      <replica exchange output file>]\n";
  usage += "\t\t[repdat      <replica exchange data file>]\n";
  usage += "\t[@template     <template filename, absolute or relative to @dir>]\n";
  usage += "\t[@queue        <queue flags>]\n";
  usage += "\t[@cmd          <overwrite last command>]\n";
  usage += "\t[@mail         <get ERROR email for crashed jobs: EVERY (standard), LAST or NONE> [<job ID>] ]\n";
  usage += "\t[@force        (write script regardless of errors)]\n";
  usage += "\t[@putdev       (puts the @develop flag to the simulation script files)\n";

  try {

    Arguments args(argc, argv, knowns, usage);
    
    // check if the develop flag was set
    bool putdevelop = false;
    if((args.count("putdev") >= 0 )) {
      putdevelop = true;
    }
    // check if the debug flag was set
    bool putdebug = false;
    if((args.count("debug") >= 0 )) {
      putdebug = true;
    }
    
    // error emails?
    string mail = "EVERY";
    string jobID = "";
    if (args.count("mail") >= 0) {
      if (args.count("mail") > 0 && args.count("mail") < 3) {
        Arguments::const_iterator it = args.lower_bound("mail");
        mail = it->second;
        if(args.count("mail") == 2) {
          it++;
          jobID = ": " + it->second;
        }
      }
      if((mail != "EVERY" && mail != "LAST" && mail != "NONE") || (args.count("mail") != 1 && args.count("mail") != 2)) {
        throw gromos::Exception(argv[0], "check arguments for @mail");
      }
    }
    
    // first get some input parameters
    int scriptNumber = 1, numScripts = 1;
    string simuldir;
    {
      Arguments::const_iterator iter = args.lower_bound("dir");
      if (iter != args.upper_bound("dir")) {
        simuldir = word_expansion(iter->second);
        if (simuldir[0] != '/')
          throw gromos::Exception("mk_script",
                "Specified directory (@dir) should be an absolute path: " + simuldir);
        if (chdir(simuldir.c_str()) != 0)
          throw gromos::Exception("mk_script",
                "Specified directory (@dir) does not exist");
      } else
        simuldir = "`pwd`";

      iter = args.lower_bound("script");
      if (iter != args.upper_bound("script")) {
        scriptNumber = atoi(iter->second.c_str());
        ++iter;
      }
      if (iter != args.upper_bound("script"))
        numScripts = atoi(iter->second.c_str());

      if (numScripts < 0)
        throw Arguments::Exception("Can't deal with negative number of scripts in @script argument");
    }
    string systemname = args["sys"];

    // read in the library
    string libraryfile;
    if (args.lower_bound("template") == args.upper_bound("template")) {
      // try to get it from environment
      if (getenv("MK_SCRIPT_TEMPLATE") != NULL) {
        libraryfile = getenv("MK_SCRIPT_TEMPLATE");
        cout << "MK_SCRIPT_TEMPLATE = " << libraryfile << endl;
      } else {
        throw gromos::Exception("mk_script", "Please give @template or set the "
                "MK_SCRIPT_TEMPLATE environment variable.");
      }
    } else { // supplied by @template
      libraryfile = word_expansion(args["template"]);
    }

    // parse the files
    int l_coord = 0, l_topo = 0, l_input = 0, l_refpos = 0, l_posresspec = 0, l_xray = 0, l_anatrx;
    int l_disres = 0, l_dihres = 0, l_angres = 0, l_jvalue = 0, l_order = 0, l_sym = 0, l_ledih = 0;
    int l_friction=0, l_leumb = 0, l_bsleus = 0, l_qmmm = 0, l_pttopo = 0, l_gamd = 0;
    int l_repout=0, l_repdat=0;
    int l_jin = 0;
    int l_colvarres = 0;
    string s_coord, s_topo, s_input, s_refpos, s_posresspec, s_xray, s_anatrx;
    string s_disres, s_dihres, s_angres, s_jvalue, s_order, s_sym, s_ledih, s_leumb, s_bsleus, s_qmmm;
    string s_colvarres;
    string s_friction, s_pttopo, s_jin, s_gamd;
    string s_repout, s_repdat;
    for (Arguments::const_iterator iter = args.lower_bound("files"),
            to = args.upper_bound("files"); iter != to; ++iter) {
      switch (FILETYPE[iter->second]) {
        case anatrxfile: ++iter; //NO check if file exists
          s_anatrx = iter->second;
          l_anatrx = 1;
          break;
        case coordfile: ++iter; //check if file exists comes later
          s_coord = iter->second;
          l_coord = 1;
          break;
        case inputfile: ++iter;
          s_input = iter->second;
          if(file_exists(s_input))
            l_input = 1;
          else
            printError("File " + s_input + " does not exist!");
          break;
        case topofile: ++iter;
          s_topo = iter->second;
          if(file_exists(s_topo))
            l_topo = 1;
          else
            printError("File " + s_topo + " does not exist!");
          break;
        case refposfile: ++iter;
          s_refpos = iter->second;
          if(file_exists(s_refpos))
            l_refpos = 1;
          else
            printError("File " + s_refpos + " does not exist!");
          break;
        case posresspecfile:++iter;
          s_posresspec = iter->second;
          if(file_exists(s_posresspec))
            l_posresspec = 1;
          else
            printError("File " + s_posresspec + " does not exist!");
          break;
        case xrayfile:++iter;
          s_xray = iter->second;
          if(file_exists(s_xray))
            l_xray = 1;
          else
            printError("File " + s_xray + " does not exist!");
          break;
        case disresfile: ++iter;
          s_disres = iter->second;
          if(file_exists(s_disres))
            l_disres = 1;
          else
            printError("File " + s_disres + " does not exist!");
          break;
        case colvarresfile: ++iter;
          s_colvarres = iter->second;
          if(file_exists(s_colvarres))
            l_colvarres = 1;
          else
            printError("File " + s_colvarres + " does not exist!");
          break;
        case dihresfile: ++iter;
          s_dihres = iter->second;
          if(file_exists(s_dihres))
            l_dihres = 1;
          else
            printError("File " + s_dihres + " does not exist!");
          break;
        case angresfile: ++iter;
          s_angres = iter->second;
          if(file_exists(s_angres))
            l_angres = 1;
          else
            printError("File " + s_angres + " does not exist!");
          break;
        case jvaluefile: ++iter;
          s_jvalue = iter->second;
          if(file_exists(s_jvalue))
            l_jvalue = 1;
          else
            printError("File " + s_jvalue + " does not exist!");
          break;
        case orderfile: ++iter;
          s_order = iter->second;
          if(file_exists(s_order))
            l_order = 1;
          else
            printError("File " + s_order + " does not exist!");
          break;
        case symfile: ++iter;
          s_sym = iter->second;
          if(file_exists(s_sym))
            l_sym = 1;
          else
            printError("File " + s_sym + " does not exist!");
          break;
        case frictionfile: ++iter;
          s_friction = iter->second;
          if(file_exists(s_friction))
            l_friction = 1;
          else
            printError("File " + s_friction + " does not exist!");
          break;
        case ledihfile: ++iter;
          s_ledih = iter->second;
          if(file_exists(s_ledih))
            l_ledih = 1;
          else
            printError("File " + s_ledih + " does not exist!");
          break;
        case leumbfile: ++iter;
          s_leumb = iter->second;
          if(file_exists(s_leumb))
            l_leumb = 1;
          else
            printError("File " + s_leumb + " does not exist!");
          break;
        case bsleusfile: ++iter;
          s_bsleus = iter->second;
          if(file_exists(s_bsleus))
            l_bsleus = 1;
          else
            printError("File " + s_bsleus + " does not exist!");
          break;
        case qmmmfile: ++iter;
          s_qmmm = iter->second;
          if(file_exists(s_qmmm))
            l_qmmm = 1;
          else
            printError("File " + s_qmmm + " does not exist!");
          break;
        case jinfile: ++iter;
          s_jin = iter->second;
          if(file_exists(s_jin))
            l_jin = 1;
          else
            printError("File " + s_jin + " does not exist!");
          break;
        case pttopofile: ++iter;
          s_pttopo = iter->second;
          if(file_exists(s_pttopo))
            l_pttopo = 1;
          else
            printError("File " + s_pttopo + " does not exist!");
          break;
        case gamdfile: ++iter;
          s_gamd = iter->second;
          if(file_exists(s_gamd))
            l_gamd = 1;
          else
            printError("File " + s_gamd + " does not exist!");
          break;
        case repoutfile: ++iter; //NO check if file exists: output files
          s_repout = iter->second;
          l_repout=1;
          break;
        case repdatfile: ++iter; //NO check if file exists: output files
          s_repdat = iter->second;
          l_repdat=1;
          break;
        case outputfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrxfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrvfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrffile: ++iter;
          printWarning(iter->second + " not used");
        case outtrefile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrgfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outbaefile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outbagfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case outtrsfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case joutfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case jtrjfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case scriptfile: ++iter;
          printWarning(iter->second + " not used");
          break;
        case unknownfile: printError("Don't know how to handle file "
                  + iter->second);
      }
    }

    // check which outformat we want (gromos96 or gromosXX)
    // in g08 only relevant for job scripts!
    bool gromosXX = false;
    if (args["version"] == "md++")
      gromosXX = true;
    else if (args["version"] == "promd")
      gromosXX = false;
    else
      throw gromos::Exception("mk_script", "Please specify version (md++ or promd).");


    // read topology
    if (!l_topo) {
      throw gromos::Exception("mk_script", "You have to specify a topology\n" + usage);
    }
    InTopology it(s_topo);
    System sys = it.system();
    fileInfo topo;
    {
      Ginstream itopo(s_topo);
      itopo >> topo;
      itopo.close();
    }


    // read the input file
    int numErrorsBefore = numErrors;
    if (!l_input) {
      throw gromos::Exception("mk_script", "You have to specify an input file\n" + usage);
    }
    Ginstream imd(s_input);

    input gin;
    cout << "\nReading the input file (" << s_input << "):\n";
    cout << "----------------------------------------\n\n";
    imd >> gin;
    if (numErrors - numErrorsBefore != 0) { // no forcing possible, if the input file was wrongly read
                          // you really should not go on writing the scrips
      cout << "\n--------------------------------------------------------------------------------\n\n";
      if (numErrors  - numErrorsBefore > 1) {
        cout << "THERE WERE " << numErrors << " ERRORS WHEN READING THE INPUT FILE\n";
      } else {
        cout << "THERE WAS 1 ERROR WHEN READING THE INPUT FILE\n";
      }
      cout << "No scripts have been written\n";
      exit(1);
    } else {
      cout << "DONE\n";
      cout << "\n--------------------------------------------------------------------------------\n\n";
    }

    imd.close();
    
    // adapt time if replica exchange -- this seems to assume we do an integer number of ps?
    int steps = gin.step.nstlim * gin.step.dt;
    if (gin.replica.found)
      steps *= (gin.replica.nretrial + gin.replica.nrequil);
    else if (gin.reeds.found)
      steps *= (gin.reeds.nretrial + gin.reeds.nrequil);

    // read the jobinfo file
    if (args.count("joblist") > 0 && args.count("script") > 0)
      throw gromos::Exception("mk_script", "You can only specify @script OR "
            "@joblist");

    map<int, jobinfo> joblist;
    if (args.count("joblist") > 0) {
      scriptNumber = 0;
      readJobinfo(args["joblist"], joblist);
      map<int, jobinfo>::iterator iter = joblist.begin(), to = joblist.end();
      {
        ostringstream os;
        os << gin.step.t;
        iter->second.param["T"] = os.str();
      }
      {
        if(iter->second.param.find("NSTLIM")!=iter->second.param.end()){
          int nstlimnew = atoi(iter->second.param["NSTLIM"].c_str());
          // because steps is assumed to be integer, this would work in most cases
          cerr << "steps" << steps << endl;
          //steps *= nstlimnew;
          //steps = steps / gin.step.nstlim;
          steps -= gin.step.nstlim * gin.step.dt;
          steps += nstlimnew * gin.step.dt; 
          cerr << "nstlimnew " << nstlimnew << endl;
          cerr << "steps" << steps << endl;
          cerr << "gin.step.nstlim" << gin.step.nstlim << endl;
        }
        ostringstream os;
        os << gin.step.t + steps;
        iter->second.param["ENDTIME"] = os.str();
      }
      for (++iter; iter != to; ++iter) {
        if (joblist.find(iter->second.prev_id) != joblist.end() &&
                iter->first != iter->second.prev_id) {
          iter->second.param["T"] = "-1";
          iter->second.param["ENDTIME"] = "-1";
        }
      }
    } else { // no joblist
      for (int i = 0; i < numScripts; i++) {
        jobinfo job;
        {
          ostringstream os;
          os << gin.step.t + i * steps;
          job.param["T"] = os.str();
        }
        {
          ostringstream os;
          os << gin.step.t + (i + 1) * steps;
          job.param["ENDTIME"] = os.str();
        }
        job.dir = ".";
        job.prev_id = i + scriptNumber - 1;
        joblist[i + scriptNumber] = job;
      }
    }

    // replace the %queue% variable?
    string queue = "";
    {
      ostringstream os;
      Arguments::const_iterator iter = args.lower_bound("queue"),
              to = args.upper_bound("queue");
      for (; iter != to; ++iter) {
        std::string s = iter->second;
        if (s.find("\\n") != string::npos)
          s.replace(s.find("\\n"), 2, "\n");
        else s += " ";

        // os << iter->second << " ";
        os << s;
      }
      queue = os.str();
    }
    
    // put debug @verb flag
    string debug = "";
    {
      ostringstream os;
      Arguments::const_iterator iter = args.lower_bound("debug"),
              to = args.upper_bound("debug");
      for (; iter != to; ++iter) {
        std::string s = iter->second;
        if (s.find("\\n") != string::npos)
          s.replace(s.find("\\n"), 2, " ");
        else s += " ";

        // os << iter->second << " ";
        os << s;
      }
      debug = os.str();
    }

    vector<directive> directives;
    // create names for automated file names
    vector<filename> filenames;
    vector<filename> misc;
    vector<int> linkadditions;
    vector<string> linknames;

    for (int i = 0; i < numFiletypes; i++) {
      filename newname(systemname, gin.step.t, steps,
              scriptNumber, queue);
      filenames.push_back(newname);
    }
    // workdir lastcommand firstcommand mpicommand stopcommand
    for (int i = 0; i < 5; i++) {
      filename newname(systemname, gin.step.t, steps,
              scriptNumber, queue);
      misc.push_back(newname);
    }

    // set the standard templates
    filenames[FILETYPE["script"]].setTemplate("%system%_%number%.run");
    filenames[FILETYPE["input"]].setTemplate("%system%_%number%.imd");
    filenames[FILETYPE["topo"]].setTemplate("%system%.top");
    filenames[FILETYPE["refpos"]].setTemplate("%system%_%number%.rpr");
    filenames[FILETYPE["posresspec"]].setTemplate("%system%_%number%.por");
    filenames[FILETYPE["anatrj"]].setTemplate("%system%_%number%.trc.gz");
    filenames[FILETYPE["xray"]].setTemplate("%system%_%number%.xrs");
    filenames[FILETYPE["disres"]].setTemplate("%system%_%number%.dsr");
    filenames[FILETYPE["colvarres"]].setTemplate("%system%_%number%.cvr");
    filenames[FILETYPE["pttopo"]].setTemplate("%system%.ptp");
    filenames[FILETYPE["gamd"]].setTemplate("%system%.gamd");
    filenames[FILETYPE["dihres"]].setTemplate("%system%_%number%.dhr");
    filenames[FILETYPE["angres"]].setTemplate("%system%_%number%.bar");
    filenames[FILETYPE["jvalue"]].setTemplate("%system%_%number%.jvr");
    filenames[FILETYPE["order"]].setTemplate("%system%_%number%.ord");
    filenames[FILETYPE["sym"]].setTemplate("%system%_%number%.sym");
    filenames[FILETYPE["ledih"]].setTemplate("%system%_%number%.led");
    filenames[FILETYPE["leumb"]].setTemplate("%system%_%number%.lud");
    filenames[FILETYPE["bsleus"]].setTemplate("%system%.bsleus");
    filenames[FILETYPE["jin"]].setTemplate("%system%_%number%.jin");
    filenames[FILETYPE["jout"]].setTemplate("%system%_%number%.jin");
    filenames[FILETYPE["jtrj"]].setTemplate("%system%_%number%.jtrj");
    filenames[FILETYPE["friction"]].setTemplate("%system%_%number%.frc");
    filenames[FILETYPE["coord"]].setTemplate("%system%_%number%.cnf");
    filenames[FILETYPE["qmmm"]].setTemplate("%system%_%number%.qmm");
    filenames[FILETYPE["output"]].setTemplate("%system%_%number%.omd");
    filenames[FILETYPE["outtrx"]].setTemplate("%system%_%number%.trc");
    filenames[FILETYPE["outtrv"]].setTemplate("%system%_%number%.trv");
    filenames[FILETYPE["outtrf"]].setTemplate("%system%_%number%.trf");
    filenames[FILETYPE["outtre"]].setTemplate("%system%_%number%.tre");
    filenames[FILETYPE["outtrg"]].setTemplate("%system%_%number%.trg");
    filenames[FILETYPE["outbae"]].setTemplate("%system%_%number%.bae");
    filenames[FILETYPE["outbag"]].setTemplate("%system%_%number%.bag");
    filenames[FILETYPE["outtrs"]].setTemplate("%system%_%number%.trs");
    if(!l_repout) filenames[FILETYPE["repout"]].setTemplate("repout_%system%_%number%.dat");
    else filenames[FILETYPE["repout"]].setTemplate(s_repout);
    if(!l_repdat) filenames[FILETYPE["repdat"]].setTemplate("repdat_%system%_%number%.dat");
    else filenames[FILETYPE["repdat"]].setTemplate(s_repdat);

    cout << "Reading the library...";
    // And here is a gromos-like function call!
    readLibrary(libraryfile, directives, filenames, misc,
            linknames, linkadditions,
            systemname, queue, gin.step.t,
            steps,
            scriptNumber);
    cout << "done!\n";

    // overwrite last command if given as argument
    if (args.count("cmd") > 0) {
      ostringstream os;
      Arguments::const_iterator iter = args.lower_bound("cmd"),
              to = args.upper_bound("cmd");
      for (; iter != to; ++iter) {
        std::string s = iter->second;
        if (s.find("\\n") != string::npos)
          s.replace(s.find("\\n"), 2, "\n");
        else s += " ";

        // os << iter->second << " ";
        os << s;
      }

      misc[1].setTemplate(os.str());
    }

    // read what is in the coordinate file
    if (!l_coord) { //no coordinate file was specified
      // try to open it from the template
      ifstream fin(filenames[FILETYPE["coord"]].name(-1).c_str());
      if (fin) {
        ostringstream os;
        os << "No coordinate file specified, but I found "
                << filenames[FILETYPE["coord"]].name(-1)
                << " which I will use\n";
        printWarning(os.str());
        s_coord = filenames[FILETYPE["coord"]].name(-1);
        l_coord = 1;
      } else {
        ostringstream os;
        os << "No coordinate file is specified, checks involving this file are not performed\n";
        os << "Assuming it does not exist yet, I will use "
                << filenames[FILETYPE["coord"]].name(-1)
                << " in the script\n";
        s_coord = filenames[FILETYPE["coord"]].name(-1);

        printWarning(os.str());
      }

    }
    // read the coordinate file now; don't read if replica exchange continuation run
    if (l_coord && gin.replica.cont != 1) {
      if(!file_exists(s_coord)){
        string msg = "Coordinate file " + s_coord + " not found! No checks involving this file will be performed.\n";
        printWarning(msg); //only print warning, not error
      }
      else{ //the file exists: make some checks
        InG96 ic;
        ic.open(s_coord);
        ic.select("ALL");
        int count = 0;
        while (!ic.eof()) {
          ic >> sys;
          count++;
        }
        ic.close();
        if (count > 1) {
          stringstream msg;
          msg << count << "frames were read from the coordinate file " << s_coord;
          msg << "The last frame was used as initial configuration.";
          printWarning(msg.str());
        }
      }
    }

    // calculate some standard numbers
    int numSoluteAtoms = 0;
    for (int i = 0; i < sys.numMolecules(); i++)
      numSoluteAtoms += sys.mol(i).topology().numAtoms();
    int numSolventAtoms = sys.sol(0).topology().numAtoms();
    int numVirtualAtoms = 0;
    if(gin.virtualatom.found && gin.virtualatom.virt==1) numVirtualAtoms+=gin.virtualatom.numvirt;
    int numTotalAtoms = gin.system.npm * numSoluteAtoms +
            gin.system.nsm * numSolventAtoms;

    // carry out a thousand tests:

    // Does the binary exist?
    std::string gromosbin = word_expansion(args["bin"]);
    {
      ifstream fin(gromosbin.c_str());
      if (!fin)
        printWarning("Specified binary not found! "
              + gromosbin + "\n");
      else
        fin.close();
    }

    // ***** LOOP OVER JOBS *******
    unsigned int writtenScripts = 0;
    map<int, jobinfo>::iterator iter = joblist.begin(), to = joblist.end();
    for (; iter != to; ++iter) {

      //make sure we start in the right directory
      if(chdir(simuldir.c_str()) != 0) {
        throw gromos::Exception(argv[0], "could not change to the simuldir directory");
      }

      l_coord = l_coord && iter == joblist.begin();

      // update the input parameters
      setParam(gin, iter->second);
      {
        ostringstream os;
        os << steps;
        iter->second.param["DELTAT"] = os.str();
      }
      if (iter != joblist.begin()) {
        double time = atof(joblist[iter->second.prev_id].param["ENDTIME"].c_str());
        // some of the parameters used to compute steps may have changed, so we determine it newly
        steps = gin.step.nstlim * gin.step.dt;
        if (gin.replica.found)
          steps *= (gin.replica.nretrial + gin.replica.nrequil);

        double endtime = time + steps;
        ostringstream os;
        os << endtime;
        iter->second.param["T"] = joblist[iter->second.prev_id].param["ENDTIME"];
        iter->second.param["ENDTIME"] = os.str();
        gin.step.t = time;

      }

      for (unsigned int i = 0; i < directives.size(); i++) {
        directives[i].setInfo(systemname, gin.step.t, steps,
                iter->first, queue);
      }
      for (unsigned int i = 0; i < filenames.size(); i++) {
        filenames[i].setInfo(systemname, gin.step.t, steps,
                iter->first, queue);
      }
      for (unsigned int i = 0; i < misc.size(); i++) {
        misc[i].setInfo(systemname, gin.step.t, steps,
                iter->first, queue);
      }

      // do the checks for every script (some variables may change due to a
      // joblist file
      cout << "Performing checks for script " << iter->first << endl;
      cout << "----------------------------------------" << endl;
      numWarnings = 0;
      numErrors = 0;

      // Ignore md++ specific blocks if promd input is to be written
      // and vice versa (start)
      if (!gromosXX) { // Ignore md++ specific blocks
        if (gin.multibath.found) {
          if (gin.thermostat.found) {
            printWarning("Ignored md++ specific block MULTIBATH\n");
            gin.multibath.found = 0;
          } else {
            printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                    "block MULTIBATH. Maybe you want to specify the THERMOSTAT block instead?\n");
          }
        }
        if (gin.pressurescale.found) {
          if (gin.barostat.found) {
            printWarning("Ignored md++ specific block PRESSURESCALE\n");
            gin.pressurescale.found = 0;
          } else {
            printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                    "block PRESSURESCALE. Maybe you want to specify the BAROSTAT block instead?\n");
          }
        }
        if (gin.comtransrot.found) {
          if (gin.overalltransrot.found) {
            printWarning("Ignored md++ specific block COMTRANSROT\n");
            gin.comtransrot.found = 0;
          } else {
            printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                    "block COMTRANSROT. Maybe you want to specify the OVERALLTRANSROT block instead?\n");
          }
        }
        if (gin.ewarn.found) {
          printWarning("Ignored md++ specific block EWARN\n");
          gin.ewarn.found = 0;
        }
        if (gin.constraint.found) {
          if (gin.geomconstraints.found) {
            printWarning("Ignored md++ specific block CONSTRAINT\n");
            gin.constraint.found = 0;
          } else {
            printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                    "block CONSTRAINT. Maybe you want to specify the GEOMCONSTRAINTS block instead?\n");
          }
        }
        if (gin.pairlist.found) {
          if (gin.neighbourlist.found) {
            printWarning("Ignored md++ specific block PAIRLIST\n");
            gin.pairlist.found = 0;
          } else {
            printError("You want to perform a PROMD run and you have specified the md++ specific\n"
                    "block PAIRLIST. Maybe you want to specify the NEIGHBOURLIST block instead?\n");
          }
        }
        if (gin.cgrain.found) {
          printWarning("Ignored md++ specific block CGRAIN\n");
          gin.cgrain.found = 0;
        }
        if (gin.rottrans.found) {
          printWarning("Ignored md++ specific block ROTTRANS\n");
          gin.rottrans.found = 0;
        }
        if (gin.lambdas.found) {
          printWarning("Ignored md++ specific block LAMBDAS\n");
          gin.lambdas.found = 0;
        }
        if (gin.perscale.found) {
          printWarning("Ignored md++ specific block PERSCALE\n");
          gin.perscale.found = 0;
        }
        if (gin.polarise.found) {
          printWarning("Ignored md++ specific block POLARISE\n");
          gin.polarise.found = 0;
        }
        if (gin.electric.found) {
          printWarning("Ignored md++ specific block ELECTRIC\n");
          gin.electric.found = 0;
        }
        if (gin.replica.found) {
          printWarning("Ignored md++ specific block REPLICA\n");
          gin.replica.found = 0;
        }
        if (gin.innerloop.found) {
          printWarning("Ignored md++ specific block INNERLOOP\n");
          gin.innerloop.found = 0;
        }
        if (gin.integrate.found) {
          printWarning("Ignored md++ specific block INTEGRATE\n");
          gin.integrate.found = 0;
        }
        if (gin.randomnumbers.found) {
          printWarning("Ignored md++ specific block RANDOMNUMBERS\n");
          gin.randomnumbers.found = 0;
        }
        if (gin.multistep.found) {
          printWarning("Ignored md++ specific block MULTISTEP\n");
          gin.multistep.found = 0;
        }
        if(gin.eds.found) {
          printWarning("Ignored md++ specific block EDS\n");
          gin.eds.found = 0;
        }
        if(gin.aeds.found) {
          printWarning("Ignored md++ specific block AEDS\n");
          gin.eds.found = 0;
        }
        if(gin.gamd.found) {
          printWarning("Ignored md++ specific block GAMD\n");
          gin.gamd.found = 0;
        }
        if(gin.bsleus.found) {
          printWarning("Ignored md++ specific block BSLEUS\n");
          gin.bsleus.found = 0;
        }
        if(gin.sasa.found) {
          printWarning("Ignored md++ specific block SASA\n");
          gin.sasa.found = 0;
        }
      } else { // Ignore promd specific blocks
        if (gin.consistencycheck.found) {
          printWarning("Ignored promd specific block CONSISTENCYCHECK\n");
          gin.consistencycheck.found = 0;
        }
        if (gin.thermostat.found) {
          if (gin.multibath.found) {
            printWarning("Ignored promd specific block THERMOSTAT\n");
            gin.thermostat.found = 0;
          } else {
            printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                    "block THERMOSTAT. Maybe you want to specify the MULTIBATH block instead?\n");
          }
        }
        if (gin.barostat.found) {
          if (gin.pressurescale.found) {
            printWarning("Ignored promd specific block BAROSTAT\n");
            gin.barostat.found = 0;
          } else {
            printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                    "block BAROSTAT. Maybe you want to specify the PRESSURESCALE block instead?\n");
          }
        }
        if (gin.virial.found) {
          printWarning("Ignored promd specific block VIRIAL\n");
          gin.virial.found = 0;
        }
        if (gin.overalltransrot.found) {
          if (gin.comtransrot.found) {
            printWarning("Ignored promd specific block OVERALLTRANSROT\n");
            gin.overalltransrot.found = 0;
          } else {
            printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                    "block OVERALLTRANSROT. Maybe you want to specify the COMTRANSROT block instead?\n");
          }
        }
        if (gin.debug.found) {
          printWarning("Ignored promd specific block DEBUG\n");
          gin.debug.found = 0;
        }
        if (gin.geomconstraints.found) {
          if (gin.constraint.found) {
            printWarning("Ignored promd specific block GEOMCONSTRAINTS\n");
            gin.geomconstraints.found = 0;
          } else {
            printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                    "block GEOMCONSTRAINTS. Maybe you want to specify the CONSTRAINT block instead?\n");
          }
        }
        if (gin.neighbourlist.found) {
          if (gin.pairlist.found) {
            printWarning("Ignored promd specific block NEIGHBOURLIST\n");
            gin.neighbourlist.found = 0;
          } else {
            printError("You want to perform a md++ run and you have specified the PROMD specific\n"
                    "block NEIGHBOURLIST. Maybe you want to specify the PAIRLIST block instead?\n");
          }
        }
        if (gin.umbrella.found) {
          printWarning("Ignored promd specific block UMBRELLA\n");
          gin.umbrella.found = 0;
        }
        if (gin.pathint.found) {
          printWarning("Ignored promd specific block PATHINT\n");
          gin.pathint.found = 0;
        }
      }

      // And check if all compulsory blocks have been specified
      if (!gin.system.found) {
        printError("Could not find SYSTEM block\n");
      }
      if (!gin.step.found) {
        printError("Could not find STEP block\n");
      }
      if (!gin.boundcond.found) {
        printError("Could not find BOUNDCOND block\n");
      }
      if (!gin.force.found) {
        printError("Could not find FORCE block\n");
      }
      if (gromosXX && !gin.constraint.found) {
        printError("Could not find CONSTRAINT block\n");
      }
      if (!gromosXX && !gin.geomconstraints.found) {
        printError("Could not find GEOMCONSTRAINTS block\n");
      }
      if (!gromosXX && !gin.neighbourlist.found) {
        printError("Could not find NEIGHBOURLIST block\n");
      }
      if (gromosXX && !gin.pairlist.found) {
        printError("Could not find PAIRLIST block\n");
      }
      if (!gin.nonbonded.found) {
        printError("Could not find NONBONDED block\n");
      }
      if (!gin.initialise.found) {
        printError("Could not find INITIALISE block\n");
      }

      // Check all variables for allowed values (no cross checks are done here)
      if (gin.barostat.found) {
        if (gin.barostat.ntp < 0 || gin.barostat.ntp > 3) {
          stringstream read;
          read << gin.barostat.ntp;
          printIO("BAROSTAT", "NPT", read.str(), "0..3");
        }
        if (gin.barostat.npvar < 0) {
          stringstream read;
          read << gin.barostat.npvar;
          printIO("BAROSTAT", "NPVAR", read.str(), ">=0");
        }
        // NPBTH is corrected when reading the input file
        // [see mk_script.h: istringstream & operator>>(istringstream &is, ibarostat &s)]
        if (gin.barostat.comp <= 0.0) {
          stringstream read;
          read << gin.barostat.comp;
          printIO("BAROSTAT", "COMP", read.str(), ">=0.0");
        }
        for (unsigned int i = 0; i < gin.barostat.pbaths.size(); i++) {
          if (gin.barostat.pbaths[i].prsbth < 0.0) {
            stringstream read, blockName;
            read << gin.barostat.pbaths[i].prsbth;
            blockName << "PRSBTH[" << i + 1 << "]";
            printIO("BAROSTAT", blockName.str(), read.str(), ">=0.0");
          }
          for (unsigned int j = 0; j < gin.barostat.pbaths[0].taubba.size(); j++) {
            if (gin.barostat.pbaths[0].taubba[j] <= 0.0) {
              stringstream read, blockName;
              read << gin.barostat.pbaths[i].taubba[j];
              blockName << "TAUBBA[" << i + 1 << "," << j + 1 << "]";
              printIO("BAROSTAT", blockName.str(), read.str(), ">0.0");
            }
          }
        }
        for (int i = 0; i < 6; i++) {
          if (gin.barostat.npcpl[i] < 0 || gin.barostat.npcpl[i] > (int) gin.barostat.pbaths.size()) {
            stringstream read, blockName;
            read << gin.barostat.npcpl[i];
            blockName << "NPCPL[" << i + 1 << "]";
            printIO("BAROSTAT", blockName.str(), read.str(), "0..NPBTH");
          }
        }
      }
      if (gin.boundcond.found) {
        if (gin.boundcond.ntb < -1 || gin.boundcond.ntb > 2) {
          stringstream read;
          read << gin.boundcond.ntb;
          printIO("BOUNDCOND", "NTB", read.str(), "-1..2");
        }
        if (gin.boundcond.ndfmin < 0) {
          stringstream read;
          read << gin.boundcond.ndfmin;
          printIO("BOUNDCOND", "NDFMIN", read.str(), ">=0");
        }
      }
      if (gin.bsleus.found) {
        if (gin.bsleus.write < 0) {
          stringstream read;
          read << gin.bsleus.write;
          printIO("BSLEUS", "WRITE", read.str(), ">= 0");
        }
        if (gin.bsleus.bsleus < 0 || gin.bsleus.bsleus > 1) {
          stringstream read;
          read << gin.bsleus.bsleus;
          printIO("BSLEUS", "BSLEUS", read.str(), "0..1");
        }
        if (gin.bsleus.build < 0 || gin.bsleus.build > 1) {
          stringstream read;
          read << gin.bsleus.build;
          printIO("BSLEUS", "BUILD", read.str(), "0..1");
        }
        if (! l_bsleus && ! l_jin){
          printError("No B&S-LEUS file given!\n");
        }
      }
      if (gin.cgrain.found) {
        if (gin.cgrain.ntcgran < 0 || gin.cgrain.ntcgran > 4) {
          stringstream read;
          read << gin.cgrain.ntcgran;
          printIO("CGRAIN", "NTCGRAN", read.str(), "0..4");
        }
        if (gin.cgrain.eps < 0.0) {
          stringstream read;
          read << gin.cgrain.eps;
          printIO("CGRAIN", "EPS", read.str(), ">=0.0");
        }
        if (gin.cgrain.epsm < 0.0) {
          stringstream read;
          read << gin.cgrain.epsm;
          printIO("CGRAIN", "EPSM", read.str(), ">=0.0");
        }
      }
      // no checks necessary for COMTRANSROT block
      if (gin.consistencycheck.found) {
        if (gin.consistencycheck.ntchk < 0 || gin.consistencycheck.ntchk > 1) {
          stringstream read;
          read << gin.consistencycheck.ntchk;
          printIO("CONSISTENCYCHECK", "NTCHK", read.str(), "0,1");
        }
        if (gin.consistencycheck.ntckf < 0 || gin.consistencycheck.ntckf > 1) {
          stringstream read;
          read << gin.consistencycheck.ntckf;
          printIO("CONSISTENCYCHECK", "NTCKF", read.str(), "0,1");
        }
        if (gin.consistencycheck.fdckf <= 0.0) {
          stringstream read;
          read << gin.consistencycheck.fdckf;
          printIO("CONSISTENCYCHECK", "FDCKF", read.str(), ">0.0");
        }
        if (gin.consistencycheck.ntckv < 0 || gin.consistencycheck.ntckv > 1) {
          stringstream read;
          read << gin.consistencycheck.ntckv;
          printIO("CONSISTENCYCHECK", "NTCKV", read.str(), "0,1");
        }
        if (gin.consistencycheck.fdckv <= 0.0) {
          stringstream read;
          read << gin.consistencycheck.fdckv;
          printIO("CONSISTENCYCHECK", "FDCKV", read.str(), ">0.0");
        }
        if (gin.consistencycheck.ntckt < 0 || gin.consistencycheck.ntckt > 1) {
          stringstream read;
          read << gin.consistencycheck.ntckt;
          printIO("CONSISTENCYCHECK", "NTCKT", read.str(), "0,1");
        }
        if (gin.consistencycheck.ntcke < 0 || gin.consistencycheck.ntcke > 1) {
          stringstream read;
          read << gin.consistencycheck.ntcke;
          printIO("CONSISTENCYCHECK", "NTCKE", read.str(), "0,1");
        }
        if (gin.consistencycheck.ntckr < 0 || gin.consistencycheck.ntckr > 1) {
          stringstream read;
          read << gin.consistencycheck.ntckr;
          printIO("CONSISTENCYCHECK", "NTCKR", read.str(), "0,1");
        }
        if (gin.consistencycheck.ntckl < 0 || gin.consistencycheck.ntckl > 1) {
          stringstream read;
          read << gin.consistencycheck.ntckl;
          printIO("CONSISTENCYCHECK", "NTCKL", read.str(), "0,1");
        }
        if (gin.consistencycheck.fdckl <= 0.0) {
          stringstream read;
          read << gin.consistencycheck.fdckl;
          printIO("CONSISTENCYCHECK", "FDCKL", read.str(), ">0.0");
        }
        // no need to check NACKF
        for (unsigned int i = 0; i < gin.consistencycheck.nckf.size(); i++) {
          if (gin.consistencycheck.nckf[i] < 1) {
            stringstream read, blockName;
            read << gin.consistencycheck.fdckl;
            blockName << "NCKF[" << i + 1 << "]";
            printIO("CONSISTENCYCHECK", blockName.str(), read.str(), ">=1");
          }
        }
      }
      if (gin.constraint.found) {
        if (gin.constraint.ntc < 0 || gin.constraint.ntc > 4) {
          stringstream read;
          read << gin.constraint.ntc;
          printIO("CONSTRAINT", "NTC", read.str(), "0..4");
        }
        if (gin.constraint.ntcp < 1 || gin.constraint.ntcp > 3) {
          stringstream read;
          read << gin.constraint.ntcp;
          printIO("CONSTRAINT", "NTCP", read.str(), "1..3");
        }
        if (gin.constraint.ntcp0[0] < 0) {
          stringstream read;
          read << gin.constraint.ntcp0[0];
          printIO("CONSTRAINT", "NTCP0[1]", read.str(), ">=0");
        }
        if (gin.constraint.ntcp == 3) {
          if (gin.constraint.ntcp0[1] < 0) {
            stringstream read;
            read << gin.constraint.ntcp0[1];
            printIO("CONSTRAINT", "NTCP0[2]", read.str(), ">=0");
          }
          if (gin.constraint.ntcp0[2] < 0) {
            stringstream read;
            read << gin.constraint.ntcp0[2];
            printIO("CONSTRAINT", "NTCP0[3]", read.str(), ">=0");
          }
        }
        if (gin.constraint.ntcs < 0 || gin.constraint.ntcs > 6) {
          stringstream read;
          read << gin.constraint.ntcs;
          printIO("CONSTRAINT", "NTCS", read.str(), "1..6");
        }
        if (gin.constraint.ntcs0[0] < 0) {
          stringstream read;
          read << gin.constraint.ntcs0[0];
          printIO("CONSTRAINT", "NTCS0[1]", read.str(), ">=0");
        }
        if (gin.constraint.ntcs == 6) {
          if (gin.constraint.ntcg < 1) {
            stringstream read;
            read << gin.constraint.ntcg;
            printIO("CONSTRAINT", "NTCG", read.str(), ">0");
          }
          for (int g = 0; g < gin.constraint.ntcg; g++) {
            if (gin.constraint.ntcd[g] < -1) {
              stringstream read;
              read << gin.constraint.ntcd[g];
              printIO("CONSTRAINT", "NTCD", read.str(), ">=-1");
            }
          }
        }
        if (gin.constraint.ntcs == 3) {
          if (gin.constraint.ntcs0[1] < 0) {
            stringstream read;
            read << gin.constraint.ntcs0[1];
            printIO("CONSTRAINT", "NTCS0[2]", read.str(), ">=0");
          }
          if (gin.constraint.ntcs0[2] < 0) {
            stringstream read;
            read << gin.constraint.ntcs0[2];
            printIO("CONSTRAINT", "NTCS0[3]", read.str(), ">=0");
          }
        }
      }
      if (gin.covalentform.found) {
        if (gin.covalentform.ntbbh < 0 || gin.covalentform.ntbbh > 1) {
          stringstream read;
          read << gin.covalentform.ntbbh;
          printIO("COVALENTFORM", "NTBBH", read.str(), "0,1");
        }
        if (gin.covalentform.ntbah < 0 || gin.covalentform.ntbah > 1) {
          stringstream read;
          read << gin.covalentform.ntbah;
          printIO("COVALENTFORM", "NTBAH", read.str(), "0,1");
        }
        if (gin.covalentform.ntbdn < 0 || gin.covalentform.ntbdn > 1) {
          stringstream read;
          read << gin.covalentform.ntbdn;
          printIO("COVALENTFORM", "NTBDN", read.str(), "0,1");
        }
      }
      // DEBUG has to be added
      if (gin.dihedralres.found) {
        if (gin.dihedralres.ntdlr < 0 || gin.dihedralres.ntdlr > 3) {
          stringstream read;
          read << gin.dihedralres.ntdlr;
          printIO("DIHEDRALRES", "NTDLR", read.str(), "0..3");
        }
        if (gin.dihedralres.cdlr < 0.0) {
          stringstream read;
          read << gin.dihedralres.cdlr;
          printIO("DIHEDRALRES", "CDLR", read.str(), ">=0.0");
        }
        if (gin.dihedralres.philin < 0.0 || gin.dihedralres.philin > 180.0) {
          stringstream read;
          read << gin.dihedralres.philin;
          printIO("DIHEDRALRES", "PHILIN", read.str(), "0..180");
        }
        if (gin.dihedralres.vdih < 0 || gin.dihedralres.vdih > 1) {
          stringstream read;
          read << gin.dihedralres.vdih;
          printIO("DIHEDRALRES", "VDIH", read.str(), "0,1");
        }
        if (gin.dihedralres.ntwdlr < 0) {
          stringstream read;
          read << gin.dihedralres.ntwdlr;
          printIO("DIHEDRALRES", "NTWDLR", read.str(), ">=0");
        }
        if (gin.dihedralres.toldac < 0) {
          stringstream read;
          read << gin.dihedralres.toldac;
          printIO("DIHEDRALRES", "TOLDAC", read.str(), ">=0");
        }
      }
      if (gin.angleres.found) {
        if (gin.angleres.ntalr < 0 || gin.angleres.ntalr > 6) {
          stringstream read;
          read << gin.angleres.ntalr;
          printIO("ANGLERES", "NTALR", read.str(), "0..6");
        }
        if (gin.angleres.calr < 0.0) {
          stringstream read;
          read << gin.angleres.calr;
          printIO("ANGLERES", "CALR", read.str(), ">=0.0");
        }
        if (gin.angleres.vares < 0 || gin.angleres.vares > 1) {
          stringstream read;
          read << gin.angleres.vares;
          printIO("ANGLERES", "VARES", read.str(), "0,1");
        }
        if (gin.angleres.ntwalr < 0) {
          stringstream read;
          read << gin.angleres.ntwalr;
          printIO("ANGLERES", "NTWALR", read.str(), ">=0");
        }
        if (gin.angleres.tolbac < 0) {
          stringstream read;
          read << gin.angleres.tolbac;
          printIO("ANGLERES", "TOLBAC", read.str(), ">=0");
        }
      }
      if(gin.distancefield.found) {
	if(gin.distancefield.ntdfr !=0 && gin.distancefield.ntdfr !=1){
	  stringstream read;
	  read<< gin.distancefield.ntdfr;
	  printIO("DISTANCEFIELD", "NTDFR", read.str(), "0,1");
	}
	if(gin.distancefield.grid <= 0.0){
	  stringstream read;
	  read<< gin.distancefield.grid;
	  printIO("DISTANCEFIELD", "GRID", read.str(), ">0.0");
	}
	if(gin.distancefield.proteinoffset <= 0.0){
	  stringstream read;
	  read<< gin.distancefield.grid;
	  printIO("DISTANCEFIELD", "PROTEINOFFSET", read.str(), ">0.0");
	}
	if(gin.distancefield.proteincutoff <= 0.0){
	  stringstream read;
	  read<< gin.distancefield.grid;
	  printIO("DISTANCEFIELD", "PROTEINCUTOFF", read.str(), ">0.0");
	}
	if(gin.distancefield.update <= 0){
	  stringstream read;
	  read<< gin.distancefield.grid;
	  printIO("DISTANCEFIELD", "UPDATE", read.str(), ">0");
	}
	if(gin.distancefield.smooth <= 0){
	  stringstream read;
	  read<< gin.distancefield.grid;
	  printIO("DISTANCEFIELD", "SMOOTH", read.str(), ">0");
	}
	if(gin.distancefield.rl <= 0.0){
	  stringstream read;
	  read<< gin.distancefield.grid;
	  printIO("DISTANCEFIELD", "RL", read.str(), ">0.0");
	}
	if(gin.distancefield.ntwdf < 0){
	  stringstream read;
	  read<< gin.distancefield.ntwdf;
	  printIO("DISTANCEFIELD", "NTWDF", read.str(), ">=0");
	}
	if(gin.distancefield.printgrid != 0 && gin.distancefield.printgrid != 1){
	  stringstream read;
	  read << gin.distancefield.printgrid;
	  printIO("DISTANCEFIELD", "PRINTGRID", read.str(), "0,1");
	}
        if(gin.distancefield.protect < 0){
	  stringstream read;
	  read << gin.distancefield.protect;
	  printIO("DISTANCEFIELD", "PROTECT", read.str(), "<0");
	}
      }
      
      if (gin.distanceres.found) {
        if (gin.distanceres.ntdir < -2 || gin.distanceres.ntdir > 3) {
          stringstream read;
          read << gin.distanceres.ntdir;
          printIO("DISTANCERES", "NTDIR", read.str(), "-2..3");
        }
        if (gin.distanceres.ntdira < 0 || gin.distanceres.ntdira > 1) {
          stringstream read;
          read << gin.distanceres.ntdira;
          printIO("DISTANCERES", "NTDIRA", read.str(), "0,1");
        }
        if (gin.distanceres.cdir < 0.0) {
          stringstream read;
          read << gin.distanceres.cdir;
          printIO("DISTANCERES", "CDIR", read.str(), ">=0.0");
        }
        if (gin.distanceres.dir0 < 0.0) {
          stringstream read;
          read << gin.distanceres.dir0;
          printIO("DISTANCERES", "DIR0", read.str(), ">=0.0");
        }
        if (gin.distanceres.taudir < 0.0) {
          stringstream read;
          read << gin.distanceres.taudir;
          printIO("DISTANCERES", "TAUDIR", read.str(), ">=0.0");
        }
        if (gin.distanceres.ntwdir < 0) {
          stringstream read;
          read << gin.distanceres.ntwdir;
          printIO("DISTANCERES", "NTWDIR", read.str(), ">=0");
        }
	if(gin.distanceres.vdir != 0 && gin.distanceres.vdir != 1){
	  stringstream read;
	  read << gin.distanceres.vdir;
	  printIO("DISTANCERES", "VDIR", read.str(), "0,1");
	}
	if(gin.distanceres.forcescale < 0 || gin.distanceres.forcescale > 2){
	  stringstream read;
	  read << gin.distanceres.forcescale;
	  printIO("DISTANCERES", "FORCESCALE", read.str(), "0..2");
	}
	
      }
      if (gin.colvarres.found) {
        if (gin.colvarres.cvr < 0 || gin.colvarres.cvr > 1) {
          stringstream read;
          read << gin.colvarres.cvr;
          printIO("COLVARRES", "CVR", read.str(), "0..1");
        }
        if (gin.colvarres.cvk < 0.0) {
          stringstream read;
          read << gin.colvarres.cvk;
          printIO("COLVARRES", "CVK", read.str(), ">=0.0");
        }
        if (gin.colvarres.taucvr < 0) {
          stringstream read;
          read << gin.colvarres.taucvr;
          printIO("COLVARRES", "TAUCVR", read.str(), ">=0.0");
        }
        if (gin.colvarres.vcvr < 0.0) {
          stringstream read;
          read << gin.colvarres.vcvr;
          printIO("COLVARRES", "VCVR", read.str(), "0..1");
        }
        if (gin.colvarres.ntwcv < 0 || gin.colvarres.cvr > 1) {
          stringstream read;
          read << gin.colvarres.ntwcv;
          printIO("COLVARRES", "NTWCV", read.str(), ">=0");
        }	
      }
      if (gin.eds.found) {
        if (gin.eds.eds < 0 || gin.eds.eds > 1) {
          stringstream read;
          read << gin.eds.eds;
          printIO("EDS", "EDS", read.str(), "0,1");
        }
        if (gin.eds.form < 1 || gin.eds.form > 3) {
          stringstream read;
          read << gin.eds.form;
          printIO("EDS", "FORM", read.str(), "1..3");
        }
        if (gin.eds.numstates <= 1) {
          stringstream read;
          read << gin.eds.numstates;
          printIO("EDS", "NUMSTATES", read.str(), ">1");
        }
        switch (gin.eds.form) {
          case 1:
            if (gin.eds.smooth[0] <= 0.0) {
              stringstream read;
              read << gin.eds.smooth[0];
              printIO("EDS", "S", read.str(), ">1");
            }
            break;
          case 2:
            for (int i = 0; i < gin.eds.numstates * (gin.eds.numstates - 1) / 2; i++) {
              if (gin.eds.smooth[i] <= 0.0) {
                stringstream read, blockName;
                read << gin.eds.smooth[i];
                blockName << "S[" << i + 1 << "]";
                printIO("EDS", blockName.str(), read.str(), ">0.0");
              }
            }
            break;
          case 3:
            for (int N = 0; N < gin.eds.numstates - 1; N++) {
              if (gin.eds.tree[N][0] < 1) {
                stringstream read, blockName;
                read << gin.eds.tree[N][0];
                blockName << "i[" << N + 1 << "]";
                printIO("EDS", blockName.str(), read.str(), ">1");
              }
              if (gin.eds.tree[N][1] < 1) {
                stringstream read, blockName;
                read << gin.eds.tree[N][1];
                blockName << "j[" << N + 1 << "]";
                printIO("EDS", blockName.str(), read.str(), ">1");
              }
              if (gin.eds.smooth[N] <= 0.0) {
                stringstream read, blockName;
                read << gin.eds.smooth[N];
                blockName << "S[" << N + 1 << "]";
                printIO("EDS", blockName.str(), read.str(), ">0.0");
              }
            }
            break;
          default:
            break;
        }
      }
      if (gin.aeds.found) {
        if (gin.aeds.aeds < 0 || gin.aeds.aeds > 1) {
          stringstream read;
          read << gin.aeds.aeds;
          printIO("AEDS", "AEDS", read.str(), "0,1");
        }
        if (gin.aeds.form < 1 || gin.aeds.form > 6) {
          stringstream read;
          read << gin.aeds.form;
          printIO("AEDS", "FORM", read.str(), "1..6");
        }
        if (gin.aeds.numstates <= 1) {
          stringstream read;
          read << gin.aeds.numstates;
          printIO("AEDS", "NUMSTATES", read.str(), ">1");
        }
        if (gin.aeds.ntiaedss < 0 || gin.aeds.ntiaedss > 1) {
          stringstream read;
          read << gin.aeds.ntiaedss;
          printIO("AEDS", "NTIAEDSS", read.str(), "0,1");
        }
        if (gin.aeds.restremin < 0 || gin.aeds.restremin > 1) {
          stringstream read;
          read << gin.aeds.restremin;
          printIO("AEDS", "RESTREMIN", read.str(), "0,1");
        }
        if (gin.aeds.bmaxtype < 1 || gin.aeds.bmaxtype > 2) {
          stringstream read;
          read << gin.aeds.restremin;
          printIO("AEDS", "BMAXTYPE", read.str(), "1,2");
        }
        if (gin.aeds.bmax <= 0.0) {
          stringstream read;
          read << gin.aeds.bmax;
          printIO("AEDS", "BMAX", read.str(), ">0.0");
        }
        if (gin.aeds.asteps <= 0) {
          stringstream read;
          read << gin.aeds.asteps;
          printIO("AEDS", "ASTEPS", read.str(), ">0");
        }
        if (gin.aeds.bsteps <= 0) {
          stringstream read;
          read << gin.aeds.bsteps;
          printIO("AEDS", "BSTEPS", read.str(), ">0");
        }
      }
      if (gin.energymin.found) {
        if (gin.energymin.ntem < 0 || gin.energymin.ntem > 3) {
          stringstream read;
          read << gin.energymin.ntem;
          printIO("ENERGYMIN", "NTEM", read.str(), "0..3");
        }
        if (gin.energymin.ncyc < 0 && gin.energymin.ntem > 1) {
          stringstream read;
          read << gin.energymin.ncyc;
          printIO("ENERGYMIN", "NCYC", read.str(), ">=0");
        }
        if (gin.energymin.dele == 0) {
          stringstream read;
          read << gin.energymin.dele;
          printIO("ENERGYMIN", "DELE", read.str(), ">0 or <0");
        }
        if (gin.energymin.dx0 <= 0.0) {
          stringstream read;
          read << gin.energymin.dx0;
          printIO("ENERGYMIN", "DX0", read.str(), ">0.0");
        }
        if (gin.energymin.dxm <= 0.0) {
          stringstream read;
          read << gin.energymin.dxm;
          printIO("ENERGYMIN", "DXM", read.str(), ">0.0");
        }
        if (gin.energymin.nmin <= 0) {
          stringstream read;
          read << gin.energymin.nmin;
          printIO("ENERGYMIN", "NMIN", read.str(), ">0");
        }
        if (gin.energymin.flim < 0) {
          stringstream read;
          read << gin.energymin.flim;
          printIO("ENERGYMIN", "FLIM", read.str(), ">=0.0");
        }
        if (gin.energymin.ntem > 1 && gin.energymin.cgim < 0) {
          stringstream read;
          read << gin.energymin.flim;
          printIO("ENERGYMIN", "CGIM", read.str(), ">0");
        }
        if (gin.energymin.ntem > 1 && gin.energymin.cgic < 0.0) {
          stringstream read;
          read << gin.energymin.flim;
          printIO("ENERGYMIN", "CGIM", read.str(), ">=0.0");
        }
      }
      // EWARN does not need any checking here
      if (gin.force.found) {
        for (int i = 0; i < 6; i++) {
          if (gin.force.ntf[i] < 0 || gin.force.ntf[i] > 1) {
            stringstream read, blockName;
            read << gin.force.ntf[i];
            blockName << "NTF[" << i + 1 << "]";
            printIO("FORCE", blockName.str(), read.str(), "0,1");
          }
        }
        for (unsigned int i = 0; i < gin.force.nre.size(); i++) {
          if (gin.force.nre[i] < 1) {
            stringstream read, blockName;
            read << gin.force.nre[i];
            blockName << "NRE[" << i + 1 << "]";
            printIO("FORCE", blockName.str(), read.str(), ">=1");
          }
        }
      }
      if (gin.gamd.found) {
        if (gin.gamd.gamd < 0 || gin.gamd.gamd > 1) {
          stringstream read;
          read << gin.gamd.gamd;
          printIO("GAMD", "GAMD", read.str(), "0,1");
        }
        if (gin.gamd.search < 0 || gin.gamd.search > 2) {
          stringstream read;
          read << gin.gamd.search;
          printIO("GAMD", "SEARCH", read.str(), "0..2");
        }
        if (gin.gamd.form < 1 || gin.gamd.form > 3) {
          stringstream read;
          read << gin.gamd.form;
          printIO("GAMD", "FORM", read.str(), "1..3");
        }
        if (gin.gamd.thresh < 1 || gin.gamd.thresh > 2) {
          stringstream read;
          read << gin.gamd.thresh;
          printIO("GAMD", "THRESH", read.str(), "1,2");
        }
        if (gin.gamd.ntigamds < 0 || gin.gamd.ntigamds > 1) {
          stringstream read;
          read << gin.gamd.ntigamds;
          printIO("GAMD", "NTIGAMDS", read.str(), "0,1");
        }
        if (gin.gamd.agroups <= 0) {
          stringstream read;
          read << gin.gamd.agroups;
          printIO("GAMD", "AGROUPS", read.str(), ">0");
        }
        if (gin.gamd.igroups <= 0) {
          stringstream read;
          read << gin.gamd.igroups;
          printIO("GAMD", "IGROUPS", read.str(), ">0");
        }
        if (gin.gamd.dihstd <= 0.0) {
          stringstream read;
          read << gin.gamd.dihstd;
          printIO("GAMD", "DIHSTD", read.str(), ">0.0");
        }
        if (gin.gamd.totstd <= 0.0) {
          stringstream read;
          read << gin.gamd.totstd;
          printIO("GAMD", "TOTSTD", read.str(), ">0.0");
        }
        if (gin.gamd.eqsteps < 0) {
          stringstream read;
          read << gin.gamd.eqsteps;
          printIO("GAMD", "EQSTEPS", read.str(), ">=0");
        }
        if (gin.gamd.window < 0) {
          stringstream read;
          read << gin.gamd.window;
          printIO("GAMD", "WINDOW", read.str(), ">0");
        }
      }
      if (gin.geomconstraints.found) {
        if (gin.geomconstraints.ntcph < 0 || gin.geomconstraints.ntcph > 1) {
          stringstream read;
          read << gin.geomconstraints.ntcph;
          printIO("GEOMCONSTRAINTS", "NTCPH", read.str(), "0,1");
        }
        if (gin.geomconstraints.ntcpn < 0 || gin.geomconstraints.ntcpn > 1) {
          stringstream read;
          read << gin.geomconstraints.ntcpn;
          printIO("GEOMCONSTRAINTS", "NTCPN", read.str(), "0,1");
        }
        if (gin.geomconstraints.ntcs < 0 || gin.geomconstraints.ntcs > 1) {
          stringstream read;
          read << gin.geomconstraints.ntcs;
          printIO("GEOMCONSTRAINTS", "NTCS", read.str(), "0,1");
        }
        if (gin.geomconstraints.shktol <= 0.0) {
          stringstream read;
          read << gin.geomconstraints.shktol;
          printIO("GEOMCONSTRAINTS", "SHKTOL", read.str(), ">0.0");
        }
      }
      if (gin.gromos96compat.found) {
        if (gin.gromos96compat.ntnb96 < 0 || gin.gromos96compat.ntnb96 > 1) {
          stringstream read;
          read << gin.gromos96compat.ntnb96;
          printIO("CROMOS96COMPAT", "NTNB96", read.str(), "0,1");
        }
        if (gin.gromos96compat.ntr96 < 0 || gin.gromos96compat.ntr96 > 1) {
          stringstream read;
          read << gin.gromos96compat.ntr96;
          printIO("CROMOS96COMPAT", "NTR96", read.str(), "0,1");
        }
        if (gin.gromos96compat.ntp96 < 0 || gin.gromos96compat.ntp96 > 1) {
          stringstream read;
          read << gin.gromos96compat.ntp96;
          printIO("CROMOS96COMPAT", "NTP96", read.str(), "0,1");
        }
        if (gin.gromos96compat.ntg96 < 0 || gin.gromos96compat.ntg96 > 1) {
          stringstream read;
          read << gin.gromos96compat.ntg96;
          printIO("CROMOS96COMPAT", "NTG96", read.str(), "0,1");
        }
      }
      if (gin.initialise.found) {
        if (gin.initialise.ntivel < 0 || gin.initialise.ntivel > 1) {
          stringstream read;
          read << gin.initialise.ntivel;
          printIO("INITIALISE", "NTIVEL", read.str(), "0,1");
        }
        if (gin.initialise.ntishk < 0 || gin.initialise.ntishk > 3) {
          stringstream read;
          read << gin.initialise.ntishk;
          printIO("INITIALISE", "NTISHK", read.str(), "0..3");
        }
        if (gin.initialise.ntinht < 0 || gin.initialise.ntinht > 1) {
          stringstream read;
          read << gin.initialise.ntinht;
          printIO("INITIALISE", "NTINHT", read.str(), "0,1");
        }
        if (gin.initialise.ntinhb < 0 || gin.initialise.ntinhb > 1) {
          stringstream read;
          read << gin.initialise.ntinhb;
          printIO("INITIALISE", "NTINHB", read.str(), "0,1");
        }
        if (gin.initialise.ntishi < 0 || gin.initialise.ntishi > 1) {
          stringstream read;
          read << gin.initialise.ntishi;
          printIO("INITIALISE", "NTISHI", read.str(), "0,1");
        }
        if (gin.initialise.ntirtc < 0 || gin.initialise.ntirtc > 1) {
          stringstream read;
          read << gin.initialise.ntirtc;
          printIO("INITIALISE", "NTIRTC", read.str(), "0,1");
        }
        if (gin.initialise.nticom < 0 || gin.initialise.nticom > 3) {
          stringstream read;
          read << gin.initialise.nticom;
          printIO("INITIALISE", "NTICOM", read.str(), "0..3");
        }
        if (gin.initialise.ntisti < 0 || gin.initialise.ntisti > 1) {
          stringstream read;
          read << gin.initialise.ntisti;
          printIO("INITIALISE", "NTISTI", read.str(), "0,1");
        }
        if (gin.initialise.ig <= 0) {
          stringstream read;
          read << gin.initialise.ig;
          printIO("INITIALISE", "IG", read.str(), ">0");
        }
        if (gin.initialise.tempi < 0.0) {
          stringstream read;
          read << gin.initialise.tempi;
          printIO("INITIALISE", "TEMPI", read.str(), ">=0.0");
        }
      }
      if (gin.innerloop.found) {
        if (gin.innerloop.ntilm < 0 || gin.innerloop.ntilm > 4) {
          stringstream read;
          read << gin.innerloop.ntilm;
          printIO("INNERLOOP", "NTILM", read.str(), "0..4");
        }
        if (gin.innerloop.ntils < 0 || gin.innerloop.ntils > 1) {
          stringstream read;
          read << gin.innerloop.ntils;
          printIO("INNERLOOP", "NTILS", read.str(), "0,1");
        }
        // if method 4
        // if a gpu should be used also check if the device numbers are valid
        if (gin.innerloop.ntilm == 4) {
            if (gin.innerloop.ngpus < 1) {
                stringstream read;
                read << gin.innerloop.ngpus;
                printIO("INNERLOOP", "NGPUS", read.str(), ">0");
            }
            for (int g = 0; g < gin.innerloop.ngpus; g++) {
              if (gin.innerloop.ndevg[g] < -1) {
                stringstream read;
                read << gin.innerloop.ndevg[g];
                printIO("INNERLOOP", "NDEVG", read.str(), ">=-1");
              }
            }
        }
      } // INNERLOOP END
      if (gin.integrate.found) {
        if (gin.integrate.nint < 0 || gin.integrate.nint > 1) {
          stringstream read;
          read << gin.integrate.nint;
          printIO("INTEGRATE", "NINT", read.str(), "0,1");
        }
      }
      if (gin.jvalueres.found) {
        if (gin.jvalueres.ntjvr < -3 || gin.jvalueres.ntjvr > 2) {
          stringstream read;
          read << gin.jvalueres.ntjvr;
          printIO("JVALUERES", "NTJVR", read.str(), "-3..2");
        }
        if (gin.jvalueres.ntjvra < 0 || gin.jvalueres.ntjvra > 1) {
          stringstream read;
          read << gin.jvalueres.ntjvra;
          printIO("JVALUERES", "NTJVRA", read.str(), "0,1");
        }
        if (gin.jvalueres.cjvr < 0.0) {
          stringstream read;
          read << gin.jvalueres.cjvr;
          printIO("JVALUERES", "CJVR", read.str(), ">=0.0");
        }
        if (gin.jvalueres.taujvr < 0.0) {
          stringstream read;
          read << gin.jvalueres.taujvr;
          printIO("JVALUERES", "TAUJVR", read.str(), ">=0.0");
        }
        if (gin.jvalueres.njvrtars < 0 || gin.jvalueres.njvrtars > 1){
          stringstream read;
          read << gin.jvalueres.njvrtars;
          printIO("JVALUERES", "NJVRTARS", read.str(), "0,1");
        }
        if (gin.jvalueres.njvrbiqw < 0 || gin.jvalueres.njvrbiqw > 2){
          stringstream read;
          read << gin.jvalueres.njvrbiqw;
          printIO("JVALUERES", "NJVRBIQW", read.str(), "0..2");
        }
        if (gin.jvalueres.le < 0 || gin.jvalueres.le > 1) {
          stringstream read;
          read << gin.jvalueres.le;
          printIO("JVALUERES", "LE", read.str(), "0,1");
        }
        if (gin.jvalueres.ngrid <= 0) {
          stringstream read;
          read << gin.jvalueres.ngrid;
          printIO("JVALUERES", "NGRID", read.str(), ">=0");
        }
        if (gin.jvalueres.delta < 0.0) {
          stringstream read;
          read << gin.jvalueres.delta;
          printIO("JVALUERES", "DELTA", read.str(), ">=0.0");
        }
        if (gin.jvalueres.write < 0) {
          stringstream read;
          read << gin.jvalueres.write;
          printIO("JVALUERES", "NTWJV", read.str(), ">=0.0");
        }
      }
      if (gin.orderparamres.found) {
        if (gin.orderparamres.ntopr < -2 || gin.orderparamres.ntopr > 2) {
          stringstream read;
          read << gin.orderparamres.ntopr;
          printIO("ORDERPARAMRES", "NTOPR", read.str(), "-2..2");
        }
        if (gin.orderparamres.ntopra < 0 || gin.orderparamres.ntopra > 1) {
          stringstream read;
          read << gin.orderparamres.ntopra;
          printIO("ORDERPARAMRES", "NTOPRA", read.str(), "0,1");
        }
        if (gin.orderparamres.copr < 0.0) {
          stringstream read;
          read << gin.orderparamres.copr;
          printIO("ORDERPARAMRES", "COPR", read.str(), ">=0.0");
        }
        if (gin.orderparamres.tauopr < 0.0) {
          stringstream read;
          read << gin.orderparamres.tauopr;
          printIO("ORDERPARAMRES", "TAUOPR", read.str(), ">=0.0");
        }
        if (gin.orderparamres.updopr < 1) {
          stringstream read;
          read << gin.orderparamres.updopr;
          printIO("ORDERPARAMRES", "UPDOPR", read.str(), ">0");
        }
        if (gin.orderparamres.ntwop < 0) {
          stringstream read;
          read << gin.orderparamres.ntwop;
          printIO("ORDERPARAMRES", "NTWOP", read.str(), ">=0");
        }
      }
      if (gin.symres.found) {
        if (gin.symres.ntsym < 0 || gin.symres.ntsym > 2) {
          stringstream read;
          read << gin.symres.ntsym;
          printIO("SYMRES", "NTSYM", read.str(), "0..2");
        }
        if (gin.symres.csym < 0.0) {
          stringstream read;
          read << gin.symres.ntsym;
          printIO("SYMRES", "CSYM", read.str(), ">=0.0");
        }
      }
      if (gin.lambdas.found) {
        if (gin.lambdas.ntil < 0 || gin.lambdas.ntil > 1) {
          stringstream read;
          read << gin.lambdas.ntil;
          printIO("LAMBDAS", "NTIL", read.str(), "0,1");
        }
        for (unsigned int i = 0; i < gin.lambdas.lambints.size(); i++) {
          if (gin.lambdas.lambints[i].ntli < 1 || gin.lambdas.lambints[i].ntli > 13) {
            stringstream read, blockName;
            read << gin.lambdas.lambints[i].ntli;
            blockName << "NTLI[" << i + 1 << "]";
            printIO("LAMBDAS", blockName.str(), read.str(), "1..13");
          }
          if (gin.lambdas.lambints[i].nilg1 < 0) {
            stringstream read, blockName;
            read << gin.lambdas.lambints[i].nilg1;
            blockName << "NILG1[" << i + 1 << "]";
            printIO("LAMBDAS", blockName.str(), read.str(), ">=0");
          }
          if (gin.lambdas.lambints[i].nilg2 < 0) {
            stringstream read, blockName;
            read << gin.lambdas.lambints[i].nilg2;
            blockName << "NILG2[" << i + 1 << "]";
            printIO("LAMBDAS", blockName.str(), read.str(), ">=0");
          }
          // ALI, BLI, CLI, DLI and ELI does not need to be ckecked, can be
          // any double
        }
      }
      if (gin.localelev.found) {
        if (gin.localelev.ntles < 0 || gin.localelev.ntles > 5) {
          stringstream read;
          read << gin.localelev.ntles;
          printIO("LOCALELEV", "NTLES", read.str(), "0,5");
        }
        if (gin.localelev.nlepot != (int) gin.localelev.nlepid_ntlepfr.size()) {
          stringstream read, msg;
          read << gin.localelev.nlepot;
          msg << "NLEPOT is " << read.str() << " but "
                  << gin.localelev.nlepid_ntlepfr.size() << " potential(s) "
                  "(with different ID) are listed.";
          printErrMsg("LOCALELEV", "NLEPOT", msg.str());
        }
        if (gin.localelev.ntlesa < 0 || gin.localelev.ntlesa > 2) {
          stringstream read;
          read << gin.localelev.ntlesa;
          printIO("LOCALELEV", "NTLESA", read.str(), "0..2");
        }
        // IDs of potentials does not have to be checked
        int i = 0;
        for (map<int, int>::iterator it = gin.localelev.nlepid_ntlepfr.begin();
                it != gin.localelev.nlepid_ntlepfr.end(); it++) {
          i++;
          if (it->second < 0 || it->second > 1) {
            stringstream read, blockName;
            read << it->second;
            blockName << "NTLEFR[" << i << "]";
            printIO("LOCALELEV", blockName.str(), read.str(), "0..1");
          }
        }
      }
      if (gin.electric.found) {
        if (gin.electric.field < 0 || gin.electric.field > 1) {
          stringstream read;
          read << gin.electric.field;
          printIO("ELECTRIC", "FIELD", read.str(), "0..1");
        }
        if (gin.electric.dipole < 0 || gin.electric.dipole > 1) {
          stringstream read;
          read << gin.electric.dipole;
          printIO("ELECTRIC", "DIPOLE", read.str(), "0..1");
        }
        if (gin.electric.current < 0 || gin.electric.current > 1) {
          stringstream read;
          read << gin.electric.current;
          printIO("ELECTRIC", "CURRENT", read.str(), "0..1");
        }
        if (gin.electric.dipole > 0){
          if (gin.electric.dipgrp < 0 || gin.electric.dipgrp > 2) {
          stringstream read;
          read << gin.electric.dipgrp;
          printIO("ELECTRIC", "DIPGRP", read.str(), "0..2");
          }
          if (gin.electric.ntwdip < 1) {
          stringstream read;
          read << gin.electric.ntwdip;
          printIO("ELECTRIC", "NTWDIP", read.str(), ">=1");
          }
        }
        if (gin.electric.current > 0){
          if (gin.electric.ntwcur < 1){
            stringstream read;
            read << gin.electric.ntwcur;
            printIO("ELECTRIC", "NTWCUR", read.str(), ">=1");
          }
          if (gin.electric.ncurgrp < 1){
            stringstream read;
            read << gin.electric.ncurgrp;
            printIO("ELECTRIC", "NCURGRP", read.str(), ">=1");
          }
          for (int i=0; i< gin.electric.ncurgrp; ++i){
            if (gin.electric.curgrp[i] < 0){
              stringstream read;
              read << gin.electric.curgrp[i];
              printIO("ELECTRIC", "CURGRP", read.str(), ">=0");
            }
          }
        }
        
      }
      if (gin.multibath.found) {
        if (gin.multibath.ntbtyp < 0 || gin.multibath.ntbtyp > 2) {
          stringstream read;
          read << gin.multibath.ntbtyp;
          printIO("MULTIBATH", "NTBTYP", read.str(), "0..2");
        }
        if (gin.multibath.ntbtyp == 2) {
          if (gin.multibath.num < 0) {
            stringstream read;
            read << gin.multibath.num;
            printIO("MULTIBATH", "NUM", read.str(), ">=0");
          }
        }
        if (gin.multibath.nbaths < 0) {
          stringstream read;
          read << gin.multibath.nbaths;
          printIO("MULTIBATH", "NBATHS", read.str(), ">=0");
        }
        for (int i = 0; i < gin.multibath.nbaths; i++) {
          if (gin.multibath.temp0[i] < 0.0) {
            stringstream read, blockName;
            read << gin.multibath.temp0[i];
            blockName << "TEMP0[" << i + 1 << "]";
            printIO("MULTIBATH", blockName.str(), read.str(), ">=0.0");
          }
          if (gin.multibath.tau[i] < 0.0) {
            stringstream read, blockName;
            read << gin.multibath.tau[i];
            blockName << "TAU[" << i + 1 << "]";
            printIO("MULTIBATH", blockName.str(), read.str(), ">=0.0");
          }
        }
        if (gin.multibath.dofset < 0) {
          stringstream read;
          read << gin.multibath.dofset;
          printIO("MULTIBATH", "DOFSET", read.str(), ">=0");
        }
        for (int i = 0; i < gin.multibath.dofset; i++) {
          if (gin.multibath.combath[i] < 1) {
            stringstream read, blockName;
            read << gin.multibath.combath[i];
            blockName << "COM-BATH[" << i + 1 << "]";
            printIO("MULTIBATH", blockName.str(), read.str(), ">=1");
          }
          if (gin.multibath.irbath[i] < 1) {
            stringstream read, blockName;
            read << gin.multibath.irbath[i];
            blockName << "IR-BATH[" << i + 1 << "]";
            printIO("MULTIBATH", blockName.str(), read.str(), ">=1");
          }
        }
      }
      if (gin.multicell.found) {
        if (gin.multicell.ntm < 0 || gin.multicell.ntm > 1) {
          stringstream read;
          read << gin.multicell.ntm;
          printIO("MULTICELL", "NTM", read.str(), "0..1");
        }
        if (gin.multicell.ncella < 1) {
          stringstream read;
          read << gin.multicell.ncella;
          printIO("MULTICELL", "NCELLA", read.str(), ">=1");
        }
        if (gin.multicell.ncellb < 1) {
          stringstream read;
          read << gin.multicell.ncellb;
          printIO("MULTICELL", "NCELLB", read.str(), ">=1");
        }
        if (gin.multicell.ncellc < 1) {
          stringstream read;
          read << gin.multicell.ncellc;
          printIO("MULTICELL", "NCELLC", read.str(), ">=1");
        }
        if (gin.multicell.tolpx < 0.0) {
          stringstream read;
          read << gin.multicell.tolpx;
          printIO("MULTICELL", "TOLPX", read.str(), ">=0.0");
        }
        if (gin.multicell.tolpv < 0.0) {
          stringstream read;
          read << gin.multicell.tolpv;
          printIO("MULTICELL", "TOLPV", read.str(), ">=0.0");
        }
        if (gin.multicell.tolpf < 0.0) {
          stringstream read;
          read << gin.multicell.tolpf;
          printIO("MULTICELL", "TOLPF", read.str(), ">=0.0");
        }
        if (gin.multicell.tolpfw < 0.0) {
          stringstream read;
          read << gin.multicell.tolpfw;
          printIO("MULTICELL", "TOLPFW", read.str(), ">=0.0");
        }
      }
      if (gin.multigradient.found) {
        if (gin.multigradient.ntmgre < 0 || gin.multigradient.ntmgre > 1) {
          stringstream read;
          read << gin.multigradient.ntmgre;
          printIO("MULTIGRADIENT", "NTMGRE", read.str(), "0,1");
        }
        if (gin.multigradient.ntmgrp < 0 || gin.multigradient.ntmgrp > 3) {
          stringstream read;
          read << gin.multigradient.ntmgrp;
          printIO("MULTIGRADIENT", "NTMGRP", read.str(), "0-3");
        }
        if (gin.multigradient.ntmgrn < 0) {
          stringstream read;
          read << gin.multigradient.ntmgrn;
          printIO("MULTIGRADIENT", "NTMGRN", read.str(), ">=0");
        }
      }
      if (gin.multistep.found) {
        if (gin.multistep.steps < 1) {
          stringstream read;
          read << gin.multistep.steps;
          printIO("MULTISTEP", "STEPS", read.str(), ">=1");
        }
        if (gin.multistep.boost < 0 || gin.multistep.boost > 1) {
          stringstream read;
          read << gin.multistep.boost;
          printIO("MULTISTEP", "BOOST", read.str(), "0..1");
        }
      }
      if (gin.neighbourlist.found) {
        if (gin.neighbourlist.plalgo < 0 || gin.neighbourlist.plalgo > 2) {
          stringstream read;
          read << gin.neighbourlist.plalgo;
          printIO("NEIGHBOURLIST", "PLALGO", read.str(), "0..2");
        }
        if (gin.neighbourlist.nupdpl <= 0) {
          stringstream read;
          read << gin.neighbourlist.nupdpl;
          printIO("NEIGHBOURLIST", "NUPDPL", read.str(), ">0");
        }
        if (gin.neighbourlist.nupdis <= 0) {
          stringstream read;
          read << gin.neighbourlist.nupdis;
          printIO("NEIGHBOURLIST", "NUPDIS", read.str(), ">0");
        }
        if (gin.neighbourlist.nupdii <= 0) {
          stringstream read;
          read << gin.neighbourlist.nupdii;
          printIO("NEIGHBOURLIST", "NUPDII", read.str(), ">0");
        }
        if (gin.neighbourlist.rcuts < 0.0) {
          stringstream read;
          read << gin.neighbourlist.rcuts;
          printIO("NEIGHBOURLIST", "RCUTS", read.str(), ">= 0.0");
        }
        if (gin.neighbourlist.gridszx < 0.0) {
          stringstream read;
          read << gin.neighbourlist.gridszx;
          printIO("NEIGHBOURLIST", "GRIDSZX", read.str(), ">= 0.0");
        }
        if (gin.neighbourlist.gridszy < 0.0) {
          stringstream read;
          read << gin.neighbourlist.gridszy;
          printIO("NEIGHBOURLIST", "GRIDSZY", read.str(), ">= 0.0");
        }
        if (gin.neighbourlist.gridszz < 0.0) {
          stringstream read;
          read << gin.neighbourlist.gridszz;
          printIO("NEIGHBOURLIST", "GRIDSZZ", read.str(), ">= 0.0");
        }
        if (gin.neighbourlist.type < 0 || gin.neighbourlist.type > 1) {
          stringstream read;
          read << gin.neighbourlist.type;
          printIO("NEIGHBOURLIST", "TYPE", read.str(), "0..1");
        }
        if (gin.neighbourlist.ncgcen < -2) {
          stringstream read;
          read << gin.neighbourlist.ncgcen;
          printIO("NEIGHBOURLIST", "NCGCEN", read.str(), ">=-2");
        }
      }
      if (gin.nonbonded.found) {
        if (gin.nonbonded.nlrele < -4 || gin.nonbonded.nlrele > 4) {
          stringstream read;
          read << gin.nonbonded.nlrele;
          printIO("NONBONDED", "NLRELE", read.str(), "-4..4");
        }
        if (gin.nonbonded.appak < 0.0) {
          stringstream read;
          read << gin.nonbonded.appak;
          printIO("NONBONDED", "APPAK", read.str(), ">=0.0");
        }
        if (gin.nonbonded.rcrf < 0.0) {
          stringstream read;
          read << gin.nonbonded.rcrf;
          printIO("NONBONDED", "RCRF", read.str(), ">=0.0");
        }
        if (gin.nonbonded.epsrf != 0.0 && gin.nonbonded.epsrf < 1.0) {
          stringstream read;
          read << gin.nonbonded.epsrf;
          printIO("NONBONDED", "EPSRF", read.str(), "0.0 or >=1.0");
        }
        if (!(gin.nonbonded.nslfexcl == 0 || gin.nonbonded.nslfexcl == 1)){
          stringstream read;
          read << gin.nonbonded.nslfexcl;
          printIO("NONBONDED", "NSLFEXCL", read.str(), "0 or 1 (default)");
        }
        if (gin.nonbonded.nshape < -1 || gin.nonbonded.nshape > 10) {
          stringstream read;
          read << gin.nonbonded.nshape;
          printIO("NONBONDED", "NSHAPE", read.str(), "-1..10");
        }
        if (gin.nonbonded.ashape <= 0.0) {
          stringstream read;
          read << gin.nonbonded.ashape;
          printIO("NONBONDED", "ASHAPE", read.str(), ">0.0");
        }
        if (gin.nonbonded.na2clc < 0 || gin.nonbonded.na2clc > 4) {
          stringstream read;
          read << gin.nonbonded.na2clc;
          printIO("NONBONDED", "NA2CLC", read.str(), "0..4");
        }
        if (gin.nonbonded.tola2 <= 0.0) {
          stringstream read;
          read << gin.nonbonded.tola2;
          printIO("NONBONDED", "TOLA2", read.str(), ">0.0");
        }
        if (gin.nonbonded.epsls != 0.0 && gin.nonbonded.epsls < 1.0) {
          stringstream read;
          read << gin.nonbonded.epsls;
          printIO("NONBONDED", "EPSLS", read.str(), "0.0 or >=1.0");
        }
        if (gin.nonbonded.nkx <= 0) {
          stringstream read;
          read << gin.nonbonded.nkx;
          printIO("NONBONDED", "NKX", read.str(), ">0");
        }
        if (gin.nonbonded.nky <= 0) {
          stringstream read;
          read << gin.nonbonded.nky;
          printIO("NONBONDED", "NKY", read.str(), ">0");
        }
        if (gin.nonbonded.nkz <= 0) {
          stringstream read;
          read << gin.nonbonded.nkz;
          printIO("NONBONDED", "NKZ", read.str(), ">0");
        }
        if (gin.nonbonded.kcut <= 0.0) {
          stringstream read;
          read << gin.nonbonded.kcut;
          printIO("NONBONDED", "KCUT", read.str(), ">0.0");
        }
        if (gin.nonbonded.ngx <= 0) {
          stringstream read;
          read << gin.nonbonded.ngx;
          printIO("NONBONDED", "NGX", read.str(), ">0");
        }
        if (gin.nonbonded.ngy <= 0) {
          stringstream read;
          read << gin.nonbonded.ngy;
          printIO("NONBONDED", "NGY", read.str(), ">0");
        }
        if (gin.nonbonded.ngz <= 0) {
          stringstream read;
          read << gin.nonbonded.ngz;
          printIO("NONBONDED", "NGZ", read.str(), ">0");
        }
        if (gin.nonbonded.nasord < 1 || gin.nonbonded.nasord > 5) {
          stringstream read;
          read << gin.nonbonded.nasord;
          printIO("NONBONDED", "NASORD", read.str(), "1..5");
        }
        if (gin.nonbonded.nfdord < 0 || gin.nonbonded.nfdord > 5) {
          stringstream read;
          read << gin.nonbonded.nfdord;
          printIO("NONBONDED", "NFDORD", read.str(), "0..5");
        }
        if (gin.nonbonded.nspord < 0.0) {
          stringstream read;
          read << gin.nonbonded.nspord;
          printIO("NONBONDED", "NSPORD", read.str(), ">=0.0");
        }
        if (gin.nonbonded.nqeval < 0) {
          stringstream read;
          read << gin.nonbonded.nqeval;
          printIO("NONBONDED", "NQEVAL", read.str(), ">=0");
        }
        if (gin.nonbonded.faccur <= 0.0) {
          stringstream read;
          read << gin.nonbonded.faccur;
          printIO("NONBONDED", "FACCUR", read.str(), ">0.0");
        }
        if (gin.nonbonded.nrdgrd < 0 || gin.nonbonded.nrdgrd > 1) {
          stringstream read;
          read << gin.nonbonded.nrdgrd;
          printIO("NONBONDED", "NRDGRD", read.str(), "0..1");
        }
        if (gin.nonbonded.nwrgrd < 0 || gin.nonbonded.nwrgrd > 1) {
          stringstream read;
          read << gin.nonbonded.nwrgrd;
          printIO("NONBONDED", "NWRGRD", read.str(), "0..1");
        }
        if (gin.nonbonded.nlrlj < 0 || gin.nonbonded.nlrlj > 1) {
          stringstream read;
          read << gin.nonbonded.nlrlj;
          printIO("NONBONDED", "NLRLJ", read.str(), "0..1");
        }
        if (gin.nonbonded.slvdns <= 0.0) {
          stringstream read;
          read << gin.nonbonded.slvdns;
          printIO("NONBONDED", "SLVDNS", read.str(), ">0.0");
        }
      }
      if (gin.overalltransrot.found) {
        if (gin.overalltransrot.ncmtr < 0 || gin.overalltransrot.ncmtr > 1) {
          stringstream read;
          read << gin.overalltransrot.ncmtr;
          printIO("OVERALLTRANSROT", "NCMTR", read.str(), "0,1");
        }
        if (gin.overalltransrot.ncmro < 0 || gin.overalltransrot.ncmro > 1) {
          stringstream read;
          read << gin.overalltransrot.ncmro;
          printIO("OVERALLTRANSROT", "NCMRO", read.str(), "0,1");
        }
        // CMAMAX, CMAMAY and CMAMY have still to be checked here
      }
      if (gin.pairlist.found) {
        if (gin.pairlist.algorithm < 0 || gin.pairlist.algorithm > 2) {
          stringstream read;
          read << gin.pairlist.algorithm;
          printIO("PAIRLIST", "ALGORITHM", read.str(), "0,1,2");
        }
        if (gin.pairlist.nsnb <= 0) {
          stringstream read;
          read << gin.pairlist.nsnb;
          printIO("PAIRLIST", "NSNB", read.str(), ">0");
        }
        if (gin.pairlist.rcutp <= 0.0) {
          stringstream read;
          read << gin.pairlist.rcutp;
          printIO("PAIRLIST", "RCUTP", read.str(), ">0.0");
        }
        if (gin.pairlist.rcutl <= 0.0) {
          stringstream read;
          read << gin.pairlist.rcutl;
          printIO("PAIRLIST", "RCUTL", read.str(), ">0.0");
        }
        if (gin.pairlist.size <= 0.0) {
          stringstream read;
          read << gin.pairlist.size;
          printIO("PAIRLIST", "SIZE", read.str(), ">0.0");
        }
        if (gin.pairlist.type < 0 || gin.pairlist.type > 2) {
          stringstream read;
          read << gin.pairlist.type;
          printIO("PAIRLIST", "TYPE", read.str(), "0-2");
        }
      }
      if (gin.pathint.found) {
        if (gin.pathint.ntpi < 0 || gin.pathint.ntpi > 1) {
          stringstream read;
          read << gin.pathint.ntpi;
          printIO("PATHINT", "NTPI", read.str(), "0,1");
        }
      }
      if (gin.perscale.found) {
        if (gin.perscale.restype < 0 || gin.perscale.restype > 1) {
          stringstream read;
          read << gin.perscale.restype;
          printIO("PERSCALE", "RESTYPE", read.str(), "0,1");
        }
        if (gin.perscale.kdih < 0.0) {
          stringstream read;
          read << gin.perscale.kdih;
          printIO("PERSCALE", "KDIH", read.str(), ">=0.0");
        }
        if (gin.perscale.kj < 0.0) {
          stringstream read;
          read << gin.perscale.kj;
          printIO("PERSCALE", "KJ", read.str(), ">=0.0");
        }
        if (gin.perscale.t <= 0.0) {
          stringstream read;
          read << gin.perscale.t;
          printIO("PERSCALE", "T", read.str(), ">0.0");
        }
        if (gin.perscale.diff < 0.0) {
          stringstream read;
          read << gin.perscale.diff;
          printIO("PERSCALE", "DIFF", read.str(), ">=0.0");
        }
        if (gin.perscale.ratio <= 0.0) {
          stringstream read;
          read << gin.perscale.ratio;
          printIO("PERSCALE", "RATIO", read.str(), ">0.0");
        }
        if (gin.perscale.read < 0 || gin.perscale.read > 1) {
          stringstream read;
          read << gin.perscale.read;
          printIO("PERSCALE", "READ", read.str(), "0,1");
        }
      }
      if (gin.perturbation.found) {
        if (gin.perturbation.ntg < 0 || gin.perturbation.ntg > 1) {
          stringstream read;
          read << gin.perturbation.ntg;
          printIO("PERTURBATION", "NTG", read.str(), "0,1");
        }
        if (gin.perturbation.nrdgl < 0 || gin.perturbation.nrdgl > 1) {
          stringstream read;
          read << gin.perturbation.nrdgl;
          printIO("PERTURBATION", "NRDGL", read.str(), "0,1");
        }
        if (gin.perturbation.rlam < 0.0 || gin.perturbation.rlam > 1.0) {
          stringstream read;
          read << gin.perturbation.rlam;
          printIO("PERTURBATION", "RLAM", read.str(), "0.0..1.0");
        }
        if (gin.perturbation.dlamt < 0.0) {
          stringstream read;
          read << gin.perturbation.dlamt;
          printIO("PERTURBATION", "DLAMT", read.str(), ">=0.0");
        }
        if (gin.perturbation.alphlj < 0.0) {
          stringstream read;
          read << gin.perturbation.alphlj;
          printIO("PERTURBATION", "ALPHLJ", read.str(), ">=0.0");
        }
        if (gin.perturbation.alphc < 0.0) {
          stringstream read;
          read << gin.perturbation.alphc;
          printIO("PERTURBATION", "ALPHC", read.str(), ">=0.0");
        }
        if (gin.perturbation.nlam <= 0) {
          stringstream read;
          read << gin.perturbation.nlam;
          printIO("PERTURBATION", "NLAM", read.str(), ">0");
        }
        if (gin.perturbation.nscale < 0 || gin.perturbation.nscale > 2) {
          stringstream read;
          read << gin.perturbation.nscale;
          printIO("PERTURBATION", "NSCALE", read.str(), "0..2");
        }
      }
      if (gin.polarise.found) {
        if (gin.polarise.cos < 0 || gin.polarise.cos > 2) {
          stringstream read;
          read << gin.polarise.cos;
          printIO("POLARISE", "COS", read.str(), "0,2");
        }
        if (gin.polarise.efield < 0 || gin.polarise.efield > 1) {
          stringstream read;
          read << gin.polarise.efield;
          printIO("POLARISE", "EFIELD", read.str(), "0,1");
        }
        if (gin.polarise.minfield <= 0.0) {
          stringstream read;
          read << gin.polarise.minfield;
          printIO("POLARISE", "MINFIELD", read.str(), ">0.0");
        }
        if (gin.polarise.damp < 0 || gin.polarise.damp > 1) {
          stringstream read;
          read << gin.polarise.damp;
          printIO("POLARISE", "DAMP", read.str(), "0,1");
        }
        if (gin.polarise.write < 0) {
          stringstream read;
          read << gin.polarise.write;
          printIO("POLARISE", "WRITE", read.str(), ">=0");
        }
      }
      if (gin.positionres.found) {
        if (gin.positionres.ntpor < 0 || gin.positionres.ntpor > 3) {
          stringstream read;
          read << gin.positionres.ntpor;
          printIO("POSITIONRES", "NTPOR", read.str(), "0..3");
        }
        if (gin.positionres.ntporb < 0 || gin.positionres.ntporb > 1) {
          stringstream read;
          read << gin.positionres.ntporb;
          printIO("POSITIONRES", "NTPORB", read.str(), "0,1");
        }
        if (gin.positionres.ntpors < 0 || gin.positionres.ntpors > 1) {
          stringstream read;
          read << gin.positionres.ntpors;
          printIO("POSITIONRES", "NTPORS", read.str(), "0,1");
        }
        if (gin.positionres.cpor < 0.0) {
          stringstream read;
          read << gin.positionres.cpor;
          printIO("POSITIONRES", "CPOR", read.str(), ">=0.0");
        }
      }
      if (gin.precalclam.found) {
        if (gin.precalclam.nrlam < 0) {
          stringstream read;
          read << gin.precalclam.nrlam;
          printIO("PRECALCLAM", "NRLAM", read.str(), ">=0");
        }
        if (gin.precalclam.minlam < 0 || gin.precalclam.minlam > 1) {
          stringstream read;
          read << gin.precalclam.minlam;
          printIO("PRECALCLAM", "MINLAM", read.str(), "0..1");
        }
        if (gin.precalclam.maxlam < 0 || gin.precalclam.maxlam > 1) {
          stringstream read;
          read << gin.precalclam.maxlam;
          printIO("PRECALCLAM", "MAXLAM", read.str(), "0..1");
        }
        if (gin.precalclam.minlam >= gin.precalclam.maxlam ) {
          stringstream msg;
          msg << "MINLAM >= MAXLAM is not allowed.";
          printErrMsg("PRECALCLAM", "MINLAM", msg.str());
        }
      }
      if (gin.pressurescale.found) {
        if (gin.pressurescale.couple < 0 || gin.pressurescale.couple > 2) {
          stringstream read;
          read << gin.pressurescale.couple;
          printIO("PRESSURESCALE", "COUPLE", read.str(), "0..2");
        }
        if (gin.pressurescale.scale < 0 || gin.pressurescale.scale > 4) {
          stringstream read;
          read << gin.pressurescale.scale;
          printIO("PRESSURESCALE", "SCALE", read.str(), "0..4");
        }
        if (gin.pressurescale.comp <= 0.0) {
          stringstream read;
          read << gin.pressurescale.comp;
          printIO("PRESSURESCALE", "COMP", read.str(), ">0.0");
        }
        if (gin.pressurescale.taup <= 0.0) {
          stringstream read;
          read << gin.pressurescale.taup;
          printIO("PRESSURESCALE", "TAUP", read.str(), ">0.0");
        }
        if (gin.pressurescale.virial < 0 || gin.pressurescale.virial > 2) {
          stringstream read;
          read << gin.pressurescale.virial;
          printIO("PRESSURESCALE", "VIRIAL", read.str(), "0..2");
        }
        if (gin.pressurescale.x_semi < 0 || gin.pressurescale.x_semi > 2) {
          stringstream read;
          read << gin.pressurescale.x_semi;
          printIO("PRESSURESCALE", "SEMIX", read.str(), "0..2");
        }
        if (gin.pressurescale.y_semi < 0 || gin.pressurescale.y_semi > 2) {
          stringstream read;
          read << gin.pressurescale.y_semi;
          printIO("PRESSURESCALE", "SEMIY", read.str(), "0..2");
        }
        if (gin.pressurescale.z_semi < 0 || gin.pressurescale.z_semi > 2) {
          stringstream read;
          read << gin.pressurescale.z_semi;
          printIO("PRESSURESCALE", "SEMIZ", read.str(), "0..2");
        }
        // what to do to check the virial (PRES0[][])?
      }
      if (gin.printout.found) {
        if (gin.printout.ntpr < 0) {
          stringstream read;
          read << gin.printout.ntpr;
          printIO("PRINTOUT", "NTPR", read.str(), ">=0");
        }
        if (gin.printout.ntpp < 0 || gin.printout.ntpp > 1) {
          stringstream read;
          read << gin.printout.ntpp;
          printIO("PRINTOUT", "NTPP", read.str(), "0,1");
        }
      }
      if (gin.randomnumbers.found) {
        if (gin.randomnumbers.ntrng < 0 || gin.randomnumbers.ntrng > 1) {
          stringstream read;
          read << gin.randomnumbers.ntrng;
          printIO("RANDOMNUMBERS", "NTRNG", read.str(), "0,1");
        }
        if (gin.randomnumbers.ntgsl < -1) {
          stringstream read;
          read << gin.randomnumbers.ntgsl;
          printIO("RANDOMNUMBERS", "NTGSL", read.str(), ">=-1");
        }
      }
      if (gin.readtraj.found) {
        if (gin.readtraj.ntrd < 0 || gin.readtraj.ntrd > 1) {
          stringstream read;
          read << gin.readtraj.ntrd;
          printIO("READTRAJ", "NTRD", read.str(), "0,1");
        }
        if (gin.readtraj.ntstr < 1 ) {
          stringstream read;
          read << gin.readtraj.ntstr;
          printIO("READTRAJ", "NTSTR", read.str(), ">0");
        }
        if (gin.readtraj.ntrb < 0 || gin.readtraj.ntrb > 1) {
          stringstream read;
          read << gin.readtraj.ntrb;
          printIO("READTRAJ", "NTRB", read.str(), "0,1");
        }
        if (gin.readtraj.ntshk < 0 || gin.readtraj.ntshk > 2) {
          stringstream read;
          read << gin.readtraj.ntshk;
          printIO("READTRAJ", "NTSHK", read.str(), "0,1,2");
        }
      }
      if (gin.replica.found) {
        // NRET is checked when reading it from the input file
        for (unsigned int i = 0; i < gin.replica.ret.size(); i++) {
          if (gin.replica.ret[i] < 0.0) {
            stringstream read, blockName;
            read << gin.replica.ret[i];
            blockName << "RET[" << i + 1 << "]";
            printIO("REPLICA", blockName.str(), read.str(), ">=0.0");
          }
        }
        if (gin.replica.lrescale < 0 || gin.replica.lrescale > 1) {
          stringstream read;
          read << gin.replica.lrescale;
          printIO("REPLICA", "LRESCALE", read.str(), "0,1");
        }
        for (unsigned int i = 0; i < gin.replica.relam.size(); i++) {
          if (gin.replica.relam[i] < 0.0) {
            stringstream read, blockName;
            read << gin.replica.relam[i];
            blockName << "RELAM[" << i + 1 << "]";
            printIO("REPLICA", blockName.str(), read.str(), ">=0.0");
          }
          if (gin.replica.rets[i] < 0.0) {
            stringstream read, blockName;
            read << gin.replica.rets[i];
            blockName << "RETS[" << i + 1 << "]";
            printIO("REPLICA", blockName.str(), read.str(), ">=0.0");
          }
        }
        if (gin.replica.nretrial < 0) {
          stringstream read;
          read << gin.replica.nretrial;
          printIO("REPLICA", "NRETRIAL", read.str(), ">=0");
        }
        if (gin.replica.nrequil < 0) {
          stringstream read;
          read << gin.replica.nrequil;
          printIO("REPLICA", "NREQUIL", read.str(), ">=0");
        }
        if (gin.replica.cont < 0 || gin.replica.cont > 1) {
          stringstream read;
          read << gin.replica.cont;
          printIO("REPLICA", "CONT", read.str(), "0,1");
        }
      }
      if (gin.reeds.found) {
        if (gin.reeds.reeds < 0 || gin.reeds.reeds > 1) {
          stringstream read;
          read << gin.reeds.reeds;
          printIO("REPLICA_EDS", "REEDS", read.str(), "0,1");
        }
        if (gin.reeds.cont < 0 || gin.reeds.cont > 1) {
          stringstream read;
          read << gin.reeds.cont;
          printIO("REPLICA_EDS", "CONT", read.str(), "0,1");
        }
        if (gin.reeds.nretrial < 0) {
          stringstream read;
          read << gin.reeds.nretrial;
          printIO("REPLICA_EDS", "NRETRIAL", read.str(), ">=0");
        }
        if (gin.reeds.nrequil < 0) {
          stringstream read;
          read << gin.reeds.nrequil;
          printIO("REPLICA_EDS", "NREQUIL", read.str(), ">=0");
        }
      }
      if (gin.rottrans.found) {
        if (gin.rottrans.rtc < 0 || gin.rottrans.rtc > 1) {
          stringstream read;
          read << gin.rottrans.rtc;
          printIO("ROTTRANS", "RTC", read.str(), "0,1");
        }
        if (gin.rottrans.rtclast <= 0) {
          stringstream read;
          read << gin.rottrans.rtclast;
          printIO("ROTTRANS", "RTCLAST", read.str(), ">0");
        }
      }
      if (gin.sasa.found) {
        if(gin.sasa.ntsasa < 0 || gin.sasa.ntsasa >1) {
          stringstream read;
          read << gin.sasa.ntsasa;
          printIO("SASA", "NTSASA", read.str(), "0,1");
        }
        if(gin.sasa.ntvol < 0 || gin.sasa.ntvol >1) {
          stringstream read;
          read << gin.sasa.ntvol;
          printIO("SASA", "NTVOL", read.str(), "0,1");
        }
        if(gin.sasa.p12 <= 0 || gin.sasa.p12 >=1) {
          stringstream read;
          read << gin.sasa.p12;
          printIO("SASA", "P_12", read.str(), ">0 and <1");
        }
        if(gin.sasa.p13 <= 0 || gin.sasa.p13 >=1) {
          stringstream read;
          read << gin.sasa.p13;
          printIO("SASA", "P_13", read.str(), ">0 and <1");
        }
        if(gin.sasa.p1x <= 0 || gin.sasa.p1x >=1) {
          stringstream read;
          read << gin.sasa.p1x;
          printIO("SASA", "P_1X", read.str(), ">0 and <1");
        }
        if(gin.sasa.rsolv < 0) {
          stringstream read;
          read << gin.sasa.rsolv;
          printIO("SASA", "RSOLV", read.str(), ">0");
        }
        if(gin.sasa.as1 < 0) {
          stringstream read;
          read << gin.sasa.as1;
          printIO("SASA", "AS1", read.str(), ">0");
        }
        if(gin.sasa.as2 < 0) {
          stringstream read;
          read << gin.sasa.as2;
          printIO("SASA", "AS2", read.str(), ">0");
        }
      }
      if (gin.step.found) {
        if (gin.step.nstlim < 0) {
          stringstream read;
          read << gin.step.nstlim;
          printIO("STEP", "NSTLIM", read.str(), ">=0");
        }
        if (gin.step.t < 0.0 && gin.step.t != -1) {
          stringstream read;
          read << gin.step.t;
          printIO("STEP", "T", read.str(), ">=0.0 or -1");
        }
        if (gin.step.dt <= 0.0) {
          stringstream read;
          read << gin.step.dt;
          printIO("STEP", "DT", read.str(), ">0.0");
        }
      }
      if (gin.stochdyn.found) {
        if (gin.stochdyn.ntsd < 0 || gin.stochdyn.ntsd > 1) {
          stringstream read;
          read << gin.stochdyn.ntsd;
          printIO("STOCHDYN", "NTSD", read.str(), "0,1");
        }
        if (gin.stochdyn.ntfr < 0 || gin.stochdyn.ntfr > 3) {
          stringstream read;
          read << gin.stochdyn.ntfr;
          printIO("STOCHDYN", "NTFR", read.str(), "0..3");
        }
        if (gin.stochdyn.nsfr <= 0) {
          stringstream read;
          read << gin.stochdyn.nsfr;
          printIO("STOCHDYN", "NSFR", read.str(), ">0");
        }
        if (gin.stochdyn.nbref <= 0) {
          stringstream read;
          read << gin.stochdyn.nbref;
          printIO("STOCHDYN", "NBREF", read.str(), ">0");
        }
        if (gin.stochdyn.rcutf < 0.0) {
          stringstream read;
          read << gin.stochdyn.rcutf;
          printIO("STOCHDYN", "RCUTF", read.str(), ">=0.0");
        }
        if (gin.stochdyn.cfric < 0.0) {
          stringstream read;
          read << gin.stochdyn.cfric;
          printIO("STOCHDYN", "CFRIC", read.str(), ">=0.0");
        }
        if (gin.stochdyn.tempsd < 0.0) {
          stringstream read;
          read << gin.stochdyn.tempsd;
          printIO("STOCHDYN", "TEMPSD", read.str(), ">=0.0");
        }
      }
      if (gin.system.found) {
        if ((gin.system.npm < 0 || gin.system.npm > 1) && gromosXX) {
          stringstream read;
          read << gin.system.npm;
          printIO("SYSTEM", "NPM", read.str(), "0,1");
        } else if (gin.system.npm < 0) {
          stringstream read;
          read << gin.system.npm;
          printIO("SYSTEM", "NPM", read.str(), ">=0");
        }
        if (gin.system.nsm < 0) {
          stringstream read;
          read << gin.system.nsm;
          printIO("SYSTEM", "NSM", read.str(), ">=0");
        }
      }
      if (gin.thermostat.found) {
        if (gin.thermostat.ntt < 0 || gin.thermostat.ntt > 1) {
          stringstream read;
          read << gin.thermostat.ntt;
          printIO("THERMOSTAT", "NTT", read.str(), "0,1");
        }
        if (gin.thermostat.ntbth < 0) {
          stringstream read;
          read << gin.thermostat.ntbth;
          printIO("THERMOSTAT", "NTBTH", read.str(), ">=0");
        }
        if (gin.thermostat.ntset < 0) {
          stringstream read;
          read << gin.thermostat.ntset;
          printIO("THERMOSTAT", "NTSET", read.str(), ">=0");
        }
        for (int i = 0; i < gin.thermostat.ntbth; i++) {
          if (gin.thermostat.baths[i].index != i + 1) {
            stringstream read, blockName, allowed;
            read << gin.thermostat.baths[i].index;
            blockName << "I[" << i + 1 << "]";
            allowed << i + 1;
            printIO("THERMOSTAT", blockName.str(), read.str(), allowed.str());
          }
          if (gin.thermostat.baths[i].ntbtyp < 0 || gin.thermostat.baths[i].ntbtyp > 3) {
            stringstream read, blockName;
            read << gin.thermostat.baths[i].ntbtyp;
            blockName << "NTBTYP[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), "0..3");
          }
          if (gin.thermostat.baths[i].tembth < 0.0) {
            stringstream read, blockName;
            read << gin.thermostat.baths[i].tembth;
            blockName << "TEMBTH[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), ">=0.0");
          }
          if (gin.thermostat.baths[i].ntbvar < 0) {
            stringstream read, blockName;
            read << gin.thermostat.baths[i].ntbvar;
            blockName << "NTBVAR[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), ">=0");
          }
          for (int j = 0; j < gin.thermostat.baths[i].ntbvar; j++) {
            if (gin.thermostat.baths[i].taubth[j] < 0.0) {
              stringstream read, blockName;
              read << gin.thermostat.baths[i].taubth[j];
              blockName << "TAUBTH[" << i + 1 << "]";// << "[" << j << "]";
              printIO("THERMOSTAT", blockName.str(), read.str(), ">=0.0");
            }
          }
        }
        for (int i = 0; i < gin.thermostat.ntset; i++) {
          if (gin.thermostat.dofgroups[i].ntscpl < 0) {
            stringstream read, blockName;
            read << gin.thermostat.dofgroups[i].ntscpl;
            blockName << "NTSCPL[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), ">=0");
          }
          if (gin.thermostat.dofgroups[i].ntstyp < 0 || gin.thermostat.dofgroups[i].ntstyp > 2) {
            stringstream read, blockName;
            read << gin.thermostat.dofgroups[i].ntstyp;
            blockName << "NTSTYP[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), "0..2");
          }
          if (gin.thermostat.dofgroups[i].ntscns < 0 || gin.thermostat.dofgroups[i].ntscns > 1) {
            stringstream read, blockName;
            read << gin.thermostat.dofgroups[i].ntscns;
            blockName << "NTSCNS[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), "0..1");
          }
          if (gin.thermostat.dofgroups[i].ntsgt < -2) {
            stringstream read, blockName;
            read << gin.thermostat.dofgroups[i].ntsgt;
            blockName << "NTSGT[" << i + 1 << "]";
            printIO("THERMOSTAT", blockName.str(), read.str(), ">=-2");
          }
        }
      }
      if (gin.umbrella.found) {
        if (gin.umbrella.ntus < 0 || gin.umbrella.ntus > 1) {
          stringstream read;
          read << gin.umbrella.ntus;
          printIO("UMBRELLA", "NTUS", read.str(), "0,1");
        }
        if (gin.umbrella.uscst1 < 0.0) {
          stringstream read;
          read << gin.umbrella.uscst1;
          printIO("UMBRELLA", "USCST1", read.str(), ">=0.0");
        }
        if (gin.umbrella.uscst2 < 0.0) {
          stringstream read;
          read << gin.umbrella.uscst2;
          printIO("UMBRELLA", "USCST2", read.str(), ">=0.0");
        }
        if (gin.umbrella.usref1 < 0.0) {
          stringstream read;
          read << gin.umbrella.usref1;
          printIO("UMBRELLA", "USREF1", read.str(), ">=0.0");
        }
        if (gin.umbrella.usref2 < 0.0) {
          stringstream read;
          read << gin.umbrella.usref2;
          printIO("UMBRELLA", "USREF2", read.str(), ">=0.0");
        }
      }
      if (gin.virial.found) {
        if (gin.virial.ntv < 0 || gin.virial.ntv > 1) {
          stringstream read;
          read << gin.virial.ntv;
          printIO("VIRIAL", "NTV", read.str(), "0,1");
        }
        if (gin.virial.ntvg < 0 || gin.virial.ntvg > 3) {
          stringstream read;
          read << gin.virial.ntvg;
          printIO("VIRIAL", "NTVG", read.str(), "0,3");
        }
      }
      if (gin.virtualatom.found) {
        if (gin.virtualatom.virt < 0 || gin.virtualatom.virt > 1) {
          stringstream read;
          read << gin.virtualatom.virt;
          printIO("VIRTUALATOM", "VIRT", read.str(), "0,1");
        }
        if (gin.virtualatom.numvirt < 0) {
          stringstream read;
          read << gin.virtualatom.numvirt;
          printIO("VIRTUALATOM", "NUMVIRT", read.str(), ">=0");
        }
        if (gin.virtualatom.lastvirt < 0 || 
            gin.virtualatom.lastvirt < gin.virtualatom.numvirt) {
          stringstream read;
          read << gin.virtualatom.lastvirt;
          printIO("VIRTUALATOM", "LASTVIRT", read.str(), "0..NUMVIRT");
        } 
      }
      if (gin.writetraj.found) {
        // no checks needed for NTWX
        if (gin.writetraj.ntwse < 0) {
          stringstream read;
          read << gin.writetraj.ntwse;
          printIO("WRITETRAJ", "NTWSE", read.str(), ">=0");
        }
        // no checks needed for NTWV nad NTWF
        if (gin.writetraj.ntwe < 0) {
          stringstream read;
          read << gin.writetraj.ntwe;
          printIO("WRITETRAJ", "NTWE", read.str(), ">=0");
        }
        if (gin.writetraj.ntwg < 0) {
          stringstream read;
          read << gin.writetraj.ntwg;
          printIO("WRITETRAJ", "NTWG", read.str(), ">=0");
        }
        if (gin.writetraj.ntwb < 0) {
          stringstream read;
          read << gin.writetraj.ntwb;
          printIO("WRITETRAJ", "NTWB", read.str(), ">=0");
        }
      }

      // here some more complicated checks follow
      //
      // compare the LAST atom number from MULTIBATH block
      // with the total number of atoms
      if (gin.multibath.found) {
        int mxlast = 0;
        for (unsigned int i = 0; i < gin.multibath.last.size(); i++) {
          if (gin.multibath.last[i] > mxlast)
            mxlast = gin.multibath.last[i];
        }
        if (mxlast != numTotalAtoms+numVirtualAtoms) {
          std::stringstream ss;
          ss << "Highest occuring LAST atom in MULTIBATH ("
                  << mxlast << ") should be equal to the total\n"
                  "number of atoms (" << numTotalAtoms << ")";
          if(gin.virtualatom.found && gin.virtualatom.virt)
	    ss << " + the number of virtual atoms (" << numVirtualAtoms << ")";
          printWarning(ss.str());
        }
      }
      // check the stepsize with or without a CONSTRAINT or GEOMCONSTRAINTS block
      if (gin.system.found) {
        if (gin.constraint.found) {
          if (gin.system.npm == 0 && gin.constraint.ntc > 1)
            printError("No solute molecules (NPM=0 in SYSTEM block), solvent only simulation does not work with SHAKE for solute (NTC>1 in CONSTRAINT block)");

          if ((gin.system.npm != 0 && gin.constraint.ntc == 1 && gin.step.dt > 0.0005) ||
                  (gin.constraint.ntc == 2 && gin.step.dt > 0.001) ||
                  (gin.constraint.ntc == 3 && gin.step.dt > 0.002) ||
                  (gin.constraint.ntc == 4 && gin.step.dt > 0.0005)) {
            ostringstream os;
            string comment;
            double suggest = 0.0005;
            if (gin.constraint.ntc == 1) {
              comment = "no constraints on solute";
              suggest = 0.0005;
            } else if (gin.constraint.ntc == 2) {
              comment = "constraints on bonds involving H";
              suggest = 0.001;
            } else if (gin.constraint.ntc == 3) {
              comment = "constraints on all bonds";
              suggest = 0.002;
            } else if (gin.constraint.ntc == 4) {
              comment = "constraints on some bonds, and no constraints on other";
              suggest = 0.0005;
            }
            os << "DT in STEP block is set to " << gin.step.dt << ", which is "
                    << "considered to be too large \n";
            os << "if NTC = " << gin.constraint.ntc
                    << " in CONSTRAINT block.\n";
            os << "For NTC = " << gin.constraint.ntc << " (" << comment
                    << ") rather "
                    << "use DT = " << suggest << ".\n";
            printWarning(os.str());
          }
        }
        if (gin.geomconstraints.found) {
          if (gin.system.npm == 0 && gin.geomconstraints.ntcph == 1)
            printError("No solute molecules (NPM=0 in SYSTEM block), what do you want to constrain (NTCPH=1 in GEOMCONSTRAINTS block)");
          if (gin.system.npm == 0 && gin.geomconstraints.ntcpn == 1)
            printError("No solute molecules (NPM=0 in SYSTEM block), what do you want to constrain (NTCPN=1 in GEOMCONSTRAINTS block)");

          if ((gin.geomconstraints.ntcs == 0 && gin.step.dt > 0.0005)
                  ||
                  (gin.geomconstraints.ntcph == 0 && gin.step.dt > 0.0005)
                  ||
                  (gin.geomconstraints.ntcpn == 0 && gin.step.dt > 0.001)
                  || gin.step.dt > 0.002) {
            ostringstream os;
            double suggest = 0.002;
            if (gin.geomconstraints.ntcpn == 0) {
              suggest = 0.001;
            }
            if (gin.geomconstraints.ntcph == 0) {
              suggest = 0.0005;
            }
            if (gin.geomconstraints.ntcs == 0) {
              suggest = 0.0005;
            }
            os << "DT in STEP block is set to " << gin.step.dt << ", which is "
                    << "considered to be too large \n";
            os << "for NTCPH,NTCPN,NTCS = " << gin.geomconstraints.ntcph << ","
                    << gin.geomconstraints.ntcpn << "," << gin.geomconstraints.ntcs
                    << " in GEOMCONSTRAINTS block.\n";
            os << "Suggestion: use DT = " << suggest;
            printWarning(os.str());
          }
        }
      }
      // are FORCEs calculated for bonds that are shaken?
      if (gin.force.found) {
        if (gin.constraint.found) {
          if (gin.constraint.found && gin.constraint.ntc == 2 && gin.force.ntf[0] == 1)
            printWarning("NTF[1]=1 in FORCE block, but bond lengths are constraint");
          //if (gin.constraint.found && gin.constraint.ntc == 2 && gin.force.ntf[1] == 0)
          //  printWarning("NTF[2]=0 in FORCE block, and bond lengths are not constraint");
          if (gin.constraint.found && gin.constraint.ntc == 3 && gin.force.ntf[0] == 1)
            printWarning("NTF[1]=1 in FORCE block, but bond lengths are constraint");
          //if (gin.constraint.found && gin.constraint.ntc == 3 && gin.force.ntf[1] == 1)
          //  printWarning("NTF[2]=1 in FORCE block, but bond lengths are constraint");
          if (gin.constraint.found && gin.constraint.ntc < 2 && gin.force.ntf[0] == 0 && gin.system.npm != 0)
            printWarning("NTF[1]=0 in FORCE block, and bond lengths are not constraint");
          //if (gin.constraint.found && gin.constraint.ntc < 2 && gin.force.ntf[1] == 0)
          //  printWarning("NTF[2]=0 in FORCE block, and bond lengths are not constraint");
          if (gin.constraint.found && gin.constraint.ntc == 4 && gin.force.ntf[0] == 0)
            printWarning("NTF[1]=0 in FORCE block, and bond lengths may not be constraint");
          //if (gin.constraint.found && gin.constraint.ntc == 4 && gin.force.ntf[1] == 0)
          //  printWarning("NTF[2]=0 in FORCE block, and bond lengths may not be constraint");
          if (gin.constraint.found && gin.constraint.ntc == 4 && gin.force.ntf[0] == 1)
            printWarning("NTF[1]=0 in FORCE block, but bond lengths may be constraint");
          //if (gin.constraint.found && gin.constraint.ntc == 4 && gin.force.ntf[1] == 1)
          //  printWarning("NTF[2]=0 in FORCE block, but bond lengths may be constraint");
        }
        if (gin.geomconstraints.found) {
          if (gin.geomconstraints.ntcph == 1 && gin.force.ntf[0] == 1)
            printWarning("NTF[1]=1 in FORCE block, but bond lengths are constraint");
          if (gin.geomconstraints.ntcph == 0 && gin.force.ntf[0] == 0)
            printWarning("NTF[1]=0 in FORCE block, and bond lengths are not constraint");
          //if (gin.geomconstraints.ntcpn == 1 && gin.force.ntf[1] == 1)
          //  printWarning("NTF[2]=1 in FORCE block, but bond lengths are constraint");
          //if (gin.geomconstraints.ntcpn == 0 && gin.force.ntf[1] == 0)
          //  printWarning("NTF[2]=0 in FORCE block, and bond lengths are not constraint");
        }
      }

      // is the box big enough for the RCUTL-value?
      if (gin.boundcond.ntb != 0) { // no vacuum simulation
        // find the minimum box dimension
        Vec K = sys.box().K();
        Vec L = sys.box().L();
        Vec M = sys.box().M();
        if (gin.multicell.ntm) {
          K *= gin.multicell.ncella;
          L *= gin.multicell.ncellb;
          M *= gin.multicell.ncellc;
        }
        double h_K = (K - (K.dot(L.normalize()) + K.dot(M.normalize()))).abs();
        double h_L = (L - (L.dot(K.normalize()) + L.dot(M.normalize()))).abs();
        double h_M = (M - (M.dot(K.normalize()) + M.dot(L.normalize()))).abs();
        // take the minimum
        double minDim = h_K;
        if (h_L < minDim)
          minDim = h_L;
        if (h_M < minDim)
          minDim = h_M;

        // check the longrange cut off for gromosXX and promd
        if ((gromosXX && minDim / 2 < gin.pairlist.rcutl) || (!gromosXX && minDim / 2 < gin.neighbourlist.rcuti)) {
          stringstream msg;
          string var = "RCUTL";
          if (!gromosXX) {
            var = "RCUTI";
          }
          msg.precision(9);
          msg << "RCUTL in the PAIRLIST block is larger than half of the smallest "
                  "box dimension\n";
          msg << "Suggestion: use " << var << " no bigger than " << minDim / 2.0;
          printWarning(msg.str());
        }
      }
      // the reaction field cut-off distance is not equal to the RCUTL/RCUTI vailue
      if ((gromosXX && gin.nonbonded.rcrf != gin.pairlist.rcutl) ||
              ((!gromosXX && gin.nonbonded.rcrf != gin.neighbourlist.rcuti))) {
        string var = "RCUTL", blockName = "PAIRLIST";
        if (!gromosXX) {
          var = "RCUTI";
          blockName = "NEIGHBOURLIST";
        }
        stringstream msg;
        msg << var << " in the " << blockName << " block is not equal to the\n";
        msg << "reaction-field radius RCRF in the NONBONDED block";
        printWarning(msg.str());
      }
      // perturbation topology was given but no perturbation is requested from the input file
      if (l_pttopo > 0) {
        if (gin.perturbation.found == 0 && gin.eds.found == 0 && gin.aeds.found == 0) {
          stringstream msg;
          msg << "A perturbation topology was given but there is no PERTURBATION block"
                  " in the input file";
          printWarning(msg.str());
        } else if (gin.perturbation.ntg == 0 && gin.eds.eds == 0 && gin.aeds.aeds == 0) {
          stringstream msg;
          msg << "A perturbation topology was given but NTG = 0 in the PERTURBATION block";
          printWarning(msg.str());
        }
      }
      // gamd input file was given but no gamd is requested from the input file
      if (l_gamd > 0) {
        if (gin.gamd.found == 0) {
          stringstream msg;
          msg << "A gamd input file was given but there is no GAMD block"
                  " in the input file";
          printWarning(msg.str());
      }
      }
      // lambda values in a perturbation run get bigger than 1
      if (gin.perturbation.found && gin.perturbation.ntg == 1) {
        if ((gin.perturbation.rlam + gin.step.nstlim * gin.step.dt * gin.perturbation.dlamt) > 1.0) {
          stringstream msg;
          msg << "The combination of RLAM,DLAMT,NSTLIM = "
                  << gin.perturbation.rlam << "," << gin.perturbation.dlamt << ","
                  << gin.step.nstlim << " in the PERTURBATION\n"
                  "and STEP block leads to lambda > 1.0";
          printWarning(msg.str());
        }
      }
      // we want to precalculate free energy derivs but perturbation is off
      if (gin.precalclam.found && gin.precalclam.nrlam > 0) {
        if (gin.perturbation.found == 0 || gin.perturbation.ntg == 0) {
          stringstream msg;
          msg << "PRECALCLAM cannot be on without perturbation";
          printError(msg.str());
        }
      }
      // we want to precalculate free energy derivs but no energy trajectories are written
      if (gin.precalclam.found && gin.precalclam.nrlam > 0) {
        if (gin.writetraj.ntwe == 0 ||  gin.writetraj.ntwg == 0) {
          stringstream msg;
          msg << "NTWE and NTWG can not be 0 in WRITETRAJ block "
              << "when we want to precalculate free energy derivatives "
              << "(PRECALCLAM block)";
          printError(msg.str());
        }
      }
      // is the indicated and the read box in the same box format?
      if (sys.hasBox && gin.boundcond.ntb != sys.box().ntb()) {
        string iboxShape = "truncated-octahedral";
        string cboxShape = "truncated-octahedral";
        if (gin.boundcond.ntb == 0) {
          iboxShape = "vacuum";
        } else if (gin.boundcond.ntb == 1) {
          iboxShape = "rectangular";
        } else if (gin.boundcond.ntb == 2) {
          iboxShape = "triclinic";
        }
        if (sys.box().ntb() == 0) {
          cboxShape = "vacuum";
        } else if (sys.box().ntb() == 1) {
          cboxShape = "rectangular";
        } else if (sys.box().ntb() == 2) {
          cboxShape = "triclinic";
        }
        stringstream msg;
        msg << "NTB = " << gin.boundcond.ntb << " (" << iboxShape << ") does not match the box shape\n"
                "in " << s_coord << " (" << cboxShape << ")";
        printWarning(msg.str());
      }

      // Here goes the ERRORS
      if (gin.replica.cont != 1 && file_exists(s_coord)) { // not done if replica exchange continuation run and if the coord file does not exist
        if (iter == joblist.begin() && gin.initialise.ntivel == 0 && sys.hasVel == false && !(gin.energymin.found && gin.energymin.ntem != 0)) {
          stringstream msg;
          msg << "NTIVEL = 0 in INITIALISE block but no VELOCITY block in the initial coordinate file";
          printError(msg.str());
        }
        if (gin.initialise.ntishi == 0 && gin.boundcond.ntb == 0) {
          stringstream msg;
          msg << "NTISHI = 0 (INITIALISE block) is not allowed for a vacuum simulation\n";
          printError(msg.str());
        }
        if (iter == joblist.begin() &&  gin.initialise.ntishi != 1 && !sys.hasLatticeshifts && gin.boundcond.ntb != 0) {
          stringstream msg;
          msg << "There is no LATTICESHIFTS block in " << s_coord << " but NTISHI = 0 (INITIALISE block)\n"
                "indicates to read the lattice-shift vectors from there";
          printError(msg.str());
        }
        if (iter == joblist.begin() && gin.boundcond.ntb != 0 && !sys.hasBox) { //only if the file exists and it has no box, then complain
          stringstream msg;
          msg << "NTB != 0 in BOUNDARY block needs a box in " << s_coord;
          printError(msg.str());
        }
      }
      // number of atom in topology and force block
      if (gin.force.nre.size() && gin.force.nre[gin.force.nre.size() - 1] != numTotalAtoms + numVirtualAtoms) {
        stringstream msg;
        msg << "NRE[" << gin.force.nre.size() << "] = " << gin.force.nre[gin.force.nre.size() - 1]
                << " in FORCE block is not equal to the total number\n"
                << "of atoms calculated from the topology and SYSTEM block ("
                << numTotalAtoms << ")";
        if(numVirtualAtoms) 
          msg << " and the VIRTUALATOM block (" << numVirtualAtoms << ")"; 
        printError(msg.str());
      }
      // number of atoms from topology vs number of atoms from coordinate file
      if (gin.replica.cont != 1 && file_exists(s_coord)) { // not done for a replica exchange continuation run and if coord file does not exist
        int numAtoms = 0;
        if (sys.numSolvents() > 0) {
          numAtoms += sys.numSolvents() * sys.sol(0).numAtoms();
        }
        for (int i = 0; i < sys.numMolecules(); i++) {
          numAtoms += sys.mol(i).numAtoms();
        }
        if (numAtoms != numTotalAtoms) {
          stringstream msg;
          msg << "The number of positions in " << s_coord << " (" << numAtoms << ") does not fit the number of\n"
                  "atoms calculated from the topology an SYSTEM block (" << numTotalAtoms << ")";
          printError(msg.str());
        }
      }
      // number of virtual atoms in input and topology
      if(gin.virtualatom.found && gin.virtualatom.virt){
        if(sys.vas().numVirtualAtoms() != gin.virtualatom.numvirt){
          stringstream msg;
          msg << "The number of virtual atoms in VIRTUALATOM block (" 
              << numVirtualAtoms 
              << ") does not match the number of virtual atoms in topology (" 
              << sys.vas().numVirtualAtoms() << ")";
          printError(msg.str());
        }
        if(gin.virtualatom.lastvirt != gin.system.npm * numSoluteAtoms + sys.vas().numVirtualAtoms()){
          stringstream msg;
          msg << "The index of the last virtual atom in VIRTUALATOM block (" 
              << gin.virtualatom.lastvirt 
              << ") does not match the value expected from the topology ("
              << gin.system.npm * numSoluteAtoms + sys.vas().numVirtualAtoms() << ")";
          printError(msg.str()); 
        }
      }
      // RCUTP must be <= RCUTL
      if ((gromosXX && gin.pairlist.rcutl < gin.pairlist.rcutp) ||
              (!gromosXX && gin.neighbourlist.rcuts > gin.neighbourlist.rcuti)) {
        string sr = "RCUTP", lr = "RCUTL", blockName = "PAIRLIST";
        if (!gromosXX) {
          sr = "RCUTS";
          lr = "RCUTI";
          blockName = "NEIGHBOURLIST";
        }
        stringstream msg;
        msg << sr << " > " << lr << " (" << blockName << " block) is not allowed";
        printError(msg.str());
      }
      // position restraining demanded but no position restraints specification file
      if (gin.positionres.found && gin.positionres.ntpor > 0 && l_posresspec == 0) {
        stringstream msg;
        msg << "NTPOR > 0 in POSITIONRES block indicates a position restraining but no\n"
                "posresspec file has been found";
        printError(msg.str());
      }
      // distancefield specified but no disres file
      if(gin.distancefield.found && gin.distancefield.ntdfr != 0 && l_disres == 0){
	stringstream msg;
	msg << "Distancefield restraint is turned on (NTDFR = " << gin.distancefield.ntdfr << " in DISTANCEFIELD block)\n"
	    << "but no distance restraint file has been found";
	printError(msg.str());
      }
      // colvarres specified but no colvarres file
      if(gin.colvarres.found && gin.colvarres.cvr != 0 && l_colvarres == 0){
	stringstream msg;
	msg << "Collective variable restraint is turned on (CVR = " << gin.colvarres.cvr << " in COLVARRES block)\n"
	    << "but no colvar restraint file has been found";
	printError(msg.str());
      }
      
      // no perturbation topology specified but perturbation is switched on
      if (gin.perturbation.found && gin.perturbation.ntg != 0 && l_pttopo == 0) {
        stringstream msg;
        msg << "Perturbation is turned on (NTG = " << gin.perturbation.ntg << " in PERTURBATION block)\n"
                "but no perturbation topology has been found - might be ok if you only have perturbed restraints";
        printWarning(msg.str());
      }
      // no gamd input file specified but perturbation is switched on
      if (gin.gamd.found && gin.gamd.gamd != 0 && l_gamd == 0) {
        stringstream msg;
        msg << "GAMD is turned on but no gamd input file has been found";
        printError(msg.str());
      }
      // electric block specified with more atoms than SYSTEM
      if (gin.electric.found && gin.electric.curgrp.size() > 0 && gin.electric.current > 0){
        if(gin.electric.curgrp[gin.electric.curgrp.size() -1] > numTotalAtoms){
          stringstream msg;
          msg << "CURGRP[" << gin.electric.curgrp.size()
                << " in ELECTRIC block is bigger than the total number\n"
                << "of atoms calculated from the topology and SYSTEM block ("
                << numTotalAtoms << ")";
        printError(msg.str());
        }
      }
      // Hamiltonian replica exchange but no perturbation
      if (gin.replica.relam.size() > 1 && (gin.perturbation.ntg != 1 || !gin.perturbation.found)) {
        stringstream msg;
        msg << "Hamiltonian replica exchange but no PERTURBATION block";
        printWarning(msg.str());
      }


      // Now, write the script
      if (numErrors + numWarnings > 0) {
        if (numErrors + numWarnings == 1) {
          cout << "\nTHERE WAS ";
        } else {
          cout << "\nTHERE WERE ";
        }
        if (numErrors > 1) {
          cout << numErrors << " ERRORS";
        } else if(numErrors == 1) {
          cout << "1 ERROR";
        }
        if (numErrors > 0 && numWarnings > 0) {
          cout << " AND ";
        }
        if (numWarnings > 1) {
          cout << numWarnings << " WARNINGS";
        } else if(numWarnings == 1) {
          cout << "1 WARNING";
        }
        cout << endl;
        // are the script forced to be written?
        if (numErrors > 0 && args.count("force") == -1) {
            cout << "No script has been written\n";
            cout << "\n--------------------------------------------------------------------------------\n\n";
            continue;
        }
      } else {
        cout << "OK (no warnings, no errors)\n";
      }


      //first check whether we should be in a different directory
      string subdir = simuldir;
      if (iter->second.dir != ".") {
        subdir += "/" + iter->second.dir;
      }
      mkdir(subdir.c_str(), 00755);
      if (chdir(subdir.c_str()) != 0) {
        throw gromos::Exception(argv[0], "could not change to the subdir directory");
      }
      if(numErrors == 0) {
        cout << "Writing script: " << filenames[FILETYPE["script"]].name(0) << endl;
      } else {
        cout << "Forcing script: " << filenames[FILETYPE["script"]].name(0) << endl;
      }

      ofstream fout(filenames[FILETYPE["script"]].name(0).c_str());
      fout.setf(ios::left, ios::adjustfield);
      fout << "#!/bin/sh" << endl;
      for (unsigned i = 0; i < directives.size(); ++i)
        fout << directives[i].name(0) << endl;
      fout << "\n# first we set some variables\n";
      fout << "NAME=`whoami`\n";
      fout << "PROGRAM=" << gromosbin << endl;
      fout << "SIMULDIR=" << simuldir << endl;
      fout << "\n# create temporary directory\n";
      fout << "WORKDIR=" << misc[0].name(0) << endl;
      fout << "mkdir -p ${WORKDIR}\n";
      fout << "cd       ${WORKDIR}\n";
      fout << "\n# set the input files\n";
      fout << "TOPO=${SIMULDIR}/" << s_topo << endl;
      fout << "IUNIT=${SIMULDIR}/";
      if (iter->second.dir != ".") fout << iter->second.dir << "/";
      fout << filenames[FILETYPE["input"]].name(0) << endl;

      if (iter != joblist.begin() || iter->second.dir != "." ||
              filenames[FILETYPE["input"]].name(0) != s_input) {
        // write the new input files
        cout << "     and input: " << filenames[FILETYPE["input"]].name(0) << endl;
        printInput(filenames[FILETYPE["input"]].name(0), gin);
      }

      cout << "\n--------------------------------------------------------------------------------\n\n";

      fout << "INPUTCRD=${SIMULDIR}/";
      if (iter == joblist.begin()) {
        fout << s_coord << endl;
      } else {
        if (iter->second.prev_id == -1) {
          //if(iter->second.dir!=".") fout << iter->second.dir << "/";
          fout << s_coord << endl;
        } else {
          jobinfo prevjob = joblist[iter->second.prev_id];
          if (prevjob.dir != ".") fout << prevjob.dir << "/";
          filenames[FILETYPE["coord"]].setInfo(systemname,
                  atof(prevjob.param["T"].c_str()),
                  atof(prevjob.param["DELTAT"].c_str()),
                  iter->second.prev_id,
                  queue);
          fout << filenames[FILETYPE["coord"]].name(0) << endl;
          filenames[FILETYPE["coord"]].setInfo(systemname,
                  atof(iter->second.param["T"].c_str()),
                  atof(iter->second.param["DELTAT"].c_str()),
                  iter->first,
                  queue);
        }
      }
      // re-analyzing?
      if (gin.readtraj.found && gin.readtraj.ntrd == 1) {
        fout << "ANATRX=${SIMULDIR}/";
        if (l_anatrx) fout  << s_anatrx << endl;
        else fout << filenames[FILETYPE["anatrj"]].name(0) << endl;
      }
      // EVTRL
      if (l_refpos) fout << "REFPOS=${SIMULDIR}/" << s_refpos << endl;
      if (l_posresspec) fout << "POSRESSPEC=${SIMULDIR}/"
              << s_posresspec << endl;
      if (l_xray) fout << "XRAY=${SIMULDIR}/" << s_xray << endl;
      if (l_disres) fout << "DISRES=${SIMULDIR}/" << s_disres << endl;
      if (l_dihres) fout << "DIHRES=${SIMULDIR}/" << s_dihres << endl;
      if (l_angres) fout << "ANGRES=${SIMULDIR}/" << s_angres << endl;
      if (l_colvarres) fout << "COLVARRES=${SIMULDIR}/" << s_colvarres << endl;
      if (l_jvalue) fout << "JVALUE=${SIMULDIR}/" << s_jvalue << endl;
      if (l_order) fout << "ORDER=${SIMULDIR}/" << s_order << endl;
      if (l_sym) fout << "SYM=${SIMULDIR}/" << s_sym << endl;
      if (l_ledih) fout << "LEDIH=${SIMULDIR}/" << s_ledih << endl;
      if (l_friction) fout << "FRICTION=${SIMULDIR}/" << s_friction << endl;
      if (l_leumb) fout << "LEUMB=${SIMULDIR}/" << s_leumb << endl;
      if (l_bsleus) fout << "BSLEUS=${SIMULDIR}/" << s_bsleus << endl;
      if (l_pttopo) fout << "PTTOPO=${SIMULDIR}/" << s_pttopo << endl;
      if (l_gamd) fout << "GAMD=${SIMULDIR}/" << s_gamd << endl;
      if (l_qmmm) fout << "QMMM=${SIMULDIR}/" << s_qmmm << endl;
      
      // any additional links?
      for (unsigned int k = 0; k < linkadditions.size(); k++)
        if (linkadditions[k] < 0)
          fout << linknames[k] << "=${SIMULDIR}/"
                << filenames[numFiletypes + k].name(0)
          << endl;

      fout << "\n#set the output files\n";
      fout << "OUNIT=" << filenames[FILETYPE["output"]].name(0) << endl;
      fout << "OUTPUTCRD=" << filenames[FILETYPE["coord"]].name(0) << endl;
      if (gin.writetraj.ntwx)
        fout << "OUTPUTTRX="
              << filenames[FILETYPE["outtrx"]].name(0)
        << endl;
      if (gin.writetraj.ntwv)
        fout << "OUTPUTTRV="
              << filenames[FILETYPE["outtrv"]].name(0)
        << endl;
      if (gin.writetraj.ntwf)
        fout << "OUTPUTTRF="
              << filenames[FILETYPE["outtrf"]].name(0)
        << endl;
      if (gin.writetraj.ntwe)
        fout << "OUTPUTTRE="
              << filenames[FILETYPE["outtre"]].name(0)
        << endl;
      if (gin.writetraj.ntwg)
        fout << "OUTPUTTRG="
              << filenames[FILETYPE["outtrg"]].name(0)
        << endl;
      if (gin.writetraj.ntwb)
        fout << "OUTPUTBAE="
              << filenames[FILETYPE["outbae"]].name(0)
        << endl;
      
      if (l_jin) {
        fout << "JIN=";
        if (iter == joblist.begin()) {
          fout << s_jin << endl;
        } else {
          if (iter->second.prev_id == -1) {
            //if(iter->second.dir!=".") fout << iter->second.dir << "/";
            fout << s_coord << endl;
          } else {
            jobinfo prevjob = joblist[iter->second.prev_id];
            if (prevjob.dir != ".") fout << prevjob.dir << "/";
            filenames[FILETYPE["jin"]].setInfo(systemname,
                    atof(prevjob.param["T"].c_str()),
                    atof(prevjob.param["DELTAT"].c_str()),
                    iter->second.prev_id,
                    queue);
            fout << filenames[FILETYPE["jin"]].name(0) << endl;
            filenames[FILETYPE["jin"]].setInfo(systemname,
                    atof(iter->second.param["T"].c_str()),
                    atof(iter->second.param["DELTAT"].c_str()),
                    iter->first,
                    queue);
          }
        }
        fout << "JOUT="
             << filenames[FILETYPE["jin"]].name(0) << endl;
        fout << "JTRJ="
             << filenames[FILETYPE["jtrj"]].name(0) << endl;
      }
      
      bool write_trs = gin.polarise.write || gin.jvalueres.write || gin.orderparamres.ntwop|| gin.xrayres.ntwxr ||
              gin.localelev.ntwle || gin.bsleus.write || gin.addecouple.write || gin.nemd.write|| gin.printout.ntpp == 1
              || gin.electric.dipole == 1 || gin.electric.current == 1 || gin.distanceres.ntwdir > 0 
              || gin.distancefield.ntwdf > 0 || gin.dihedralres.ntwdlr > 0 || gin.angleres.ntwalr > 0 || gin.colvarres.ntwcv > 0;
      if (write_trs) {
        fout << "OUTPUTTRS="
	     << filenames[FILETYPE["outtrs"]].name(0)
	     << endl;
      }

      if (gin.writetraj.ntwb &&
              (gin.perturbation.found && gin.perturbation.ntg > 0))
        fout << "OUTPUTBAG="
              << filenames[FILETYPE["outbag"]].name(0)
        << endl;
      if (gin.replica.found) {
        //Write repout and repdat if not specified in @files
        fout << "REPOUT="
             << filenames[FILETYPE["repout"]].name(0)
             << endl;
        fout << "REPDAT="
             << filenames[FILETYPE["repdat"]].name(0)
             << endl;
      }

      // any additional links?
      for (unsigned int k = 0; k < linkadditions.size(); k++)
        if (linkadditions[k] > 0)
          fout << linknames[k] << "="
                << filenames[numFiletypes + k].name(0) << endl;

      if (misc[2].name(0) != "") {
        fout << "\n# first command\n";
        fout << misc[2].name(0) << "\n";
      }


      fout << "\n\n";

      fout << "MDOK=1\n\n";
      fout << misc[3].name(0) << "${PROGRAM}";

      fout << " \\\n\t" << setw(12) << "@topo" << " ${TOPO}";
      fout << " \\\n\t" << setw(12) << "@conf" << " ${INPUTCRD}";
      fout << " \\\n\t" << setw(12) << "@input" << " ${IUNIT}";
      if (gin.readtraj.found && gin.readtraj.ntrd == 1) fout << " \\\n\t"
              << setw(12) << "@anatrj" << " ${ANATRX}";
      if (l_pttopo) fout << " \\\n\t"
              << setw(12) << "@pttopo" << " ${PTTOPO}";
      if (l_gamd) fout << " \\\n\t"
              << setw(12) << "@gamd" << " ${GAMD}";
      if (l_posresspec) fout << " \\\n\t"
              << setw(12) << "@posresspec" << " ${POSRESSPEC}";
      if (l_refpos) fout << " \\\n\t"
              << setw(12) << "@refpos" << " ${REFPOS}";
      if (l_xray) fout << " \\\n\t"
              << setw(12) << "@xray" << " ${XRAY}";
      if (l_disres) fout << " \\\n\t"
              << setw(12) << "@distrest" << " ${DISRES}";
      if (l_dihres) fout << " \\\n\t"
              << setw(12) << "@dihrest" << " ${DIHRES}";
      if (l_angres) fout << " \\\n\t"
              << setw(12) << "@angrest" << " ${ANGRES}";
      if (l_colvarres) fout << " \\\n\t"
              << setw(12) << "@colvarres" << " ${COLVARRES}";
      if (l_bsleus) fout << " \\\n\t"
              << setw(12) << "@bsleus" << " ${BSLEUS}";
      if (l_qmmm) fout << " \\\n\t"
              << setw(12) << "@qmmm" << " ${QMMM}";
      if (l_jin){
        fout << " \\\n\t"
             << setw(12) << "@jin" << " ${JIN}"
             << " \\\n\t"
             << setw(12) << "@jout" << " ${JOUT}"
             << " \\\n\t"
             << setw(12) << "@jtraj" << " ${JTRJ}";
      }
      if (l_friction) fout << " \\\n\t"
              << setw(12) << "@friction" << " ${FRICTION}";
      if (l_jvalue) fout << " \\\n\t"
              << setw(12) << "@jval" << " ${JVALUE}";
      if (l_order) fout << "\\\n\t"
              << setw(12) << "@order" << " ${ORDER}";
      if (l_sym) fout << " \\\n\t"
              << setw(12) << "@sym" << " ${SYM}";
      if (l_ledih) fout << " \\\n\t"
              << setw(12) << "@led" << " ${LEDIH}";
      if (l_leumb) fout << " \\\n\t"
              << setw(12) << "@lud" << " ${LEUMB}";

      fout << " \\\n\t" << setw(12) << "@fin" << " ${OUTPUTCRD}";
      if (gin.writetraj.ntwx) fout << " \\\n\t" << setw(12) << "@trc"
        << " ${OUTPUTTRX}";
      if (gin.writetraj.ntwv) fout << " \\\n\t" << setw(12) << "@trv"
        << " ${OUTPUTTRV}";
      if (gin.writetraj.ntwf) fout << " \\\n\t" << setw(12) << "@trf"
        << " ${OUTPUTTRF}";
      if (gin.writetraj.ntwe) fout << " \\\n\t" << setw(12) << "@tre"
        << " ${OUTPUTTRE}";
      if (gin.writetraj.ntwg) fout << " \\\n\t" << setw(12) << "@trg"
        << " ${OUTPUTTRG}";
      if (gin.writetraj.ntwb) fout << " \\\n\t" << setw(12) << "@bae"
        << " ${OUTPUTBAE}";
      if (write_trs)
        fout << " \\\n\t" << setw(12) << setw(12) << "@trs" << " ${OUTPUTTRS}";

      if (gin.writetraj.ntwb > 0 &&
              gin.perturbation.found && gin.perturbation.ntg > 0)
        fout << " \\\n\t" << setw(12) << "@bag"
        << " ${OUTPUTBAG}";

      // any additional links
      for (unsigned int k = 0; k < linkadditions.size(); k++)
        fout << " \\\n\t@" << setw(11) << linknames[k]
              << " ${" << linknames[k] << "}";

      if (gin.replica.found) {
        fout << "\\\n\t" << setw(12) << "@repout" << " ${REPOUT}";
        fout << "\\\n\t" << setw(12) << "@repdat" << " ${REPDAT}";
      }

      if (gromosXX) {
        if(putdevelop) {
          fout << " \\\n\t" << setw(12) << "@develop";
        }
        if (putdebug) {
          fout << " \\\n\t" << setw(12) << "@verb " << debug;
        }
        fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT}\n";
        fout << "grep \"finished successfully\" ${OUNIT} > /dev/null || MDOK=0";
      } else {
        fout << "\\\n\t" << setw(12) << ">" << " ${OUNIT}     || MDOK=0";
      }
      fout << "\n\n";


      fout << "uname -a >> ${OUNIT}\n";

      if (gin.writetraj.ntwx || gin.writetraj.ntwv || gin.writetraj.ntwf ||
              gin.writetraj.ntwe || gin.writetraj.ntwg || gin.printout.ntpp == 1)
        fout << "\n# compress some files\n";
      // replica exchange
      if (gin.replica.relam.size() > 1 || gin.replica.ret.size() > 1) {
        // adapt file names
        std::stringstream tmp;
        tmp << "_*";
        string outtrx_name = filenames[FILETYPE["outtrx"]].name(0);
        int pos = outtrx_name.find_last_of(".");
        outtrx_name = outtrx_name.insert(pos, tmp.str());
        string outtrv_name = filenames[FILETYPE["outtrv"]].name(0);
        pos = outtrv_name.find_last_of(".");
        outtrv_name = outtrv_name.insert(pos, tmp.str());
        string outtrf_name = filenames[FILETYPE["outtrf"]].name(0);
        pos = outtrf_name.find_last_of(".");
        outtrf_name = outtrf_name.insert(pos, tmp.str());
        string outtre_name = filenames[FILETYPE["outtre"]].name(0);
        pos = outtre_name.find_last_of(".");
        outtre_name = outtre_name.insert(pos, tmp.str());
        string outtrg_name = filenames[FILETYPE["outtrg"]].name(0);
        pos = outtrg_name.find_last_of(".");
        outtrg_name = outtrg_name.insert(pos, tmp.str());
        string outbae_name = filenames[FILETYPE["outbae"]].name(0);
        pos = outbae_name.find_last_of(".");
        outbae_name = outbae_name.insert(pos, tmp.str());
        string outbag_name = filenames[FILETYPE["outbag"]].name(0);
        pos = outbag_name.find_last_of(".");
        outbag_name = outbag_name.insert(pos, tmp.str());
        string outtrs_name = filenames[FILETYPE["outtrs"]].name(0);
        pos = outtrs_name.find_last_of(".");
        outtrs_name = outtrs_name.insert(pos, tmp.str());
        string repout_name = filenames[FILETYPE["repout"]].name(0);
        pos = repout_name.find_last_of(".");
        repout_name = repout_name.insert(pos, tmp.str());
        string coord_name = filenames[FILETYPE["coord"]].name(0);
        pos = coord_name.find_last_of(".");
        coord_name = coord_name.insert(pos, tmp.str());

        if (gin.writetraj.ntwx) fout << "gzip " << outtrx_name << "\n";
        if (gin.writetraj.ntwv) fout << "gzip " << outtrv_name << "\n";
        if (gin.writetraj.ntwf) fout << "gzip " << outtrf_name << "\n";
        if (gin.writetraj.ntwe) fout << "gzip " << outtre_name << "\n";
        if (gin.writetraj.ntwg) fout << "gzip " << outtrg_name << "\n";
        if (gin.writetraj.ntwb) fout << "gzip " << outbae_name << "\n";
        if (gin.writetraj.ntwb &&
                gin.perturbation.found && gin.perturbation.ntg > 0)
          fout << "gzip " << outbag_name << "\n";
        if (write_trs)
          fout << "gzip " << outtrs_name << "\n";

        fout << "\n# copy the files back\n";
        fout << "OK=1\n";
        fout << setw(25) << "cp ${OUNIT}" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        fout << setw(25) << "cp " << coord_name << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        if (gin.writetraj.ntwx) {
          fout << setw(25) << "cp " << outtrx_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwv) {
          fout << setw(25) << "cp " << outtrv_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwf) {
          fout << setw(25) << "cp " << outtrf_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwe) {
          fout << setw(25) << "cp " << outtre_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwg) {
          fout << setw(25) << "cp " << outtrg_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwb) {
          fout << setw(25) << "cp " << outbae_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";

          if (gin.perturbation.found && gin.perturbation.ntg > 0) {
            fout << setw(25) << "cp " << outbag_name << ".gz ${SIMULDIR}";
            if (iter->second.dir != ".") fout << "/" << iter->second.dir;
            fout << " || OK=0\n";
          }
        }
        if (write_trs) {
          fout << setw(25) << "cp " << outtrs_name << ".gz ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        fout << setw(25) << "cp " << repout_name << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        fout << setw(25) << "cp ${REPDAT}" << "${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        
      } else { // no replica exchange
        if (gin.writetraj.ntwx) fout << "gzip ${OUTPUTTRX}\n";
        if (gin.writetraj.ntwv) fout << "gzip ${OUTPUTTRV}\n";
        if (gin.writetraj.ntwf) fout << "gzip ${OUTPUTTRF}\n";
        if (gin.writetraj.ntwe) fout << "gzip ${OUTPUTTRE}\n";
        if (gin.writetraj.ntwg) fout << "gzip ${OUTPUTTRG}\n";
        if (gin.writetraj.ntwb) fout << "gzip ${OUTPUTBAE}\n";
        if (gin.writetraj.ntwb &&
                gin.perturbation.found && gin.perturbation.ntg > 0)
          fout << "gzip ${OUTPUTBAG}\n";
        if (write_trs)
          fout << "gzip ${OUTPUTTRS}\n";
        if (l_jin)
          fout << "gzip ${JTRJ}\n";

        fout << "\n# copy the files back\n";
        fout << "OK=1\n";
        fout << setw(25) << "cp ${OUNIT}" << " ${SIMULDIR}";
        if (iter->second.dir != ".") fout << "/" << iter->second.dir;
        fout << " || OK=0\n";
        if (!gin.readtraj.found || (gin.readtraj.found && gin.readtraj.ntrd == 0)) {
          fout << setw(25) << "cp ${OUTPUTCRD}" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwx) {
          fout << setw(25) << "cp ${OUTPUTTRX}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwv) {
          fout << setw(25) << "cp ${OUTPUTTRV}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwf) {
          fout << setw(25) << "cp ${OUTPUTTRF}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwe) {
          fout << setw(25) << "cp ${OUTPUTTRE}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwg) {
          fout << setw(25) << "cp ${OUTPUTTRG}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (gin.writetraj.ntwb) {
          fout << setw(25) << "cp ${OUTPUTBAE}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";

          if (gin.perturbation.found && gin.perturbation.ntg > 0) {
            fout << setw(25) << "cp ${OUTPUTBAG}.gz" << " ${SIMULDIR}";
            if (iter->second.dir != ".") fout << "/" << iter->second.dir;
            fout << " || OK=0\n";
          }
        }
        if (write_trs) {
          fout << setw(25) << "cp ${OUTPUTTRS}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
        if (l_jin) {
          fout << setw(25) << "cp ${JTRJ}.gz" << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
      }
      // any additional links
      for (unsigned int k = 0; k < linkadditions.size(); k++) {
        if (linkadditions[k] > 0) {
          string s("cp ${" + linknames[k] + "}");
          string sgz("cp ${" + linknames[k] + "}.gz");
          fout << setw(25) << s << " ${SIMULDIR}";
          if (iter->second.dir != ".") fout << "/" << iter->second.dir;
          fout << " || OK=0\n";
        }
      }

      fout << "\n# clean up after successful run\n";
      map<int, jobinfo>::iterator testiter = iter;
      testiter++;
      if (mail == "EVERY" || (mail == "LAST" && testiter == joblist.end())) {
        fout << "if `test ${OK} -eq 0`; then\n";
        fout << "  uname -a > mess;\n";
        fout << "  echo 'cp failed for " << systemname << ", run "
                << iter->first << "' >> mess;\n";
        string subject = "ERROR" + jobID;
        fout << "  Mail -s \"" + subject + "\" ${NAME} < mess;\n";
        fout << "  cd ${SIMULDIR};\n";
        fout << "else\n";
        fout << "  cd ${SIMULDIR};\n";
        fout << "  rm ${WORKDIR}/*;\n";
        fout << "  rmdir ${WORKDIR};\n";
        fout << "fi\n";
      } else {
        fout << "if `test ${OK} -eq 0`; then\n";
        fout << "  cd ${SIMULDIR};\n";
        fout << "else\n";
        fout << "  cd ${SIMULDIR};\n";
        fout << "  rm ${WORKDIR}/*;\n";
        fout << "  rmdir ${WORKDIR};\n";
        fout << "fi\n";
      }

      fout << "\n# stop if MD was not succesfull\n";
      fout << "if `test ${MDOK} -eq 0`; then\n";
      if (misc[4].name(0) != "") {
        fout << "  " << misc[4].name(0) << "\n";
      }
      fout << "  exit\n";
      fout << "fi\n";

      fout << "\n# perform last command (usually submit next job)\n";
      // which job do we have to submit (also check in the earlier ones?)
      map<int, jobinfo>::const_iterator it = joblist.begin();
      while (it != to) {
        if (it->second.prev_id == iter->first) {
          setParam(gin, it->second);
          misc[1].setInfo(systemname,
                  atof(iter->second.param["ENDTIME"].c_str()),
                  gin.step.dt * gin.step.nstlim, it->first, queue);
          if (it->first != iter->first) {

            fout << "cd ${SIMULDIR}";
            if (it->second.dir != ".") fout << "/" << it->second.dir;
            fout << "\n";
            fout << misc[1].name(0) << endl;
          }

        }
        ++it;
      }

      fout.close();
      chmod(filenames[FILETYPE["script"]].name(0).c_str(), 00755);

      writtenScripts++;

    }

    if (numTotErrors + numTotWarnings > 0) {
      if (numTotErrors + numTotWarnings == 1) {
        cout << "THERE WAS ";
      } else {
        cout << "THERE WERE ";
      }
      if (numTotErrors > 1) {
        cout << numTotErrors << " ERRORS";
      } else if (numTotErrors == 1) {
        cout << "1 ERROR";
      }
      if (numTotErrors > 0 && numTotWarnings > 0) {
        cout << " AND ";
      }
      if (numTotWarnings > 1) {
        cout << numTotWarnings << " WARNINGS";
      } else if (numTotWarnings == 1) {
        cout << "1 WARNING";
      }
      cout << " IN TOTAL\n";
      string scriptFile = "script file";
      string toBe = "has";
      if(writtenScripts > 1) {
        toBe = "have";
      }
      if(joblist.size() > 1) {
        scriptFile = "script files";
      }
     
      cout << writtenScripts << " of " << joblist.size() << " " << scriptFile << " "
              << toBe << " been ";
      if (args.count("force") >= 0) {
        cout << "(forced to be) ";
      }
      cout << "written\n";

      if (args.count("force") == -1 && writtenScripts < joblist.size()){
        cout << "Use @force to force the writing of the scripts\n";
      }
    } else if(writtenScripts == joblist.size()) {
      cout << "NO ERRORS OR WARNINGS\n";
      cout << "All " << writtenScripts << " job scripts have been written\n";
    } else {
      cout << "UNKNOWN NUMBER OF ERRORS OR WARNINGS\n";
      cout << "There is a bug in this program...\n";
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void printInput(string ofile, input gin) {
  ofstream fout(ofile.c_str());
  const time_t t = time(0);
  fout << "TITLE\n";
  fout << "\tAutomatically generated input file\n\t";
  fout << getenv("USER") << " " << ctime(&t);
  fout << "END\n";
  fout << gin;
}

void readJobinfo(string file, map<int, jobinfo> &ji) {
  Ginstream gin(file);
  vector<string> buffer;
  gin.getblock(buffer);
  if (buffer[0] != "JOBSCRIPTS")
    throw gromos::Exception("mk_script", "Reading of jobscript file failed. "
          "No JOBSCRIPTS block");
  istringstream iss(buffer[1]);
  vector<string> head;
  string b;
  while (!iss.eof()) {
      iss >> b;
      head.push_back(b);
  }
  if (head[0] != "job_id" || head.back() != "run_after"
          || head[head.size() - 2] != "subdir")
    throw gromos::Exception("mk_script", "Reading of jobscript file failed.\n"
          "First line syntax:\n"
          "job_id PARAM PARAM ... subdir run_after");
  if (buffer.back().find("END") != 0)
    throw gromos::Exception("mk_script", "Jobscript file "
          + gin.name() +
          " is corrupted. No END in JOBSCRIPTS"
          " block. Got\n"
          + buffer.back());
  int id = 0;
  for (unsigned int i = 2; i < buffer.size() - 1; i++) {
    vector<string> tmp(head.size());
    iss.clear();
    iss.str(buffer[i]);

    for (unsigned int j = 0; j < head.size(); j++) iss >> tmp[j];
    jobinfo job;
    id = atoi(tmp[0].c_str());
    for (unsigned int j = 1; j < head.size() - 2; j++)
      job.param[head[j]] = tmp[j];
    job.dir = tmp[head.size() - 2];
    job.prev_id = atoi(tmp.back().c_str());
    // check restrictions for job id's
    if ((i > 2 && id <= job.prev_id)||(i==2 && id == job.prev_id))
      throw gromos::Exception("mk_script", "Jobscript file "
            + gin.name() +
            " is corrupted. For the first job, job_id must not be equal"
            " to run_after. For all subsequent jobs,"
            " job_id has to be larger than run_after.");
    ji[id] = job;
  }
}

void readLibrary(string file, vector<directive> &directives,
        vector<filename> &names, vector<filename> &misc,
        vector<string> &linknames, vector<int> &linkadditions,
        string system, string queue, double t,
        double dt, int ns) {
  // Open the file
  Ginstream templates(file);
  string sdum, temp, first;
  templates.getline(first);

  while (!templates.stream().eof()) {
    vector<string> buffer;
    templates.getblock(buffer);
    if (buffer.size() && first == "DIRECTIVES") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file "
              + templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);
      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        stringstream ss;
        ss << "#" << buffer[j];
        directive d(system, t, dt, ns, queue);
        d.setTemplate(ss.str());
        directives.push_back(d);
      }
    }
    if (buffer.size() && first == "FILENAMES") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file "
              + templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);
      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        istringstream iss(buffer[j]);
        iss >> sdum >> temp;
        switch (FILETYPE[sdum]) {
          case inputfile: names[inputfile].setTemplate(temp);
            break;
          case topofile: names[topofile].setTemplate(temp);
            break;
          case coordfile: names[coordfile].setTemplate(temp);
            break;
          case refposfile: names[refposfile].setTemplate(temp);
            break;
          case anatrxfile: names[anatrxfile].setTemplate(temp);
            break;
          case posresspecfile: names[posresspecfile].setTemplate(temp);
            break;
          case xrayfile: names[xrayfile].setTemplate(temp);
            break;
          case disresfile: names[disresfile].setTemplate(temp);
            break;
          case colvarresfile: names[colvarresfile].setTemplate(temp);
            break;
          case pttopofile: names[pttopofile].setTemplate(temp);
            break;
          case gamdfile: names[gamdfile].setTemplate(temp);
            break;
          case dihresfile: names[dihresfile].setTemplate(temp);
            break;
          case angresfile: names[angresfile].setTemplate(temp);
            break;
          case jvaluefile: names[jvaluefile].setTemplate(temp);
            break;
          case orderfile: names[orderfile].setTemplate(temp);
            break;
          case symfile: names[symfile].setTemplate(temp);
            break;
          case ledihfile: names[ledihfile].setTemplate(temp);
            break;
          case frictionfile: names[frictionfile].setTemplate(temp);
            break;
          case leumbfile: names[leumbfile].setTemplate(temp);
            break;
          case bsleusfile: names[bsleusfile].setTemplate(temp);
            break;
          case qmmmfile: names[qmmmfile].setTemplate(temp);
            break;
          case jinfile: names[jinfile].setTemplate(temp);
            break;
          case joutfile: names[joutfile].setTemplate(temp);
            break;
          case jtrjfile: names[jtrjfile].setTemplate(temp);
            break;
          case outputfile: names[outputfile].setTemplate(temp);
            break;
          case outtrxfile: names[outtrxfile].setTemplate(temp);
            break;
          case outtrvfile: names[outtrvfile].setTemplate(temp);
            break;
          case outtrffile: names[outtrffile].setTemplate(temp);
            break;
          case outtrefile: names[outtrefile].setTemplate(temp);
            break;
          case outtrgfile: names[outtrgfile].setTemplate(temp);
            break;
          case outbaefile: names[outbaefile].setTemplate(temp);
            break;
          case outbagfile: names[outbagfile].setTemplate(temp);
            break;
          case scriptfile: names[scriptfile].setTemplate(temp);
            break;
          case outtrsfile: names[outtrsfile].setTemplate(temp);
            break;
          case repoutfile: names[repoutfile].setTemplate(temp);
            break;
          case repdatfile: names[repdatfile].setTemplate(temp);
            break;
          case unknownfile:
            printWarning("Don't know how to handle template for " + sdum
                    + ". Ingoring");
        }
      }
    }
    if (buffer.size() && first == "MISCELLANEOUS") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file " +
              templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        istringstream iss(buffer[j]);
        iss >> sdum;
        if (sdum == "workdir") {
          iss >> temp;
          temp = word_expansion(temp);
          misc[0].setTemplate(temp);
        }
        if (sdum == "lastcommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[1].setTemplate(os.str());
        }
        if (sdum == "firstcommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[2].setTemplate(os.str());
        }
        if (sdum == "mpicommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[3].setTemplate(os.str());
        }
        if (sdum == "stopcommand") {
          ostringstream os;
          while (!iss.eof()) {
            iss >> sdum;
            os << sdum << " ";
          }
          misc[4].setTemplate(os.str());
        }
      }
    }
    if (buffer.size() && first == "LINKADDITION") {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("mk_script", "Template file " +
              templates.name() +
              " is corrupted. No END in " + first +
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 0; j < buffer.size() - 1; j++) {
        istringstream iss(buffer[j]);
        int k;
        string varname;

        iss >> sdum >> varname >> temp >> k;
        filename newlink(system, t, dt, ns, queue);
        newlink.setTemplate(temp);
        names.push_back(newlink);
        if (sdum == "input") k *= -1;
        linkadditions.push_back(k);
        linknames.push_back(varname);

      }
    }
    templates.getline(first);
  }
}

void setParam(input &gin, jobinfo const &job) {
  map<string, string>::const_iterator iter = job.param.begin(),
          to = job.param.end();
  for (; iter != to; ++iter) {
    // std::cout << "iter = " << iter->first << std::endl;
    if (iter->first == "ENDTIME")
      ; // do nothing but avoid warning

      // AEDS
    else if (iter->first.substr(0, 4) == "EIR[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.aeds.numstates)
        gin.aeds.eir[i - 1] = atof(iter->second.c_str());
      else {
        std::stringstream ss;
        ss << iter->second.c_str() << " in joblistfile out of range";
        printError(ss.str());
      }
    }
    else if (iter->first == "NTIAEDSS")
      gin.aeds.ntiaedss = atoi(iter->second.c_str());
    else if (iter->first == "ASTEPS")
      gin.aeds.asteps = atoi(iter->second.c_str());
    else if (iter->first == "BSTEPS")
      gin.aeds.bsteps = atoi(iter->second.c_str());

    // GAMD ToDo: Allow more variables to be changed
    else if (iter->first == "GAMD")
      gin.gamd.gamd = atoi(iter->second.c_str());
    else if (iter->first == "SEARCH")
      gin.gamd.search = atoi(iter->second.c_str());
    else if (iter->first == "FORM")
      gin.gamd.form = atoi(iter->second.c_str());
    else if (iter->first == "NTIGAMDS")
      gin.gamd.ntigamds = atoi(iter->second.c_str());
    else if (iter->first == "EQSTEPS")
      gin.gamd.eqsteps = atoi(iter->second.c_str());
    else if (iter->first == "WINDOW")
      gin.gamd.window = atoi(iter->second.c_str());


      // BAROSTAT
    else if (iter->first == "NTP")
      gin.barostat.ntp = atoi(iter->second.c_str());
    else if (gin.barostat.found && iter->first == "COMP")
      gin.barostat.comp = atof(iter->second.c_str());
    else if (iter->first.substr(0, 7) == "PRSBTH[") {
      unsigned int i = atoi(iter->first.substr(7, iter->first.find("]")).c_str());
      if (i <= gin.barostat.pbaths.size())
        gin.barostat.pbaths[i - 1].prsbth = atof(iter->second.c_str());
      else {
        std::stringstream ss;
        ss << iter->second.c_str() << " in joblistfile out of range";
        printError(ss.str());
      }
    }      // BOUNDCOND
    else if (iter->first == "NTB")
      gin.boundcond.ntb = atoi(iter->second.c_str());
    else if (iter->first == "NDFMIN")
      gin.boundcond.ndfmin = atoi(iter->second.c_str());

      // CGRAIN
    else if (iter->first == "NTCGRAN")
      gin.cgrain.ntcgran = atoi(iter->second.c_str());
    else if (iter->first == "EPS")
      gin.cgrain.eps = atof(iter->second.c_str());
    else if (iter->first == "EPSM")
      gin.cgrain.epsm = atof(iter->second.c_str());

      // COMTRANSROT
    else if (iter->first == "NSCM")
      gin.comtransrot.nscm = atoi(iter->second.c_str());

      // CONSISTENCYCHECK
    else if (iter->first == "NTCHK")
      gin.consistencycheck.ntchk = atoi(iter->second.c_str());
    else if (iter->first == "NTCKF")
      gin.consistencycheck.ntckf = atoi(iter->second.c_str());
    else if (iter->first == "FDCKF")
      gin.consistencycheck.fdckf = atof(iter->second.c_str());
    else if (iter->first == "NTCKV")
      gin.consistencycheck.ntckv = atoi(iter->second.c_str());
    else if (iter->first == "FDCKV")
      gin.consistencycheck.fdckv = atof(iter->second.c_str());
    else if (iter->first == "NTCKT")
      gin.consistencycheck.ntckt = atoi(iter->second.c_str());
    else if (iter->first == "NTCKE")
      gin.consistencycheck.ntcke = atoi(iter->second.c_str());
    else if (iter->first == "NTCKR")
      gin.consistencycheck.ntckr = atoi(iter->second.c_str());
    else if (iter->first == "NTCKL")
      gin.consistencycheck.ntckl = atoi(iter->second.c_str());
    else if (iter->first == "FDCKL")
      gin.consistencycheck.fdckl = atof(iter->second.c_str());

      // CONSTRAINT
    else if (iter->first == "NTC")
      gin.constraint.ntc = atoi(iter->second.c_str());
    else if (iter->first == "NTCP")
      gin.constraint.ntcp = atoi(iter->second.c_str());
    else if (iter->first == "NTCS")
      gin.constraint.ntcs = atoi(iter->second.c_str());
    else if (iter->first == "NTCG")
      gin.constraint.ntcg = atoi(iter->second.c_str());
    else if (iter->first == "NTCD")
      gin.constraint.ntcd[0] = atoi(iter->second.c_str());

      // COVALENTFORM
    else if (iter->first == "NTBBH")
      gin.covalentform.ntbbh = atoi(iter->second.c_str());
    else if (iter->first == "NTBAH")
      gin.covalentform.ntbah = atoi(iter->second.c_str());
    else if (iter->first == "NTBDN")
      gin.covalentform.ntbdn = atoi(iter->second.c_str());

      // DEBUG
      // DIHEDRALRES
    else if (iter->first == "NTDLR")
      gin.dihedralres.ntdlr = atoi(iter->second.c_str());
    else if (iter->first == "CDLR")
      gin.dihedralres.cdlr = atof(iter->second.c_str());
    else if (iter->first == "PHILIN")
      gin.dihedralres.philin = atof(iter->second.c_str());
    else if (iter->first == "VDIH")
      gin.dihedralres.vdih = atof(iter->second.c_str());
    else if (iter->first == "NTWDLR")
      gin.dihedralres.ntwdlr = atof(iter->second.c_str());
    else if (iter->first == "TOLDAC")
      gin.dihedralres.toldac = atof(iter->second.c_str());
      
      // ANGLERES
    else if (iter->first == "NTALR")
      gin.angleres.ntalr = atoi(iter->second.c_str());
    else if (iter->first == "CALR")
      gin.angleres.calr = atof(iter->second.c_str());
    else if (iter->first == "VARES")
      gin.angleres.vares = atof(iter->second.c_str());
    else if (iter->first == "NTWALR")
      gin.angleres.ntwalr = atof(iter->second.c_str());
    else if (iter->first == "TOLBAC")
      gin.angleres.tolbac = atof(iter->second.c_str());

    // DISTANCEFIELD
    else if(iter->first == "NTDFR")
      gin.distancefield.ntdfr = atoi(iter->second.c_str());
    else if(iter->first == "GRID")
      gin.distancefield.grid = atof(iter->second.c_str());
    else if(iter->first == "PROTEINOFFSET")
      gin.distancefield.proteinoffset = atof(iter->second.c_str());
    else if(iter->first == "PROTEINCUTOFF")
      gin.distancefield.proteincutoff = atof(iter->second.c_str());
    else if(iter->first == "UPDATE")
      gin.distancefield.update = atoi(iter->second.c_str());
    else if(iter->first == "SMOOTH")
      gin.distancefield.smooth = atoi(iter->second.c_str());
    else if(iter->first == "RL")
      gin.distancefield.rl = atof(iter->second.c_str());
    else if(iter->first == "NTWDF")
      gin.distancefield.ntwdf = atoi(iter->second.c_str());
    else if(iter->first == "PRINTGRID")
      gin.distancefield.printgrid = atoi(iter->second.c_str());
    else if(iter->first == "PROTECT")
      gin.distancefield.protect = atof(iter->second.c_str());
    
      //DISTANCERES
    else if (iter->first == "NTDIR")
      gin.distanceres.ntdir = atoi(iter->second.c_str());
    else if (iter->first == "NTDIRA")
      gin.distanceres.ntdira = atoi(iter->second.c_str());
    else if (iter->first == "CDIR")
      gin.distanceres.cdir = atof(iter->second.c_str());
    else if (iter->first == "DIR0")
      gin.distanceres.dir0 = atof(iter->second.c_str());
    else if (iter->first == "TAUDIR")
      gin.distanceres.taudir = atoi(iter->second.c_str());
    else if (iter->first == "FORCESCALE")
      gin.distanceres.forcescale = atoi(iter->second.c_str());
    else if(iter->first == "VDIR")
      gin.distanceres.vdir = atoi(iter->second.c_str());
      
      //COLVARRES
    else if (iter->first == "CVR")
      gin.colvarres.cvr = atoi(iter->second.c_str());
    else if (iter->first == "CVK")
      gin.colvarres.cvk = atof(iter->second.c_str());
    else if(iter->first == "TAUCVR")
      gin.colvarres.taucvr = atoi(iter->second.c_str());
    else if(iter->first == "VCVR")
      gin.colvarres.vcvr = atoi(iter->second.c_str());
    else if(iter->first == "NTWCV")
      gin.colvarres.ntwcv = atoi(iter->second.c_str());

      // ELECTRIC
     else if(iter->first == "EF_x")
      gin.electric.ef_x = atof(iter->second.c_str());
     else if(iter->first == "EF_y")
      gin.electric.ef_y = atof(iter->second.c_str());
     else if(iter->first == "EF_z")
      gin.electric.ef_z = atof(iter->second.c_str());

      // ENERGYMIN
    else if (iter->first == "NTEM")
      gin.energymin.ntem = atoi(iter->second.c_str());
    else if (iter->first == "NCYC")
      gin.energymin.ncyc = atoi(iter->second.c_str());
    else if (iter->first == "DELE")
      gin.energymin.dele = atof(iter->second.c_str());
    else if (iter->first == "DX0")
      gin.energymin.dx0 = atof(iter->second.c_str());
    else if (iter->first == "DXM")
      gin.energymin.dxm = atof(iter->second.c_str());
    else if (iter->first == "NMIN")
      gin.energymin.nmin = atoi(iter->second.c_str());
    else if (iter->first == "FLIM")
      gin.energymin.flim = atof(iter->second.c_str());
    else if (iter->first == "CGIM")
      gin.energymin.cgim = atof(iter->second.c_str());
    else if (iter->first == "CGIC")
      gin.energymin.cgic = atof(iter->second.c_str());

      // EWARN
    else if (iter->first == "MAXENER")
      gin.ewarn.maxener = atof(iter->second.c_str());

      // FORCE
    else if (iter->first == "NTF[1]")
      gin.force.ntf[0] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[2]")
      gin.force.ntf[1] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[3]")
      gin.force.ntf[2] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[4]")
      gin.force.ntf[3] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[5]")
      gin.force.ntf[4] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[6]")
      gin.force.ntf[5] = atoi(iter->second.c_str());
    /*else if (iter->first == "NTF[7]")
      gin.force.ntf[6] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[8]")
      gin.force.ntf[7] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[9]")
      gin.force.ntf[8] = atoi(iter->second.c_str());
    else if (iter->first == "NTF[10]")
      gin.force.ntf[9] = atoi(iter->second.c_str());*/

      // GEOMCONSTRAINTS
    else if (iter->first == "NTCPH")
      gin.geomconstraints.ntcph = atoi(iter->second.c_str());
    else if (iter->first == "NTCPN")
      gin.geomconstraints.ntcpn = atoi(iter->second.c_str());
    else if (iter->first == "NTCS")
      gin.geomconstraints.ntcs = atoi(iter->second.c_str());
    else if (iter->first == "SHKTOL")
      gin.geomconstraints.shktol = atof(iter->second.c_str());

      // GROMOS96COMPAT
    else if (iter->first == "NTNB96")
      gin.gromos96compat.ntnb96 = atoi(iter->second.c_str());
    else if (iter->first == "NTR96")
      gin.gromos96compat.ntr96 = atoi(iter->second.c_str());
    else if (iter->first == "NTP96")
      gin.gromos96compat.ntp96 = atoi(iter->second.c_str());
    else if (iter->first == "NTG96")
      gin.gromos96compat.ntg96 = atoi(iter->second.c_str());

      // INITIALISE
    else if (iter->first == "NTIVEL")
      gin.initialise.ntivel = atoi(iter->second.c_str());
    else if (iter->first == "NTISHK")
      gin.initialise.ntishk = atoi(iter->second.c_str());
    else if (iter->first == "NTINHT")
      gin.initialise.ntinht = atoi(iter->second.c_str());
    else if (iter->first == "NTINHB")
      gin.initialise.ntinhb = atoi(iter->second.c_str());
    else if (iter->first == "NTISHI")
      gin.initialise.ntishi = atoi(iter->second.c_str());
    else if (iter->first == "NTIRTC")
      gin.initialise.ntirtc = atoi(iter->second.c_str());
    else if (iter->first == "NTICOM")
      gin.initialise.nticom = atoi(iter->second.c_str());
    else if (iter->first == "NTISTI")
      gin.initialise.ntisti = atoi(iter->second.c_str());
    else if (iter->first == "IG")
      gin.initialise.ig = atoi(iter->second.c_str());
    else if (iter->first == "TEMPI")
      gin.initialise.tempi = atof(iter->second.c_str());

      // INNERLOOP
    else if (iter->first == "NTILM")
      gin.innerloop.ntilm = atoi(iter->second.c_str());
    else if (iter->first == "NTILS")
      gin.innerloop.ntils = atoi(iter->second.c_str());

      // INTEGRATE
    else if (iter->first == "NINT")
      gin.integrate.nint = atoi(iter->second.c_str());

      // JVALUERES
    else if (iter->first == "NTJVR")
      gin.jvalueres.ntjvr = atoi(iter->second.c_str());
    else if (iter->first == "NTJVRA")
      gin.jvalueres.ntjvra = atoi(iter->second.c_str());
    else if (iter->first == "CJVR")
      gin.jvalueres.cjvr = atof(iter->second.c_str());
    else if (iter->first == "TAUJVR")
      gin.jvalueres.taujvr = atof(iter->second.c_str());
    else if (iter->first == "NJVRTARS")
      gin.jvalueres.njvrtars = atoi(iter->second.c_str());
    else if (iter->first == "NJVRBIQW")
      gin.jvalueres.njvrbiqw = atoi(iter->second.c_str());
    else if (iter->first == "LE")
      gin.jvalueres.le = atoi(iter->second.c_str());
    else if (iter->first == "NGRID")
      gin.jvalueres.ngrid = atoi(iter->second.c_str());
    else if (iter->first == "DELTA")
      gin.jvalueres.delta = atof(iter->second.c_str());
    else if (iter->first == "NTWJV")
      gin.jvalueres.write = atoi(iter->second.c_str());
    
      // ORDERPARAMRES
    else if (iter->first == "NTOPR")
      gin.orderparamres.ntopr = atoi(iter->second.c_str());
    else if (iter->first == "NTOPRA")
      gin.orderparamres.ntopra = atoi(iter->second.c_str());
    else if (iter->first == "COPR")
      gin.orderparamres.copr = atof(iter->second.c_str());
    else if (iter->first == "TAUOPR")
      gin.orderparamres.tauopr = atof(iter->second.c_str());
    else if (iter->first == "UPDOPR")
      gin.orderparamres.updopr = atoi(iter->second.c_str());
    else if (iter->first == "NTWOP")
      gin.orderparamres.ntwop = atoi(iter->second.c_str());
    
      // SYMRES
    else if (iter->first == "NTSYM")
      gin.symres.ntsym = atoi(iter->second.c_str());
    else if (iter->first == "CSYM")
      gin.symres.csym = atof(iter->second.c_str());

      //LAMBDAS
    else if (iter->first.substr(0, 4) == "ALI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].ali = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "BLI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].bli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "CLI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].cli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "DLI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].dli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "ELI[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.lambdas.lambints.size())
        gin.lambdas.lambints[i - 1].eli = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    }

    // LOCALELEV
    else if (iter->first == "NTLES")
      gin.localelev.ntles = atoi(iter->second.c_str());
    else if (iter->first == "NLEPOT")
      gin.localelev.nlepot = atoi(iter->second.c_str());
    else if (iter->first == "NTLESA")
      gin.localelev.ntlesa = atoi(iter->second.c_str());
    else if (iter->first.substr(0,6) == "NLEPID[") {
      int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if(i > gin.localelev.nlepot) {
        stringstream msg;
        msg << "NLEPID[" << i << "] cannot be set when NLEPOT = " << gin.localelev.nlepot
                << " (LOCALELEV block)";
        printError(msg.str());
      }
      int count = 0;
      map<int, int> tmp;
      for(map<int,int>::iterator it = gin.localelev.nlepid_ntlepfr.begin();
      it != gin.localelev.nlepid_ntlepfr.end(); ++it){
        if(count == i) {
          tmp.insert( pair<int, int>(atoi(iter->second.c_str()), it->second));
        } else {
          tmp.insert( pair<int, int>(it->first, it->second));
        }
        count++;
      }
      gin.localelev.nlepid_ntlepfr = tmp;
      if(int(gin.localelev.nlepid_ntlepfr.size()) != gin.localelev.nlepot) {
        printError("NLEPID in LOCALELEV block is ambiguous");
      }
    }
    else if (iter->first.substr(0,6) == "NTLEFR[") {
      int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      int count = 0;
      for(map<int,int>::iterator it = gin.localelev.nlepid_ntlepfr.begin();
      it != gin.localelev.nlepid_ntlepfr.end(); ++it){
        if(count == i) {
          it->second = atoi(iter->second.c_str());
          if (it->second < 0 || it->second >1) {
            printError("NTLEFR in LOCALELEV must be 0 or 1");
          }
          break;
        }
        count++;
      }
    }
    

      // MULTIBATH
    else if (iter->first == "NTBTYP")
      gin.multibath.ntbtyp = atoi(iter->second.c_str());
    else if (iter->first.substr(0, 6) == "TEMP0[") {
      int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if (i <= gin.multibath.nbaths)
        gin.multibath.temp0[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 4) == "TAU[") {
      int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.multibath.nbaths)
        gin.multibath.tau[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    }
    else if (iter->first == "NUM")
      gin.multibath.num = atoi(iter->second.c_str());


    // MULTICELL

    // MULTIGRADIENT
    else if (iter->first == "NTMGRE")
      gin.multigradient.ntmgre = atoi(iter->second.c_str());
    else if (iter->first.substr(0, 7) == "MGRCPT[") {
      size_t firstBracket = iter->first.find("[");
      size_t secondBracket = iter->first.find("[", firstBracket+1);
      int i = int(atof(iter->first.substr(firstBracket+1, iter->first.find("]", firstBracket)).c_str()) - 1);
      int j = int(atof(iter->first.substr(secondBracket+1, iter->first.find("]", secondBracket)).c_str()) - 1);
      double d = gin.multigradient.curves[i][j].second;
      pair<double, double> p (atof(iter->second.c_str()), d);
      gin.multigradient.curves[i][j] = p;
    }
    else if  (iter->first.substr(0, 7) == "MGRCPV[") {
      size_t firstBracket = iter->first.find("[");
      size_t secondBracket = iter->first.find("[", firstBracket+1);
      int i = int(atof(iter->first.substr(firstBracket+1, iter->first.find("]", firstBracket)).c_str()) - 1);
      int j = int(atof(iter->first.substr(secondBracket+1, iter->first.find("]", secondBracket)).c_str()) - 1);
      double d = gin.multigradient.curves[i][j].first;
      pair<double, double> p (d, atof(iter->second.c_str()));
      gin.multigradient.curves[i][j] = p;
    }

    // MULTISTEP
    else if (iter->first == "STEPS")
      gin.multistep.steps = atoi(iter->second.c_str());
    else if (iter->first == "BOOST")
      gin.multistep.boost = atoi(iter->second.c_str());

      // NEIGHBOURLIST
    else if (iter->first == "PLALGO")
      gin.neighbourlist.plalgo = atoi(iter->second.c_str());
    else if (iter->first == "NUPDPL")
      gin.neighbourlist.nupdpl = atoi(iter->second.c_str());
    else if (iter->first == "NUPDIS")
      gin.neighbourlist.nupdis = atoi(iter->second.c_str());
    else if (iter->first == "NUPDII")
      gin.neighbourlist.nupdii = atoi(iter->second.c_str());
    else if (iter->first == "RCUTS")
      gin.neighbourlist.rcuts = atof(iter->second.c_str());
    else if (iter->first == "RCUTI")
      gin.neighbourlist.rcuti = atof(iter->second.c_str());
    else if (iter->first == "GRIDSZX")
      gin.neighbourlist.gridszx = atof(iter->second.c_str());
    else if (iter->first == "GRIDSZY")
      gin.neighbourlist.gridszy = atof(iter->second.c_str());
    else if (iter->first == "GRIDSZZ")
      gin.neighbourlist.gridszz = atof(iter->second.c_str());
    else if (iter->first == "TYPE")
      gin.neighbourlist.type = atoi(iter->second.c_str());
    else if (iter->first == "NCGCEN")
      gin.neighbourlist.ncgcen = atoi(iter->second.c_str());

      // NONBONDED
    else if (iter->first == "NLRELE")
      gin.nonbonded.nlrele = atoi(iter->second.c_str());
    else if (iter->first == "APPAK")
      gin.nonbonded.appak = atof(iter->second.c_str());
    else if (iter->first == "RCRF")
      gin.nonbonded.rcrf = atof(iter->second.c_str());
    else if (iter->first == "EPSRF")
      gin.nonbonded.epsrf = atof(iter->second.c_str());
    else if (iter->first == "NSHAPE")
      gin.nonbonded.nshape = atoi(iter->second.c_str());
    else if (iter->first == "NSLFEXCL")
      gin.nonbonded.nslfexcl = atof(iter->second.c_str());
    else if (iter->first == "ASHAPE")
      gin.nonbonded.ashape = atof(iter->second.c_str());
    else if (iter->first == "NA2CLC")
      gin.nonbonded.na2clc = atoi(iter->second.c_str());
    else if (iter->first == "TOLA2")
      gin.nonbonded.tola2 = atof(iter->second.c_str());
    else if (iter->first == "EPSLS")
      gin.nonbonded.epsls = atof(iter->second.c_str());
    else if (iter->first == "NKX")
      gin.nonbonded.nkx = atoi(iter->second.c_str());
    else if (iter->first == "NKY")
      gin.nonbonded.nky = atoi(iter->second.c_str());
    else if (iter->first == "NKZ")
      gin.nonbonded.nkz = atoi(iter->second.c_str());
    else if (iter->first == "KCUT")
      gin.nonbonded.kcut = atof(iter->second.c_str());
    else if (iter->first == "NGX")
      gin.nonbonded.ngx = atoi(iter->second.c_str());
    else if (iter->first == "NGY")
      gin.nonbonded.ngy = atoi(iter->second.c_str());
    else if (iter->first == "NGZ")
      gin.nonbonded.ngz = atoi(iter->second.c_str());
    else if (iter->first == "NASORD")
      gin.nonbonded.nasord = atoi(iter->second.c_str());
    else if (iter->first == "NFDORD")
      gin.nonbonded.nfdord = atoi(iter->second.c_str());
    else if (iter->first == "NALIAS")
      gin.nonbonded.nalias = atoi(iter->second.c_str());
    else if (iter->first == "NSPORD")
      gin.nonbonded.nspord = atoi(iter->second.c_str());
    else if (iter->first == "NQEVAL")
      gin.nonbonded.nqeval = atoi(iter->second.c_str());
    else if (iter->first == "FACCUR")
      gin.nonbonded.faccur = atof(iter->second.c_str());
    else if (iter->first == "NRDGRD")
      gin.nonbonded.nrdgrd = atoi(iter->second.c_str());
    else if (iter->first == "NWRGRD")
      gin.nonbonded.nwrgrd = atoi(iter->second.c_str());
    else if (iter->first == "NLRLJ")
      gin.nonbonded.nlrlj = atoi(iter->second.c_str());
    else if (iter->first == "SLVDNS")
      gin.nonbonded.slvdns = atoi(iter->second.c_str());

      // OVERALLTRANSROT
    else if (iter->first == "NCMTR")
      gin.overalltransrot.ncmtr = atoi(iter->second.c_str());
    else if (iter->first == "NCMRO")
      gin.overalltransrot.ncmro = atoi(iter->second.c_str());
    else if (iter->first == "CMAMX")
      gin.overalltransrot.cmamx = atoi(iter->second.c_str());
    else if (iter->first == "CMAMY")
      gin.overalltransrot.cmamy = atoi(iter->second.c_str());
    else if (iter->first == "CMAMZ")
      gin.overalltransrot.cmamz = atoi(iter->second.c_str());

      // PAIRLIST
    else if (iter->first == "NSNB")
      gin.pairlist.nsnb = atoi(iter->second.c_str());
    else if (iter->first == "RCUTP")
      gin.pairlist.rcutp = atof(iter->second.c_str());
    else if (iter->first == "RCUTL")
      gin.pairlist.rcutl = atof(iter->second.c_str());
    else if (iter->first == "TYPE")
      gin.pairlist.type = atoi(iter->second.c_str());

      // PATHINT
    else if (iter->first == "NTPI")
      gin.pathint.ntpi = atoi(iter->second.c_str());

      // PERSCALE
    else if (iter->first == "RESTYPE")
      gin.perscale.restype = atoi(iter->second.c_str());
    else if (iter->first == "KDIH")
      gin.perscale.kdih = atof(iter->second.c_str());
    else if (iter->first == "KJ")
      gin.perscale.kj = atof(iter->second.c_str());
    else if (iter->first == "DIFF")
      gin.perscale.diff = atof(iter->second.c_str());
    else if (iter->first == "RATIO")
      gin.perscale.ratio = atof(iter->second.c_str());
    else if (iter->first == "READ")
      gin.perscale.read = atoi(iter->second.c_str());

      // PERTURBATION
    else if (iter->first == "NTG")
      gin.perturbation.ntg = atoi(iter->second.c_str());
    else if (iter->first == "NRDGL")
      gin.perturbation.nrdgl = atoi(iter->second.c_str());
    else if (iter->first == "RLAM")
      gin.perturbation.rlam = atof(iter->second.c_str());
    else if (iter->first == "DLAMT")
      gin.perturbation.dlamt = atof(iter->second.c_str());
    else if (iter->first == "ALPHLJ")
      gin.perturbation.alphlj = atof(iter->second.c_str());
    else if (iter->first == "ALPHC")
      gin.perturbation.alphc = atof(iter->second.c_str());
    else if (iter->first == "NLAM")
      gin.perturbation.nlam = atoi(iter->second.c_str());
    else if (iter->first == "NSCALE")
      gin.perturbation.nscale = atoi(iter->second.c_str());

      // POLARISE
    else if (iter->first == "COS")
      gin.polarise.cos = atoi(iter->second.c_str());
    else if (iter->first == "EFIELD")
      gin.polarise.efield = atoi(iter->second.c_str());
    else if (iter->first == "MINFIELD")
      gin.polarise.minfield = atof(iter->second.c_str());
    else if (iter->first == "DAMP")
      gin.polarise.damp = atoi(iter->second.c_str());
    else if (iter->first == "WRITE")
      gin.polarise.write = atoi(iter->second.c_str());

      // POSITIONRES
    else if (iter->first == "NTPOR")
      gin.positionres.ntpor = atoi(iter->second.c_str());
    else if (iter->first == "NTPORB")
      gin.positionres.ntporb = atoi(iter->second.c_str());
    else if (iter->first == "NTPORS")
      gin.positionres.ntpors = atoi(iter->second.c_str());
    else if (iter->first == "CPOR")
      gin.positionres.cpor = atof(iter->second.c_str());
      
      // PRECALCLAM
    else if (iter->first == "NRLAM")
      gin.precalclam.nrlam = atoi(iter->second.c_str());
    else if (iter->first == "MINLAM")
      gin.precalclam.nrlam = atoi(iter->second.c_str());
    else if (iter->first == "MAXLAM")
      gin.precalclam.nrlam = atoi(iter->second.c_str());

      // PRESSURESCALE
    else if (iter->first == "COUPLE")
      gin.pressurescale.couple = atoi(iter->second.c_str());
    else if (iter->first == "SCALE")
      gin.pressurescale.scale = atoi(iter->second.c_str());
    else if (gin.pressurescale.found && iter->first == "COMP")
      gin.pressurescale.comp = atof(iter->second.c_str());
    else if (iter->first == "TAUP")
      gin.pressurescale.taup = atof(iter->second.c_str());
    else if (iter->first == "VIRIAL")
      gin.pressurescale.virial = atoi(iter->second.c_str());
    else if (iter->first == "X_SEMI")
      gin.pressurescale.x_semi = atoi(iter->second.c_str());
    else if (iter->first == "Y_SEMI")
      gin.pressurescale.y_semi = atoi(iter->second.c_str());
    else if (iter->first == "Z_SEMI")
      gin.pressurescale.z_semi = atoi(iter->second.c_str());

      // PRINTOUT
    else if (iter->first == "NTPR")
      gin.printout.ntpr = atoi(iter->second.c_str());
    else if (iter->first == "NTPP")
      gin.printout.ntpp = atoi(iter->second.c_str());

      // QMMM
    else if (iter->first == "NTQMMM")
      gin.qmmm.ntqmmm = atoi(iter->second.c_str());
    else if (iter->first == "NTQMSW")
      gin.qmmm.ntqmsw = atoi(iter->second.c_str());
    else if (iter->first == "RCUTQM")
      gin.qmmm.rcutqm = atof(iter->second.c_str());
    else if (iter->first == "NTWQMMM")
      gin.qmmm.ntwqmmm = atoi(iter->second.c_str());
    else if (iter->first == "QMLJ")
      gin.qmmm.qmlj = atoi(iter->second.c_str());
    else if (iter->first == "QMCON")
      gin.qmmm.qmcon = atoi(iter->second.c_str());
    else if (iter->first == "MMSCALE")
      gin.qmmm.mmscale = atof(iter->second.c_str());

      // RANDOMNUMBERS
    else if (iter->first == "NTRNG")
      gin.randomnumbers.ntrng = atoi(iter->second.c_str());
    else if (iter->first == "NTGSL")
      gin.randomnumbers.ntgsl = atoi(iter->second.c_str());

      // READTRAJ
    else if (iter->first == "NTRD")
      gin.readtraj.ntrd = atoi(iter->second.c_str());
    else if (iter->first == "NTSTR")
      gin.readtraj.ntstr = atoi(iter->second.c_str());
    else if (iter->first == "NTRB")
      gin.readtraj.ntrb = atoi(iter->second.c_str());
    else if (iter->first == "NTSHK")
      gin.readtraj.ntshk = atoi(iter->second.c_str());

      // REPLICA
    else if (iter->first.substr(0, 4) == "RET[") {
      unsigned int i = atoi(iter->first.substr(4, iter->first.find("]")).c_str());
      if (i <= gin.replica.ret.size())
        gin.replica.ret[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 6) == "RELAM[") {
      unsigned int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if (i <= gin.replica.relam.size())
        gin.replica.relam[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first.substr(0, 5) == "RETS[") {
      unsigned int i = atoi(iter->first.substr(6, iter->first.find("]")).c_str());
      if (i <= gin.replica.rets.size())
        gin.replica.rets[i - 1] = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist out of range");
    } else if (iter->first == "LRESCALE")
      gin.replica.lrescale = atoi(iter->second.c_str());
    else if (iter->first == "NRETRIAL" and gin.replica.found)
      gin.replica.nretrial = atoi(iter->second.c_str());
    else if (iter->first == "NREQUIL" and gin.replica.found)
      gin.replica.nrequil = atoi(iter->second.c_str());
    else if (iter->first == "CONT" and gin.replica.found)
      gin.replica.cont = atoi(iter->second.c_str());

    // REEDS
    else if (iter->first == "REEDS")
      gin.reeds.reeds = atoi(iter->second.c_str());
    else if (iter->first == "NRETRIAL" and gin.reeds.found)
      gin.reeds.nretrial = atoi(iter->second.c_str());
    else if (iter->first == "NREQUIL" and gin.reeds.found)
      gin.reeds.nrequil = atoi(iter->second.c_str());

      // ROTTRANS
    else if (iter->first == "RTC")
      gin.rottrans.rtc = atoi(iter->second.c_str());

      // STEP
    else if (iter->first == "NSTLIM")
      gin.step.nstlim = atoi(iter->second.c_str());
    else if (iter->first == "T")
      gin.step.t = atof(iter->second.c_str());
    else if (iter->first == "DT")
      gin.step.dt = atof(iter->second.c_str());

      // STOCHDYN
    else if (iter->first == "NTSD")
      gin.stochdyn.ntsd = atoi(iter->second.c_str());
    else if (iter->first == "NTFR")
      gin.stochdyn.ntfr = atoi(iter->second.c_str());
    else if (iter->first == "NSFR")
      gin.stochdyn.nsfr = atoi(iter->second.c_str());
    else if (iter->first == "NBREF")
      gin.stochdyn.nbref = atoi(iter->second.c_str());
    else if (iter->first == "RCUTF")
      gin.stochdyn.rcutf = atof(iter->second.c_str());
    else if (iter->first == "CFRIC")
      gin.stochdyn.cfric = atof(iter->second.c_str());
    else if (iter->first == "TEMPSD")
      gin.stochdyn.tempsd = atof(iter->second.c_str());


      // SYSTEM
    else if (iter->first == "NPM")
      gin.system.npm = atoi(iter->second.c_str());
    else if (iter->first == "NSM")
      gin.system.nsm = atoi(iter->second.c_str());

      // THERMOSTAT
    else if (iter->first == "NTT")
      gin.thermostat.ntt = atoi(iter->second.c_str());
    else if (iter->first.substr(0, 7) == "TEMBTH[") {
      unsigned int i = atoi(iter->first.substr(7, iter->first.find("]")).c_str());
      if (i <= gin.thermostat.baths.size())
        gin.thermostat.baths[i - 1].tembth = atof(iter->second.c_str());
      else
        printError(iter->first + " in joblist is out of range");
    }// I don't really dare to let the use change other things...

      // UMBRELLA
    else if (iter->first == "NTUS")
      gin.umbrella.ntus = atoi(iter->second.c_str());
    else if (iter->first == "USCST1")
      gin.umbrella.uscst1 = atof(iter->second.c_str());
    else if (iter->first == "USCST2")
      gin.umbrella.uscst2 = atof(iter->second.c_str());
    else if (iter->first == "USREF1")
      gin.umbrella.usref1 = atof(iter->second.c_str());
    else if (iter->first == "USREF2")
      gin.umbrella.usref2 = atof(iter->second.c_str());

      // VIRIAL
    else if (iter->first == "NTV")
      gin.virial.ntv = atoi(iter->second.c_str());
    else if (iter->first == "NTVG")
      gin.virial.ntvg = atoi(iter->second.c_str());

      //WRITETRAJ
    else if (iter->first == "NTWX")
      gin.writetraj.ntwx = atoi(iter->second.c_str());
    else if (iter->first == "NTWSE")
      gin.writetraj.ntwse = atoi(iter->second.c_str());
    else if (iter->first == "NTWV")
      gin.writetraj.ntwv = atoi(iter->second.c_str());
    else if (iter->first == "NTWF")
      gin.writetraj.ntwf = atoi(iter->second.c_str());
    else if (iter->first == "NTWE")
      gin.writetraj.ntwe = atoi(iter->second.c_str());
    else if (iter->first == "NTWG")
      gin.writetraj.ntwe = atoi(iter->second.c_str());
    else if (iter->first == "NTWB")
      gin.writetraj.ntwb = atoi(iter->second.c_str());


      // XRAYRES
    else if (iter->first == "NTXR")
      gin.xrayres.ntxr = atoi(iter->second.c_str());
    else if (iter->first == "NTXLE")
      gin.xrayres.ntxle = atoi(iter->second.c_str());
    else if (iter->first == "CXR")
      gin.xrayres.cxr = atof(iter->second.c_str());
    else if (iter->first == "NTWXR")
      gin.xrayres.ntwxr = atoi(iter->second.c_str());
    else if (iter->first == "NTWDE")
      gin.xrayres.ntwde = atoi(iter->second.c_str());
    else if (iter->first == "NTWX<")
      gin.xrayres.ntwxm = atoi(iter->second.c_str());
    else if (iter->first == "CXTAU")
      gin.xrayres.cxtau = atof(iter->second.c_str());
    else if (iter->first == "RDAVG")
      gin.xrayres.rdavg = atoi(iter->second.c_str());

    else
      throw gromos::Exception("mk_script", "Cannot automatically change "
            + iter->first + " in input file");
  }
}

bool file_exists(string file){
    ifstream infile(file.c_str()); //destructor closes the file
    return infile.is_open();
}

#ifdef HAVE_WORDEXP_H
#include <wordexp.h>
#endif

string word_expansion(const string & arg) {
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE) && defined(HAVE_WORDEXP_H)
  wordexp_t p;
  if (wordexp(arg.c_str(), &p, WRDE_SHOWERR | WRDE_UNDEF) != 0 || p.we_wordc == 0) {
    throw gromos::Exception("wordexp", "Cannot expand \"" + arg + "\"\n");
  }
  string str("");
  for(size_t i = 0; i < p.we_wordc; ++i)
    str += string(p.we_wordv[i]);
  wordfree(&p);
  return str;
#else
  return arg;
#endif
}

