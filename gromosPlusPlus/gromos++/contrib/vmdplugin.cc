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
 * @file vmdplugin.cc
 * A plugin for VMD to read GROMOS files.
 *
 * @page programs Program Documentation
 *
 * @anchor vmdplugin
 * @section vmdplugin Read GROMOS files in VMD
 * @author @ref ns
 * @date 08.11.2010
 *
 * This program is not an individual program but a plugin library which
 * runs in the VMD (Visual Molecular Dynamics) program. It is used to
 * open GROMOS configuration files directly in VMD.
 *
 * Once a file is opened using one of the GROMOS plugins VMD will prompt for the
 * arguments in the VMD console. Because the topological data in GROMOS is
 * seperated from the configurational data the plugin can only roughly guess
 * the data from a configuration file. For this reason it is recomended to give
 * information about the topology and periodic boundary conditions using the arguments.
 * The arguments can also be read from a file (\@f).
 * Gathering and rotational fitting (to a reference structure or the first structure
 * in the trajectory) can be carried out directly in the plugin.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> [\@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; ]</td></tr>
 * <tr><td> [\@include</td><td>&lt;SOLUTE (default), SOLVENT or ALL&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference structure to fit to&gt;] </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to fit to&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;time and dt&gt;] </td></tr>
 * <tr><td> [\@factor</td><td>&lt;factor to convert length unit to Angstrom, 10.0&gt;]</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
    @topo        ex.top
    @pbc         r
 @endverbatim
 *
 * <hr>
 */

#include "../config.h"
#ifdef HAVE_VMD
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <cassert>
#include <cstring>
#include <algorithm>

#include "../gcore/System.h"
#include "../gcore/Box.h"
#include "../args/Arguments.h"

#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Constraint.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/utils/groTime.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/CommandLine.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/utils/AtomSpecifier.h"
#include "molfile_plugin.h"

//#define DEBUG(x) std::cerr << x << std::endl
#define DEBUG(X)

using namespace std;
using namespace gio;
using namespace gcore;
using namespace utils;
using namespace gmath;
using namespace args;
using namespace bound;
using namespace fit;

#define COLOR_RED "\033[1;31m"
#define COLOR_RESET "\033[22;0m"

// this enum is used to determine the plugin type
enum plugin_type {
  ptUnknown, // Don't know what type
  ptGuess, // the system was guessed based on POSITION block
  ptTop // the system was created from a topology
};

// this struct is passed by VMD to all the functions
struct plugin_state {
  // create an empty set of pointers
  plugin_state() : pt(ptUnknown), sys(NULL), refSys(NULL), incrd(NULL), 
  time(NULL), pbc(NULL), factor(10.0), fitatoms(NULL), rf(NULL), bond_from(NULL), bond_to(NULL) {
  }

  // free memory if this is required.
  ~plugin_state() {
    if (sys != NULL)
      delete sys;
    if (refSys != NULL)
      delete refSys;
    if (incrd != NULL) {
      incrd->close();
      delete incrd;
    }
    if (time != NULL)
      delete time;
    if (pbc != NULL)
      delete pbc;
    if (fitatoms != NULL)
      delete fitatoms;
    if (rf != NULL) {
      delete rf->getReference();
      delete rf;
    }
    if (bond_from != NULL)
      delete[] bond_from;
    if (bond_to != NULL)
      delete[] bond_to;
  }

  plugin_type pt;                   // type of the plugin
  System * sys;                     // the system, either guessed or read
  System * refSys;                  // the reference system
  InG96 * incrd;                    // the coordinate file
  Time * time;                      // the time, read, arged or guessed
  Boundary *pbc;                    // boundary type
  Boundary::MemPtr gathmethod;      // gathering method
  double factor;                    // factor to convert nm to A.
  AtomSpecifier * fitatoms;         // atoms for fitting
  RotationalFit * rf;               // the fit
  int * bond_from;                  // the bonds array (from)
  int * bond_to;                    // the bonds array (end)
};

#ifdef  __cplusplus
extern "C" {
#endif

  // this function is called by VMD at the very beginning when the file is
  // opened. We prompt for the arguments and create all the objects
  // we need and save them in the plugin state.
  static void * open_read(const char *filename, const char *, int *natoms) {
    plugin_state * s = new plugin_state;
    *natoms = 0;

    string select = "SOLUTE";

    // this is to avoid exception to be thrown at VMD and crashing the program.
    try {
      // try to do it GROMOS++ style. Guessing we do if this fails...
      try {
        Argument_List knowns;
        knowns << "topo" << "atomsfit" << "pbc" << "ref" << "time" << "include" << "factor";
 
        cerr << COLOR_RED << "The GROMOS VMD plugin requires several command line arguments to proceed because" << endl
             << "the configuration file do not contain enough topological data." << endl
             << "Just hit enter to let the plugin guess the parameters from the configuration file." << COLOR_RESET << endl;

        string usage = "Options for the GROMOS VMD plugin:";
        usage += "\n\t@topo       <molecular topology file>\n";
        usage += "\t[@pbc        <boundary type> [<gathermethod>]]\n";
        usage += "\t[@time       <time and dt>]\n";
        usage += "\t[@atomsfit   <atoms to consider for fit>]\n";
        usage += "\t[@ref        <reference coordinates (if absent, the first frame of @traj is reference)>]\n";
        usage += "\t[@include    <SOLUTE, SOLVENT, ALL, default: SOLUTE)>]\n";
        usage += "\t[@factor     <factor to convert length unit, default 10.0>]\n";
        cerr << endl << usage;

        // this is nasty. We have to convert the string of arguments to an
        // argv C-style array.
        // This is done via a vector of strings.
        string cmdargs = CommandLine::getLine("Arguments: ", cerr);
        istringstream is(cmdargs);
        vector<string> strargs;
        strargs.push_back("vmdplugin");
        while(!is.eof()) {
          string s;
          is >> s;
          strargs.push_back(s);
        }

        // no convert vector<string> to char[][]
        char** argv = new char*[strargs.size()];
        DEBUG("creating argv");
        for(unsigned int i = 0; i < strargs.size(); ++i) {
          argv[i] = new char[strargs[i].size()+1];
          strcpy(argv[i], strargs[i].c_str());
          DEBUG("arg " << i << ": " << argv[i]);
        }

        DEBUG("parsing args");
        Arguments args(strargs.size(), (char**)argv, knowns, usage);

        // free argv again
        for(unsigned int i = 0; i < strargs.size(); ++i) {
          delete[] argv[i];
        }
        delete[] argv;

        s->pt = ptTop;
        DEBUG("parsing topo");
        InTopology it(args["topo"]);
        s->sys = new System(it.system());
        s->refSys = new System(it.system());
        InG96 ic(filename);

        if (args.count("include") > 0) {
          string sel = args["include"];
          transform(sel.begin(), sel.end(), sel.begin(), static_cast<int (*)(int)> (std::toupper));
          if (sel != "SOLUTE" && sel != "SOLVENT" && sel != "ALL")
            throw gromos::Exception("vmdplugin", "@include should be SOLUTE, SOLVENT or ALL.");
          select = sel;
        }

        // read the frame. This is just done to get the solvent. We need to know
        // the actual number of atoms here.
        ic.select(select);
        ic >> *(s->sys);
        ic.close();

        if (args.count("ref") > 0) {
          ic.open(args["ref"]);
        } else {
          ic.open(filename);
        }
        ic.select(select);
        ic >> *(s->refSys);
        ic.close();

        DEBUG("parsing gathering");
        s->pbc = BoundaryParser::boundary(*(s->sys), args);
        // GatherParser
        s->gathmethod = args::GatherParser::parse(*(s->sys), *(s->refSys), args);

        // factor
        if (args.count("factor") > 0) {
          istringstream is(args["factor"]);
          if (!(is >> s->factor)) {
            s->factor = 10.0;
            throw gromos::Exception("vmdplugin", "The factor has to be numeric");
          }
        }

        // fitatoms
        s->fitatoms = new AtomSpecifier(*(s->refSys));
        if(args.count("atomsfit") > 0){
          Arguments::const_iterator iter = args.lower_bound("atomsfit");
          Arguments::const_iterator to = args.upper_bound("atomsfit");

          for(;iter!=to;iter++){
            string spec=iter->second.c_str();
            s->fitatoms->addSpecifier(spec);
          }
        }
        if (s->fitatoms->size()) {
          Reference * ref = new Reference(s->refSys);
          ref->addAtomSpecifier(*(s->fitatoms));
          s->rf = new RotationalFit(ref);
        }

        // time
        s->time = new Time(args);

        // compute the total number of atoms.
        for (int m = 0; m < s->sys->numMolecules(); ++m)
          *natoms += s->sys->mol(m).numAtoms();

        for (int sol = 0; sol < s->sys->numSolvents(); ++sol)
          *natoms += s->sys->sol(sol).numPos();
      } catch (const gromos::Exception & exp) {
        cerr << "Error: " << exp.what() << endl;
        cerr << COLOR_RED << "No topology provided. Trying to guess..."  << COLOR_RESET << endl;
        s->pt = ptGuess;

        // just open the file and look for a POSITION block.
        DEBUG("opening file" << filename);
        Ginstream file(filename);
        s->sys = new System;

        vector<string> block;
        bool found = false;
        while (file.getblock(block)) {
          if (block[0] == "POSITION") {
            if (block[block.size() - 1] != "END") {
              throw gromos::Exception("vmdplugin", "POSITION block is corrupt. No END.");
            }

            // the number of atoms is the size of the content of the position block.
            *natoms = block.size() - 2;
            DEBUG("Number of lines:" << *natoms);

            // create the topology based on the information in the block.
            // this is very rough. There is no way to really identify
            // new molecules. So the whole solute and solvent will be one
            // molecule each.
            MoleculeTopology mt;
            int currentRes = -1;
            for (int i = 1; i <= *natoms; ++i) {
              istringstream line(block[i]);
              int resnum, atomnr;
              string atomname;
              string resname;
              double dummy;
              line >> resnum >> resname >> atomname >> atomnr >> dummy >> dummy >> dummy;
              if (line.fail()) {
                throw gromos::Exception("vmdplugin", "Bad line in POSITION block.");
              }

              if (currentRes == -1) {
                currentRes = resnum;
              }

              if (currentRes != resnum) {
                // a new molecule started. This means that the residue number
                // starts from 1 again.
                if (resnum < currentRes) {
                  s->sys->addMolecule(Molecule(mt));
                  mt = MoleculeTopology();
                }
                mt.setResName(resnum - 1, resname);
                currentRes = resnum;
              }

              // add the atom
              AtomTopology atom;
              atom.setName(atomname);
              mt.addAtom(atom);
              mt.setResNum(mt.numAtoms() - 1, resnum - 1);
            }
            s->sys->addMolecule(Molecule(mt));
            found = true;
            break;
          }
        }

        // if there was no POSITION block (e.g. only a POSITIONRED block) we abort.
        if (!found)
          throw gromos::Exception("vmdplugin", "Cannot read system data from file.");
        s->time = new Time;
      }

      // Here the file is finally opened as a coordinate stream. The frames
      // are not read yet.
      s->incrd = new InG96(filename);
      s->incrd->select(select);
    } catch (const gromos::Exception &e) {
      // if everything fails...
      cerr << COLOR_RED << e.what() << COLOR_RESET  << endl;
      return NULL;
    }
    return (void*) s;
  }

  // this is called by VMD at the very end to clean up. The main work
  // is carried out by the destructor of plugin_state.
  static void close_read(void * data) {
    DEBUG("finishing");
    if (data == NULL)
      return;
    plugin_state * s = static_cast<plugin_state *> (data);
    delete s;
  }

  // this is called by VMD to read the structure
  // the name of the function is slighlty confusing. In fact the
  // topology is read. We can get the toplogocal data from the system.
  //
  // resid: GROMOS residue number (per molecule)
  // chain: A-Z for solvent molecules, 0 for solvent.
  //
  // so far we keep bfactor and occ empty. This could be used to display
  // other atomic properties (i.e. RMSF read from a file)
  static int read_structure(void * data, int *optflags, molfile_atom_t *atoms) {
    DEBUG("reading structure");
    plugin_state * s = static_cast<plugin_state *> (data);
    System & sys = *(s->sys);

    molfile_atom_t * atom = atoms;

    *optflags = MOLFILE_NOOPTIONS;

    if (s->pt == ptTop) {
      // if a toplogy was read we get some extra nice stuff.
      *optflags = MOLFILE_MASS | MOLFILE_CHARGE;
    }

    // copy the solute toplogy to VMD
    for (int m = 0; m < sys.numMolecules(); ++m) {
      for (int a = 0; a < sys.mol(m).numAtoms(); ++a, ++atom) {
        //DEBUG("mol: " << m << " atom: " << a);
        strcpy(atom->name, sys.mol(m).topology().atom(a).name().c_str());
        //DEBUG("name: " << atom->name);
        strcpy(atom->type, sys.mol(m).topology().atom(a).name().c_str());
        atom->resid = sys.mol(m).topology().resNum(a) + 1;
        strcpy(atom->resname, sys.mol(m).topology().resName(atom->resid-1).c_str());
        //DEBUG("resname: " << atom->resname);

        atom->chain[0] = 'A' + char(m);
        if (atom->chain[0] < 'A' || atom->chain[0] > 'Z')
          atom->chain[0] = 'Z';

        atom->chain[1] = '\0';
        atom->segid[0] = '\0';

        if (s->pt == ptTop) {
          atom->mass = sys.mol(m).topology().atom(a).mass();
          atom->charge = sys.mol(m).topology().atom(a).charge();
        }
      }
    }

    // copy the solvent to VMD
    int resid = 0;
    for (int sol = 0; sol < sys.numSolvents(); ++sol) {
      for(int a = 0; a < sys.sol(sol).numPos(); ++a, ++atom) {
        int index = a % sys.sol(sol).topology().numAtoms();
        if (index == 0) resid++;
        strcpy(atom->name, sys.sol(sol).topology().atom(index).name().c_str());
        //DEBUG("name: " << atom->name);
        strcpy(atom->type, sys.sol(sol).topology().atom(index).name().c_str());
        atom->resid = resid;
        strcpy(atom->resname, "SOLV");
        //DEBUG("resname: " << atom->resname);
        atom->chain[0] = '0';
        atom->chain[1] = '\0';
        atom->segid[0] = '\0';

        if (s->pt == ptTop) {
          atom->mass = sys.sol(sol).topology().atom(index).mass();
          atom->charge = sys.sol(sol).topology().atom(index).charge();
        }
      }
    }
    return MOLFILE_SUCCESS;
  }

  // This is called by VMD as soon as the structure is loaded.
  // it is used to determine the bonds. We do a simplistic approach and
  // just give the bonds but nothing else.
  static int read_bonds(void *data, int *nbonds, int **fromptr, int **toptr,
                        float ** bondorder,int **bondtype,
                        int *nbondtypes, char ***bondtypename) {

    *nbonds = 0;
    *fromptr = NULL;
    *toptr = NULL;
    *bondorder = NULL; 
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    plugin_state * s = static_cast<plugin_state *> (data);
    if (s->pt != ptTop)
      return MOLFILE_SUCCESS;

    const System & sys = *(s->sys);

    // first get the bonds for solute and solvent and save them in this vector.
    vector<pair<int,int> > bonds;
    int offatom = 1;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      for (BondIterator bit(sys.mol(i).topology()); bit; ++bit) {
        bonds.push_back(pair<int,int>(bit()[0] + offatom, bit()[1] + offatom));
      }
      offatom += sys.mol(i).numAtoms();
    }
    // solvent
    for(int i = 0; i < sys.numSolvents(); ++i) {  
      int num_atoms = sys.sol(i).topology().numAtoms();
      int num_solv_mol = sys.sol(i).numPos() / num_atoms;
      for(int s = 0; s < num_solv_mol; ++s) {
        for (ConstraintIterator bit(sys.sol(i).topology()); bit; ++bit) {
          bonds.push_back(pair<int,int>(bit()[0] + offatom, bit()[1] + offatom));
        }
        offatom += num_atoms;
      }
    }

    // allocate space for the bond arrays
    s->bond_from = new int[bonds.size()];
    s->bond_to = new int[bonds.size()];
    // copy the data
    for(unsigned int i = 0; i < bonds.size(); ++i) {
      s->bond_from[i] = bonds[i].first;
      s->bond_to[i] = bonds[i].second;
    }

    // tell VMD about the bonds.
    *nbonds = bonds.size();
    *fromptr = s->bond_from;
    *toptr = s->bond_to;

    return MOLFILE_SUCCESS;
  }

  // this is called by VMD to read a structure.
  // so it does not only read the TIMESTEP block but the whole configuration
  // at a given step.
  static int read_timestep(void * data, int natoms, molfile_timestep_t *ts) {
    DEBUG("reading a frame");
    plugin_state * s = static_cast<plugin_state *> (data);

    InG96 & incrd = *(s->incrd);
    System & sys = *(s->sys);
    Time & time = *(s->time);

    // don't throw exceptions at VMD!
    try {
      if (incrd.eof()) {
        DEBUG("finished. EOF.");
        return MOLFILE_EOF;
      }

      // read a frame and a time.
      DEBUG("reading system.");
      incrd >> sys;
      try {
        DEBUG("reading time.");
        incrd >> time;
      } catch (const gromos::Exception & e) {
        time.time() += time.dt();
      }

      // VMD signals to skip the frame by setting ts to 0.
      if (ts) {
        if (s->pt == ptTop) {
          // try to gather and fit
          DEBUG("Gathering and fitting");
          try {
            (*(s->pbc).*(s->gathmethod))();
            if (s->rf != NULL)
              s->rf->fit(s->sys);
          } catch (const gromos::Exception & exp) {
            cout << exp.what() << endl;
          }
        }

        DEBUG("Setting parameters");
        // set the time and the box
        ts->physical_time = time.time();
        if (sys.hasBox && sys.box().ntb() != Box::vacuum) {
          DEBUG("Box");
          ts->A = sys.box().K().abs()*s->factor;
          ts->B = sys.box().L().abs()*s->factor;
          ts->C = sys.box().M().abs()*s->factor;
          ts->alpha = sys.box().alpha();
          ts->beta = sys.box().beta();
          ts->gamma = sys.box().gamma();
        }

        // copy the positions to VMDs arrays.
        DEBUG("positions");
        float * crd = ts->coords;
        int num_pos = 0;
        for (int m = 0; m < sys.numMolecules(); ++m) {
          for (int a = 0; a < sys.mol(m).numAtoms(); ++a, ++num_pos) {
            if (sys.hasPos) {
              *(crd++) = float(sys.mol(m).pos(a)[0] * s->factor);
              *(crd++) = float(sys.mol(m).pos(a)[1] * s->factor);
              *(crd++) = float(sys.mol(m).pos(a)[2] * s->factor);
            }
          } // for atoms
        } // for mol
        DEBUG("sol: " << sys.sol(0).numPos());
        for (int m = 0; m < sys.numSolvents(); ++m) {
          for (int a = 0; a < sys.sol(m).numPos(); ++a, ++num_pos) {
            // here we have to make sure that we don't read more atoms than
            // contained in the first frame as VMD only has space
            // for that many ayomd;
            if (num_pos >= natoms) 
              break;
            
            if (sys.hasPos) {
              *(crd++) = float(sys.sol(m).pos(a)[0] * s->factor);
              *(crd++) = float(sys.sol(m).pos(a)[1] * s->factor);
              *(crd++) = float(sys.sol(m).pos(a)[2] * s->factor);
            }
          } // for atoms
        } // for mol
      } // if ts

    } catch (const gromos::Exception &e) {
      cout << e.what() << endl;
      return MOLFILE_ERROR;
    }
    return MOLFILE_SUCCESS;

  }
#ifdef GROMOS_MAJOR_VERSION
  static const int maj_vers = GROMOS_MAJOR_VERSION;
  static const int min_vers = GROMOS_MINOR_VERSION;
#else
  static const int maj_vers = 0;
  static const int min_vers = 0;
#endif


  // the GROMOS plugins are all identical but
  // different in the file extensions they register for.
  static molfile_plugin_t gromos_cnf_plugin = {
    vmdplugin_ABIVERSION, // ABI version
    MOLFILE_PLUGIN_TYPE, // type of plugin
    "grocnf", // short name of plugin
    "GROMOS configuration", // pretty name of plugin
    "Nathan Schmid",
    (maj_vers), // major version
    (min_vers), // minor version
    VMDPLUGIN_THREADUNSAFE, // is not reentrant
    "cnf", // filename extension
    open_read,
    read_structure,
    read_bonds,
    read_timestep,
    close_read,
    0, // open_write
    0, // write_structure
    0, // write_timestep
    0, // close_write
    0, // read_volumetric_metadata
    0, // read_volumetric_data
    0 // read_rawgraphics
  };
  static molfile_plugin_t gromos_trc_plugin = {
    vmdplugin_ABIVERSION, // ABI version
    MOLFILE_PLUGIN_TYPE, // type of plugin
    "grotrc", // short name of plugin
    "GROMOS trajectory", // pretty name of plugin
    "Nathan Schmid",
    (maj_vers), // major version
    (min_vers), // minor version
    VMDPLUGIN_THREADUNSAFE, // is not reentrant
    "trc", // filename extension
    open_read,
    read_structure,
    read_bonds,
    read_timestep,
    close_read,
    0, // open_write
    0, // write_structure
    0, // write_timestep
    0, // close_write
    0, // read_volumetric_metadata
    0, // read_volumetric_data
    0 // read_rawgraphics
  };
  static molfile_plugin_t gromos_gz_plugin = {
    vmdplugin_ABIVERSION, // ABI version
    MOLFILE_PLUGIN_TYPE, // type of plugin
    "grogz", // short name of plugin
    "GROMOS compressed data", // pretty name of plugin
    "Nathan Schmid",
    (maj_vers), // major version
    (min_vers), // minor version
    VMDPLUGIN_THREADUNSAFE, // is not reentrant
    "gz", // filename extension
    open_read,
    read_structure,
    read_bonds,
    read_timestep,
    close_read,
    0, // open_write
    0, // write_structure
    0, // write_timestep
    0, // close_write
    0, // read_volumetric_metadata
    0, // read_volumetric_data
    0 // read_rawgraphics
  };

  VMDPLUGIN_API int VMDPLUGIN_init() {
    return 0;
  }

  VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
    DEBUG("registering plugin...");
    // register the three plugins.
    (*cb)(v, (vmdplugin_t *) & gromos_cnf_plugin);
    (*cb)(v, (vmdplugin_t *) & gromos_trc_plugin);
    (*cb)(v, (vmdplugin_t *) & gromos_gz_plugin);
    return 0;
  }

  VMDPLUGIN_API int VMDPLUGIN_fini() {
    return 0;
  }

#ifdef  __cplusplus
}
#endif
#endif
