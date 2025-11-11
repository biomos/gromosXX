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
 * @file cos_dipole.cc
 * Calculate molecular dipoles including COS charges
 */

/**
 * @page programs Program Documentation
 *
 * @anchor cos_dipole
 * @section cos_dipole Calculate molecular dipoles including COS charges
 * @author @ref sb @ref mp
 * @date 11-01-2018
 *
 * Program cos_dipole calculates the average dipole moments over a selected set of molecules, taking into account also polarizable sites.
 * Standardly it outputs the magnitude of the average total, fixed and induced molecular dipoles, but if needed, the x-, y- and z-components can be written to a file by specifying the \@xyz flag.
 * 
 * Note that the total dipole moment is only well-defined for 
 * systems consisting of neutral molecules. If the system carries a net-charge,
 * the dipole moment will depend on the position of the origin. In cases where
 * the overall system is neutral but contains ions, the dipole moment will
 * depend on which periodic copy of the ions is taken. In these cases, the program
 * issues a warning that results will critically depend on the choice of
 * gathering method.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@fac</td><td>&lt;conversion factor for the unit of the dipole, default: 1; use 48.032045 to convert from e*nm to Debye&gt;] </td></tr>
 * <tr><td> [\@molecules</td><td>&lt;solute molecules to average over, e.g. 1-5,37,100-101&gt;] </td></tr>
 * <tr><td> [\@solv</td><td>&lt;include solvent&gt;]</td></tr>
 * <tr><td> [\@xyz</td><td>&lt;filename for writing out dipole x-,y-,z-components, Default: Mxyz.out&gt;]</td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@trs</td><td>&lt;special trajectory files with COS displacements&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  cos_dipole
    @topo  ex.top
    @pbc   r
    @fac 48.032045
    @molecules 1-5,8-10
    @solv
    [@time  0 0.2]
    @traj  ex.trc
    @trs   ex.trs
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/utils/parse.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

void get_atom_dipole(const AtomTopology &current_atom, Vec &atom_pos,
                     Vec &offsite_i, Vec &offsite_j, Vec &cosdispl,
                     int &pol_count,
                     Vec &mol_dip, Vec &fix_dip, Vec &ind_dip,
                     double &nchargepa)
{
  Vec pol_site_pos(0, 0, 0);
  double pol_site_q = 0;
  if (current_atom.isPolarisable())
  {
    pol_count++;
    if (current_atom.poloffsiteGamma() > 0.000001)
    {
      double off_set_gamma;
      off_set_gamma = current_atom.poloffsiteGamma();

      //position of the offset atom
      pol_site_pos = atom_pos + off_set_gamma *
                                    (offsite_i + offsite_j - 2 * atom_pos);
    }
    else
    {
      pol_site_pos = atom_pos;
    }
    pol_site_q = current_atom.charge() + (current_atom.cosCharge() * -1);
    Vec cos_pos = pol_site_pos + cosdispl;

    // contribution of the cos charges with cos position
    mol_dip += cos_pos * current_atom.cosCharge();
    ind_dip += cosdispl * current_atom.cosCharge();
  }
  else
  {
    // case no polarisation nor offset
    pol_site_pos = atom_pos;
    pol_site_q = current_atom.charge();
  }
  mol_dip += (pol_site_q - nchargepa) * pol_site_pos;
  fix_dip += (current_atom.charge() - nchargepa) * pol_site_pos;
}

int main(int argc, char **argv)
{

  Argument_List knowns;
  knowns << "topo"
         << "pbc"
         << "time"
         << "molecules"
         << "fac"
         << "solv"
         << "xyz"
         << "traj"
         << "trs";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t[@fac   <conversion factor for the unit of the dipole, default: 1; use 48.032045 to convert from e*nm to Debye>]\n";
  usage += "\t[@molecules <solute molecules to average over, e.g. 1-5,37,100-101>]\n";
  usage += "\t[@solv <include solvent>]\n";
  usage += "\t[@xyz <filename for writing out dipole x-,y-,z-components, Default: Mxyz.out>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@trs    <special traj with cosdisplacement>\n";

  try
  {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time(args), time_trs(args);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // include solvent
    bool include_solvent = false;
    if (args.count("solv") > -1)
      include_solvent = true;

    // write x-, y- z-components?
    std::string fname_xyz;
    bool write_xyz = false;
    if (args.count("xyz") == 0)
    {
      fname_xyz = "Mxyz.out";
      write_xyz = true;
    }
    else if (args.count("xyz") > 0)
    {
      fname_xyz = args.getValue<std::string>("xyz", false, "Mxyz.out");
      write_xyz = true;
    }

    // create an atomspecifier that contains all atoms (including the solvent)
    // or just the ones that got specified by @molecules
    utils::AtomSpecifier atoms(sys);

    for (int m = 0; m < sys.numMolecules(); ++m)
      for (int a = 0; a < sys.mol(m).numAtoms(); ++a)
        atoms.addAtom(m, a);

    //determine net charge
    double ncharge = 0;
    double nchargepa = 0;
    for (unsigned int i = 0; i < atoms.size(); i++)
    {
      ncharge += atoms.charge(i);
    }
    // atoms.size() does not include the solvent, so for solvent-only systems:
    int num_solv_mol = 0;
    bool charged_solvent = false;
    if (include_solvent)
    {
      // read first frame to get number of solvent atoms
      InG96 ic;
      if (args.count("traj") > 0)
      {
        ic.open(args.lower_bound("traj")->second);

        ic.select("ALL");
        ic >> sys;
        ic.close();
      }
      else
      {
        throw gromos::Exception("dipole", "no trajectory specified (@traj)");
      }
      atoms.addSpecifier("s:a");
      for (int i = 0; i < sys.numSolvents(); i++)
      {
        double solv_charge_mol = 0;
        for (int a = 0; a < sys.sol(i).topology().numAtoms(); a++)
          solv_charge_mol += sys.sol(i).topology().atom(a).charge();
        if (solv_charge_mol)
          charged_solvent = true;

        int _num_solv_mol = sys.sol(i).numAtoms() / sys.sol(i).topology().numAtoms();
        num_solv_mol += _num_solv_mol;
        ncharge += num_solv_mol * solv_charge_mol;
      }
    }

    if (atoms.size() + atoms.numSolventAtoms())
      nchargepa = ncharge / (atoms.size() + atoms.numSolventAtoms());
    else
      nchargepa = 0;

    int num_mol = 0;
    std::vector<int> molecules;
    std::string mol_parse_string;
    if (args.count("molecules") == -1)
    {
      mol_parse_string = "all";
      for (int i = 0; i < sys.numMolecules(); i++)
      {
        molecules.push_back(i);
      }
      num_mol = sys.numMolecules();
      cout << "# average dipole and quadrupole moments over all " << num_mol << " solute molecules" << endl;
    }
    else
    {

      mol_parse_string = args.getValue<std::string>("molecules", false, "1-2");
      utils::parse_range(mol_parse_string, molecules);
      num_mol = molecules.size();

      for (unsigned int i = 0; i < molecules.size(); i++)
      {
        if (molecules[i] > sys.numMolecules())
        {
          std::ostringstream oss;
          oss << "flag @molecules: " << mol_parse_string << " , but there are only " << sys.numMolecules() << " molecules in the system";
          throw gromos::Exception("cos_dipole", oss.str());
        }
      }
      cout << "# average dipole and quadrupole moments of molecules " << mol_parse_string << endl;
    }

    if (include_solvent)
      cout << "# including " << num_solv_mol << " solvent molecules " << std::endl;

    // read conversion factor
    double conversion_fac = args.getValue<double>("fac", false, 1);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // write a title
    cout << "#\n";
    if (ncharge != 0.0)
    {
      cout << "# WARNING the system carries a net charge ( "
           << ncharge << ")\n"
           << "#         this means that the dipole depends on the origin!\n"
           << "#\n";
    }
    else
    {
      // we also want to know if there are any ions
      int ion_count = 0;
      for (int m = 0; m < sys.numMolecules(); ++m)
      {
        double tc = 0.0;
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a)
          tc += sys.mol(m).topology().atom(a).charge();
        if (tc != 0.0)
          ion_count++;
      }
      if (ion_count != 0 || charged_solvent)
      {
        cout << "# WARNING although the system is neutral overall,\n";
        if (ion_count != 0)
          cout << "#        - there are" << ion_count << " molecules that carry a charge\n";
        if (charged_solvent)
          cout << "#        - there is charged solvent\n";
        cout << "#\n"
             << "# -> the system-dipole will depend on their positions in "
             << "the periodic box\n"
             << "#         this is likely to give random results\n";
      }
    }

    cout << "#" << setw(14) << " time"
         << setw(15) << " dipole_tot"
         << setw(15) << " dipole_fix"
         << setw(15) << " dipole_induced"
         << setw(15) << " Q_xx"
         << setw(15) << " Q_yy"
         << setw(15) << " Q_zz\n";

    ofstream os(fname_xyz.c_str());
    if (write_xyz)
    {
      os << "# x-, y- and z- components of the total, fixed and reduced dipole moments\n";
      os << "# averages over the selected solute molecules: " << mol_parse_string << "\n";
      if (include_solvent)
        os << "# including " << num_solv_mol << " solvent molecules " << std::endl;
      os << "#\n";
      os << "#" << setw(9) << "Time"
         << setw(12) << "Mtot_x" << setw(12) << "Mtot_y" << setw(12) << "Mtot_z"
         << setw(12) << "Mfix_x" << setw(12) << "Mfix_y" << setw(12) << "Mfix_z"
         << setw(12) << "Mind_x" << setw(12) << "Mind_y" << setw(12) << "Mind_z"
         << std::endl;
    }

    int numFrames = 0;

    // sum over the timeseries
    double tot_sum_mol_dip = 0, tot_sum_fix_dip = 0, tot_sum_ind_dip_frame = 0;

    vector<Vec> dipole_vector;

    // define input coordinate
    InG96 ic, is;

    // loop over all trajectories
    for (Arguments::const_iterator
             iter = args.lower_bound("traj"),
             to = args.upper_bound("traj"),
             siter = args.lower_bound("trs"),
             sto = args.upper_bound("trs");
         iter != to || siter != sto; ++iter, ++siter)
    {

      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");

      //open special traj
      is.open((siter->second).c_str());
      is.select("ALL");

      // loop over single trajectory
      while (!ic.eof() && !is.eof())
      {
        is >> sys >> time_trs;

        ic >> sys >> time;

        // check if the times are the same
        if (time_trs.time() != time.time())
        {
          cerr << "Time is not corresponding between trajectories\n"
               << "Special traj:\t" << time_trs.time() << "\n"
               << "Position traj:\t" << time.time() << endl;
          return -1;
        }

        // we have to reconnect the molecules
        // in the case of neutral molecules, this will be sufficient
        // non-neutral molecules will behave weirdly.
        (*pbc.*gathmethod)();

        Vec tot_dipole(0.0, 0.0, 0.0), fix_dipole(0.0, 0.0, 0.0), ind_dipole(0.0, 0.0, 0.0);
        double sum_mol_dip_frame = 0;
        double sum_fix_dip_frame = 0;
        double sum_ind_dip_frame = 0;

        int pol_count = 0;

//calculate molecular dipoles for selected solute molecules
        for (unsigned int m = 0; m < molecules.size(); m++)
        {
          Molecule molecule = sys.mol(molecules[m] - 1);

          Vec mol_com(0, 0, 0);
          double mass = 0;
          for (int a = 0; a < molecule.numAtoms(); a++)
          {
            mol_com += molecule.topology().atom(a).mass() * molecule.pos(a);
            mass += molecule.topology().atom(a).mass();
          }
          mol_com /= mass;
          //cerr << "mass " << mass << " com " <<  v2s(mol_com) << endl;

          //cout << "id" << omp_get_thread_num() << " mol" << m << " " <<std::endl;
          Vec mol_dip(0, 0, 0);
          Vec fix_dip(0, 0, 0);
          Vec ind_dip(0, 0, 0);

          Vec offsite_i(0, 0, 0), offsite_j(0, 0, 0);

          //cerr << "#Amount of atoms in molecule:\t" << sys.mol(m).numAtoms() << endl;
          for (int a = 0; a < molecule.numAtoms(); a++)
          {
            const AtomTopology &current_atom = molecule.topology().atom(a);

            Vec atom_pos = molecule.pos(a) - mol_com;
            if (current_atom.isPolarisable() && current_atom.poloffsiteGamma() > 0.000001)
            {
              offsite_i = molecule.pos(current_atom.poloffsiteI()) - mol_com;
              offsite_j = molecule.pos(current_atom.poloffsiteJ()) - mol_com;
            }
            Vec cosdispl = molecule.cosDisplacement(a);

            get_atom_dipole(current_atom, atom_pos, offsite_i, offsite_j, cosdispl, pol_count, mol_dip, fix_dip, ind_dip, nchargepa);
          }
          //cerr << m << ":\t" << mol_dip.abs() << endl;
          tot_dipole += mol_dip;
          fix_dipole += fix_dip;
          ind_dipole += ind_dip;
          sum_mol_dip_frame += mol_dip.abs();
          sum_fix_dip_frame += fix_dip.abs();
          sum_ind_dip_frame += ind_dip.abs();
        }

        // calculate molecular dipoles for solvent molecules
        if (include_solvent)
        {
          for (int s = 0; s < sys.numSolvents(); s++)
          {
            int num_solvent_atoms = sys.sol(s).topology().numAtoms();
            for (int i = 0; i < sys.sol(s).numAtoms(); i += num_solvent_atoms)
            {
              Vec mol_dip(0, 0, 0);
              Vec fix_dip(0, 0, 0);
              Vec ind_dip(0, 0, 0);

              Vec offsite_i(0, 0, 0), offsite_j(0, 0, 0);

              //cerr << "#Amount of atoms in molecule:\t" << sys.mol(m).numAtoms() << endl;
              for (int a = 0; a < num_solvent_atoms; a++)
              {
                const AtomTopology &current_atom = sys.sol(s).topology().atom(a);

                Vec atom_pos = sys.sol(s).pos(i + a);

                if (current_atom.isPolarisable() && current_atom.poloffsiteGamma() > 0.000001)
                {
                  offsite_i = sys.sol(s).pos(i + a + current_atom.poloffsiteI());
                  offsite_j = sys.sol(s).pos(i + a + current_atom.poloffsiteJ());
                }
                Vec cosdispl = sys.sol(s).cosDisplacement(i + a);

                get_atom_dipole(current_atom, atom_pos, offsite_i, offsite_j, cosdispl, pol_count, mol_dip, fix_dip, ind_dip, nchargepa);
              }
              //cerr << m << ":\t" << mol_dip.abs() << endl;
              tot_dipole += mol_dip;
              fix_dipole += fix_dip;
              ind_dipole += ind_dip;
              sum_mol_dip_frame += mol_dip.abs();
              sum_fix_dip_frame += fix_dip.abs();
              sum_ind_dip_frame += ind_dip.abs();
            }
          }
        }

        int tot_num_mol = num_mol + num_solv_mol;
        tot_sum_mol_dip += sum_mol_dip_frame / tot_num_mol;
        tot_sum_fix_dip += sum_fix_dip_frame / tot_num_mol;

        double ave_ind_dip_frame;
        if (pol_count != 0)
        {
          tot_sum_ind_dip_frame += sum_ind_dip_frame / pol_count;
          ave_ind_dip_frame = sum_ind_dip_frame * conversion_fac / pol_count;
        }
        else
        {
          ave_ind_dip_frame = 0;
        }
        // do some bookkeeping
        numFrames++;

        cout << time
             << setw(15) << setprecision(8) << sum_mol_dip_frame * conversion_fac / tot_num_mol
             << setw(15) << setprecision(8) << sum_fix_dip_frame * conversion_fac / tot_num_mol
             << setw(15) << setprecision(8) << ave_ind_dip_frame
             << endl;
        if (write_xyz)
        {
          os << setw(10) << setprecision(5) << time.time()
             << setw(12) << tot_dipole[0] * conversion_fac / tot_num_mol << setw(12) << tot_dipole[1] * conversion_fac / tot_num_mol << setw(12) << tot_dipole[2] * conversion_fac / tot_num_mol
             << setw(12) << fix_dipole[0] * conversion_fac / tot_num_mol << setw(12) << fix_dipole[1] * conversion_fac / tot_num_mol << setw(12) << fix_dipole[2] * conversion_fac / tot_num_mol
             << setw(12) << ind_dipole[0] * conversion_fac / tot_num_mol << setw(12) << ind_dipole[1] * conversion_fac / tot_num_mol << setw(12) << ind_dipole[2] * conversion_fac / tot_num_mol
             << std::endl;
        }
      }
      ic.close();
      is.close();
    }

    delete pbc;
    os.close();

    //now print average over the systems for all three values
    cout << "#\n#" << setw(14) << "Averages:"
         << setw(15) << setprecision(8) << tot_sum_mol_dip * conversion_fac / numFrames
         << setw(15) << setprecision(8) << tot_sum_fix_dip * conversion_fac / numFrames
         << setw(15) << setprecision(8) << tot_sum_ind_dip_frame * conversion_fac / numFrames
         << endl;
  }

  catch (const gromos::Exception &e)
  {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
