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
 * @file cos_epsilon.cc
 * Calculate relative permittivity and box dipole moment autocorrelations
 */

/**
 * @page programs Program Documentation
 *
 * @anchor cos_epsilon
 * @section cos_epsilon Calculate relative permittivity and box dipole moment autocorrelations
 * @author @ref sb @ref mp
 * @date 11-01-2018
 *
 * Program cos\_epsilon estimates the relative dielectric permittivity, 
 * @f$\epsilon(0)@f$, of a simulation box from a Kirkwood - Fr&ouml;hlich type of
 * equation, as derived by Neumann [Mol. Phys. 50, 841 (1983)], taking into account also polarizable centers.
 *
 * @f[ (\epsilon(0) - 1) \frac{2\epsilon_{RF}+1}{2\epsilon_{RF}+\epsilon(0)} = \frac{<\vec{M}^2> - <\vec{M}>^2}{3\epsilon_0 Vk_BT} @f]
 *
 * where @f$\vec{M}@f$ is the total dipole moment of the system, 
 * @f$\epsilon_0@f$ is the dielectric permittivity of vacuum,
 * @f$\epsilon_{RF}@f$ is a reaction-field epsilon value, @f$V@f$ is the volume
 * and @f$k_BT@f$ is the absolute temperature multiplied by the Boltzmann 
 * constant. 
 * 
 * If the \@autocorr flag is given, it also writes out the normalized autocorrelation function @f$\Theta(\tau)@f$ of the total dipole moment of the system: 
 * @f[ \Theta(\tau) = \frac{<\vec{M}(t)\cdot\vec{M}(t+\tau)>_t}{<M(t)^2>_t} = exp(-\frac{t}{\tau_\phi})@f]
 * from which the Debye relaxation time @f$\tau_D@f$ can be calculated [J. Chem. Phys. 82, 5663 (1985)]:
 * @f[ \tau_D = \frac{2\epsilon_{rf} + \epsilon(0)}{2\epsilon_{rf} + 1} \tau_{\phi}@f] 
 * Using the \@truncate flag, the autocorrelation output can be truncated when the number of independent contributing frames drops below a given threshold.
 *
 * Note that the total dipole moment of the system is only well-defined for 
 * systems consisting of neutral molecules. If the system carries a net-charge,
 * the dipole moment will depend on the position of the origin. In cases where
 * the overall system is neutral but contains ions, the dipole moment will
 * depend on which periodic copy of the ions is taken. In these cases, cos\_epsilon
 * issues a warning that results will critically depend on the choice of
 * gathering method.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> [\@e_rf</td><td>&lt;reaction field epsilon&gt;] </td></tr>
 * <tr><td> [\@fac</td><td>&lt;conversion factor for the unit of the dipole, default: 1; use 48.032045 to convert from e*nm to Debye&gt;] </td></tr>
 * <tr><td> \@temp</td><td>&lt;temperature&gt; </td></tr>
 * <tr><td> [\@autocorr</td><td>&lt;filename for storing time autocorrelation&gt;] </td></tr>
 * <tr><td> [\@truncate</td><td>&lt;minimum number of independent contributing frames after which to truncate the correlation function&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@trs</td><td>&lt;special trajectory files with COS displacements&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  cos_epsilon
    @topo  ex.top
    @pbc   r
    [@time  0 0.2]
    @e_rf  61
    @fac 48.032045
    @temp  300
    @traj  ex.trc
    @trs   ex.trs
 @endverbatim
 *
 * <hr>
 */

#include <cstdlib>
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
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"
#include "../src/gmath/Correlation.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
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
 Vec &mol_dip,
 double &nchargepa)
{
  Vec pol_site_pos(0, 0, 0);
  double pol_site_q = 0;
  if (current_atom.isPolarisable())
  {
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
  }
  else
  {
    // case no polarisation nor offset
    pol_site_pos = atom_pos;
    pol_site_q = current_atom.charge();
  }
  mol_dip += (pol_site_q - nchargepa) * pol_site_pos;
}

int main(int argc, char **argv)
{

  Argument_List knowns;
  knowns << "topo"
         << "pbc"
         << "time"
         << "e_rf"
         << "fac"
         << "temp"
         << "autocorr"
         << "truncate"
         << "traj"
         << "trs";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pbc    <boundary type> [<gather method>]\n";
  usage += "\t[@time   <time and dt>]\n";
  usage += "\t@e_rf   <reaction field epsilon>\n";
  usage += "\t[@fac   <conversion factor for the unit of the dipole, default: 1; use 48.032045 to convert from e*nm to Debye>]\n";
  usage += "\t@temp   <temperature>\n";
  usage += "\t[@autocorr    <filename for storing time autocorrelation, if no name given: Mxyz.out>]\n";
  usage += "\t[@truncate    <minimum number of independent contributing frames after which to truncate the correlation function>]\n";
  usage += "\t@traj   <trajectory files>\n";
  usage += "\t@trs    <special traj with cosdisplacement>\n";

  try
  {
    Arguments args(argc, argv, knowns, usage);

    // get the @time argument
    utils::Time time_trj(args), time_trs(args);
    std::vector<double> time;

    // read the temperature
    double temp = args.getValue<double>("temp");

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    unsigned int truncate_at = args.getValue<unsigned int>("truncate", false, 0);

    bool do_autocorr = false;
    if (args.count("autocorr") >= 0)
      do_autocorr = true;

    std::string autocorr_fname;
    if (args.count("autocorr") == 0)
    {
      autocorr_fname = "Mcorr.out";
    }
    else if (args.count("autocorr") > 0)
    {
      autocorr_fname = args.getValue<std::string>("autocorr", false, "Mcorr.out");
    }

    // create an atomspecifier that contains all atoms (including the solvent)
    utils::AtomSpecifier atoms(sys);
    for (int m = 0; m < sys.numMolecules(); ++m)
      for (int a = 0; a < sys.mol(m).numAtoms(); ++a)
        atoms.addAtom(m, a);
    atoms.addSpecifier("s:a");

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
      throw gromos::Exception("cos_epsilon", "no trajectory specified (@traj)");
    }
    
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

    if (atoms.size() + atoms.numSolventAtoms())
      nchargepa = ncharge / (atoms.size() + atoms.numSolventAtoms());
    else
      nchargepa = 0;

    // read e_rf
    double e_rf = args.getValue<double>("e_rf", true);
    // read conversion factor
    double conversion_fac = args.getValue<double>("fac", false, 1);

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // write a title
    if (ncharge != 0.0)
    {
      cout << "# WARNING the system carries a net charge ( "
           << ncharge << ")\n"
           << "#         this means that the dipole depends on the origin\n"
           << "#         and your estimate of eps is probably wrong\n"
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
         << setw(15) << " boxdipole"
         << setw(15) << " <mol.dipole>"
         << setw(15) << " epsilon\n";

    // prepare the calculation of the average volume
    double vol = 0, sum_vol = 0, vcorr = 1.0;
    Arguments::const_iterator iter = args.lower_bound("pbc");
    if (iter != args.upper_bound("pbc"))
      if (iter->second[0] == 't')
        vcorr = 0.5;

    // and values to calculate the dipole fluctuation
    Vec sum_dip(0.0, 0.0, 0.0);
    vector<Vec> dipole_vector;
    double sum_dip2 = 0.0, fluc = 0.0;
    int numFrames = 0;

    // define input coordinate
    InG96 is;

    // we also need these factors
    double fac, a, b, eps;
    double f = 3.0 * gmath::physConst.get_eps0() * gmath::physConst.get_boltzmann() * temp * (2 * e_rf + 1.0);

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

        ic >> sys >> time_trj;

        // check if the times are the same
        if (time_trs.time() != time_trj.time())
        {
          cerr << "Time is not corresponding between trajectories\n"
               << "Special traj:\t" << time_trs.time() << "\n"
               << "Position traj:\t" << time_trj.time() << endl;
          return -1;
        }

        // I store it because even when no @time flag is given I want to determine a timestep
        // as the difference between the first and second time point for the correlation
        time.push_back(time_trs.time());

        (*pbc.*gathmethod)();

        // calculate the volume
        sys.box().update_triclinic();
        vol = vcorr * sys.box().K_L_M();
        sum_vol += vol;

        Vec tot_dipole(0.0, 0.0, 0.0);
        double sum_mol_dip_frame = 0;
        for (int m = 0; m < sys.numMolecules(); m++)
        {
          Vec mol_dip(0, 0, 0);

          Vec offsite_i(0, 0, 0), offsite_j(0, 0, 0);

          Molecule molecule = sys.mol(m);

          for (int a = 0; a < molecule.numAtoms(); a++)
          {
            const AtomTopology &current_atom = molecule.topology().atom(a);

            Vec atom_pos = molecule.pos(a);
            if (current_atom.isPolarisable() && current_atom.poloffsiteGamma() > 0.000001)
            {
              offsite_i = molecule.pos(current_atom.poloffsiteI());
              offsite_j = molecule.pos(current_atom.poloffsiteJ());
            }
            Vec cosdispl = molecule.cosDisplacement(a);

            get_atom_dipole(current_atom, atom_pos, offsite_i, offsite_j, cosdispl, mol_dip, nchargepa);
          }
          tot_dipole += mol_dip;
          sum_mol_dip_frame += mol_dip.abs();
        }

        unsigned int num_solv_molecules = 0;
        for (int s = 0; s < sys.numSolvents(); s++)
        {
          int num_solvent_atoms = sys.sol(s).topology().numAtoms();
          //for (unsigned int i=0; i < sys.sol(s).numAtoms(); i++)
          //unsigned int i=0;

          for (int i = 0; i < sys.sol(s).numAtoms(); i += num_solvent_atoms)
          {
            num_solv_molecules++;
            Vec mol_dip(0, 0, 0);

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

              get_atom_dipole(current_atom, atom_pos, offsite_i, offsite_j, cosdispl, mol_dip, nchargepa);
            }

            tot_dipole += mol_dip;
            sum_mol_dip_frame += mol_dip.abs();
          }
        }
        sum_mol_dip_frame = sum_mol_dip_frame / (sys.numMolecules() + num_solv_molecules);

        // store for autocorrelation
        dipole_vector.push_back(tot_dipole);

        // do some bookkeeping
        numFrames++;

        sum_dip2 += tot_dipole.abs2();
        sum_dip += tot_dipole;

        // calculate the fluctuations of the dipole moment
        fluc = sum_dip2 / numFrames - (sum_dip / numFrames).abs2();

        // calculate the current estimate of eps
        fac = f * sum_vol / numFrames;
        a = 2 * e_rf * fluc + fac;
        b = -fluc + fac;
        eps = a / b;

        cout << time_trj
             << setw(15) << setprecision(8) << tot_dipole.abs() * conversion_fac
             << setw(15) << setprecision(8) << sum_mol_dip_frame * conversion_fac
             << setw(15) << setprecision(8) << eps
             << endl;
      }
      ic.close();
      is.close();
    }

    sum_dip2 /= numFrames;

    double dt;
    if (time.size() > 1)
      dt = time[1] - time[0];
    else
      dt = 1;

    if (do_autocorr)
    {
      //now doing the auto correlation calculation

      ofstream os(autocorr_fname.c_str());
      os << "#\n";
      os << "# normalized autocorrelation function C(t) of the box dipole\n";
      os << "#" << setw(11) << "time" << setw(20) << "C(t)" << setw(15) << "indep.frames\n";

      gmath::Correlation *corr;
      // several cases
      corr = new gmath::Correlation(dipole_vector, dipole_vector);
      corr->calc_direct();

      double tau = 0;
      for (unsigned int i = 0; i < corr->size(); i++, tau += dt)
      {
        unsigned int independent_contrib = 0;
        if (i != 0)
        {
          independent_contrib = corr->size() / i;
          if (independent_contrib < truncate_at)
            break;
        }

        os << setw(12) << tau
           << setw(20) << (*corr)[i] / sum_dip2
           << setw(15) << independent_contrib
           << std::endl;
      }

      delete pbc;
      delete corr;
      os.close();
    }
  }

  catch (const gromos::Exception &e)
  {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
