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
 * @file ener.cc
 * Recalculates interaction energies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor ener
 * @section ener Recalculates interaction energies
 * @author @ref co
 * @date 22-11-2004
 *
 * Program ener can recalculate interaction energies over molecular trajectory 
 * files using the interaction parameters specified in the molecular topology 
 * file.
 *
 * Nonbonded interactions are calculated for all selected atoms with all other
 * atoms in the system. Some atoms can be specified as being soft, indicating 
 * that interactions involving any of these atoms have a specified softness
 * parameter, for all other atoms in the system, the softness parameter
 * @f$\alpha = 0@f$. Vanderwaals interactions between particles i and j are
 * calculated as
 * 
 * @f[ V^{vdw}_{ij}=\left[\frac{C_{12}(i,j)}{ ( r_{ij}^6 + \alpha_{LJ} \lambda ^2 C_{126})}-C_6(i,j)\right] \frac{1}{(r_{ij}^6 + \alpha_{LJ} \lambda ^2 C_{126})} @f]
 *
 * with @f$C_{126} = C_{12}/C_6 @f$ for @f$C_{12}@f$ and @f$C_6@f$ unequal 0,
 * @f$C_{126} = 0@f$ otherwise. @f$C_{12}@f$ and @f$C_6@f$ are the interaction
 * parameters taken from the topology, @f$\lambda@f$ and @f$\alpha_{LJ}@f$ are
 * specified by the user. Similarly, the electrostatic interaction, including
 * reaction field contribution for a homogeneous medium outside the cutoff
 * sphere is calculated as 
 *
 * @f[ V^{ele}_{ij}=\frac{q_iq_j}{4\pi\epsilon_0}\left[\frac{1}{(r^2_{ij}+\alpha_{C}\lambda^2)^{1/2}} - \frac{\frac{1}{2}C_{RF}r_{ij}^2}{(R_{RF}^2+\alpha_{C}\lambda^2)^{3/2}} - \frac{(1-\frac{1}{2}C_{RF})}{R_{RF}}\right] @f]
 *
 * where @f$\epsilon_0@f$ is the dielectric permittivity of vacuum and 
 * @f$q_i@f$ and @f$q_j@f$ are the atomic partial charges. @f$R_{RF}@f$ is the
 * reaction field cutoff distance, here assumed to be the same as the
 * interaction cutoff. @f$\alpha_{C}@f$ and @f$\lambda@f$ are again user 
 * specified. @f$C_{RF}@f$ is calculated from the reaction field dielectric
 * constant @f$\epsilon_{RF}@f$ and @f$\kappa_{RF}^{-1}@f$ (user specified) as
 *
 * @f[ C_{RF} = \frac{ (2 - 2 \epsilon_{RF}) (1 + \kappa_{RF}^{-1} R_{RF}) - \epsilon_{RF} (\kappa_{RF}^{-1} R_{RF})^2 }{ (1 + 2 \epsilon_{RF}) (1 + \kappa_{RF}^{-1} R_{RF}) + \epsilon_{RF} (\kappa_{RF}^{-1} R_{RF})^2 } @f]
 *
 * The bonded interactiona are calculated for all specified properties using 
 * the following interaction functions. For bonds we use:
 *
 * @f[ V^{bond}=\frac{1}{4}K_{b_n}\left[b_n^2 - b_{0_n}^2\right]^2 @f]
 *
 * with @f$b_n@f$ the actual bond length, @f$K_{b_n}@f$ and @f$b_{0_n}@f$ the 
 * force constant and optimal bond length, respectively. For angles we use:
 *
 * @f[ V^{angle}=\frac{1}{2}K_{\theta_n}\left[\cos{\theta_n} - \cos{\theta_{0_n}}\right]^2 @f]
 *
 * with @f$\theta_n@f$ the actual bond angle, @f$K_{\theta_n}@f$ and 
 * @f$\theta_{0_n}@f$ the force constant and optimal bond angle respectively.
 * For proper torsional dihedral angle terms we use:
 *
 * @f[ V^{trig}=K_{\phi_n}\left[1+\cos(\delta_n)\cos(m_n\phi_n)\right] @f]
 *
 * with @f$\phi_n@f$ the actual dihedral angle value, @f$K_{\phi_n}@f$ the
 * force constant and @f$\delta_n@f$ and @f$m_n@f$ the phase shift and
 * multiplicity, respectively. Improper dihedral energy contributions are 
 * calculated from:
 * @f[ V^{har}=\frac{1}{2}K_{\xi_n}\left[\xi_n - \xi_{0_n}\right]^2 @f]
 *
 * with @f$\xi_n@f$ the actual dihedral angle value, @f$K_{\xi_n}@f$ and
 * @f$\xi_{0_n}@f$ are the force constant and optimal improper dihedral angle 
 * value.
 *
 * The programme prints out the total bonded and nonbonded energies separately,
 * as well as the overall total energy. It is easily modified to print out more
 * detailed energy contributions as well.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" for nonbonded interaction&gt; </td></tr>
 * <tr><td> \@props</td><td>&lt;@ref PropertySpecifier "properties" to be calculated&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field contribution&gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa_{RF}^{-1} for reaction field contribution&gt; </td></tr>
 * <tr><td> [\@RFex</td><td>&lt;calculate RF contribution for excluded atoms: on/off&gt;] </td></tr>
 * <tr><td> \@soft</td><td>&lt;soft @ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@softpar</td><td>&lt;lam&gt; &lt;a_lj&gt; &lt;a_crf&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ener
    @topo    ex.top
    @pbc     r
    @atoms   1:3-13
    @props   d%1:1,2 a%1:1,2,3 t%1:1,2,4,6 t%1:4,2,5,6
    [@time    0 0.2]
    @cut     1.4
    @eps     61
    @kap     0.0
    @RFex    on
    @soft    1:4
    @softpar 0.5 1.51 0.5
    @traj    ex.tr
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/groTime.h"
#include "../src/utils/IntegerInputParser.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "atoms" << "energies" << "props" << "time" << "cut"
         << "eps" << "kap" << "soft" << "softpar" << "RFex" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type> [<gather method>]\n";
  usage += "\t@atoms    <atoms for nonbonded>\n";
  usage += "\t@energies <integer identifier for the energies to write out:\n";
  usage += "\t\t   can take comma-separated identifiers, ranges, or a combination (e.g. 1,4,7-9)\n";
  usage += "\t\t   1..covalent (use @props to specify properties)\n";
  usage += "\t\t   2..elec with solute\n";
  usage += "\t\t   3..elec with solvent\n";
  usage += "\t\t   4..elec total\n";
  usage += "\t\t   5..vdw with solute\n";
  usage += "\t\t   6..vdw with solvent\n";
  usage += "\t\t   7..vdw total\n";
  usage += "\t\t   8..nonbonded total\n";
  usage += "\t\t   9..total\n";
  usage += "\t\t   all..write out 1-9>\n";
  usage += "\t[@props   <propertyspecifier for covalent energies>]\n";
  usage += "\t@cut      <cut-off distance in nm>\n";
  usage += "\t@eps      <epsilon for reaction field correction>\n";
  usage += "\t[@kap      <kappa_{RF}^{-1} for reaction field correction (default: 0)]>\n";
  usage += "\t@RFex     <calculate RF for excluded atoms: on/off (default: on)>\n";
  usage += "\t@coulomb_scaling <use scaling of 1/1.2 for electrostatic 1-4 interactions (default: off)\n>";
  usage += "\t[@soft    <soft atoms>]\n";
  usage += "\t[@softpar <lam> <a_lj> <a_c>]\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@traj     <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  // get the @time argument
    utils::Time time(args);

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // declare the energy class
  Energy en(sys, gff, *pbc);

  //  set atoms
  AtomSpecifier atoms(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  en.setAtoms(atoms);

  // RF for excluded atoms?
  {
    std::string s=args.getValue<string>("RFex",false,"on");
    if (s!="off" && s!="on")
      throw gromos::Exception("ener","unknown value for @RFex flag - use on/off");
    if(s=="off")
      en.setRFexclusions(false);
    else
      en.setRFexclusions(true);
  }
  
  // coulomb scaling for electrostatic 1-4 interactions?
  {
    std::string s=args.getValue<string>("coulomb_scaling",false,"off");
    if(s=="off")
      en.setCoulombScaling(false);
    else
      en.setCoulombScaling(true);
  }

  // set properties
  PropertyContainer props(sys, pbc);
  {
    Arguments::const_iterator iter=args.lower_bound("props");
    Arguments::const_iterator to=args.upper_bound("props");
    for(;iter!=to;iter++){
      string p=iter->second.c_str();
      props.addSpecifier(p);
    }
  }
  en.setProperties(props);
  // set non-bonded parameters
  //   get cut-off distance
  {
    double cut = args.getValue<double>("cut", true, 1.4);
    en.setCutOff(cut);
  }
  //  get epsilon and kappa_{RF}^{-1}
  {
    double eps = args.getValue<double>("eps", true, 61.0);
    double kap = args.getValue<double>("kap", false, 0.0);
    en.setRF(eps, kap);
  }
  // get soft atom list
  AtomSpecifier soft(sys);
  {
    bool lsoft=false;
    Arguments::const_iterator iter=args.lower_bound("soft");
    Arguments::const_iterator to=args.upper_bound("soft");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      soft.addSpecifier(spec);
      lsoft=true;
    }
    //  get softness parameter
    std::vector<double> softpar = args.getValues<double>("softpar", 3, lsoft, 
            Arguments::Default<double>() << 0.0 << 0.0 << 0.0);
    if (lsoft)
      en.setSoft(soft, softpar[0], softpar[1], softpar[2]);
  }

  // get energies that user requested
  IntegerInputParser iip;
  iip.addSpecifier(args["energies"], 9); // we have 9 different energy types
  
  // define input coordinate
  InG96 ic;

  // check user requests for covalent properties
  if ( (props.size() > 0) && (*iip.begin() > 1) ) {
      throw gromos::Exception("ener","At least one property for covalent interactions was requested but covalent energies are not written out; use '@energies 1'");
    } else if ( (props.size() == 0) && (*iip.begin() == 1) ) {
      throw gromos::Exception("ener","Covalent energies were requested but no properties specified; use '@props' specifier");
    }
    
  // print titles
  cout << setw(15) << "# Time";
  for (IntegerInputParser::const_iterator iter = iip.begin(), to = iip.end(); iter != to; ++iter) {
    switch(*iter) {
      //    case 1: cout << setw(17) << "bonds"; break;
      //    case 2: cout << setw(17) << "angles"; break;
      //    case 3: cout << setw(17) << "impdihedrals"; break;
      //    case 4: cout << setw(17) << "torsdihedrals"; break;
    case 1: cout << setw(17) << "covalent"; break;
    case 2: cout << setw(17) << "elec_solute"; break;
    case 3: cout << setw(17) << "elec_solvent"; break;
    case 4: cout << setw(17) << "elec_tot"; break;
    case 5: cout << setw(17) << "vdw_solute"; break;
    case 6: cout << setw(17) << "vdw_solvent"; break;
    case 7: cout << setw(17) << "vdw_tot"; break;
    case 8: cout << setw(17) << "nonbonded"; break;
    case 9: cout << setw(17) << "total"; break;
    default: throw gromos::Exception("ener: ","wrong energy specifier");
    }
  }
  cout << endl;
  
  // declare some variables for averaging
  int num_frames=0;
  double cov=0.0;
  double el_m=0.0;
  double el_s=0.0;
  double el_t=0.0;
  double vdw_m=0.0;
  double vdw_s=0.0;
  double vdw_t=0.0;
  double nb=0.0;
  double tot=0.0;
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");

    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys >> time;
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected
      pbc->gathergr();

      // calculate the energies
      en.calc();

      
      // declare a temporary variable for energies
      double e = 0.0;
      // set cout variables
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      // print user requested output
      cout << setw(15) << time;
      for (IntegerInputParser::const_iterator iter = iip.begin(), to = iip.end(); iter != to; ++iter) {
	switch(*iter) {
	  //	case 1: cout << setw(17) << en.dist(); break;
	  //	case 2: cout << setw(17) << en.angle(); break;
	  //	case 3: cout << setw(17) << en.impdihed(); break;
	  //	case 4: cout << setw(17) << en.torsdihed(); break;
	case 1:
	  e = en.cov();
	  cov += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 2:
	  e = en.el_m();
	  el_m += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 3:
	  e = en.el_s();
	  el_s += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 4:
	  e = en.el();
	  el_t += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 5:
	  e = en.vdw_m();
	  vdw_m += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 6:
	  e = en.vdw_s();
	  vdw_s += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 7:
	  e = en.vdw();
	  vdw_t += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 8:
	  e = en.el()+en.vdw();
	  nb += e;
	  cout << ' ' << setw(16) << e;
	  break;
	case 9:
	  e = en.tot();
	  tot += e;
	  cout << ' ' << setw(16) << e;
	  break;
	}
      }
      cout << endl;
      
      //increase the frames counter
      num_frames++;
    
    }
  }
  // print out averages
  if(num_frames>1){
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << setw(15) << "# ave.";
    for (IntegerInputParser::const_iterator iter = iip.begin(), to = iip.end(); iter != to; ++iter) {
      switch(*iter) {
	case 1: cout << ' ' << setw(16) << cov/num_frames; break;
	case 2: cout << ' ' << setw(16) << el_m/num_frames; break;
	case 3: cout << ' ' << setw(16) << el_s/num_frames; break;
	case 4: cout << ' ' << setw(16) << el_t/num_frames; break;
	case 5: cout << ' ' << setw(16) << vdw_m/num_frames; break;
	case 6: cout << ' ' << setw(16) << vdw_s/num_frames; break;
	case 7: cout << ' ' << setw(16) << vdw_t/num_frames; break;
	case 8: cout << ' ' << setw(16) << nb/num_frames; break;
	case 9: cout << ' ' << setw(16) << tot/num_frames; break;
      }
    }
    cout << endl;
  }
 }
  catch (const gromos::Exception &e){
    cerr << endl << e.what() << endl;
    exit(1);
  }
 return 0;
}







