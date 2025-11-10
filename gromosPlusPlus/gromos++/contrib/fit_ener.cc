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
 * @file fit_ener.cc
 * Recalculates interaction energies for a solute molecule superposed on a
 * specific location
 */

/**
 * @page programs Program Documentation
 *
 * @anchor fit_ener
 * @section ener Recalculates interaction energies for a solute molecule 
 * superposed on a specific location
 * @author @ref co
 * @date 16-07-2015
 *
 * THIS PROGRAM IS ONLY INTENDED FOR VERY SPECIFIC USE:
 * IT WILL TAKE THE COORDINATES OF THE SOLUTE SPECIFIED BY @FITTOPO AND ADD
 * THE SOLVENT FROM THE TRAJECTORY. THEN IT CALCULATES THE ENERGY OF THE 
 * SELECTED ATOMS
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
 * <tr><td> \@props</td><td>&lt;@ref PropertyContainer "properties" to be calculated&gt; </td></tr>
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
    @fittopo  solute.top
    @fitcoord solute.cnf
    @fitatoms <atomspecifier> <atomspecifier>
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

#include "../src/gmath/Vec.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/utils/groTime.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gromos/Exception.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "atoms" << "props" << "time" << "cut"
         << "eps" << "kap" << "soft" << "softpar" << "RFex" << "traj"
	 << "fittopo" << "fitcoord" << "fitatoms";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t@pbc      <boundary type> [<gather method>]\n";
  usage += "\t@atoms    <atoms for nonbonded>\n";
  usage += "\t@props    <propertyspecifier>\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@cut      <cut-off distance>\n";
  usage += "\t@eps      <epsilon for reaction field correction>\n";
  usage += "\t@kap      <kappa_{RF}^{-1} for reaction field correction>\n";
  usage += "\t[@RFex    <calculate RF for excluded atoms: on/off>]\n";
  usage += "\t[@soft    <soft atoms>]\n";
  usage += "\t[@softpar <lam> <a_lj> <a_c>]\n";
  usage += "\t[@fittopo  <topology for molecule to superpose>]\n";
  usage += "\t[@fitcoord <coordinate file for molecule to superpose>]\n";
  usage += "\t[@fitatoms <atomspecifier for topo> <atomspecifier for fittopo>]\n";
  usage += "\t@traj    <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  cerr << "THIS PROGRAM IS ONLY INTENDED FOR VERY SPECIFIC USE:\n"
       << "  IT WILL TAKE THE COORDINATES OF THE SOLUTE SPECIFIED BY @FITTOPO AND ADD\n"
       << "  THE SOLVENT FROM THE TRAJECTORY. THEN IT CALCULATES THE ENERGY OF THE\n"
       << "  SELECTED ATOMS\n";
  
  // get the @time argument
  utils::Time time(args);

  //  read topologies
  InTopology it(args["topo"]);
  InTopology fitit(args["fittopo"]);
  
  System sys(it.system());
  System fitsys(fitit.system());
  System calcsys(fitit.system());
  
  if(it.forceField().ForceField() != fitit.forceField().ForceField()){
    throw gromos::Exception("fit_ener", 
			    "force fields for topo and fit_topo are not the same\n");
  }
    
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(calcsys, args);

  // declare the energy class
  Energy en(calcsys, gff, *pbc);

  //  set atoms
  AtomSpecifier atoms(calcsys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  en.setAtoms(atoms);

  // get fitatoms
  AtomSpecifier fit_atoms_sys(sys);
  AtomSpecifier fit_atoms_fit(fitsys);
  {
    Arguments::const_iterator iter=args.lower_bound("fitatoms");
    Arguments::const_iterator to=args.upper_bound("fitatoms");
    string spec=iter->second.c_str();
    fit_atoms_sys.addSpecifier(spec);
    iter++;
    spec=iter->second.c_str();
    fit_atoms_fit.addSpecifier(spec);
    iter++;
    if(iter!=to){
      throw gromos::Exception("fit_ener", "you have to give two arguments for fitatoms");
    }
    if(fit_atoms_sys.size() != 1){
      throw gromos::Exception("fit_ener", "give one atom each for fitatoms (trajectory)");
      
    }    
    if(fit_atoms_fit.size() != 1){
      throw gromos::Exception("fit_ener", "give one atom each for fitatoms (fittopo)");
    }
  }

  // read in coordinates 
  {
    InG96 fitic;
    fitic.open(args["fitcoord"].c_str());
    fitic >> fitsys;
  }
  
    
  // RF for excluded atoms?
  {
    std::string s=args.getValue<string>("RFex",1,"on");
    if(s=="off")
      en.setRFexclusions(false);
    else
      en.setRFexclusions(true);
  }

  // set properties
  PropertyContainer props(calcsys, pbc);
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
    double cut = args.getValue<double>("cut", false, 1.4);
    en.setCutOff(cut);
  }
  //  get epsilon and kappa_{RF}^{-1}
  {
    double eps = args.getValue<double>("eps", false, 1.0);
    double kap = args.getValue<double>("kap", false, 0.0);
    en.setRF(eps, kap);
  }
  // get soft atom list
  AtomSpecifier soft(calcsys);
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
 
  // define input coordinate
  InG96 ic;
  
  
  // print titles
  cout << "# Time"
       << "              covalent"
       << "             nonbonded"
       << "                 Total"
       << endl;

  // declare some variables for averaging
  int num_frames=0;
  double cov=0.0;
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
      
      // move all atoms of fitsys to the correct position
      fit::PositionUtils::translate(&fitsys, 
				    fit_atoms_sys.pos(0)-fit_atoms_fit.pos(0));

      // calcsys gets the solute from fitsys
      for(int m=0; m< calcsys.numMolecules(); m++)
	for(int a=0; a < calcsys.mol(m).numAtoms(); a++)
	  calcsys.mol(m).pos(a) = fitsys.mol(m).pos(a);
      // calcsys gets the solvent and the box from sys
      calcsys.sol(0) = sys.sol(0);
      calcsys.box() = sys.box();
      calcsys.hasBox = true;
      
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected
      pbc->gathergr();

      // calculate the energies
      en.calc();

      // print any ouput you like
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << time
	   << ' ' << setw(21) << en.cov()
           << ' ' << setw(21) << en.nb()
           << ' ' << setw(21) << en.tot()
	   << endl;

      //store some averages
      cov+=en.cov();
      nb+=en.nb();
      tot+=en.tot();
      
      num_frames++;
    }
  }
  // print out averages
  if(num_frames>1){
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << endl << "# ave."
         << ' ' << setw(21) << cov/num_frames 
         << ' ' << setw(21) << nb/num_frames
         << ' ' << setw(21) << tot/num_frames
         << endl;
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







