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
 * @file pert_top.cc
 * Creates a perturbation topology to remove interactions for specified atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pert_top
 * @section pert_top Creates a perturbation topology to alter interactions for 
 * specified atoms
 * @author @ref co @ref ns
 * @date 7-6-07
 *
 * Creates a perturbation topology to perturb specified atoms. A perturbation 
 * topology is written that defines a perturbation to alter the specified atoms
 * to the specified atom types, charges and masses. Each of the arguments
 * \@types, \@masses and \@charges can be omitted. In this case the values from
 * the topology are taken. If not sufficient values are given, the last given 
 * value is taken for all the remaining atoms.

 * Use program @ref pt_top to convert the resulting perturbation topology to a
 * different format or to a regular molecular topology.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to be modified </td></tr>
 * <tr><td> \@types</td><td>&lt;IACS of the perturbed atoms&gt; </td></tr>
 * <tr><td> \@charges</td><td>&lt;charges of the perturbed atoms&gt; </td></tr>
 * <tr><td> \@masses</td><td>&lt;masses of the perturbed atoms&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pert_top
    @topo    ex.top
    @atoms   1:34-51
    @types   13 19
    @charges 0.1 0.0
    @masses  1.008
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/PtTopology.h"
#include "../src/gio/OutPtTopology.h"
#include "../src/gcore/System.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char *argv[]){

  Argument_List knowns;
  knowns << "topo" << "atoms" << "types" << "charges" << "masses";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@atoms   <atoms to be modified\n";
  usage += "\t@types   <IACS of the perturbed atoms>\n";
  usage += "\t@charges <charges of the perturbed atoms>\n";
  usage += "\t@masses  <masses of the perturbed atoms>\n";
   
  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    
    // get the atoms to perturb
    utils::AtomSpecifier at(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      for(; iter!=args.upper_bound("atoms"); ++iter){
	at.addSpecifier(iter->second);
      }
    }
    
    if (at.empty())
      throw gromos::Exception(argv[0], "please give the atoms you want to perturb.");

    // get the new IACs
    vector<int> iacs;
    {
      Arguments::const_iterator it = args.lower_bound("types"),
          to = args.upper_bound("types");
      for (; it != to; ++it) {
        istringstream is(it->second);
        int type;
        if (!(is >> type))
          throw gromos::Exception(argv[0], "value in @types is not numeric.");
        iacs.push_back(type-1);
      }
    }

    // get the new charges
    vector<double> charges;
    {
      Arguments::const_iterator it = args.lower_bound("charges"),
          to = args.upper_bound("charges");
      for (; it != to; ++it) {
        istringstream is(it->second);
        double charge;
        if (!(is >> charge))
          throw gromos::Exception(argv[0], "value in @charges is not numeric.");
        charges.push_back(charge);
      }
    }
    
    // get the new masses
    vector<double> masses;
    {
      Arguments::const_iterator it = args.lower_bound("masses"),
          to = args.upper_bound("masses");
      for (; it != to; ++it) {
        istringstream is(it->second);
        double mass;
        if (!(is >> mass))
          throw gromos::Exception(argv[0], "value in @masses is not numeric.");
        masses.push_back(mass);
      }
    }
    
    // this just makes no sense
    if (iacs.empty() && charges.empty() && masses.empty())
      throw gromos::Exception(argv[0], "please give at least one IAC, charge or mass to perturb.");

    // create the perturbation topology
    gcore::PtTopology pttop;
    pttop.setSize(at.size(), 2);
    
    OutPtTopology out_pt(cout);

    // set the title
    ostringstream title;
    title << "Perturbation of atoms: ";
    for (unsigned int i = 0; i < at.size(); i++) title << at.toString(i) << " ";
    title << endl << "From topology: " << args["topo"];
    out_pt.setTitle(title.str());
    
    for (unsigned int i = 0; i < at.size(); i++) {
      // set state A
      pttop.setAtomNum(i, at.gromosAtom(i));
      pttop.setAtomName(i, at.name(i));
      pttop.setResidueNum(i, at.resnum(i));
      pttop.setResidueName(i, at.resname(i));
      pttop.setIac(i, 0, at.iac(i));
      pttop.setMass(i, 0, at.mass(i));
      pttop.setCharge(i, 0, at.charge(i));
      
      // set state B
      // perturbation in iac?
      int iac;
      if (iacs.empty()) { // no perturbation
        iac = at.iac(i);
      } else {
        // take the right value or the last one
        iac = i < (unsigned int)(iacs.size()) ? iacs[i] : iacs.back();
      }
      pttop.setIac(i, 1, iac);
      
      // perturbation in mass?
      double mass;
      if (masses.empty()) { // no perturbation
        mass = at.mass(i);
      } else {
        // take the right value or the last one
        mass = i < (unsigned int)(masses.size()) ? masses[i] : masses.back();
      }
      pttop.setMass(i, 1, mass);
      
      // perturbation in charge?
      double charge;
      if (charges.empty()) { // no perturbation
        charge = at.charge(i);
      } else {
        // take the right value or the last one
        charge = i < (unsigned int)(charges.size()) ? charges[i] : charges.back();
      }
      pttop.setCharge(i, 1, charge);
    }
    
    // write out the resultung topology
    out_pt.write(pttop);
    return 0;
  } catch (gromos::Exception e) {
    cerr << e.what() << endl;
    return 1;
  }
}
