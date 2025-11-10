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
 * @file lb_top.cc
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor lb_top
 * @section lb_top converts the Lennard-Jones parameters into the parameters of the Lorentz-Berthelot combination rules
 * @author @ref co
 *
 * The program lb_top converts the Lennard-Jones parameters from an existing 
 * gromos96 topology into the parameters of the Lorentz-Berthelot combination 
 * rules for non-equal atomtypes
 *        GROMOS uses:
 *                    @f[C6(i,j) = \sqrt{C6(i,i) C6(j,j)} @f]
 *                    @f[C12(i,j) = \sqrt{C12(i,i) C12(j,j)} @f]
 *        Lorentz-Berthelot:
 *                    @f[\epsilon(i,j) = \sqrt{\epsilon(i,i)\epsilon(j,j)} @f]
 *                    @f[\sigma(i,j) = \frac{\sigma(i,i) + \sigma(j,j)}{2} @f]
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 */

#include <cassert>
#include <iostream>
#include <cmath>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InParameter.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/System.h"
#include "../src/gromos/Exception.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace std;


int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "topo";
  
  string usage = "#" + string(argv[0]);
  usage += "\n\t@topo <topology>\n";
  
  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);


    System sys(it.system());

    
    OutTopology ot(cout);
    string addtitle;
    addtitle+="\nConverted to the Lorentz-Berthelot combination rules";
    ot.setTitle(it.title()+addtitle);

    GromosForceField gff(it.forceField());

    //loop over all non-equal atom pairs
    for(int i=0; i< gff.numAtomTypeNames(); i++){
      for(int j=0; j<i; j++){
        //get parameters for ii and for jj
	AtomPair ii(i,i);
	AtomPair jj(j,j);
	LJType ljii(gff.ljType(ii));
	LJType ljjj(gff.ljType(jj));
	
	// calculate epsilon and sigma for ii and jj 
        double eii=0, esii=0, ejj=0, esjj=0;
	
	if(ljii.c12() !=0) eii =ljii.c6()  * ljii.c6()  / ljii.c12()  / 4;
	if(ljii.cs12()!=0) esii=ljii.cs6() * ljii.cs6() / ljii.cs12() / 4;
        if(ljjj.c12() !=0) ejj =ljjj.c6()  * ljjj.c6()  / ljjj.c12()  / 4;
	if(ljjj.cs12()!=0) esjj=ljjj.cs6() * ljjj.cs6() / ljjj.cs12() / 4;
	
	double sii=0, ssii=0, sjj=0, ssjj=0;
	
	if(ljii.c6() !=0) sii =pow(ljii.c12() /ljii.c6() , 1.0/6.0);
	if(ljii.cs6()!=0) ssii=pow(ljii.cs12()/ljii.cs6(), 1.0/6.0);
	if(ljjj.c6() !=0) sjj =pow(ljjj.c12() /ljjj.c6() , 1.0/6.0);
	if(ljjj.cs6()!=0) ssjj=pow(ljjj.cs12()/ljjj.cs6(), 1.0/6.0);
	
	// now calculate epsilon and sigma for ij
	AtomPair ij(i,j);
        double eij = sqrt(eii  * ejj );
	double esij= sqrt(esii * esjj);
	double sij = (sii  + sjj )/2;
	double ssij= (ssii + ssjj)/2;

	// and calculate them back to C6 and C12
	double s6;
	s6=sij*sij*sij*sij*sij*sij;
	double c6ij = 4*eij*s6;
	double c12ij = 4*eij*s6*s6;
	s6=ssij*ssij*ssij*ssij*ssij*ssij;
	double cs6ij = 4*esij*s6;
	double cs12ij= 4*esij*s6*s6;
	//if(i==3&&j==0){
	//  cout << "eii " << eii << endl;
	//  cout << "sii " << sii << endl;
	//  cout << "ejj " << ejj << endl;
	//  cout << "sjj " << sjj << endl;
	//  cout << "eij " << eij << endl;
	//  cout << "sij " << sij << endl;
	//  cout << "c6ij " << c6ij << endl;
	//  cout << "c12ij " << c12ij << endl;
	//}
	
	
        LJType ljij(c12ij, c6ij, cs12ij, cs6ij);
	gff.setLJType(ij, ljij);
      }
    }
    
    ot.write(sys,gff);
    
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}





