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
#include "FfExpert.h"

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "FfExpertGraph.h"
#include "../gcore/BuildingBlock.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/BbSolute.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Angle.h"
#include "../gcore/Dihedral.h"
#include "../gcore/Improper.h"
#include "../gcore/MoleculeTopology.h"

void utils::FfExpert::learn(gcore::BuildingBlock const & mtb, const utils::FfExpertGraphMapper * mapper)
{
  // loop over all buildingblock (solute only)
  for(int b=0; b<mtb.numBbSolutes(); b++){
    // loop over all atoms
    for(int a=0; a<mtb.bb(b).numAtoms(); a++){
      std::string s=mtb.bb(b).atom(a).name().substr(0,1);
      int iac=mtb.bb(b).atom(a).iac();
      bool found=false;
      for(std::multimap<std::string,counter>::iterator it=d_name2iac.lower_bound(s),
	    to=d_name2iac.upper_bound(s); it!=to; ++it){
	if (it->second.type == iac) {
	  it->second.occurence++;
	  found=true;
	}
      }
      
      if(!found) 
	d_name2iac.insert(std::multimap<std::string,counter>::value_type(s,counter(iac,1)));
      found=false;
      int mass=int(rint(mtb.bb(b).atom(a).mass()));


      for(std::multimap<int,counter>::iterator it=d_iac2mass.lower_bound(iac),
	    to=d_iac2mass.upper_bound(iac); it!=to; ++it){
	if (it->second.type == mass) {
	  it->second.occurence++;
	  found=true;
	}
      }
      
      if(!found) 
	d_iac2mass.insert(std::multimap<int,counter>::value_type(iac,counter(mass,1))); 
      found=false;
      double charge=mtb.bb(b).atom(a).charge();
      int icharge=-1;
      for(unsigned int ii=0; ii<d_chargeType.size(); ii++){
	if(d_chargeType[ii]==charge) {
	  icharge=ii;
	  break;
	}
      }
      if(icharge==-1){
	icharge=d_chargeType.size();
	d_chargeType.push_back(charge);
      }
      for(std::multimap<int,counter>::iterator it=d_iac2charge.lower_bound(iac),
	    to=d_iac2charge.upper_bound(iac); it!=to; ++it){
	if (it->second.type == icharge) {
	  it->second.occurence++;
	  found=true;
	}
      }
      
      if(!found) 
	d_iac2charge.insert(std::multimap<int,counter>::value_type(iac,counter(icharge,1)));
    }
    {
      // loop over the bonds
      gcore::BondIterator bi(mtb.bb(b));
      for(;bi; ++bi){
	int ii=bi()[0];
	if(ii<0) ii+=mtb.bb(b).numAtoms();
	if(ii>=mtb.bb(b).numAtoms()) ii-=mtb.bb(b).numAtoms();
	int ij=bi()[1];
	if(ij<0) ij+=mtb.bb(b).numAtoms();
	if(ij>=mtb.bb(b).numAtoms()) ij-=mtb.bb(b).numAtoms();

	gcore::Bond iacbond(mtb.bb(b).atom(ii).iac(),
		            mtb.bb(b).atom(ij).iac(), 0);

	int type=bi().type();
	bool found=false;
	for(std::multimap<gcore::Bond,counter>::iterator 
	      it=d_iac2bond.lower_bound(iacbond),
	      to=d_iac2bond.upper_bound(iacbond);
	    it!=to;
	    ++it){
	  if(it->second.type==type){
	    it->second.occurence++;
	    found=true;
	  }
	}
	if(!found)
	  d_iac2bond.insert(std::multimap<gcore::Bond,counter>::value_type(iacbond, counter(type,1)));
      }
    }

    // loop over the angles
    {
      gcore::AngleIterator bi(mtb.bb(b));
      for(;bi; ++bi){
	int ii=bi()[0];
	if(ii<0) ii+=mtb.bb(b).numAtoms();
	if(ii>=mtb.bb(b).numAtoms()) ii-=mtb.bb(b).numAtoms();
	int ij=bi()[1];
	if(ij<0) ij+=mtb.bb(b).numAtoms();
	if(ij>=mtb.bb(b).numAtoms()) ij-=mtb.bb(b).numAtoms();
	int ik=bi()[2];
	if(ik<0) ik+=mtb.bb(b).numAtoms();
	if(ik>=mtb.bb(b).numAtoms()) ik-=mtb.bb(b).numAtoms();

	gcore::Angle iacangle(mtb.bb(b).atom(ii).iac(),
	   	              mtb.bb(b).atom(ij).iac(),
		              mtb.bb(b).atom(ik).iac(), 0);
	int type=bi().type();
	bool found=false;
	for(std::multimap<gcore::Angle,counter>::iterator 
	      it=d_iac2angle.lower_bound(iacangle),
	      to=d_iac2angle.upper_bound(iacangle);
	    it!=to;
	    ++it){
	  if(it->second.type==type){
	    it->second.occurence++;
	    found=true;
	  }
	}
	if(!found)
	  d_iac2angle.insert(std::multimap<gcore::Angle,counter>::value_type(iacangle, counter(type,1)));
      }
    }
    {
      // loop over the impropers
      gcore::ImproperIterator bi(mtb.bb(b));
      for(;bi; ++bi){
	int ii=bi()[0];
	if(ii<0) ii+=mtb.bb(b).numAtoms();
	if(ii>=mtb.bb(b).numAtoms()) ii-=mtb.bb(b).numAtoms();
	int ij=bi()[1];
	if(ij<0) ij+=mtb.bb(b).numAtoms();
	if(ij>=mtb.bb(b).numAtoms()) ij-=mtb.bb(b).numAtoms();
	int ik=bi()[2];
	if(ik<0) ik+=mtb.bb(b).numAtoms();
	if(ik>=mtb.bb(b).numAtoms()) ik-=mtb.bb(b).numAtoms();
	int il=bi()[3];
	if(il<0) il+=mtb.bb(b).numAtoms();
	if(il>=mtb.bb(b).numAtoms()) il-=mtb.bb(b).numAtoms();
	gcore::Improper iacimproper(mtb.bb(b).atom(ii).iac(),
				    mtb.bb(b).atom(ij).iac(),
				    mtb.bb(b).atom(ik).iac(),
				    mtb.bb(b).atom(il).iac(), 0);
	int type=bi().type();
	bool found=false;
	for(std::multimap<gcore::Improper,counter>::iterator 
	      it=d_iac2improper.lower_bound(iacimproper),
	      to=d_iac2improper.upper_bound(iacimproper);
	    it!=to;
	    ++it){
	  if(it->second.type==type){
	    it->second.occurence++;
	    found=true;
	  }
	}
	if(!found)
	  d_iac2improper.insert(std::multimap<gcore::Improper,counter>::
				value_type(iacimproper, counter(type,1)));
      }
    }
    {
      // loop over the dihedrals
      gcore::DihedralIterator bi(mtb.bb(b));
      for(;bi; ++bi){
	bool stop=false;
	
	int ii=bi()[0];
	if(ii<-2) stop=true;
	if(ii<0) ii+=mtb.bb(b).numAtoms();
	if(ii>=mtb.bb(b).numAtoms()) ii-=mtb.bb(b).numAtoms();

	int ij=bi()[1];
	if(ij<-2) stop=true;
	
	if(ij<0) ij+=mtb.bb(b).numAtoms();
	if(ij>=mtb.bb(b).numAtoms()) ij-=mtb.bb(b).numAtoms();

	int ik=bi()[2];
	if(ik<-2) stop=true;
	
	if(ik<0) ik+=mtb.bb(b).numAtoms();
	if(ik>=mtb.bb(b).numAtoms()) ik-=mtb.bb(b).numAtoms();
	int il=bi()[3];
	if(il<-2) stop=true;
	
	if(il<0) il+=mtb.bb(b).numAtoms();
	if(il>=mtb.bb(b).numAtoms()) il-=mtb.bb(b).numAtoms();

	if(stop) continue;
	
	gcore::Dihedral iacdihedral(mtb.bb(b).atom(ii).iac(),
				    mtb.bb(b).atom(ij).iac(),
				    mtb.bb(b).atom(ik).iac(),
				    mtb.bb(b).atom(il).iac(),
                                    0);
	int type=bi().type();
	bool found=false;
	for(std::multimap<gcore::Dihedral,counter>::iterator 
	      it=d_iac2dihedral.lower_bound(iacdihedral),
	      to=d_iac2dihedral.upper_bound(iacdihedral);
	    it!=to;
	    ++it){
	  if(it->second.type==type){
	    it->second.occurence++;
	    found=true;
	  }
	}
	if(!found)
	  d_iac2dihedral.insert(std::multimap<gcore::Dihedral,counter>::
				value_type(iacdihedral, counter(type,1)));
      }
    }
    
    if (mapper != NULL) { // create the graphs
      d_graphs.push_back(utils::FfExpertGraph(*mapper, mtb.bb(b)));
    }
    
  }
  return;
}

void utils::FfExpert::substructure2iac(unsigned int i, const utils::FfExpertGraph & query, 
        std::vector<std::vector<utils::Vertex> > & iacs) const {
  assert(i < query.vertices().size());
  // go up to oder 4
  const unsigned int max_order = 4;
  iacs.clear();
  iacs.resize(max_order);
  for(unsigned int radius = max_order; radius != 0; --radius) {
    const unsigned int r = radius - 1;
    const utils::FfExpertGraph & substructure = query.cut_subgraph(query.vertices()[i], radius);
    // loop over the dataset
    iacs[r].clear();
    for(std::vector<utils::FfExpertGraph>::const_iterator it = d_graphs.begin(),
            to = d_graphs.end(); it != to; ++it) {
      const std::vector<utils::Vertex> & hits = it->equal_subgraph(substructure, utils::Vertex::equal_type, radius);
      iacs[r].insert(iacs[r].end(), hits.begin(), hits.end());
    }
  }
}

void utils::FfExpert::substructure2charge(unsigned int i, const utils::FfExpertGraph & query, 
        std::vector<std::vector<utils::Vertex> > & charge) {
  assert(i < query.vertices().size());
  // go up to oder 4
  const unsigned int max_order = 4;
  charge.clear();
  charge.resize(max_order);
  for(unsigned int radius = max_order; radius != 0; --radius) {
    const unsigned int r = radius - 1;
    const utils::FfExpertGraph & substructure = query.cut_subgraph(query.vertices()[i], radius);
    // loop over the dataset
    charge[r].clear();
    for(std::vector<utils::FfExpertGraph>::const_iterator it = d_graphs.begin(),
            to = d_graphs.end(); it != to; ++it) {
      const std::vector<utils::Vertex> & hits = it->equal_subgraph(substructure, utils::Vertex::equal_iac, radius);
      charge[r].insert(charge[r].end(), hits.begin(), hits.end());
    }
  }
}

int utils::sort(std::vector<FfExpert::counter> &v, bool tt)
{
  int max_occur = 0, max_index = 0;
  if (tt) {
    for (unsigned int i = 1; i < v.size(); i++) {

      FfExpert::counter t = v[i];
      int j = i - 1;
      while ((j >= 0) && t.type < v[j].type) {
        v[j + 1] = v[j];
        j--;
      }
      v[j + 1] = t;
    }
    for (unsigned int i = 0; i < v.size(); i++) {
      if (v[i].occurence > max_occur) {
        max_occur = v[i].occurence;
        max_index = i;
      }
    }
    return max_index;
  } else {
    for (unsigned int i = 1; i < v.size(); i++) {
      FfExpert::counter t = v[i];
      int j = i - 1;
      while ((j >= 0) && t.occurence > v[j].occurence) {
        v[j + 1] = v[j];
        j--;
      }
      v[j + 1] = t;
    }
    return 0;
  }
  return 0;
}
