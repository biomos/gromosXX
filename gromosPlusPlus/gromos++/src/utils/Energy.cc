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
#include "Energy.h"

#include <cassert>
#include <math.h>
#include <string>
#include <vector>
#include <vector>
#include <set>
#include <map>

#include "AtomSpecifier.h"
#include "SimplePairlist.h"
#include "PropertyContainer.h"
#include "Property.h"
#include "Value.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJType.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Exclusion.h"
#include "../gcore/Bond.h"
#include "../gcore/BondType.h"
#include "../gcore/Angle.h"
#include "../gcore/AngleType.h"
#include "../gcore/Improper.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/CrossDihedral.h"
#include "../gcore/DihedralType.h"
#include "../gcore/LJExceptionType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/GromosForceField.h"
#include "../gmath/Vec.h"
#include "../bound/Boundary.h"
#include "../gcore/MoleculeTopology.h"
#include "../gromos/Exception.h"

using namespace gcore;
using namespace std;
using namespace utils;
//using utils::Energy;
namespace utils {

  Energy::Energy(gcore::System &sys, gcore::GromosForceField &gff,
      bound::Boundary &pbc) {
    d_sys = &sys;
    d_gff = &gff;
    d_pbc = &pbc;
    d_as = new utils::AtomSpecifier(sys);
    d_pc = new utils::PropertyContainer(sys, &pbc);
    d_soft = new utils::AtomSpecifier(sys);
    d_lam = 0.0;
    d_alj = 0.0;
    d_ac = 0.0;
    d_eps = 1.0;
    d_kap = 0.0;
    d_cut = 1.4;
    d_RFex = true;
  }

  void Energy::calc() {
    calcNb();
    calcCov();
  }

  void Energy::calcNb() {
    // Make a pairlist
    calcPairlist();

    // and calculate the interactions
    calcNb_interactions();
  }

  void Energy::calcField() {
    // Make a pairlist
    calcPairlist();

    // and calculate the forces on the atoms

    // define some variables that we will need
    double qq, d1, d3, drf, el;
    double cut3 = d_cut * d_cut*d_cut;
    double crf = ((2 - 2 * d_eps)*(1 + d_kap * d_cut) - d_eps * (d_kap * d_kap * d_cut * d_cut)) /
        ((1 + 2 * d_eps)*(1 + d_kap * d_cut) + d_eps * (d_kap * d_kap * d_cut * d_cut));

    // loop over the atoms
    for (unsigned int i = 0; i < d_as->size(); i++) {

      double qi = d_as->charge(i);
      gmath::Vec vi = d_as->pos(i);

      // set the arrays for this atom to zero
      d_f_el_m[i] = Vec(0.0, 0.0, 0.0);
      d_f_el_s[i] = Vec(0.0, 0.0, 0.0);

      // now, loop over the pairlist
      for (unsigned int j = 0; j < d_pl[i].size(); j++) {

        // determine parameters

        qq = qi * d_pl[i].charge(j);

        // now, we calculate the distance between atoms
        gmath::Vec dd = d_pbc->nearestImage(vi, d_pl[i].pos(j), d_sys->box());

        dd = vi - dd;
        d1 = dd.abs();
        d3 = d1 * d1*d1;

        drf = 1 / d3 + crf / cut3;

        el = qq * drf * d_gff->fpepsi();

        // and store the energy in the correct array
        if (d_pl[i].mol(j) < 0) {
          d_f_el_s[i] += el*dd;
        } else {
          d_f_el_m[i] += el*dd;
        }
      }
    }
  }

  void Energy::calcNb_interactions() {
    // define some variables that we will need
    double cut3 = d_cut * d_cut*d_cut;
    double l2alj = d_lam * d_lam*d_alj;
    double l2ac = d_lam * d_lam*d_ac;
    const double crf = ((2 - 2 * d_eps)*(1 + d_kap * d_cut) - d_eps * (d_kap * d_kap * d_cut * d_cut)) /
        ((1 + 2 * d_eps)*(1 + d_kap * d_cut) + d_eps * (d_kap * d_kap * d_cut * d_cut));
    const double dirf = (1 - 0.5 * crf) / d_cut;

    // loop over the atoms
    double tmp_el = 0.0, tmp_vdw = 0.0;
#ifdef OMP
#pragma omp parallel for reduction(+ : tmp_el, tmp_vdw)
#endif
    for (unsigned int i = 0; i < d_as->size(); i++) {
      int mi = d_as->mol(i);
      int ai = d_as->atom(i);
      int gi = d_as->gromosAtom(i);
      int iaci = d_as->iac(i);
      const double qi = d_as->charge(i);

      gmath::Vec & vi = *(d_as->coord(i));

      // check if this atom is soft
      bool sft = d_soft->findAtom(mi, ai) != -1 ? true : false;

      // set the arrays for this atom to zero
      d_vdw_m[i] = 0.0;
      d_el_m[i] = 0.0;
      d_vdw_s[i] = 0.0;
      d_el_s[i] = 0.0;

      // now, loop over the pairlist
      for (unsigned int j = 0; j < d_pl[i].size(); j++) {
        int mj = d_pl[i].mol(j);
        int aj = d_pl[i].atom(j);
        int gj = d_pl[i].gromosAtom(j);

        // determine parameters
        gcore::LJType lj(d_gff->ljType(AtomPair(iaci, d_pl[i].iac(j))));
        double c6 = 0.0, c12 = 0.0;
        if (d_third[i].count(aj) && mj == mi) {
          c6 = lj.cs6();
          c12 = lj.cs12();
        } else {
          c6 = lj.c6();
          c12 = lj.c12();
        }
        const double qq = qi * d_pl[i].charge(j);

        // overwrite the LJ parameters in case of a LJ exception
        map<AtomPair, LJExceptionType>::const_iterator lje = d_gff->ljException().find(AtomPair(gi, gj));
        if(lje != d_gff->ljException().end()) {
          c6 = lje->second.c6();
          c12 = lje->second.c12();
        }

        // now, we calculate the distance between atoms
        gmath::Vec dd = d_pbc->nearestImage(vi, *d_pl[i].coord(j), d_sys->box());

        const double d1 = (vi - dd).abs();
        const double d2 = d1*d1;
        double d6 = d2 * d2*d2;
        double drf;

        // check if we have a soft atom
        if (sft || d_soft->findAtom(mj, aj) != -1) {
          if (c6 != 0.0 && c12 != 0.0) d6 += l2alj * c12 / c6;
          double cuts = l2ac + d_cut*d_cut;
          cuts = cuts * sqrt(cuts);
          drf = 1 / sqrt(l2ac + d2) - 0.5 * crf * d2 / cuts - dirf;
        } else if(coulomb_scaling && d_third[i].count(aj) && mj == mi){
          drf = (1.0 / 1.2) * 1 / d1 - 0.5 * crf * d2 / cut3 - dirf;
        } else
          drf = 1 / d1 - 0.5 * crf * d2 / cut3 - dirf;

        const double vdw = (c12 / d6 - c6) / d6;
        const double el = qq * drf * d_gff->fpepsi();

        // finally, check if atom a was also in d_as
        if (d_as->findAtom(mj, aj) != -1) {
          tmp_vdw += 0.5 * vdw;
          tmp_el += 0.5 * el;
        }
        // and store the energy in the correct array
        if (mj < 0) {
          d_vdw_s[i] += vdw;
          d_el_s[i] += el;
        } else {
          d_vdw_m[i] += vdw;
          d_el_m[i] += el;
        }
      }
      // now, loop over the exclusions (if requested)
      if(d_RFex){
	  
	for (set<int>::iterator iter=d_ex[i].begin(), 
	       to=d_ex[i].end(); iter!=to; ++iter){
	  utils::AtomSpecifier tmpas(*d_sys);
	  tmpas.addAtom(mi,*iter);
	  
	  //      for (int j = 0; j < d_pl[i].size(); j++) {
	  int mj = tmpas.mol(0);
	  int aj = tmpas.atom(0);
	  
	  // determine parameters
	  const double qq = qi * tmpas.charge(0);
	  
	  // now, we calculate the distance between atoms
	  gmath::Vec dd = d_pbc->nearestImage(vi, *tmpas.coord(0), d_sys->box());
	  
	  const double d1 = (vi - dd).abs();
	  const double d2 = d1*d1;
	  //double d6 = d2 * d2*d2;
	  double drf;
	  
	  // check if we have a soft atom
	  if (sft || d_soft->findAtom(mj, aj) != -1) {
	    double cuts = l2ac + d_cut*d_cut;
	    cuts = cuts * sqrt(cuts);
	    drf = - 0.5 * crf * d2 / cuts - dirf;
	  } else
	    drf = - 0.5 * crf * d2 / cut3 - dirf;
	  
	  const double el = qq * drf * d_gff->fpepsi();
	  
	  // finally, check if atom a was also in d_as
	  // this also includes the self term
	  if (d_as->findAtom(mj, aj) != -1) {
	    tmp_el += 0.5 * el;
	  }
	  // and store the energy in the correct array
	  d_el_m[i] += el;
	}
      }
    }
    d_p_vdw = tmp_vdw;
    d_p_el = tmp_el;
  }

  void Energy::calcCov() {
    // first calculate the values for all properties
    d_pc->calc();
    // loop over properties
    for (unsigned int i = 0; i < d_pc->size(); i++) {
      if((*d_pc)[i]->type() == "Distance"){
	d_cov[i] = calcBond((*d_pc)[i]->getValue().scalar(), d_covpar[i]);
      }
      else if((*d_pc)[i]->type() == "Angle"){
	d_cov[i] = calcAngle((*d_pc)[i]->getValue().scalar(), d_covpar[i]);
      }
      else if((*d_pc)[i]->type() == "Torsion" || 
	      (*d_pc)[i]->type() == "PeriodicTorsion"){
	//check if it is an improper based on the number of parameters
	if (d_covpar[i].size() == 2)
	  d_cov[i] = calcImproper((*d_pc)[i]->getValue().scalar(), d_covpar[i]);
	else
	  d_cov[i] = calcDihedral((*d_pc)[i]->getValue().scalar(), d_covpar[i]);
      }
      else if((*d_pc)[i]->type() == "CrossTorsion"){
	d_cov[i] = calcCrossDihedral((*d_pc)[i]->getValue().scalar(), d_covpar[i]);
      }
      else{
	throw Energy::Exception(
				" cannot compute covalent interaction energy for property " + (*d_pc)[i]->type());
	
          break;
      }
    }
  }

  void Energy::calcPair(int i, int j, double &vdw, double &el) {
    double qq, cuts, d, d1, d2, d6, drf, c6 = 0, c12 = 0;
    double cut3 = d_cut * d_cut*d_cut;
    double l2alj = d_lam * d_lam*d_alj;
    double l2ac = d_lam * d_lam*d_ac;
    double crf = ((2 - 2 * d_eps)*(1 + d_kap * d_cut) - d_eps * (d_kap * d_kap * d_cut * d_cut)) /
        ((1 + 2 * d_eps)*(1 + d_kap * d_cut) + d_eps * (d_kap * d_kap * d_cut * d_cut));
    double dirf = (1.0 - 0.5 * crf) / d_cut;

    int ai = d_as->atom(i);
    int gi = d_as->gromosAtom(i);
    int aj = d_as->atom(j);
    int gj = d_as->gromosAtom(j);
    int mi = d_as->mol(i);
    int mj = d_as->mol(j);
    int soft = 0;
    gmath::Vec dd;
    gmath::Vec chgrp1 = calcChgrp(i);
    gmath::Vec chgrp2 = calcChgrp(j);
    // check if one of the atoms is soft
    if (d_soft->findAtom(mi, ai) != -1 ||
        d_soft->findAtom(mj, aj) != -1) soft = 1;

    // calculate the distances between the chargegroups
    chgrp2 = d_pbc->nearestImage(chgrp1, chgrp2, d_sys->box());
    d = (chgrp2 - chgrp1).abs2();
    if (d <= d_cut * d_cut) {
      if (mi != mj || !d_ex[i].count(aj)) {
        //determine parameters
        gcore::LJType lj(d_gff->ljType(AtomPair(d_as->iac(i), d_as->iac(j))));
        qq = d_as->charge(i) * d_as->charge(j);

        // check third neighbour
        if (d_third[i].count(aj) && mj == mi) {
          c6 = lj.cs6();
          c12 = lj.cs12();
        } else {
          c6 = lj.c6();
          c12 = lj.c12();
        }
        // overwrite the LJ parameters in case of a LJ exception
        map<AtomPair, LJExceptionType>::const_iterator lje = d_gff->ljException().find(AtomPair(gi, gj));
        if(lje != d_gff->ljException().end()) {
          c6 = lje->second.c6();
          c12 = lje->second.c12();
        }

        // now, we calculate the distance between atoms
        dd = d_pbc->nearestImage(*(d_as->coord(i)),
            *(d_as->coord(j)),
            d_sys->box());
        d1 = (*d_as->coord(i) - dd).abs();
        d2 = d1*d1;
        d6 = d2 * d2*d2;
        if (soft) {
          if (c6 != 0 && c12 != 0) d6 += l2alj * c12 / c6;
          cuts = l2ac + d_cut*d_cut;
          cuts = cuts * sqrt(cuts);
          drf = 1 / sqrt(l2ac + d2) - 0.5 * crf * d2 / cuts - dirf;
        }
        // check third neighbours if coulomb_scaling is switched on
        else if(coulomb_scaling && d_third[i].count(aj) && mj == mi){
          drf = (1.0 / 1.2) * 1 / d1 - 0.5 * crf * d2 / cut3 - dirf;
        } else {
          drf = 1 / d1 - 0.5 * crf * d2 / cut3 - dirf;
        }

        vdw = (c12 / d6 - c6) / d6;
        el = qq * drf * d_gff->fpepsi();
      }
    }
  }

  int Energy::setAtoms(utils::AtomSpecifier &as) {
    d_ex.resize(0);
    d_third.resize(0);
    d_vdw_m.resize(0);
    d_vdw_s.resize(0);
    d_el_m.resize(0);
    d_el_s.resize(0);
    d_pl.resize(0);
    d_f_el_m.resize(0);
    d_f_el_s.resize(0);


    // CHRIS: this is a memory leak ???
    d_as = &as;
    // for all specified atoms, determine all excluded atoms and all third 
    // neighbours

    for (unsigned int i = 0; i < d_as->size(); i++) {
      if (d_as->atom()[i]->type() == spec_virtual)
        throw gromos::Exception("Energy", "Cannot calculate energy for a virtual atom");

      std::set<int> ex, third;
      int m = d_as->mol(i);
      int a = d_as->atom(i);
      if (m >= 0) {

        // first find atoms (<a) from which a is excluded
        for (int ai = 0; ai < a; ai++)
          for (int e = 0; e < d_sys->mol(m).topology().atom(ai).exclusion().size(); e++)
            if (a == d_sys->mol(m).topology().atom(ai).exclusion().atom(e))
              ex.insert(ai);
        // now, we add the exclusions of a itself
        for (int e = 0; e < d_sys->mol(m).topology().atom(a).exclusion().size(); e++)
          ex.insert(d_sys->mol(m).topology().atom(a).exclusion().atom(e));
        // and a is excluded of itself
        ex.insert(a);

        // first find atoms (<a) which have a as third neighbour
        for (int ai = 0; ai < a; ai++)
          for (int e = 0; e < d_sys->mol(m).topology().atom(ai).exclusion14().size();
              e++)
            if (a == d_sys->mol(m).topology().atom(ai).exclusion14().atom(e))
              third.insert(ai);
        // now, we add the third neighbours of a itself
        for (int e = 0; e < d_sys->mol(m).topology().atom(a).exclusion14().size(); e++)
          third.insert(d_sys->mol(m).topology().atom(a).exclusion14().atom(e));
      }

      // add things to the necessary vectors
      d_ex.push_back(ex);
      d_third.push_back(third);

      d_vdw_m.push_back(0.0);
      d_vdw_s.push_back(0.0);
      d_el_m.push_back(0.0);
      d_el_s.push_back(0.0);

      SimplePairlist spl(*d_sys, *d_pbc, d_cut);
      spl.setAtom(*d_as->atom()[i]);

      //  spl.setAtom(m,a);
      spl.setType("CHARGEGROUP");
      d_pl.push_back(spl);
      d_f_el_s.push_back(Vec(0.0, 0.0, 0.0));
      d_f_el_m.push_back(Vec(0.0, 0.0, 0.0));
    }
    return d_as->size();
  }

  int Energy::setProperties(utils::PropertyContainer &pc) {
    d_pc = &pc;
    for (unsigned int i = 0; i < d_pc->size(); i++) {
      std::vector<double> temp;
      if(pc[i]->type() == "Distance"){
	temp.resize(2);
	int t = findBond(*pc[i]);
	temp[0] = d_gff->bondType(t).b0();
	temp[1] = d_gff->bondType(t).fc();
      }
      else if(pc[i]->type() == "Angle") {
	temp.resize(2);
	int t = findAngle(*pc[i]);
	temp[0] = d_gff->angleType(t).t0();
	temp[1] = d_gff->angleType(t).fc();
      }
      else if(pc[i]->type() == "Torsion" ||
	      pc[i]->type() == "PeriodicTorsion") {

	// there can be more than one dihedral for one set of atoms
	std::vector<int> t = findDihedral(*pc[i]);

	// we use the Torsion property to specify both proper and improper 
	// dihedrals
	// Dirty fix to deal with both impropers and dihedrals:
	// if t < 0 this means it is an improper, with types counting
	// -1, -2, ...
	if (t[0] < 0) {
	  temp.resize(2);
	  t[0] = -1 * (t[0] + 1);
	  temp[0] = d_gff->improperType(t[0]).q0();
	  temp[1] = d_gff->improperType(t[0]).fc();
	} else {
	  temp.clear();
	  temp.push_back(t.size());
	  for (unsigned int i = 0; i < t.size(); ++i) {
	    temp.push_back(d_gff->dihedralType(t[i]).pd());
	    temp.push_back(d_gff->dihedralType(t[i]).np());
	    temp.push_back(d_gff->dihedralType(t[i]).fc());
	  }
	}
      }
      else if(pc[i]->type() =="CrossTorsion"){
	
	std::vector<int> t = findCrossDihedral(*pc[i]);
	temp.clear();
	temp.push_back(t.size());
	
	for (unsigned int i = 0; i < t.size(); ++i) {
	  temp.push_back(d_gff->dihedralType(t[i]).pd());
	  temp.push_back(d_gff->dihedralType(t[i]).np());
	  temp.push_back(d_gff->dihedralType(t[i]).fc());
	}
      }
      else {
	
	throw Energy::Exception(
				" cannot compute covalent interaction energy for property " + pc[i]->type() + " in " + pc[i]->toTitle());
      }
      d_covpar.push_back(temp);
      d_cov.push_back(0.0);

    }
    return d_pc->size();
  }

  gmath::Vec Energy::calcChgrp(int i) const{
    gmath::Vec chgrp(0.0, 0.0, 0.0);
    int mi = d_as->mol(i), ai = d_as->atom(i);
    if (mi < 0) {
      int nsa = d_sys->sol(0).topology().numAtoms();
      int solv = ai / nsa;
      solv *= nsa;
      return d_sys->sol(0).pos(solv);
    }

    int begin = ai - 1, end = ai;
    if (ai > 0)
      for (begin = ai - 1;
          begin >= 0 && d_sys->mol(mi).topology().atom(begin).chargeGroup() != 1;
          begin--);
    for (end = ai;
        d_sys->mol(mi).topology().atom(end).chargeGroup() != 1;
        end++);

    // charge group goes from begin+1 to end
    for (int k = begin + 1; k <= end; k++)
      chgrp += d_sys->mol(mi).pos(k);
    return chgrp / (end - begin);
  }

  int Energy::findBond(utils::Property &pp) const{
    int m, a, b, f = 0;
    if (pp.atoms().mol(0) == pp.atoms().mol(1))
      m = pp.atoms().mol(0);
    else
      throw Energy::Exception(
        " Covalent interactions are always within one molecule: " + pp.toTitle());
    if (pp.atoms().atom(0) < pp.atoms().atom(1)) {
      a = pp.atoms().atom(0);
      b = pp.atoms().atom(1);
    } else {
      a = pp.atoms().atom(1);
      b = pp.atoms().atom(0);
    }
    BondIterator bi(d_sys->mol(m).topology());
    while (bi && f == 0)
      if (bi()[0] == a && bi()[1] == b) f = 1;
      else ++bi;
    if (bi) return bi().type();
    else
      throw Energy::Exception(
        " Bond not found in topology: " + pp.toTitle());
  }

  int Energy::findAngle(utils::Property &pp) const{
    int m, a, b, c, f = 0;
    if (pp.atoms().mol(0) == pp.atoms().mol(1) && pp.atoms().mol(0) == pp.atoms().mol(2))
      m = pp.atoms().mol(0);
    else
      throw Energy::Exception(
        " Covalent interactions are always within one molecule: " + pp.toTitle());
    if (pp.atoms().atom(0) < pp.atoms().atom(2)) {
      a = pp.atoms().atom(0);
      b = pp.atoms().atom(1);
      c = pp.atoms().atom(2);
    } else {
      a = pp.atoms().atom(2);
      b = pp.atoms().atom(1);
      c = pp.atoms().atom(0);
    }
    AngleIterator ai(d_sys->mol(m).topology());
    while (ai && f == 0)
      if (ai()[0] == a && ai()[1] == b && ai()[2] == c) f = 1;
      else ++ai;
    if (ai) return ai().type();
    else
      throw Energy::Exception(
        " Angle not found in topology: " + pp.toTitle());
  }

  std::vector<int> Energy::findDihedral(utils::Property &pp) const{
    std::vector<int> result;
    int m, a, b, c, d;
    if (pp.atoms().mol(0) == pp.atoms().mol(1) &&
        pp.atoms().mol(0) == pp.atoms().mol(2) &&
        pp.atoms().mol(0) == pp.atoms().mol(3))
      m = pp.atoms().mol(0);
    else
      throw Energy::Exception(
        " Covalent interactions are always within one molecule: " + pp.toTitle());
    if (pp.atoms().atom(1) < pp.atoms().atom(2)) {
      a = pp.atoms().atom(0);
      b = pp.atoms().atom(1);
      c = pp.atoms().atom(2);
      d = pp.atoms().atom(3);
    } else {
      a = pp.atoms().atom(3);
      b = pp.atoms().atom(2);
      c = pp.atoms().atom(1);
      d = pp.atoms().atom(0);
    }
    DihedralIterator di(d_sys->mol(m).topology());
    for (; di; ++di) {
      if (di()[0] == a && di()[1] == b && di()[2] == c && di()[3] == d)
        result.push_back(di().type());
    }

    //Maybe we have an improper
    ImproperIterator ii(d_sys->mol(m).topology());
    for (; ii; ++ii) {
      if (ii()[0] == a && ii()[1] == b && ii()[2] == c && ii()[3] == d) {
        result.push_back(-1 * (ii().type() + 1));
      }
    }

    if (result.empty())
      throw Energy::Exception("(improper) Dihedral not found in topology: " + pp.toTitle());

    return result;
  }

  std::vector<int> Energy::findCrossDihedral(utils::Property &pp) const{
    std::vector<int> result;
    const utils::CrossTorsionProperty & p = (const utils::CrossTorsionProperty &) pp;
    int m;
    if (p.atoms().mol(0) == p.atoms().mol(1) &&
        p.atoms().mol(0) == p.atoms().mol(2) &&
        p.atoms().mol(0) == p.atoms().mol(3) &&
        p.atoms().mol(0) == p.atoms2().mol(0) &&
        p.atoms().mol(0) == p.atoms2().mol(1) &&
        p.atoms().mol(0) == p.atoms2().mol(2) &&
        p.atoms().mol(0) == p.atoms2().mol(3))
      m = p.atoms().mol(0);
    else
      throw Energy::Exception(
        " Covalent interactions are always within one molecule: " + p.toTitle());

    CrossDihedral dih(p.atoms().atom(0), p.atoms().atom(1), p.atoms().atom(2), p.atoms().atom(3),
            p.atoms2().atom(0), p.atoms2().atom(1), p.atoms2().atom(2), p.atoms2().atom(3));

    CrossDihedralIterator di(d_sys->mol(m).topology());
    for (; di; ++di) {
      if (di()[0] == dih[0] &&
              di()[1] == dih[1] &&
              di()[2] == dih[2] &&
              di()[3] == dih[3] &&
              di()[4] == dih[4] &&
              di()[5] == dih[5] &&
              di()[6] == dih[6] &&
              di()[7] == dih[7])
        result.push_back(di().type());
    }

    if (result.empty())
      throw Energy::Exception("Cross Dihedral not found in topology: " + p.toTitle());

    return result;
  }

  double Energy::calcBond(double val, const std::vector<double> & par) const{
    double diff = val * val - par[0] * par[0];
    return 0.25 * par[1] * diff*diff;
  }

  double Energy::calcAngle(double val, const std::vector<double> & par) const{
    val = val * M_PI / 180.0;
    double t0 = par[0] * M_PI / 180.0;
    double diff = cos(val) - cos(t0);
    return 0.5 * par[1] * diff*diff;
  }

  double Energy::calcDihedral(double val, const std::vector<double> & par) const{
    val = val * M_PI / 180.0;
    double e = 0.0;

    unsigned int num = int(par[0]);
    for (unsigned int i = 0; i < num; ++i)
      e += par[1 + 3*i + 2]*(1 + par[1 + 3*i] * cos(par[1 + 3*i + 1] * val));
    return e;
  }

  double Energy::calcCrossDihedral(double val, const std::vector<double> & par) const{
    val = val * M_PI / 180.0;
    double e = 0.0;

    unsigned int num = int(par[0]);
    for (unsigned int i = 0; i < num; ++i)
      e += par[1 + 3*i + 2]*(1 + par[1 + 3*i] * cos(par[1 + 3*i + 1] * val));
    return e;
  }

  double Energy::calcImproper(double val, const std::vector<double> & par) const{
    if (val > 180.0) val = val - 360;
    double diff = val - par[0];
    return 0.5 * par[1] * diff*diff;
  }

  double Energy::cov() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_pc->size(); i++)
      e += this->cov(i);
    return e;
  }

  /*
  double Energy::dist() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_pc->size(); i++)
      {
	if((*d_pc)[i]->type() == "Distance"){
	  e += this->cov(i);
	}
      }
    return e;
  }

  double Energy::angle() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_pc->size(); i++)
      {
	if((*d_pc)[i]->type() == "Angle"){
	  e += this->cov(i);
	}
      }
    return e;
  }

  double Energy::impdihed() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_pc->size(); i++)
      {
	if((*d_pc)[i]->type() == "Torsion" ||
	   (*d_pc)[i]->type() == "PeriodicTorsion"){
	  if (d_covpar[i].size() == 2) {
	    e += this->cov(i);
	  }
	}
      }
    return e;
  }
  
  double Energy::torsdihed() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_pc->size(); i++)
      {
	if((*d_pc)[i]->type() == "Torsion" ||
	   (*d_pc)[i]->type() == "PeriodicTorsion"){
	  if (d_covpar[i].size() != 2) {
	    e += this->cov(i);
	  }
	}
	else if((*d_pc)[i]->type() == "CrossTorsion"){
	  e += this->cov(i);
	}
      }
    return e;
  }
  */
  
  double Energy::vdw() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_as->size(); i++)
      e += this->vdw(i);
    return e - d_p_vdw;
  }

  double Energy::el() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_as->size(); i++)
      e += this->el(i);
    return e - d_p_el;
  }
  // this is new
  double Energy::vdw_m() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_as->size(); i++)
      e += d_vdw_m[i];
    return e - d_p_vdw; // need so subtract the couple counting of interactions
  }

  // this is new
  double Energy::vdw_s() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_as->size(); i++)
      e += d_vdw_s[i];
    return e;
  }

  // this is new
  double Energy::el_m() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_as->size(); i++)
      e += d_el_m[i];
    return e - d_p_el; // need so subtract the couple counting of interactions
  }

  // this is new
  double Energy::el_s() const{
    double e = 0.0;
    for (unsigned int i = 0; i < d_as->size(); i++)
      e += d_el_s[i];
    return e;
  }
  
  double Energy::vdw(unsigned int i) const{
    assert(i < d_as->size());
    return d_vdw_m[i] + d_vdw_s[i];
  }

  double Energy::el(unsigned int i) const{
    assert(i < d_as->size());
    return d_el_m[i] + d_el_s[i];
  }

  double Energy::cov(unsigned int i) const{
    assert(i < d_pc->size());
    return d_cov[i];
  }

  double Energy::vdw_m(unsigned int i) const{
    assert(i < d_as->size());
    return d_vdw_m[i];
  }

  double Energy::vdw_s(unsigned int i) const{
    assert(i < d_as->size());
    return d_vdw_s[i];
  }

  double Energy::el_m(unsigned int i) const{
    assert(i < d_as->size());
    return d_el_m[i];
  }

  double Energy::el_s(unsigned int i) const{
    assert(i < d_as->size());
    return d_el_s[i];
  }

  gmath::Vec Energy::f_el(unsigned int i) const{
    assert(i < d_as->size());
    return d_f_el_s[i] + d_f_el_m[i];
  }

  gmath::Vec Energy::f_el_m(unsigned int i) const{
    assert(i < d_as->size());
    return d_f_el_m[i];
  }

  gmath::Vec Energy::f_el_s(unsigned int i) const{
    assert(i < d_as->size());
    return d_f_el_s[i];
  }

  void Energy::calcPairlist() {
    if (d_pl.size() != d_as->size())
      throw Energy::Exception(
        " Cannot calculate pairlist without setting atoms first");
    const int size = d_pl.size();
#ifdef OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < size; ++i) {
      d_pl[i].setCutOff(d_cut);
      d_pl[i].clear();
      d_pl[i].calc();
      d_pl[i].removeExclusions();
    }
  }

  void Energy::setPairlist(int i, SimplePairlist &as) {
    assert(i < int(d_pl.size()));
    d_pl[i].clear();
    d_pl[i] = as;
    d_pl[i].removeExclusions();
  }

  void Energy::setPairlistType(string t) {
    if (d_pl.size() != d_as->size() || d_pl.size() == 0)
      throw Energy::Exception(
        " Cannot set pairlist type, without setting atoms first");
    for (unsigned int i = 0; i < d_pl.size(); ++i)
      d_pl[i].setType(t);
  }

  void Energy::setRFexclusions(bool p){
    d_RFex=p;
  }
  
  void Energy::setCoulombScaling(bool p){
    coulomb_scaling = p;
  }
  
}

