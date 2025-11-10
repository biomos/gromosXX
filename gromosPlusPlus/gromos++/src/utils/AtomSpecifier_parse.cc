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
#include "AtomSpecifier.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "VirtualAtom.h"
#include "parse.h"
#include "ExpressionParser.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gio/Ginstream.h"
#include "../bound/Boundary.h"


using namespace gcore;
using namespace std;
using utils::AtomSpecifier;

void AtomSpecifier::parse(std::string s, int x)
{
  std::string::size_type it = find_par(s, ';');

  if (it == std::string::npos){
    parse_single(s, x);
  }
  else{
    // std::cerr << "AS::parse\tx = " << x << std::endl;
    // std::cerr << "\tparsing " << s.substr(0, it) << std::endl;
    parse_single(s.substr(0, it), x);
    // std::cerr << "\tparsing " << s.substr(it+1, std::string::npos) << std::endl;
    parse(s.substr(it+1, std::string::npos), x);
  }

  if (d_not_atoms != NULL) {
    // remove atoms
    for (unsigned int i = 0; i < d_not_atoms->size(); ++i) {
      removeAtom(d_not_atoms->mol(i), d_not_atoms->atom(i));
    }
    // and solvent types
    vector<int> new_solv;
    for(unsigned int i = 0; i < d_solventType.size(); ++i) {
      if (std::find(d_not_atoms->d_solventType.begin(),
              d_not_atoms->d_solventType.end(),
              d_solventType[i]) == d_not_atoms->d_solventType.end())
        new_solv.push_back(d_solventType[i]);
    }
    d_solventType = new_solv;
    // check for specific atoms violation this
    for(unsigned int i = 0; i < size(); i++) {
      if (mol(i) >= 0) continue;

      if (std::find(d_not_atoms->d_solventType.begin(),
              d_not_atoms->d_solventType.end(),
          atom(i) % d_sys->sol(0).topology().numAtoms()) != d_not_atoms->d_solventType.end())
        removeAtom(i);
    }
  }
}


void AtomSpecifier::parseStrict(std::string s, int x)
{
  std::string::size_type it = find_par(s, ';');

  if (it == std::string::npos){
    parse_singleStrict(s, x);
  }
  else{
    // std::cerr << "AS::parse\tx = " << x << std::endl;
    // std::cerr << "\tparsing " << s.substr(0, it) << std::endl;
    parse_singleStrict(s.substr(0, it), x);
    // std::cerr << "\tparsing " << s.substr(it+1, std::string::npos) << std::endl;
    parseStrict(s.substr(it+1, std::string::npos), x);
  }

  /*
  if (d_not_atoms != NULL) {
    // remove atoms
    for (int i = 0; i < d_not_atoms->size(); ++i) {
      //removeAtom(d_not_atoms->mol(i), d_not_atoms->atom(i));
      std::cout << " here to remove atom " << d_not_atoms->mol(i) << ":" << d_not_atoms->atom(i)<< endl;
    }
    // and solvent types
    vector<int> new_solv;
    for(unsigned int i = 0; i < d_solventType.size(); ++i) {
      if (std::find(d_not_atoms->d_solventType.begin(),
              d_not_atoms->d_solventType.end(),
              d_solventType[i]) == d_not_atoms->d_solventType.end())
        new_solv.push_back(d_solventType[i]);
    }
    d_solventType = new_solv;
    // check for specific atoms violation this
    for(int i = 0; i < size(); i++) {
      if (mol(i) >= 0) continue;

      if (std::find(d_not_atoms->d_solventType.begin(),
              d_not_atoms->d_solventType.end(),
          atom(i) % d_sys->sol(0).topology().numAtoms()) != d_not_atoms->d_solventType.end())
        removeAtom(i);
    }
  } */
}

void AtomSpecifier::parse_singleStrict(std::string s, int x) {
  if (s.substr(0, 3) == "va(") {
      //std::cout << std::endl << s << "\t" << "va(" << std::endl;
    std::string::size_type it = find_matching_bracket(s, '(', 3);
    parse_va(s.substr(3, it - 4), x);
  } else if (s.substr(0, 5) == "file(") {
      //std::cout << std::endl << s << "\t" << "file" << std::endl;
    std::string::size_type it = find_matching_bracket(s, '(', 5);
    parse_atominfoStrict(s.substr(5, it - 6));
  } else if (s.substr(0, 4) == "not(") {
      //std::cout << std::endl << s << "\t" << "not(" << std::endl;
    std::string::size_type it = find_matching_bracket(s, '(', 4);
    if (d_not_atoms == NULL)
      d_not_atoms = new AtomSpecifier(*d_sys);
    d_not_atoms->addSpecifierStrict(s.substr(4, it - 5),x);
  } else if (s.substr(0, 6) == "minus(") {
      //std::cout << std::endl << s << "\t" << "minus(" << std::endl;
    std::string::size_type it = find_matching_bracket(s, '(', 6);
    parse_minus(s.substr(6, it - 7),x);
  } else {
      //std::cout << std::endl << s << "\t" << "else" << std::endl;
    std::string::size_type it = s.find(':');
    if (it == std::string::npos)
      throw(Exception("no atoms in AtomSpecifier"));

    std::vector<int> mol;
    parse_moleculeStrict(s.substr(0, it), mol, x);

    for(std::vector<int>::const_iterator it = mol.begin(), to = mol.end();
        it != to; ++it) {
      // check for too high and low molecule numbers
      if (*it < 0 || *it > d_sys->numMolecules()) {
        std::ostringstream os;
        os << "Molecule " << *it << " does not exist.";
        throw Exception(os.str());
      }
    }

    for(unsigned int i=0; i<mol.size(); ++i){
      parse_atomStrict(mol[i], s.substr(it+1, std::string::npos), x);
    }

  }
}


void AtomSpecifier::parse_single(std::string s, int x) {

  if (s == "no")
    return;

  if (s.substr(0, 3) == "va(") {
    std::string::size_type it = find_matching_bracket(s, '(', 3);
    parse_va(s.substr(3, it - 4), x);
  } else if (s.substr(0, 5) == "file(") {
    std::string::size_type it = find_matching_bracket(s, '(', 5);
    parse_atominfo(s.substr(5, it - 6));
  } else if (s.substr(0, 4) == "not(") {
    std::string::size_type it = find_matching_bracket(s, '(', 4);
    if (d_not_atoms == NULL)
      d_not_atoms = new AtomSpecifier(*d_sys);
    d_not_atoms->addSpecifier(s.substr(4, it - 5),x);
  } else if (s.substr(0, 6) == "minus(") {
    std::string::size_type it = find_matching_bracket(s, '(', 6);
    parse_minus(s.substr(6, it - 7),x);
  } else {
    std::string::size_type it = s.find(':');
    if (it == std::string::npos)
      throw(Exception("no atoms in AtomSpecifier"));
    
    std::vector<int> mol;
    parse_molecule(s.substr(0, it), mol, x);
    
    for(std::vector<int>::const_iterator it = mol.begin(), to = mol.end();
        it != to; ++it) {
      // check for too high and low molecule numbers
      if (*it < 0 || *it > d_sys->numMolecules()) {
        std::ostringstream os;
        os << "Molecule " << *it << " does not exist.";
        throw Exception(os.str());
      }
    }    
    
    for(unsigned int i=0; i<mol.size(); ++i)
      parse_atom(mol[i], s.substr(it+1, std::string::npos), x);
  }
}

void AtomSpecifier::parse_molecule(std::string s, std::vector<int> & mol, int x)
{
  if (s=="a"){
    for(int i=0; i < d_sys->numMolecules(); ++i)
      mol.push_back(i+1);
  }
  else if(s=="s"){
    mol.push_back(0);
  }
  else if (s=="x"){
    mol.push_back(x);
  }
  else{
    // parse_mol_range(s, mol);
    // use function from parse
    parse_range(s, mol, x);
  }
}

void AtomSpecifier::parse_moleculeStrict(std::string s, std::vector<int> & mol, int x)
{
  if (s=="a"){
    for(int i=0; i < d_sys->numMolecules(); ++i)
      mol.push_back(i+1);
  }
  else if(s=="s"){
    mol.push_back(0);
  }
  else if (s=="x"){
    mol.push_back(x);
  }
  else{
    // parse_mol_range(s, mol);
    // use function from parse
    //parse_range(s, mol, x);
      //char l=char(s);
      int m=atoi(s.c_str());
      //std::cout << "count mol " << m << "\t" << s.c_str() << "\t" << s << endl;
      mol.push_back(m);
  }
}

void AtomSpecifier::parse_atom(int mol, std::string s, int x)
{
  if (s.substr(0, 4) == "res("){
    std::string::size_type ket = find_matching_bracket(s, '(', 4);
    if (ket == std::string::npos)
      throw Exception("Residue: end bracket missing");

    std::string res = s.substr(4, ket - 5);
    
    std::string::size_type sep = res.find(':');
    if (sep == std::string::npos)
      throw Exception("No atoms for residue given");
  
    std::string resn = res.substr(0, sep);
    std::string atom = res.substr(sep+1, std::string::npos);

    parse_res(mol, resn, atom, x);
  }
  else{
    if (mol > 0) {
      parse_atom_range(mol, 0, d_sys->mol(mol-1).numAtoms(), s, x);
    } else
      parse_atom_range(mol, 0, d_sys->sol(0).numAtoms(), s, x);
  }
}

void AtomSpecifier::parse_atomStrict(int mol, std::string s, int x)
{
  if (s.substr(0, 4) == "res("){
    std::string::size_type ket = find_matching_bracket(s, '(', 4);
    if (ket == std::string::npos)
      throw Exception("Residue: end bracket missing");

    std::string res = s.substr(4, ket - 5);

    std::string::size_type sep = res.find(':');
    if (sep == std::string::npos)
      throw Exception("No atoms for residue given");

    std::string resn = res.substr(0, sep);
    std::string atom = res.substr(sep+1, std::string::npos);

    parse_resStrict(mol, resn, atom, x);
  }
  else{
    if (mol > 0) {
      parse_atom_rangeStrict(mol, 0, d_sys->mol(mol-1).numAtoms(), s, x);
    } else
      parse_atom_rangeStrict(mol, 0, d_sys->sol(0).numAtoms(), s, x);
  }
}

void AtomSpecifier::parse_va(std::string s, int x)
{
  std::string::size_type it = s.find(',');
  if (it == std::string::npos)
    throw Exception(" Virtual Atom: wrong format. Should be va(type, as)");
  
  std::string t = s.substr(0, it);
  VirtualAtom::virtual_type vt;

  // try the types
  if (t == "com") vt = VirtualAtom::COM;
  else if (t == "cog") vt = VirtualAtom::COG;
  else{
    std::istringstream is(t);
    int i;
    if (!(is >> i))
      throw Exception(" Virtual Atom: type not recognised: " + t);
    vt = VirtualAtom::virtual_type(i);
  }
  d_specatom.push_back(new VirtualSpecAtom(*d_sys, 
          s.substr(it+1, std::string::npos), x, vt));
}

void AtomSpecifier::parse_minus(std::string s, int x) {
  AtomSpecifier minus(*d_sys, s, x);
  for(unsigned int i = 0; i < minus.size(); ++i) {
    if (minus.mol(i) < 0)
      throw Exception("Solvent is not supported in the minus() atom specifier. Use not() instead.");
    removeAtom(minus.mol(i), minus.atom(i));
  }
}


void AtomSpecifier::parse_atominfo(std::string s)
{
  // std::cerr << "trying to parse atominfo file" << std::endl;
  
  std::ifstream aif(s.c_str());
  if (!aif.is_open()){
    throw Exception("could not open atominfo file");
  }
  
  gio::Ginstream ai(aif);
  
  std::vector<std::string> buffer;
  ai.getblock(buffer);

  if (!buffer.size() || buffer[0] != "ATOMS"){
    std::ostringstream os;
    os << "no ATOMS block found in " << s << "!";
    throw Exception(os.str());
  }
  
  for(unsigned int i=1; i<buffer.size()-1; ++i){
    std::string s = buffer[i];
    std::string::size_type beg = s.find_first_not_of(" ");
    std::string::size_type c = s.find(':');

    int a,m;
    std::string mol_string = s.substr(beg, c-beg);
    std::istringstream is(mol_string);

    if (mol_string == "s"){
      m = 0;
    }
    else{
      if (!(is >> m)){
	std::ostringstream os;
	os << "Could not parse line: " << buffer[i] << ": trying to read mol from " << s.substr(0,c);
	throw Exception(os.str());
      }
    }
    
    is.clear();
    is.str(s.substr(c+1, std::string::npos));
    if (!(is >> a)){
      std::ostringstream os;
      os << "Could not parse line: " << buffer[i] << ": trying to read atom from " << s.substr(c+1, std::string::npos);
      throw Exception(os.str());
    }
    
    // std::cerr << "trying to add " << m << ":" << a << std::endl;
    addAtom(m-1, a-1);
  }
}

void AtomSpecifier::parse_atominfoStrict(std::string s)
{
  // std::cerr << "trying to parse atominfo file" << std::endl;

  std::ifstream aif(s.c_str());
  if (!aif.is_open()){
    throw Exception("could not open atominfo file");
  }

  gio::Ginstream ai(aif);

  std::vector<std::string> buffer;
  ai.getblock(buffer);

  if (!buffer.size() || buffer[0] != "ATOMS"){
    std::ostringstream os;
    os << "no ATOMS block found in " << s << "!";
    throw Exception(os.str());
  }

  for(unsigned int i=1; i<buffer.size()-1; ++i){
    std::string s = buffer[i];
    std::string::size_type beg = s.find_first_not_of(" ");
    std::string::size_type c = s.find(':');

    int a,m;
    std::string mol_string = s.substr(beg, c-beg);
    std::istringstream is(mol_string);

    if (mol_string == "s"){
      m = 0;
    }
    else{
      if (!(is >> m)){
	std::ostringstream os;
	os << "Could not parse line: " << buffer[i] << ": trying to read mol from " << s.substr(0,c);
	throw Exception(os.str());
      }
    }

    is.clear();
    is.str(s.substr(c+1, std::string::npos));
    if (!(is >> a)){
      std::ostringstream os;
      os << "Could not parse line: " << buffer[i] << ": trying to read atom from " << s.substr(c+1, std::string::npos);
      throw Exception(os.str());
    }

    // std::cerr << "trying to add " << m << ":" << a << std::endl;
    //std::cout << i << "\t" << buffer[i] << "\t" << m-1 << "\t" << a-1 << endl;
    addAtomStrict(m-1, a-1);
  }
}

void AtomSpecifier::parse_atom_range(int mol, int beg, int end, std::string s, int x)
{  
  std::map<std::string, int> var;
  var["x"] = x;
  
  ExpressionParser<int> ep;
  
  if (find_par(s, ':') != std::string::npos)
      throw Exception("Unexpected ':' token in atom set/range parsing.");

  std::string::size_type it = find_par(s, ',');
  
  if (it == std::string::npos){
    
    std::string::size_type r_it = find_par(s, '-');

    if (r_it == std::string::npos || !r_it){

      if (s == "a"){
	if(mol >0)
	  for(int i=beg; i<end; ++i)
	    addAtom(mol-1, i);
	else
	  for(int i=0; i<d_sys->sol(0).topology().numAtoms(); ++i)
	    addSolventType(d_sys->sol(0).topology().atom(i).name());
      }
      else{
	std::string::size_type bra = s.find_first_not_of(" ");
	int i;
	if (s[bra] == '('){
	  i = ep.parse_expression(s, var);
	  addAtom(mol-1, beg+i-1);
	}
	else{
	  // single number (or type)
	  std::istringstream is(s);
	  if(!(is >> i)){
	    if(mol > 0)
	      addType(mol-1, s, beg, end);
	    else
	      addType(mol-1, s);
	  }
	  else{
	    if (((beg + i) > end || (beg + i) < 1)  && mol > 0)
	      throw Exception("Atom out of range");	  
	    addAtom(mol-1, beg+i-1);
	  }
	}
      }
    }
    else{
      int beg, end;
      beg = ep.parse_expression(s.substr(0, r_it), var);
      end = ep.parse_expression(s.substr(r_it+1, std::string::npos), var);
      for(int i=beg; i<=end; ++i){
	addAtom(mol-1, i-1);
      }
    }
  }
  else{
    parse_atom_range(mol, beg, end, s.substr(0, it), x);
    parse_atom_range(mol, beg, end, s.substr(it+1, std::string::npos), x);
  }

}

void AtomSpecifier::parse_atom_rangeStrict(int mol, int beg, int end, std::string s, int x)
{
  std::map<std::string, int> var;
  var["x"] = x;

  ExpressionParser<int> ep;

  if (find_par(s, ':') != std::string::npos)
      throw Exception("Unexpected ':' token in atom set/range parsing.");

  std::string::size_type it = find_par(s, ',');

  if (it == std::string::npos){

    std::string::size_type r_it = find_par(s, '-');

    if (r_it == std::string::npos || !r_it){

      if (s == "a"){
	if(mol >0)
	  for(int i=beg; i<end; ++i)
	    addAtomStrict(mol-1, i);
	else
	  for(int i=0; i<d_sys->sol(0).topology().numAtoms(); ++i)
	    addSolventType(d_sys->sol(0).topology().atom(i).name());
      }
      else{
	std::string::size_type bra = s.find_first_not_of(" ");
	int i;
	if (s[bra] == '('){
	  i = ep.parse_expression(s, var);
	  addAtomStrict(mol-1, beg+i-1);
	}
	else{
            // DW
	  std::istringstream is(s);
	  if(!(is >> i)){
	    if(mol > 0)
                for(int j=beg; j<end; ++j){
                    if(_checkName(mol-1, j, s)){
                         addAtomStrict(mol-1,j);
                    }
                }
	    else
	      addTypeStrict(mol-1, s);
	  }
	  else{
	    if (((beg + i) > end || (beg + i) < 1)  && mol > 0)
	      throw Exception("Atom out of range");
	    addAtomStrict(mol-1, beg+i-1);
	  }
	}
      }
    }
    else{
      int beg, end;
      beg = ep.parse_expression(s.substr(0, r_it), var);
      end = ep.parse_expression(s.substr(r_it+1, std::string::npos), var);
      for(int i=beg; i<=end; ++i){
	addAtomStrict(mol-1, i-1);
      }
    }
  }
  else{
    parse_atom_rangeStrict(mol, beg, end, s.substr(0, it), x);
    parse_atom_rangeStrict(mol, beg, end, s.substr(it+1, std::string::npos), x);
  }

}

void AtomSpecifier::parse_res(int mol, std::string res, std::string atom, int x)
{
  if (mol<=0) throw Exception("No residues in solvent");

  std::map<std::string, int> var;
  var["x"] = x;
  
  ExpressionParser<int> ep;

  // std::string::size_type it = res.find(',');
  std::string::size_type it = find_par(res, ',');

  if (it == std::string::npos){
    
    // std::string::size_type r_it = res.find('-');
    std::string::size_type r_it = find_par(res, '-');

    if (r_it == std::string::npos || !r_it){
      std::string::size_type bra = res.find_first_not_of(" ");
      int i;
      if (res[bra] == '('){
	i = ep.parse_expression(res, var);
	int beg, end;
	res_range(mol, i, beg, end);
	parse_atom_range(mol, beg, end, atom, x);
      }
      else{
	// single number (or type)
	std::istringstream is(res);
	int i;
	if(!(is >> i)){
	  parse_res_type(mol, res, atom, x);
	}
	else{
	  int beg, end;
	  res_range(mol, i, beg, end);
	  parse_atom_range(mol, beg, end, atom, x);
	}
      }
    }
    else{
      int beg, end;
      beg = ep.parse_expression(res.substr(0, r_it), var);
      end = ep.parse_expression
	(res.substr(r_it+1, std::string::npos),	var);

      for(int i=beg; i<=end; ++i){
	int rbeg, rend;
	res_range(mol, i, rbeg, rend);
	parse_atom_range(mol, rbeg, rend, atom, x);
      }
    }
  }
  else{
    parse_res(mol, res.substr(0, it), atom, x);
    parse_res(mol, res.substr(it+1, std::string::npos), atom, x);
  }
}

void AtomSpecifier::parse_resStrict(int mol, std::string res, std::string atom, int x)
{
  if (mol<=0) throw Exception("No residues in solvent");

  std::map<std::string, int> var;
  var["x"] = x;

  ExpressionParser<int> ep;

  // std::string::size_type it = res.find(',');
  std::string::size_type it = find_par(res, ',');

  if (it == std::string::npos){

    // std::string::size_type r_it = res.find('-');
    std::string::size_type r_it = find_par(res, '-');

    if (r_it == std::string::npos || !r_it){
      std::string::size_type bra = res.find_first_not_of(" ");
      int i;
      if (res[bra] == '('){
	i = ep.parse_expression(res, var);
	int beg, end;
	res_range(mol, i, beg, end);
	parse_atom_rangeStrict(mol, beg, end, atom, x);
      }
      else{
	// single number (or type)
	std::istringstream is(res);
	int i;
	if(!(is >> i)){
	  parse_res_type(mol, res, atom, x);
	}
	else{
	  int beg, end;
	  res_range(mol, i, beg, end);
	  parse_atom_rangeStrict(mol, beg, end, atom, x);
	}
      }
    }
    else{
      int beg, end;
      beg = ep.parse_expression(res.substr(0, r_it), var);
      end = ep.parse_expression
	(res.substr(r_it+1, std::string::npos),	var);

      for(int i=beg; i<=end; ++i){
	int rbeg, rend;
	res_range(mol, i, rbeg, rend);
	parse_atom_rangeStrict(mol, rbeg, rend, atom, x);
      }
    }
  }
  else{
    parse_resStrict(mol, res.substr(0, it), atom, x);
    parse_resStrict(mol, res.substr(it+1, std::string::npos), atom, x);
  }
}

void AtomSpecifier::res_range(int mol, int res, int &beg, int &end)
{
  beg = d_sys->mol(mol-1).numAtoms();
  end = 0;
  --res;

  for(int i=0; i<d_sys->mol(mol-1).numAtoms(); ++i){
    if (d_sys->mol(mol-1).topology().resNum(i) == res){
      if (i < beg) beg = i;
      if (i > end) end = i;
    }
  }
  ++end;
}

void AtomSpecifier::parse_res_type(int mol, std::string res, std::string s, int x)
{
  int beg = 0;
  bool match=false;
  int resn = 0;

  for(int i=0; i<d_sys->mol(mol-1).numAtoms(); ++i){
    
    if (_checkResName(mol-1, i, res)){
      if(!match){
	beg = i;
	match = true;
	resn = d_sys->mol(mol-1).topology().resNum(i);
      }
      else if(resn != d_sys->mol(mol-1).topology().resNum(i)){
	parse_atom_range(mol, beg, i, s, x);
	beg = i;
	resn = d_sys->mol(mol-1).topology().resNum(i);
      }
    }
    else{
      if (match){
	match = false;
	parse_atom_range(mol, beg, i, s, x);
      }
    }
  }
  if (match){
    match = false;
    parse_atom_range(mol, beg, d_sys->mol(mol-1).numAtoms(), s, x);
  }
}


