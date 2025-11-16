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
#include "SoluteWeightedDistance.h"

#include <cmath>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

#include "AtomSpecifier.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../gio/InTopology.h"
#include "../args/BoundaryParser.h"
#include "../bound/Boundary.h"
#include "../gcore/System.h"


static const int fgIndex = 0;
static const int cgIndex = 1;

namespace utils {
  SoluteWeightedDistance::SoluteWeightedDistance(gcore::System& sys, args::Arguments& args) :
          _solute(sys),
          _fgSolvent(sys),
          _cgSolvent(sys),
          _withEnergy(false),
          _withMeasures(false),
          _withWeights(false),
          _withMinimalDistances(false),
          _sys(sys)
  {
    if ( !(args.count("solute") > 0 &&
         ( args.count("fgsolv") > 0 ||
           args.count("cgsolv") > 0) ) ) {
      throw gromos::Exception("SoluteWeightedDistance",
              "You must specify solute and at least fg or cg solvent");
    }
    
    std::string keywords[] = {"solute", "fgsolv", "cgsolv"};
    AtomSpecifier * const as[] = {&_solute, &_fgSolvent, &_cgSolvent};

    for (unsigned int i = 0; i < 3; i++) {
      for (args::Arguments::const_iterator iter = args.lower_bound(keywords[i]),
              to = args.upper_bound(keywords[i]);
              iter != to; ++iter) {
        as[i]->addSpecifier(iter->second);
      }
    }
    _pbc = args::BoundaryParser::boundary(sys, args);
    
    // Should we also calculate the energy?
    if (args.count("energy") == 5){
      args::Arguments::const_iterator iter = args.lower_bound("energy");
      for (unsigned int i = 0; i < 2; i++){
        double params[] = {0.0, 0.0};
        for (unsigned int j = 0; j < 2; j++, iter++){
          std::istringstream is(iter->second);
          is >> params[j];
        }
        SWD_Param p(params[0], params[1]); // first force constant, then cutoff
        _params.push_back(p);
        _energyFile = iter->second;
      }
      std::cerr << "# FG: Force Constant = " << _params[fgIndex].forceConstant << std::endl
                << "#     Cut Off        = " << _params[fgIndex].cutoff << std::endl;
      std::cerr << "# CG: Force Constant = " << _params[cgIndex].forceConstant << std::endl
                << "#     Cut Off        = " << _params[cgIndex].cutoff << std::endl;
      _withEnergy = true;
    } else if (args.count("energy") != -1) {
      std::cerr << "# There is probably something wrong with the numbers of "
              << "arguments for energy!" << std::endl;
      std::cerr << "# There were " << args.count("energy") << " given,\n";
      std::cerr << "# but 5 are needed!\n";
    }
    
    // Exponent
    _exponent = 6;
    if (args.count("exponent") == 1){
      std::istringstream is(args["exponent"]);
      is >> _exponent;
    } else if (args.count("exponent") != -1) {
      std::cerr << "# There is something wrong with the numbers of "
              << "arguments for exponent!" << std::endl;
      std::cerr << "#   Should be 1!" << std::endl;
    }
    
    // Measure
    if (args.count("measure") == 1){
      std::string fname = args["measure"];
      _foutMeasures.open(fname.c_str());
      _foutMeasures << title();
      _withMeasures = true;
    } else if (args.count("measure") != -1) {
      std::cerr << "# There is something wrong with the numbers of "
              << "arguments for measure!" << std::endl;
      std::cerr << "#   Should be 1!" << std::endl;
    }
    // Weights
    if (args.count("weights") == 1){
      std::string fname = args["weights"];
      _foutWeights.open(fname.c_str());
      _withWeights = true;
    } else if (args.count("weights") != -1) {
      std::cerr << "# There is something wrong with the numbers of "
              << "arguments for weights!" << std::endl;
      std::cerr << "#   Should be 1!" << std::endl;
    }
    // Weights
    if (args.count("mindist") == 1){
      std::string fname = args["mindist"];
      _foutMinDists.open(fname.c_str());
      _withMinimalDistances = true;
    } else if (args.count("mindist") != -1) {
      std::cerr << "# There is something wrong with the numbers of "
              << "arguments for mindist!" << std::endl;
      std::cerr << "#   Should be 1!" << std::endl;
    }
    
  }

  SoluteWeightedDistance::~SoluteWeightedDistance() {
    if (_withMeasures) {
      _foutMeasures.close();
    }
    if (_withWeights) {
      _foutWeights.close();
    }
    if (_withMinimalDistances) {
      _foutMinDists.close();
    }
  }
  
  bool SoluteWeightedDistance::withEnergy() const {
    return _withEnergy;
  }
  
  std::string SoluteWeightedDistance::energyFile() const {
    return _energyFile;
  }
  
  void SoluteWeightedDistance::calculate(double time) {
    
    const double n = _exponent;
    if (_withMeasures) {
      _measures.clear();
    }
    if (_withMinimalDistances) {
      _minDists.clear();
    }
    AtomSpecifier * const as[] = {&_fgSolvent, &_cgSolvent};
    Distances * const ds[] = {&_fgDistances, &_cgDistances};
    for (unsigned int k = 0; k < 2; k++) {
      ds[k]->clear();
      for (unsigned int j = 0; j < as[k]->size(); j++) {
        double sum = 0.0;
        if (_withMeasures || _withWeights) {
          _weights.clear();
        }
        double minDist;
        for (unsigned int i = 0; i < _solute.size(); i++) {
          gmath::Vec nim_sol = _pbc->nearestImage(as[k]->pos(j), _solute.pos(i), _sys.box());
          gmath::Vec dist = as[k]->pos(j) - nim_sol;
          double r_ij = dist.abs();
          if (_withMeasures || _withWeights) {
            _weights.push_back(r_ij);
          }
          if (_withMinimalDistances){
            if (i == 0){
              minDist = r_ij;
            } else if (minDist > r_ij) {
              minDist = r_ij;
            }
          }
          sum += pow(r_ij, -n);
        }
        ds[k]->push_back(pow(sum, -1.0 / n));
        
        if (_withMeasures || _withWeights) {
          const double R_j = ds[k]->back();
          Weights::iterator w_it = _weights.begin(), 
                  w_to = _weights.end();
          for (; w_it != w_to; w_it++){
            const double r_ij = *w_it;
            *w_it = pow((r_ij / R_j), -n);
          }
        }
        if (_withMeasures){
          double measure_j = 0;
          Weights::const_iterator w_it = _weights.begin(), 
                  w_to = _weights.end();
          for (; w_it != w_to; w_it++){
            measure_j += *w_it * log(*w_it);
          }
          double num_solutes = _solute.size();
          measure_j /= log(num_solutes);
          _measures.push_back(1 + measure_j);
        }
        if (_withWeights){
          _writeWeights(ds[k]->back());
        }
        if (_withMinimalDistances){
          _minDists.push_back(minDist);
        }
      }
    }
    if (_withMeasures){
      _writeMeasures(time);
    }
    if (_withMinimalDistances) {
      _writeMinDists(time);
    }
    return;
  }
  
  void SoluteWeightedDistance::_writeMeasures(double time){
    _foutMeasures << time;
    Measures::iterator it = _measures.begin(),
            to = _measures.end();
    for (; it != to; it++){
      _foutMeasures << " " << *it;       
    }
    _foutMeasures << std::endl;
  }
  
  void SoluteWeightedDistance::_writeMinDists(double time){
    _foutMinDists << time;
    MinDists::iterator it = _minDists.begin(),
            to = _minDists.end();
    for (; it != to; it++){
      _foutMinDists << " " << *it;       
    }
    _foutMinDists << std::endl;
  }
  
  void SoluteWeightedDistance::_writeWeights(double R_j){
    _foutWeights << R_j;
    Weights::iterator it = _weights.begin(),
            to = _weights.end();
    for (; it != to; it++){
      _foutWeights << " " << *it;       
    }
    _foutWeights << std::endl;
  }
  
  std::string SoluteWeightedDistance::title() {
    std::string title("# Time   ");
    
    AtomSpecifier * const as[] = {&_fgSolvent, &_cgSolvent};
    for (unsigned int k = 0; k < 2; k++) {
      for (unsigned int i = 0; i < as[k]->size(); i++) {
        title += " ";
        title += as[k]->toString(i);
      }
    }
    title += "\n";
    return title;
  }
  
  void SoluteWeightedDistance::distances(std::ostream &os) const {
    const Distances * const ds[] = {&_fgDistances, &_cgDistances};
    for (unsigned int k = 0; k < 2; k++) {
      Distances::const_iterator it = ds[k]->begin(),
              to = ds[k]->end();
      for (; it != to; it++) {
        os << " ";
        os << *it;
      }
    }
  }

  void SoluteWeightedDistance::energies(std::ostream &os) const {
    const Distances * const ds[] = {&_fgDistances, &_cgDistances};
    for (int k = 0; k < 2; k++) {
      Distances::const_iterator it = ds[k]->begin(),
              to = ds[k]->end();
      for (; it != to; it++) {
        os << " ";
        const double dist = *it - _params[k].cutoff;
        if ((k == fgIndex && dist > 0.0) || (k == cgIndex && dist < 0.0)) {
          os << 0.5 * _params[k].forceConstant * dist * dist;
        } else {
          os << 0.0;
        }
      }
    }
  }
  
  SoluteWeightedDistance::Distances &
  SoluteWeightedDistance::fgDistances(){
    return _fgDistances;
  }
  
  SoluteWeightedDistance::Distances &
  SoluteWeightedDistance::cgDistances(){
    return _cgDistances;
  }
  
  AtomSpecifier &
  SoluteWeightedDistance::solute(){
    return _solute;
  }
  
  AtomSpecifier &
  SoluteWeightedDistance::fgSolvent(){
    return _fgSolvent;
  }
  
  AtomSpecifier &
  SoluteWeightedDistance::cgSolvent(){
    return _cgSolvent;
  }
  
  std::ostream &operator<<(std::ostream &os, SoluteWeightedDistance const & sad){
    sad.distances(os);
    return os;
  }
}
