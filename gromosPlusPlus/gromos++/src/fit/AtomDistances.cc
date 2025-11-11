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

// fit_AtomDistances.cc
//includes explicit calls to gsl now

#include "AtomDistances.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

#include "PositionUtils.h"
#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gmath/Matrix.h"
#include "../gmath/Vec.h"
#include "../utils/AtomSpecifier.h"


using gmath::Vec;
using namespace fit;
using namespace std;



double AtomDistances::dist(vector<Vec> const &ref, 
			       vector<Vec> const &sys)const
{
  double dist2=0.0;
  for(size_t i=0; i< ref.size(); ++i){
    for(size_t j=i+1; j < ref.size(); ++j){
      for(int b=0;b<3;++b){
	dist2 += 
	  ((ref[i][b] - ref[j][b]) - (sys[i][b] - sys[j][b])) *
	  ((ref[i][b] - ref[j][b]) - (sys[i][b] - sys[j][b]));
      }
    }
  }
  int N = ref.size() * (ref.size() - 1) / 2;
  
  return sqrt(dist2/N);
}
