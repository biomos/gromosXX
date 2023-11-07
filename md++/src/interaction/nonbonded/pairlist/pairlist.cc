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
 * @file pairlist.cc
 * methods of Pairlist
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include "../../../stdheader.h"
#include "../../../interaction/nonbonded/pairlist/pairlist.h"

#include "../../../util/debug.h"

namespace interaction
{
  std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl){
    
    const bool reduced = true;
    
    os << "Pairlist" << std::endl;
    
    Pairlist::const_iterator
      it = pl.begin(),
      to = pl.end();
    
    std::vector<unsigned int>::const_iterator j_it, j_to;
    std::vector<unsigned int> temp;
    
    for (unsigned int i=0; it != to; ++i, ++it) {
      
      int ind = 0;
      // if (!reduced)
      os << std::setw(5) << i << " : " << std::flush;
      
      temp = *it;
      std::sort(temp.begin(), temp.end());
      
      for(j_it = temp.begin(), j_to = temp.end(); j_it != j_to; ++j_it, ++ind){
	
	os << std::setw(5) << *j_it << " "; 

	if (!reduced)
	  if (!(++ind % 15)) os << "\n\t";
      }

      if (!reduced){
	if (ind) os << std::endl;
      }
      else
	os << std::endl;
      
    }    
    return os;
  }
  
} // namespace interaction

