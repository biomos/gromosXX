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

#include "exclusions.h"
namespace topology {

	bool Exclusions::find_and_remove(const value_type & j) {
		cont_t::iterator low_bound = std::lower_bound(indeces.begin(),
			indeces.end(), j);			
		if(low_bound == indeces.end()){
			return false;
		} else if(*low_bound == j) {
			indeces.erase(low_bound);
			return true;
		} else {
			return false;
		}
	}

	bool Exclusions::insert(const value_type & j) {
		cont_t::iterator low_bound = std::lower_bound(indeces.begin(),
			indeces.end(), j);			
		if(low_bound == indeces.end()) {
			// j greater than contained values
			indeces.push_back(j);
			return true;
		} else if(*low_bound == j) {
			// j already in list
			return false;
		} else {
			// insert j before greater value
			indeces.insert(low_bound, j);
			return true;
		}
	}
}
