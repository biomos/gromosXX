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

#ifndef GROMOSXX_EXCLUSIONS_H
#define GROMOSXX_EXCLUSIONS_H

#include <vector>
#include <algorithm>

namespace topology {
	/*
	 * @class Exclusions
	 * sorted storage for excluded atom indeces
	 */
	class Exclusions {
	public:
		/*
		 * container type in which indeces are stored in
		 */
		typedef std::vector<int> cont_t;

		/*
		 * standard types
		 */
		typedef cont_t::value_type value_type;
		typedef cont_t::reference reference;
		typedef cont_t::const_reference const_reference;
		typedef cont_t::pointer pointer;
		typedef cont_t::const_pointer const_pointer;
		typedef cont_t::size_type size_type;
		typedef cont_t::difference_type difference_type;

		/*
		 * only allow const iterators to ensure vector remains sorted
		 */
		typedef cont_t::const_iterator const_iterator;
		typedef const_iterator iterator;
		typedef cont_t::const_reverse_iterator const_reverse_iterator;
		typedef const_reverse_iterator reverse_iterator;

		/**
		 * find and delete atom j
		 * returns true if deletion was performed, false otherwise
		 */
		bool find_and_remove(const value_type & j);

		/**
		 * delete atom by iterator
		 * returns iterator to the item following the erased item
		 */
		inline cont_t::iterator erase(cont_t::const_iterator& it) {
			// erase should also work with const_iterator, but somehow doesnt
			// This is a workaround
			cont_t::iterator nit = indeces.begin() + (it - indeces.begin());
			cont_t::iterator new_it = indeces.erase(nit);
			return new_it;
		}

		/*
		 * checks if j is contained in exlusion list
		 */
		inline bool is_excluded(const value_type & j) const {
			// first check for common cases
			if(indeces.empty()){
				return false;
			} else if(indeces.back() < j) {
				return false;
			}
			
			// linear search ist faster than binary search for small sizes
			return indeces.end() != std::find(indeces.begin(), indeces.end(), j);
		}

		/*
		 * inserts j into exlusion list
		 * returns true if j was inserted,
		 * returns false if j was already in list
		 */
		bool insert(const value_type & j);

		/*
		 * provide standard insert functionality
		 */
		inline iterator insert(iterator, const value_type & j) {
			insert(j);
			return std::find(indeces.begin(), indeces.end(),j);
		}

		/*
		 * clears all exclusions
		 */
		inline void clear() {
			indeces.clear();
		}

		/*
		 * random access to constant references
		 */
		inline const_reference operator[](size_type pos) const {
			return indeces[pos];
		}

		/*
		 * number of exclusions
		 */
		inline size_type size() const {
			return indeces.size();
		}

		/*
		 * const iterator functions
		 */
		inline cont_t::const_iterator begin() const {
			return indeces.begin();
		}

		inline cont_t::const_iterator end() const {
			return indeces.end();
		}

		inline cont_t::const_reverse_iterator rbegin() const {
			return indeces.rbegin();
		}
		
		inline cont_t::const_reverse_iterator rend() const {
			return indeces.rend();
		}

	private:
		cont_t indeces;
	};


	/**
	 * container type used to store exlusions for all atoms
	 */
    typedef std::vector<Exclusions> excl_cont_t;
}


#endif //GROMOSXX_EXCLUSIONS_H
