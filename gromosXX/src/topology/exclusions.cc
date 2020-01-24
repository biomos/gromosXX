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
