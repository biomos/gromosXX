/**
 * @file pairlist.cc
 * methods of Pairlist
 */

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

#include <stdheader.h>
#include <interaction/nonbonded/pairlist/pairlist.h>

#include <util/debug.h>

namespace interaction
{
  inline std::ostream & 
  operator<<(std::ostream &os, Pairlist &pl){
    
    os << "Pairlist" << std::endl;
    
    Pairlist::const_iterator
      it = pl.begin(),
      to = pl.end();
    
    std::vector<unsigned int>::const_iterator j_it, j_to;
    
    for (unsigned int i=0; it != to; ++i) {
      
      int ind = 0;
      os << std::setw(5) << i << ": " << std::flush;
      
      for(j_it = it->begin(), j_to = it->end(); j_it != j_to; ++j_it, ++ind){
	
	os << std::setw(5) << *j_it << " "; 
	if (!(++ind % 15)) os << "\n\t";
      }
      if (ind) os << std::endl;
    }    
    return os;
  }
  
} // namespace interaction

