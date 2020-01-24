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

