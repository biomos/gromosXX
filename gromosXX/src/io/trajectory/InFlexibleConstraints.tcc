/**
 * @file InFlexibleConstraints.tcc
 * implements input for flexible constraints.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE trajectory

#include "../../debug.h"

/**
 * Constructor.
 */
io::InFlexibleConstraints::InFlexibleConstraints(std::istream &is) 
  : GInStream(is) 
{
  // read the whole file at beginning
  readStream();
};

/**
 * Read in a G96 trajectory into the system.
 */
inline void
io::InFlexibleConstraints::read_FLEXCON(std::vector<double> &vel,
					std::vector<simulation::compound::
					distance_constraint_struct> &constr)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  { // FLEXCON
    DEBUG(8, "reading in FLEXCON block");
    buffer = m_block["FLEXCON"];
    if (buffer.size()){
      
      it = buffer.begin() + 1;
      int num = 0;
      
      // flexcon.distance().clear();
      vel.clear();

      std::vector<simulation::compound::distance_constraint_struct>::iterator
	cit = constr.begin(),
	cto = constr.end();

      for( ; it != buffer.end() - 1; ++it, ++cit, ++num){
	assert(cit != cto);
	double dist, v;
	_lineStream.clear();
	_lineStream.str(*it);
	
	_lineStream >> dist >> v;

	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in FLEXCON block");
	
	cit->b0 = dist;
	vel.push_back(v);

      }
      assert(cit == cto);
    }
    else{
      DEBUG(8, "no FLEXCON block!");
      io::messages.add("no FLEXCON block in file",
		       "InFlexibleConstraints", io::message::error);
    }
    
  } // FLEXCON
  
}
