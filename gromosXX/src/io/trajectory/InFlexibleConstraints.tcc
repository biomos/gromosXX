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
 * Read in a flexible constraints file into a topology.
 */
inline void
io::InFlexibleConstraints::read_FLEXCON(std::vector<double> &vel,
					simulation::Topology &topo)
{
  std::vector<simulation::compound::distance_constraint_struct>
    & constr = topo.solute().distance_constraints();
  
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
	int i, j;
	double dist, v;
	_lineStream.clear();
	_lineStream.str(*it);
	
	_lineStream >> i >> j >> dist >> v;

	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in FLEXCON block");
	
	if (cit->i != i || cit->j != j)
	  io::messages.add("FLEXCON block: assigning velocities to "
			   "different constraints!",
			   "InFlexibleConstraints", io::message::warning);

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

/**
 * Read in a flexible constraints file into a perturbation topology.
 */
inline void
io::InFlexibleConstraints::
read_FLEXCON(std::vector<double> &vel,
	     simulation::Perturbation_Topology &topo)
{
  std::vector<simulation::compound::distance_constraint_struct>
    & constr = topo.solute().distance_constraints();
  
  std::vector<simulation::Perturbed_Solute::
    perturbed_distance_constraint_struct> 
    & pert_constr = topo.perturbed_solute().distance_constraints();

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

      for( ; cit != cto; ++it, ++cit, ++num){
	assert(cit != cto);
	int i, j;
	double dist, v;
	_lineStream.clear();
	_lineStream.str(*it);
	
	_lineStream >> i >> j >> dist >> v;

	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in FLEXCON block");
	
	if (cit->i != i || cit->j != j)
	  io::messages.add("FLEXCON block: assigning velocities to "
			   "different constraints!",
			   "InFlexibleConstraints", io::message::warning);

	cit->b0 = dist;
	vel.push_back(v);

      }

      // and the perturbed constraints
      std::vector<simulation::Perturbed_Solute::
	perturbed_distance_constraint_struct>::iterator
	pcit = pert_constr.begin(),
	pcto = pert_constr.end();

      for( ; pcit != pcto; ++it, ++pcit, ++num){

	int i, j;
	double dist, v;
	_lineStream.clear();
	_lineStream.str(*it);
	
	_lineStream >> i >> j >> dist >> v;

	if (_lineStream.fail() || ! _lineStream.eof())
	  throw std::runtime_error("bad line in FLEXCON block");
	
	if (pcit->i != i || pcit->j != j)
	  io::messages.add("FLEXCON block: assigning velocities to "
			   "different constraints! Reorder for perturbation?",
			   "InFlexibleConstraints", io::message::warning);

	pcit->b0 = dist;
	vel.push_back(v);

      }

    }
    else{
      DEBUG(8, "no FLEXCON block!");
      io::messages.add("no FLEXCON block in file",
		       "InFlexibleConstraints", io::message::error);
    }
    
  } // FLEXCON
  
}
