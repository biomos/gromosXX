/**
 * @file parse_verbosity.cc
 * implementation of the parse verbosity function
 */

#include "../stdheader.h"

#include "../util/error.h"
#include "../io/argument.h"
#include "../io/message.h"

#include "parse_verbosity.h"

int util::parse_verbosity(io::Argument &args, std::string flag, 
			   std::ostream &os)
{
#ifndef NDEBUG
  // parse the verbosity
  if (args.count(flag) == -1) return 0;

  io::Argument::const_iterator it = args.lower_bound(flag),
    to = args.upper_bound(flag);

  os << "setting debug level\n";
      
  for( ; it != to; ++it){

    int level;
    std::string module = "";
    std::string submodule = "";

    std::string s(it->second);
    std::string::size_type sep = s.find(':');
    if (sep == std::string::npos){
      // no module or submodule
      std::istringstream css(s);
      css >> level;
      // level = atoi(s.c_str());
      assert(!css.fail());
    }
    else{
      module = s.substr(0, sep);
      std::string second = s.substr(sep+1, std::string::npos);
	  
      sep = second.find(':');
      if (sep == std::string::npos){
	// no submodule
	std::istringstream css(second);
	css >> level;
	// level = atoi(second.c_str());
	assert(!css.fail());
      }
      else{
	std::istringstream css(second.substr(sep+1, std::string::npos));
	css >> level;
	// level = atoi(second.substr(sep+1, std::string::npos).c_str());
	assert(!css.fail());
	submodule = second.substr(0, sep);
      }
    }
	
    if (module == "")
      os << "\t" << std::setw(15) <<  "global";
    else
      os << "\t" << std::setw(15) <<  module;
    os << "\t" << std::setw(15) << submodule 
       << "\t" << std::setw(6) << level << "\n";

    if (module == ""){
      if (submodule == "") ::debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "configuration"){
      if (submodule == "") configuration::debug_level = level;
      else if (submodule == "configuration") configuration::configuration_debug_level = level;
      else if (submodule == "energy") configuration::energy_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "topology"){
      if (submodule == "") topology::debug_level = level;
      else if (submodule == "topology") topology::topology_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "algorithm"){
      if (submodule == "") algorithm::debug_level = level;
      else if (submodule == "algorithm") 
	algorithm::algorithm_debug_level = level;
      else if (submodule == "constraints") 
	algorithm::constraints_debug_level = level;
      else if (submodule == "integration") 
	algorithm::integration_debug_level = level;
      else if (submodule == "temperature") 
	algorithm::temperature_debug_level = level;
      else if (submodule == "pressure") 
	algorithm::pressure_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "interaction"){
      if (submodule == "") interaction::debug_level = level;
      else if (submodule == "interaction") 
	interaction::interaction_debug_level = level;
      else if (submodule == "bonded") 
	interaction::bonded_debug_level = level;
      else if (submodule == "nonbonded") 
	interaction::nonbonded_debug_level = level;
      else if (submodule == "latticesum")
	interaction::latticesum_debug_level = level;
      else if (submodule == "pairlist") 
	interaction::pairlist_debug_level = level;
      else if (submodule == "forcefield")
	interaction::forcefield_debug_level = level;
      else if (submodule == "filter") 
	interaction::filter_debug_level = level;
      else if (submodule == "special") 
	interaction::special_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "io"){
      if (submodule == "") io::debug_level = level;
      else if (submodule == "parameter") io::parameter_debug_level = level;
      else if (submodule == "topology") io::topology_debug_level = level;
      else if (submodule == "configuration") io::configuration_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "simulation"){
      if (submodule == "") simulation::debug_level = level;
      else if (submodule == "simulation")
	simulation::simulation_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "math"){
      if (submodule == "") math::debug_level = level;
      else if (submodule == "math") math::math_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else if (module == "util"){
      if (submodule == "") util::debug_level = level;
      else if (submodule == "util") util::util_debug_level = level;
      else if (submodule == "replica") util::replica_debug_level = level;
      else if (submodule == "leus") util::leus_debug_level = level;
      else if (submodule == "bs_leus") util::bs_leus_debug_level = level;
      else return E_NOT_IMPLEMENTED;
    }
    else{
      return E_NOT_IMPLEMENTED;
    }

  }
      
#endif
  return 0;
  
}
