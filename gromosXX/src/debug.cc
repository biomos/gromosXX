/**
 * This should be included in the main cc file.
 * @file debug.cc
 */

#ifndef NDEBUG
int debug_level = 0;
#endif

void parse_verbosity(io::Argument &args, std::string flag = "verb", std::ostream &os = std::cout)
{
  // parse the verbosity
  if (args.count(flag) == -1) return;

#ifdef NDEBUG
  throw std::string("@verb not supported with non-debug compilation");
#endif
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
      level = atoi(s.c_str());
    }
    else{
      module = s.substr(0, sep);
      std::string second = s.substr(sep+1, std::string::npos);
	  
      sep = second.find(':');
      if (sep == std::string::npos){
	// no submodule
	level = atoi(second.c_str());
      }
      else{
	level = atoi(second.substr(sep+1, std::string::npos).c_str());
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
      else throw std::string("submodule without module illegal");
    }
    else if (module == "simulation"){
      if (submodule == "") simulation::debug_level = level;
      else if (submodule == "simulation") simulation::simulation_debug_level = level;
      else if (submodule == "system") simulation::system_debug_level = level;
      else if (submodule == "topology") simulation::topology_debug_level = level;
      else throw std::string("unknown submodule");
    }
    else if (module == "algorithm"){
      if (submodule == "") algorithm::debug_level = level;
      else if (submodule == "constraint") algorithm::constraint_debug_level = level;
      else if (submodule == "integration") algorithm::integration_debug_level = level;
      else throw std::string("unknown submodule");	  
    }
    else if (module == "interaction"){
      if (submodule == "") interaction::debug_level = level;
      else if (submodule == "forcefield") interaction::forcefield_debug_level = level;
      else if (submodule == "interaction") interaction::interaction_debug_level = level;
      else if (submodule == "pairlist") interaction::pairlist_debug_level = level;
      else throw std::string("unknown submodule");	  
    }
    else if (module == "io"){
      if (submodule == "") io::debug_level = level;
      else if (submodule == "input") io::input_debug_level = level;
      else if (submodule == "topology") io::topology_debug_level = level;
      else if (submodule == "trajectory") io::trajectory_debug_level = level;
      else throw std::string("unknown submodule");	  
    }
    else if (module == "math"){
      if (submodule == "") math::debug_level = level;
      else throw std::string("unknown submodule");
    }
    else{
      throw std::string("unknown module");
    }

  }
      
}

