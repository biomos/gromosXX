/**
 * @file flexi.cc
 * the md program with flexible constraints.
 */

#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include "../src/debug.h"

#include <math/gmath.h>
#include <io/message.h>
#include <simulation/core.h>
#include <math/periodicity.h>
#include <simulation/simulation.h>
#include <simulation/perturbation.h>
#include <interaction/interaction.h>
#include <io/io.h>
#include <algorithm/algorithm.h>

// sppecial includes
#include <algorithm/integration/runge_kutta.h>

// global variables for debug
#include "../src/debug.cc"

#include "flexi.h"

using namespace math;

int main(int argc, char *argv[])
{
  try{
    
    char *knowns[] = 
      {
	"topo", "struct", "input", "verb", "pert",
	"trj", "fin", "trv", "trf", "tre", "trg", "print", "trp",
	"flexcon", "trc"
      };
    
    int nknowns = 15;
    
    string usage = argv[0];
    usage += "\n\t@topo     <topology>\n";
    usage += "\t[@pert    <perturbation topology>]\n";
    usage += "\t@struct   <coordinates>\n";
    usage += "\t@input    <input>\n";
    usage += "\t@trj      <trajectory>\n";
    usage += "\t@fin      <final structure>\n";
    usage += "\t[@trv     <velocity trajectory>]\n";
    usage += "\t[@trf     <force trajectory>]\n";
    usage += "\t[@tre     <energy trajectory>]\n";
    usage += "\t[@trg     <free energy trajectory>]\n";
    usage += "\t[@print   <pairlist/force>]\n";
    usage += "\t[@trp     <print file>]\n";
    usage += "\t[@flexcon <flexible constraints data file>]\n";
    usage += "\t[@trc     <flexible constraints output file>]\n";
#ifndef NDEBUG
    usage += "\t[@verb    <[module:][submodule:]level>]\n";
#endif

    io::Argument args(argc, argv, nknowns, knowns, usage);

    // parse the verbosity flag and set debug levels
    parse_verbosity(args);

    // determine whether we do perturbation
    bool perturbation = false;
    if (args.count("pert") == 1)
      perturbation = true;
    
    if (perturbation){ // leap frog + perturbation

      algorithm::Perturbation_MD<
	program::perturbed_FLEXI_spec
	>
	the_MD;
      
      if (the_MD.do_md(args)){
	return 1;
      }

      if (args.count("trc") > 0){
	
	std::ofstream fc(args["trc"].c_str());
	io::OutFlexibleConstraints ofc(fc);

	ofc.write_title(the_MD.title + 
			"\nflexible constraints positions and velocities");

	ofc.write_FLEXCON(the_MD.distance_constraint_algorithm().vel(),
			  the_MD.simulation().topology());
      }
      else{
	io::messages.add("Final flexible constraints velocities not written out",
			 "flexi", io::message::warning);
      }

    }
    else{ // leap frog, no perturbation

      algorithm::MD<
	program::FLEXI_spec
	>
	the_MD;
      
      if (the_MD.do_md(args)){
	return 1;
      }

      if (args.count("trc") > 0){
	
	std::ofstream fc(args["trc"].c_str());
	io::OutFlexibleConstraints ofc(fc);

	ofc.write_title(the_MD.title + 
			"\nflexible constraints positions and velocities");

	ofc.write_FLEXCON(the_MD.distance_constraint_algorithm().vel(),
			  the_MD.simulation().topology());
      }
      else{
	io::messages.add("Final flexible constraints velocities not written out",
			 "flexi", io::message::warning);
      }
      
    }
    
    std::cout << "\nMD finished successfully\n\n" << std::endl;
  
    std::cout << "messages (simulation)\n";
    io::messages.display(std::cout);
    std::cout << "\n\n";
  
  }
  catch(std::runtime_error e){
    // runtime problem
    std::cout << "severe error encountered:\n"
	      << e.what() << std::endl;
    std::cout << "messages (simulation)\n";
    
    io::messages.display();
    std::cout << std::endl;
    
    return 1;
  }
  catch(std::string s){
    // argument exception
    std::cout << s << std::endl;
    return 2;
  }
  
  return 0;
}

