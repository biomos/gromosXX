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

/**
 * perform an MD simulation.
 */
int do_md(io::Argument &args)
{
  // decide which code options to use
  // we need the input file!
  io::InInput input;
  DEBUG(7, "opening input");
  std::ifstream *input_file = new std::ifstream(args["input"].c_str());
  if (!input_file->good())
    io::messages.add("unable to open input file: " + args["input"], 
                     "md.tcc",
		     io::message::error);
  else
    io::messages.add("parsing input file: " + args["input"], "md.tcc",
		     io::message::notice);
  input.stream(*input_file);
  DEBUG(7, "reading input");
  input.readStream();
  input.auto_delete(true);

  // VIRIAL?
  int ntb, nrdbox;
  input.read_BOUNDARY(ntb, nrdbox);
  interaction::virial_enum do_vir = interaction::no_virial;
  if (abs(ntb) == 2)
    do_vir = interaction::molecular_virial;
  
  // PERTURBATION?
  int ntg, nlam;
  double rlam, dlamt;
  
  input.read_PERTURB(ntg, rlam, dlamt, nlam);
  bool do_perturb = false;
  if (ntg) do_perturb = true;
  if (do_perturb && (args.count("pert") != 1)){
    io::messages.add("PERTURB requested from input file but no "
		     "perturbation topology specified (@pert)!",
		     "md_global::do_md",
		     io::message::error);
  }

  //==================================================
  // create the algorithm
  //==================================================

  if (do_perturb){
    switch(do_vir){
      case interaction::no_virial:
	{
	  algorithm::Perturbation_MD<
	    program::perturbed_FLEXI_spec<interaction::no_virial>,
	    algorithm::Interaction_spec<
	    program::perturbed_FLEXI_spec<interaction::no_virial>::simulation_type,
	    // perturbation
	    true,
	    // virial
	    interaction::no_virial,
	    // atomic cutoff
	    false,
	    // scaling
	    false
	    >
	    > 
	    the_MD;

	  int res = the_MD.do_md(args, input);
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
	  return res;
	}
	
      case interaction::atomic_virial:
	{
	  algorithm::Perturbation_MD<
	    program::perturbed_FLEXI_spec<interaction::atomic_virial>,
	    algorithm::Interaction_spec<
	    program::perturbed_FLEXI_spec<interaction::atomic_virial>::simulation_type,
	    // perturbation
	    true,
	    // virial
	    interaction::atomic_virial,
	    // atomic cutoff
	    false,
	    // scaling
	    false
	    >
	    > 
	    the_MD;
	
	  int res = the_MD.do_md(args, input);
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
	  return res;
	  
	}
	
      case interaction::molecular_virial:
	{
	  algorithm::Perturbation_MD<
	    program::perturbed_FLEXI_spec<interaction::molecular_virial>,
	    algorithm::Interaction_spec<
	    program::perturbed_FLEXI_spec<interaction::molecular_virial>::simulation_type,
	    // perturbation
	    true,
	    // virial
	    interaction::molecular_virial,
	    // atomic cutoff
	    false,
	    // scaling
	    false
	    >
	    > 
	    the_MD;
	  
	  int res = the_MD.do_md(args, input);

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
	  return res;
	  
	}
	
      default:
	io::messages.add("wrong virial method specified",
			 "md_global::do_md",
			 io::message::error);
    }

  }
  else{
    switch(do_vir){
      case interaction::no_virial:
	{
	  algorithm::MD<
	    program::FLEXI_spec<interaction::no_virial>,
	    algorithm::Interaction_spec<
	    program::FLEXI_spec<interaction::no_virial>::simulation_type,
	    // perturbation
	    false,
	    // virial
	    interaction::no_virial,
	    // atomic cutoff
	    false,
	    // scaling
	    false
	    >
	    > 
	    the_MD;
	  
	  int res = the_MD.do_md(args, input);
	  
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
	  return res;
	}
	
      case interaction::atomic_virial:
	{
	  algorithm::MD<
	    program::FLEXI_spec<interaction::atomic_virial>,
	    algorithm::Interaction_spec<
	    program::FLEXI_spec<interaction::atomic_virial>::simulation_type,
	    // perturbation
	    false,
	    // virial
	    interaction::atomic_virial,
	    // atomic cutoff
	    false,
	    // scaling
	    false
	    >
	    > 
	    the_MD;
	  
	  int res = the_MD.do_md(args, input);

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
	  return res;
	}
	
      case interaction::molecular_virial:
	{
	  algorithm::MD<
	    program::FLEXI_spec<interaction::molecular_virial>,
	    algorithm::Interaction_spec<
	    program::FLEXI_spec<interaction::molecular_virial>::simulation_type,
	    // perturbation
	    false,
	    // virial
	    interaction::molecular_virial,
	    // atomic cutoff
	    false,
	    // scaling
	    false
	    >
	    > 
	    the_MD;
	  
	  int res = the_MD.do_md(args, input);
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
	  return res;
	  
	}
	
      default:
	io::messages.add("wrong virial method specified",
			 "md_global::do_md",
			 io::message::error);
    }

  }
  
  return 10;

}


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

    if (do_md(args)){
      std::cout << "\nMD failed\n\n" << std::endl;
    }
    else{
      std::cout << "\nMD finished successfully\n\n" << std::endl;
    }
    
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

