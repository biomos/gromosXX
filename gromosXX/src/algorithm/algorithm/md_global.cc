/**
 * @file md_global.tcc
 * global functions to get an md simulation started.
 */

#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include <debug.h>
#include <timing.h>

#include <math/gmath.h>
#include <io/message.h>
#include <io/argument.h>
#include <io/blockinput.h>
#include <io/GInStream.h>
#include <io/input/InInput.h>

#include <algorithm/md_global.h>

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE algorithm

#include "../../debug.h"

namespace algorithm
{
  // forward declarations
  int grid_scaled_molvir_perturbed_md(io::Argument & args, io::InInput &input);
  int grid_scaled_atvir_perturbed_md(io::Argument & args, io::InInput & input);
  int grid_scaled_molvir_md(io::Argument & args, io::InInput & input);
  int grid_scaled_atvir_md(io::Argument & args, io::InInput & input);
  int grid_scaled_perturbed_md(io::Argument & args, io::InInput & input);
  int grid_scaled_md(io::Argument & args, io::InInput & input);

  int grid_molvir_perturbed_md(io::Argument & args, io::InInput & input);
  int grid_atvir_perturbed_md(io::Argument & args, io::InInput & input);
  int grid_molvir_md(io::Argument & args, io::InInput & input);
  int grid_atvir_md(io::Argument & args, io::InInput & input);
  int grid_perturbed_md(io::Argument & args, io::InInput & input);
  int grid_md(io::Argument & args, io::InInput & input);

  int scaled_molvir_perturbed_md(io::Argument & args, io::InInput & input);
  int scaled_atvir_perturbed_md(io::Argument & args, io::InInput & input);
  int scaled_molvir_md(io::Argument & args, io::InInput & input);
  int scaled_atvir_md(io::Argument & args, io::InInput & input);
  int scaled_perturbed_md(io::Argument & args, io::InInput & input);
  int scaled_md(io::Argument & args, io::InInput & input);

  int molvir_perturbed_md(io::Argument & args, io::InInput & input);
  int atvir_perturbed_md(io::Argument & args, io::InInput & input);
  int molvir_md(io::Argument & args, io::InInput & input);
  int atvir_md(io::Argument & args, io::InInput & input);
  int perturbed_md(io::Argument & args, io::InInput & input);
  int md(io::Argument & args, io::InInput & input);
}

/**
 * perform an MD simulation.
 */
int algorithm::do_md(io::Argument &args)
{
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

  // GRID PAIRLIST?
  bool do_grid;
  int nsnb;
  double rcutp, rcutl, size;
  input.read_PLIST(do_grid, nsnb, rcutp, rcutl, size);

  // INTERACTION SCALING or PERTURBATION
  int ntg, nlam;
  double rlam, dlamt;
  bool do_scaled;
  
  input.read_PERTURB(ntg, rlam, dlamt, nlam, do_scaled);

  // VIRIAL?
  bool calc;
  int ntp;
  double comp, tau;
  math::Matrix pres0;
  interaction::virial_enum do_vir;
  
  input.read_PCOUPLE(calc, ntp, pres0, comp, tau, do_vir);

  // PERTURBATION?
  bool do_perturb = false;
  if (ntg) do_perturb = true;
  if (do_perturb && (args.count("pert") != 1)){
    io::messages.add("PERTURB requested from input file but no "
		     "perturbation topology specified (@pert)!",
		     "md_global::do_md",
		     io::message::error);
  }


  // now the decision tree...
  if (do_grid){
    if (do_scaled){
      if (do_vir == interaction::molecular_virial){
	if (do_perturb){
	  return grid_scaled_molvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return grid_scaled_molvir_md(args, input);
	}	
      }
      else if (do_vir == interaction::atomic_virial){
	if (do_perturb){
	  return grid_scaled_atvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return grid_scaled_atvir_md(args, input);
	}	
      }
      else { // no virial
	if (do_perturb){
	  return grid_scaled_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return grid_scaled_md(args, input);
	}
      }
    }
    else{ // not scaled
      if (do_vir == interaction::molecular_virial){
	if (do_perturb){
	  return grid_molvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return grid_molvir_md(args, input);
	}	
      }
      else if (do_vir == interaction::atomic_virial){
	if (do_perturb){
	  return grid_atvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return grid_atvir_md(args, input);
	}	
      }
      else{ // no virial
	if (do_perturb){
	  return grid_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return grid_md(args, input);
	}
      }      
    }
  }
  else{ // no grid
    if (do_scaled){
      if (do_vir == interaction::molecular_virial){
	if (do_perturb){
	  return scaled_molvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return scaled_molvir_md(args, input);
	}	
      }
      else if (do_vir == interaction::atomic_virial){
	if (do_perturb){
	  return scaled_atvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return scaled_atvir_md(args, input);
	}	
      }
      else{ // no virial
	if (do_perturb){
	  return scaled_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return scaled_md(args, input);
	}
      }
    }
    else{ // not scaled
      if (do_vir == interaction::molecular_virial){
	if (do_perturb){
	  return molvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return molvir_md(args, input);
	}	
      }
      if (do_vir == interaction::atomic_virial){
	if (do_perturb){
	  return atvir_perturbed_md(args, input);
	}
	else{ // not perturbed
	  return atvir_md(args, input);
	}	
      }
      else{ // no virial
	if (do_perturb){
	  return perturbed_md(args, input);
	}
	else{ // not perturbed
	  return md(args, input);
	}
      }      
    }    
  }
}


