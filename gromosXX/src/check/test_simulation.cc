/**
 * @file test_simulation.cc
 */

#include <algorithm>

#include "test_simulation.h"

#include "../io/argument.h"
#include "../io/read_input.h"

#include "../util/usage.h"

#include "test_parameters.h"

namespace testing {

Test_Simulation::Test_Simulation(const Parameter& parameter) 
  : parameter_(parameter), traj_(parameter.name_) {} // trajectory object requires a name as argument

int Test_Simulation::init_simulation() {
  // supported files
  util::Known knowns;
  for (const auto& input : parameter().input_) { knowns << input.first; }
  
  // create a C-style array of C-style strings
  std::vector<std::string> arguments(parameter().input_.size());
  const unsigned int size = 2 * parameter().input_.size() + 1; // binary_name + n * (feature name + file name)
  char* c_style_args[size]; 
  unsigned int j = 0; // iterator
  c_style_args[j] = &(parameter().binary_name_[0]); // start with binary name
  for (std::vector<std::pair<std::string, std::string>>::iterator it = parameter().input_.begin(), 
    end = parameter().input_.end(); it != end; ++it, ++j) {
    // Gromos expects a '@' before feature names
    arguments[j] = "@" + it->first;
    c_style_args[2 * j + 1] = &(arguments[j][0]); // feature name
    c_style_args[2 * j + 2] = &(it->second[0]); // file name
  }

  // parse arguments
  io::Argument args;
  if (args.parse(size, c_style_args, knowns)) {
    return Error::error_parse_arguments;
  }

  // initialize the simulation
  if (io::read_input(args, topo_, conf_, sim_, md_)) {
    return Error::error_read_input;
  }
  
  // initialize topology, configuration, and simulation objects
  if (md_.init(topo_, conf_, sim_)) {
    return Error::error_md_init;
  }
  
  // create the output files
  traj_.init(args, sim_.param());

  return Error::error_success;
}

int Test_Simulation::run_single_step() {
  int err = Error::error_md_steps; // assume there will not be enough steps provided in the imd file
  if (sim_.steps() < sim_.param().step.number_of_steps) {
    // write trajectory
    traj_.write(conf_, topo_, sim_);
    // run a step
    err = md_.run(topo_, conf_, sim_);
    if (!err) {
      err = Error::error_success; // the step ran successful
      traj_.print(topo_, conf_, sim_);
      sim_.steps() = sim_.steps() + sim_.param().analyze.stride;
      sim_.time() = sim_.param().step.t0 + sim_.steps() * sim_.time_step_size();

      // print final energies
      traj_.write(conf_, topo_, sim_, io::final);
      traj_.print_final(topo_, conf_, sim_);
    }
    else {
      err = Error::error_md_run; // there was a problem with the step
    }
  }
  return err; 
}

int Test_Simulation::run_simulation() {
  int err = Error::error_md_steps; // assume there will not be enough steps provided

  while (sim_.steps() < sim_.param().step.number_of_steps) {
    // write trajectory
    traj_.write(conf_, topo_, sim_);
    // run a step
    err = md_.run(topo_, conf_, sim_);

    // check for errors and potentially terminate
    if (!err) {
      err = Error::error_success; // the step worked
      traj_.print(topo_, conf_, sim_);
      sim_.steps() = sim_.steps() + sim_.param().analyze.stride;
      sim_.time() = sim_.param().step.t0 + sim_.steps() * sim_.time_step_size();
    }
    else {
      err = Error::error_md_run; // the step didn't work
      return err; // return already
    }
  } // main md loop
    
  // if there was no error, write final configuration
  if (err == Error::error_success) {
    this->traj_.write(conf_, topo_, sim_, io::final);
    this->traj_.print_final(topo_, conf_, sim_);
  }

  return err;
}

}