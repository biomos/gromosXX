#include <algorithm>

#include "test_simulation.h"

#include "../io/argument.h"
#include "../io/read_input.h"

#include "../util/usage.h"

int test::Test_Simulation::init_simulation() {
  // supported files
  std::vector<std::string> qualifiers = {"topo", "conf", "input", "qmmm", "trc", "tre", "fin"};
  util::Known knowns;
  for (std::string& e : qualifiers) { knowns << e; }
  
  // process arguments
  std::for_each(qualifiers.begin(), qualifiers.end(), [](std::string& s) {s.insert(0, "@"); });
  io::Argument args;

  // create a C-style array of C-style strings
  char* c_style_args[] = {
    &this->binary_name_[0],
    &qualifiers[0][0],
    &this->topology_file_[0],
    &qualifiers[1][0],
    &this->configuration_file_[0],
    &qualifiers[2][0],
    &this->simulation_file_[0],
    &qualifiers[3][0],
    &this->qmmm_file_[0],
    &qualifiers[4][0],
    &this->coordinate_trajectory_file_[0],
    &qualifiers[5][0],
    &this->energy_trajectory_file_[0],
    &qualifiers[6][0],
    &this->final_configuration_file_[0]
  };

  if (args.parse(15, c_style_args, knowns)) {
    return EXIT_FAILURE;
  }

  // initialize the simulation
  if (io::read_input(args, this->topo_, this->conf_, this->sim_, this->md_)) {
    return EXIT_FAILURE;
  }
  
  // initialize topology, configuration, and simulation objects
  this->md_.init(this->topo_, this->conf_, this->sim_);
  // create the output files
  this->traj_.init(args, this->sim_.param());

  return EXIT_SUCCESS;
}

int test::Test_Simulation::run_single_step() {
  int err;
  if (this->sim_.steps() < this->sim_.param().step.number_of_steps) {
    // write trajectory
    this->traj_.write(this->conf_, this->topo_, this->sim_);
    // run a step
    err = this->md_.run(this->topo_, this->conf_, this->sim_);
    if (!err) {
      this->traj_.print(this->topo_, this->conf_, this->sim_);
      this->sim_.steps() = this->sim_.steps() + this->sim_.param().analyze.stride;
      this->sim_.time() = this->sim_.param().step.t0 + this->sim_.steps() * this->sim_.time_step_size();

      // print final energies
      this->traj_.write(this->conf_, this->topo_, this->sim_, io::final);
      this->traj_.print_final(this->topo_, this->conf_, this->sim_);
    }
  }
  else {
    err = EXIT_FAILURE;
  }
  return err; 
}

int test::Test_Simulation::run_simulation() {
  int err;

  while (this->sim_.steps() < this->sim_.param().step.number_of_steps) {
    // write trajectory
    this->traj_.write(this->conf_, this->topo_, this->sim_);
    // run a step
    err = this->md_.run(this->topo_, this->conf_, this->sim_);

    // check for errors and potentially terminate
    if (!err) {
      this->traj_.print(this->topo_, this->conf_, this->sim_);
      this->sim_.steps() = this->sim_.steps() + this->sim_.param().analyze.stride;
      this->sim_.time() = this->sim_.param().step.t0 + this->sim_.steps() * this->sim_.time_step_size();
    }
    else {
      return err;
    }
  } // main md loop
    
  // final configuration
  this->traj_.write(this->conf_, this->topo_, this->sim_, io::final);
  this->traj_.print_final(this->topo_, this->conf_, this->sim_);

  return err;
}