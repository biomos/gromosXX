/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file nn_worker.cc
 * worker for the neural network interface
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction.h"

#include "../../../io/blockinput.h"

#include "../../../util/timing.h"
#include "../../../util/system_call.h"
#include "../../../util/debug.h"

#ifdef HAVE_PYBIND11
  #include <pybind11/stl.h>
  #include <pybind11/pybind11.h>
#endif

// special interactions
#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "nn_worker.h"

#ifdef OMP
  #include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

#ifdef HAVE_PYBIND11
  using namespace py::literals;
#endif

interaction::NN_Worker::NN_Worker() : QM_Worker("NN Worker"),
                                      param(nullptr),
                                      guard(),
                                      mlp_calculator() {};

int interaction::NN_Worker::init(const topology::Topology& topo
                               , const configuration::Configuration& conf
                               , simulation::Simulation& sim
                               , const interaction::QM_Zone& qm_zone) {
  DEBUG(15, "Initializing " << this->name());

#ifdef HAVE_PYBIND11
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.nn);
  QM_Worker::param = this->param;

  // writing of QM trajectory currently not supported as NN_Worker
  // has its own implementatino of run_QM()
  // -> calling QM_Worker::init is thus unnecessary
  // parent class initialization (trajectory files)
  // int err = QM_Worker::init(topo, conf, sim, qm_zone);
  // if (err) return err;

  // Initialize NN interface
  // Initialize pybind, schnetpack, python script

  #ifdef OMP
    // Python module modifies omp_num_threads value, we gonna restore it
    unsigned num_threads;
    #pragma omp parallel
      #pragma omp single
        num_threads = omp_get_num_threads();
  #endif

  // determine QM software
  const simulation::qm_software_enum software = sim.param().qmmm.software;

  // Model path
  py::str model_path = sim.param().qmmm.nn.model_path;
  DEBUG(11, "model_path: " << model_path.cast<std::string>());

  // Validation models paths
  py::list val_models_paths;
  if (!sim.param().qmmm.nn.val_model_paths.empty()) {
    std::vector<std::string> val_paths;
    for (std::vector<std::string>::iterator it = sim.param().qmmm.nn.val_model_paths.begin(); it != sim.param().qmmm.nn.val_model_paths.end(); ++it) {
        val_paths.push_back(*it);
        DEBUG(11, "val_model_paths " << *it);
    }
    val_models_paths = py::cast(val_paths);
  }

  // How often to run NN validation
  py::int_ write_val_step = sim.param().qmmm.nn.val_steps;

  // How often to write energy
  py::int_ write_energy_step = sim.param().write.energy;

  // To be able to import the module from the current directory
  // Get the directory of the executable
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  std::string exePath = std::string(result, (count > 0) ? count : 0);
  std::string exeDir = exePath.substr(0, exePath.find_last_of("/"));
  // Compute the correct module directory assuming binary in BUILD_X/bin
  std::string modulePath = exeDir + "/../contrib";
  py::module_ sys = py::module_::import("sys");
  // Convert C++ string to Python string and append to sys.path
  sys.attr("path").attr("append")(modulePath);

  // To be able to import the module from the current directory
  // py::module_ sys = py::module_::import("sys");
  // sys.attr("path").attr("append")("/local/gromosXX/md++/src/interaction/qmmm");

  // Initialize mlp_calculator Python object
  if (software == simulation::qm_schnetv1) {
    // Initialize schnet_v1 module
    py::module_ schnet_v1 = py::module_::import("schnet_v1");
    
    // Decide if perturbation is performed or not
    if (sim.param().perturbation.perturbation) {

      // get lambda parameter
      py::float_ lambda = py::cast(sim.param().perturbation.lambda);

      // get perturbed QM states
      py::list perturbed_qm_states = py::cast(sim.param().qmmm.nn.pertqm_state);

      // Initialize mlp_calculator Pert_SchNet_V1_Calculator Python object
      mlp_calculator = schnet_v1.attr("Pert_SchNet_V1_Calculator")(model_path, val_models_paths, write_val_step, write_energy_step, lambda, perturbed_qm_states);
    }

    else {
      // Initialize mlp_calculator SchNet_V1_Calculator Python object
      mlp_calculator = schnet_v1.attr("SchNet_V1_Calculator")(model_path, val_models_paths, write_val_step, write_energy_step);
    }
  }

  if (software == simulation::qm_schnetv2) {
    // Initialize schnet_v2 module
    py::module_ schnet_v2 = py::module_::import("schnet_v2");

    // Decide if perturbation is performed or not
    if (sim.param().perturbation.perturbation) {
      io::messages.add("Perturbation not implemented with schnet_v2 only schnet_v1", io::message::error);
    }

    else {
      // Initialize mlp_calculator SchNet_V2_Calculator Python object
      mlp_calculator = schnet_v2.attr("SchNet_V2_Calculator")(model_path, val_models_paths, write_val_step, write_energy_step);
    }
  }
  
  // Restore omp_num_threads
  #ifdef OMP
    omp_set_num_threads(num_threads);
  #endif

#endif

  DEBUG(15, "Initialized " << this->name());
  return 0;
}

interaction::NN_Worker::~NN_Worker() = default;

int interaction::NN_Worker::run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone) {
#ifdef HAVE_PYBIND11
  // run NN interface 

  // Prepare the input for mlp_calculator object
  double length_to_nn = 1 / this->param->unit_factor_length;

  std::set<interaction::QM_Atom>::const_iterator it, to;
  it = qm_zone.qm.begin();
  to = qm_zone.qm.end();

  // Atomic numbers
  std::vector<uint32_t> atom_nums;

  // Coordinates
  py::list system_coordinates;

  for (;it != to; ++it) {
    atom_nums.push_back(it->atomic_number);

    // convert positions
    math::Vec nn_pos = it->pos * length_to_nn;
    py::list atom_coordinates;
    atom_coordinates.attr("append")(nn_pos[0]);
    atom_coordinates.attr("append")(nn_pos[1]);
    atom_coordinates.attr("append")(nn_pos[2]);
    system_coordinates.attr("append")(atom_coordinates);
    DEBUG(15, "atom to NN: " << it->index << " : " << math::v2s(nn_pos));
  }

  // Write capping atoms
  DEBUG(15,"Writing capping atoms coordinates");
  for (std::set<QM_Link>::const_iterator it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    atom_nums.push_back(it->atomic_number);

    // convert positions
    math::Vec nn_pos = it->pos * length_to_nn;
    py::list atom_coordinates;
    atom_coordinates.attr("append")(nn_pos[0]);
    atom_coordinates.attr("append")(nn_pos[1]);
    atom_coordinates.attr("append")(nn_pos[2]);
    system_coordinates.attr("append")(atom_coordinates);
    DEBUG(15,"Capping atom to NN " << it->qm_index << "-" << it->mm_index << ": " << it->atomic_number << " " << math::v2s(nn_pos));
  }  

  // Convert std::vector<uint32_t> atom_nums to Python list
  py::list atomic_numbers = py::cast(atom_nums);

  // Get time step
  py::int_ step = sim.steps();

  // Run the method calculate_next_step to predict energy, forces and NN validation
  mlp_calculator.attr("calculate_next_step")(atomic_numbers, system_coordinates, step);
  
  // Store predicted energy
  const double energy = mlp_calculator.attr("get_energy")().cast<double>() * this->param->unit_factor_energy;
  DEBUG(13, "energy from NN, " << energy);
  qm_zone.QM_energy() = energy;

  // Store predicted forces
  std::vector<std::vector<double>> forces = mlp_calculator.attr("get_forces")().cast<std::vector<std::vector<double>>>();
  it = qm_zone.qm.begin();
  int QMBR_zone_size = 0;
  for (unsigned i = 0; it != to; ++it, ++i) {
    it->force[0] = forces[i][0];
    it->force[1] = forces[i][1];
    it->force[2] = forces[i][2];
    it->force *= this->param->unit_factor_force;
    DEBUG(15, "force from NN, atom " << it->index << " : " << math::v2s(it->force));
    QMBR_zone_size += 1;
  }

  // Also parse capping atoms
  unsigned i = 0;
  for(std::set<QM_Link>::iterator it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    it->force[0] = forces[QMBR_zone_size+i][0];
    it->force[1] = forces[QMBR_zone_size+i][1];
    it->force[2] = forces[QMBR_zone_size+i][2];
    it->force *= this->param->unit_factor_force;
    DEBUG(13, "force from NN, capping atom " << it->qm_index << "-" << it->mm_index << ": "  << math::v2s(it->force));
    i += 1;
  }

  // Free energy derivative if perturbation is performed
  if (sim.param().perturbation.perturbation) {
    const double energy_derivative = mlp_calculator.attr("get_derivative")().cast<double>() * this->param->unit_factor_energy;
    qm_zone.QM_energy_derivative() = energy_derivative;
  }
  
  // NN validation
  if (!sim.param().qmmm.nn.val_model_paths.empty()
      && (sim.steps() % sim.param().qmmm.nn.val_steps == 0
       || sim.steps() % sim.param().write.energy == 0)) {
    
    // Store NN valid. deviation
    const double nn_valid_ene = mlp_calculator.attr("get_nn_valid_ene")().cast<double>() * this->param->unit_factor_energy;
    conf.current().energies.nn_valid = nn_valid_ene;

    // Store NN valid. maximum force committee disagreement among all atoms if requested in .qmmm input file
    if(sim.param().qmmm.nn.nnvalid == simulation::nn_valid_maxF) {
      const double nn_valid_maxF = mlp_calculator.attr("get_nn_valid_maxF")().cast<double>() * this->param->unit_factor_force;
      conf.current().energies.nn_valid_maxF = nn_valid_maxF;
    }

    if (fabs(nn_valid_ene) > sim.param().qmmm.nn.val_thresh) {
        std::ostringstream msg;
        msg << "Deviation from validation model above threshold in step " << sim.steps() << " : " << nn_valid_ene;
        io::messages.add(msg.str(), this->name(), io::message::notice);
    }
  }
#endif
  return 0;
}
