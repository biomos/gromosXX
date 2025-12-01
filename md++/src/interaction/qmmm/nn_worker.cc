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

  // Initialize NN interface Initialize pybind, schnetpack, python script

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

  // Determine the total_charge
  int system_charge = sim.param().qmmm.qm_zone.charge + sim.param().qmmm.buffer_zone.charge;
  py::int_ total_charge = system_charge;

  // Determine the spin multiplicity
  // number of unpaired spins of the QM zone
  const int spin_qm = sim.param().qmmm.qm_zone.spin_mult - 1;
  // number of unpaired spins of the buffer zone
  const int spin_buf = sim.param().qmmm.buffer_zone.spin_mult - 1;
  // consider no spin pairing between the QM and buffer zone
  int spin_multiplicity = spin_qm + spin_buf + 1;
  py::int_ spin_mult = spin_multiplicity;

  // Get schnet_worker file paths relative to the binary assuming binary in BUILD_X/bin
  char result[PATH_MAX];
  ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
  std::string exePath = std::string(result, (count > 0) ? count : 0);
  std::string exeDir = exePath.substr(0, exePath.find_last_of("/"));
  std::string modulePath = exeDir + "/../share/contrib";
  py::module_ sys = py::module_::import("sys");
  sys.attr("path").attr("append")(modulePath);

  // Initialize mlp_calculator Python object
  if (software == simulation::qm_schnetv1) {
    // Initialize schnet_v1 module
    py::module_ schnet_v1 = py::module_::import("schnet_v1");
    
    // Check if dynamic charges are requested
    if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
      io::messages.add("Dynamic charges not implemented with schnet_v1 only schnet_v2", io::message::error);
    }

    if(sim.param().qmmm.nn.nnvalid == simulation::nn_valid_maxF) {
      io::messages.add("nn_valid_maxF not implemented with schnet_v1 only schnet_v2", io::message::error);
    }

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
      // Import schnet_v2 module
      py::module_ schnet_v2 = py::module_::import("schnet_v2");
      // Access the SchNet_V2_Calculator class
      py::object schnet_class = schnet_v2.attr("SchNet_V2_Calculator");
      // Check if model was trained on partial charges
      py::object model_trained_on_charges_func = schnet_class.attr("model_trained_on_partial_charges");
      bool model_has_charges = model_trained_on_charges_func(model_path).cast<bool>();

      // If dynamic charges requested but model not trained for charges, throw error
      if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic && !model_has_charges) {
          io::messages.add(
              "Dynamic QM charges requested but model was NOT trained to predict partial charges - check that your charge_key of your model is called charges",
              "NN_Worker", io::message::error
          );
      }

      // Check if the model has electronic embedding
      py::object model_trained_with_elecEmb = schnet_class.attr("model_trained_with_electronic_embedding");
      bool model_has_elecEmbed = model_trained_with_elecEmb(model_path).cast<bool>();

      if (model_has_elecEmbed) {
        std::ostringstream msg;
        msg << "Model is using electronic + nuclear embedding with Total Charge: " 
         << system_charge << " and Spin Multiplicity: " << spin_mult << ".";
        io::messages.add(msg.str(), "NN_Worker", io::message::notice);
      }
      else {
        io::messages.add("Model is using nuclear embedding",
              "NN_Worker", io::message::notice);
      }

      // Initialize the SchNet_V2_Calculator
      mlp_calculator = schnet_class(model_path, val_models_paths, write_val_step, write_energy_step, spin_mult, total_charge);
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

  // Assign dynamic charges for IR+BR
  if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
      std::vector<double> partial_charges =
          mlp_calculator.attr("get_charges")().cast<std::vector<double>>();

      const size_t n_qm    = qm_zone.qm.size();
      const size_t n_links = qm_zone.link.size();
      const size_t n_pred  = partial_charges.size();
      const double cha_to_mm = this->param->unit_factor_charge;

      // Assign QM atom charges
      double tot_qm_charge = 0.0;
      size_t idx = 0;

      for (auto &atom : qm_zone.qm) {
          if (idx >= n_pred) break;
          atom.qm_charge = partial_charges[idx] * cha_to_mm;
          tot_qm_charge += atom.qm_charge;
          ++idx;
      }

      // Assign LA charges
      double tot_link_charge = 0.0;

      for (auto &link : qm_zone.link) {
          if (idx >= n_pred) {
              io::messages.add(
                  "NN model returned fewer link-atom charges than expected; "
                  "remaining link atoms receive no predicted charge.",
                  "NN_Worker", io::message::warning);
              break;
          }
          link.qm_charge = partial_charges[idx] * cha_to_mm;
          tot_link_charge += link.qm_charge;
          ++idx;
      }

      // Total charges
      const double total_predicted_charge = tot_qm_charge + tot_link_charge;

      const double system_charge =
          sim.param().qmmm.qm_zone.charge +
          sim.param().qmmm.buffer_zone.charge;

      DEBUG(10, "NN charge summary: requested=" << system_charge
                << " predicted=" << total_predicted_charge
                << " (QM=" << tot_qm_charge
                << ", LA=" << tot_link_charge << ")");

      // Homogeneous background charge correction
      const double diff = system_charge - total_predicted_charge;

      if (system_charge != total_predicted_charge) {
          const size_t n_receivers = n_qm + n_links;

          const double delta = diff / static_cast<double>(n_receivers);

          DEBUG(10, "Applying homogeneous correction of " << delta
                      << " to all QM and link atoms (" << n_receivers
                      << " atoms).");

          // Apply correction to QM atoms
          for (auto &atom : qm_zone.qm) {
              atom.qm_charge += delta;
          }

          // Apply correction to link atoms
          for (auto &link : qm_zone.link) {
              link.qm_charge += delta;
          }

      } else {
          DEBUG(10, "Predicted charge matches requested charge; no correction applied.");
      }
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
