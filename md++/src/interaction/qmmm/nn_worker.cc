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

    // Decide if perturbation is performed or not
    if (sim.param().perturbation.perturbation) {
      // get lambda parameter
      py::float_ lambda = py::cast(sim.param().perturbation.lambda);

      // get perturbed QM states
      py::list perturbed_qm_states = py::cast(sim.param().qmmm.nn.pertqm_state);

      // get perturbed state total charge and spin multiplicity
      int pert_system_charge = sim.param().qmmm.qm_zone.pert_charge + sim.param().qmmm.buffer_zone.pert_charge;
      py::int_ pert_total_charge = pert_system_charge;

      // number of unpaired spins of the QM zone
      const int pert_spin_qm = sim.param().qmmm.qm_zone.pert_spin_mult - 1;
      // number of unpaired spins of the buffer zone
      const int pert_spin_buf = sim.param().qmmm.buffer_zone.pert_spin_mult - 1;
      // consider no spin pairing between the QM and buffer zone
      int pert_spin_multiplicity = pert_spin_qm + pert_spin_buf + 1;
      py::int_ pert_spin_mult = pert_spin_multiplicity;
      
      // Decide if reference vacuum energies should be used according to .qmmm specification file
      py::object ref_vacA_obj = py::none();
      py::object ref_vacB_obj = py::none();

      if (sim.param().qmmm.qm_zone.has_ref_vacA && sim.param().qmmm.qm_zone.has_ref_vacB) {
        ref_vacA_obj = py::float_(sim.param().qmmm.qm_zone.ref_vacA);
        ref_vacB_obj = py::float_(sim.param().qmmm.qm_zone.ref_vacB);
      }

      // Initialize mlp_calculator Pert_SchNet_V2_Calculator Python object
      mlp_calculator = schnet_v2.attr("Pert_SchNet_V2_Calculator")(model_path, val_models_paths, write_val_step, 
        write_energy_step, spin_mult, total_charge, lambda, perturbed_qm_states,pert_total_charge, pert_spin_mult,
        ref_vacA_obj, ref_vacB_obj);
    }

    else {
      // Initialize SchNet_V2_Calculator
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
/**
 * @brief Execute a single NN “QM” step: build NN inputs, call Python, and map outputs back.
 *
 * This routine replaces the traditional QM backend call by a call into a Python
 * calculator (SchNetPack wrappers) exposed via pybind11.
 *
 * ## High-level data flow
 * 1. **Build NN inputs**
 *    - Assemble `atomic_numbers` and `system_coordinates` for the NN model.
 *    - Convert positions from the MD length unit to the NN model unit using
 *      `unit_factor_length`.
 *
 * 2. **Call the Python calculator**
 *    - Invoke `mlp_calculator.calculate_next_step(...)` which performs the NN
 *      forward pass for the current MD step.
 *    - The Python side caches the predicted energy/forces (and optionally charges,
 *      committee disagreement, and dE/dlambda for TI).
 *
 * 3. **Read outputs and redistribute them to GROMOS objects**
 *    - Energy is stored in `qm_zone.QM_energy()`.
 *    - Forces are written into the `QM_Atom::force` fields of the QM/Buffer
 *      atoms and into the `QM_Link::force` fields of capping atoms.
 *    - If dynamic charges are enabled, predicted partial charges are mapped back
 *      to `qm_charge` fields and (optionally) corrected to match the requested
 *      integer total charge.
 *
 * ## Ordering contract (critical for correct force mapping)
 * The NN returns forces in the exact order in which atoms are passed in.
 * Historically, the ordering of `qm_zone.qm` can differ from the ordering used
 * by perturbation input blocks (QMZONE / PERTQMZONE). For perturbation, the
 * *partition* between IR (QMZONE) and BR (BUFFERZONE) must be consistent.
 *
 * We therefore build a deterministic order `qm_atoms_order`:
 *  - First: all atoms that satisfy `topo.is_qm(index)` (IR)
 *  - Second: all remaining atoms in `qm_zone.qm` (BR)
 *  - Finally: all link/capping atoms in `qm_zone.link` (caps)
 *
 * The same order is used to:
 *  - build the Python input arrays
 *  - interpret the Python output force array
 */
int interaction::NN_Worker::run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone) {
#ifdef HAVE_PYBIND11
  // run NN interface 

  // Prepare the input for mlp_calculator object
  double length_to_nn = 1 / this->param->unit_factor_length;

  // --- Build a deterministic order for qm_zone.qm that matches PERTQMZONE (QMZONE/IR) first ---
  std::vector<const interaction::QM_Atom*> qm_atoms_order;
  qm_atoms_order.reserve(qm_zone.qm.size());

  // 1) IR atoms (QMZONE): topo.is_qm(index) == true
  for (auto it = qm_zone.qm.begin(); it != qm_zone.qm.end(); ++it) {
    if (topo.is_qm(it->index)) {
      qm_atoms_order.push_back(&(*it));
    }
  }

  // 2) BR atoms (BUFFERZONE): everything else in qm_zone.qm
  for (auto it = qm_zone.qm.begin(); it != qm_zone.qm.end(); ++it) {
    if (!topo.is_qm(it->index)) {
      qm_atoms_order.push_back(&(*it));
    }
  }

  // Atomic numbers + coordinates in that order
  std::vector<uint32_t> atom_nums;
  atom_nums.reserve(qm_atoms_order.size() + qm_zone.link.size());

  py::list system_coordinates;

  for (const auto* a : qm_atoms_order) {
    atom_nums.push_back(a->atomic_number);

    math::Vec nn_pos = a->pos * length_to_nn;
    py::list atom_coordinates;
    atom_coordinates.attr("append")(nn_pos[0]);
    atom_coordinates.attr("append")(nn_pos[1]);
    atom_coordinates.attr("append")(nn_pos[2]);
    system_coordinates.attr("append")(atom_coordinates);
  }

  // Caps (link atoms) are appended last. This is important because the
  // redistribution below assumes: [QM+BR atoms..., caps...].
  for (auto it = qm_zone.link.begin(); it != qm_zone.link.end(); ++it) {
    atom_nums.push_back(it->atomic_number);

    math::Vec nn_pos = it->pos * length_to_nn;
    py::list atom_coordinates;
    atom_coordinates.attr("append")(nn_pos[0]);
    atom_coordinates.attr("append")(nn_pos[1]);
    atom_coordinates.attr("append")(nn_pos[2]);
    system_coordinates.attr("append")(atom_coordinates);
  }

  py::list atomic_numbers = py::cast(atom_nums);

  py::int_ step = sim.steps();
  py::int_ n_caps = static_cast<int>(qm_zone.link.size());

  // Run the method calculate_next_step to predict energy, forces and NN validation
  mlp_calculator.attr("calculate_next_step")(atomic_numbers, system_coordinates, step, n_caps);

  // Store predicted energy
  const double energy = mlp_calculator.attr("get_energy")().cast<double>() * this->param->unit_factor_energy;
  DEBUG(13, "energy from NN, " << energy);
  qm_zone.QM_energy() = energy;


  // Store predicted forces
  std::vector<std::vector<double>> forces = mlp_calculator.attr("get_forces")().cast<std::vector<std::vector<double>>>();
  // First QM atoms (IR+BR) in the same order we sent
  const int qmb_size = static_cast<int>(qm_atoms_order.size());
  for (int i = 0; i < qmb_size; ++i) {
    const interaction::QM_Atom* a = qm_atoms_order[i];
    // The Python model returns forces in its internal units; we convert back
    // to the MD force unit using unit_factor_force.
    //
    // NOTE: qm_atoms_order stores pointers into qm_zone.qm, so writing to
    // a->force updates the underlying QM_Atom stored in the zone.
    a->force[0] = forces[i][0];
    a->force[1] = forces[i][1];
    a->force[2] = forces[i][2];
    a->force *= this->param->unit_factor_force;
    DEBUG(15, "force from NN, atom " << a->index << " : " << math::v2s(a->force));
  }


  // Now redistribute capping atom forces. Caps start at index qmb_size.
  int i_cap = 0;
  for (auto it = qm_zone.link.begin(); it != qm_zone.link.end(); ++it, ++i_cap) {
    // Link/capping atoms are appended after all physical QM/BR atoms.
    // Their forces are stored on the QM_Link objects (used later to project
    // forces back onto the real MM atoms participating in the link).
    it->force[0] = forces[qmb_size + i_cap][0];
    it->force[1] = forces[qmb_size + i_cap][1];
    it->force[2] = forces[qmb_size + i_cap][2];
    it->force *= this->param->unit_factor_force;
    DEBUG(13, "force from NN, capping atom " << it->qm_index << "-" << it->mm_index << ": "  << math::v2s(it->force));
  }

  // Free energy derivative if perturbation is performed
  if (sim.param().perturbation.perturbation) {
    // dE/dlambda is provided by the Python TI wrapper.
    const double energy_derivative = mlp_calculator.attr("get_derivative")().cast<double>() * this->param->unit_factor_energy;
    qm_zone.QM_energy_derivative() = energy_derivative;
  }

  // Assign dynamic charges for IR+BR (partial charges predicted by the NN)
  //
  // IMPORTANT ORDERING CONTRACT:
  //   The Python calculator returns charges in the same order as coordinates were provided:
  //     [ IR atoms (deterministic) ][ BR atoms (deterministic) ][ link/cap atoms ]
  //   This is the same ordering used above for forces. Therefore we must *assign charges*
  //   via qm_atoms_order (IR+BR) and then via qm_zone.link (caps), NOT via qm_zone.qm's
  //   internal iteration order.
  if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
      std::vector<double> partial_charges = mlp_calculator.attr("get_charges")().cast<std::vector<double>>();

      // Defensive sanity check: expect one charge per (IR+BR) atom plus one per link/cap atom.
      const int n_caps = static_cast<int>(qm_zone.link.size());
      const int n_expected = qmb_size + n_caps;

      // Assign QM atom charges in the same order as sent (IR+BR).
      double tot_qm_charge = 0.0;
      for (int i = 0; i < qmb_size; ++i) {
        const interaction::QM_Atom* a = qm_atoms_order[i];
        a->qm_charge = partial_charges[i] * this->param->unit_factor_charge;
        tot_qm_charge += a->qm_charge;
        DEBUG(10, "qm_charge from NN, atom " << a->index << " : " << a->qm_charge);
      }
      // Assign link/cap charges in the same order as sent (iteration order of qm_zone.link).
      double tot_link_charge = 0.0;
      int i_cap = 0;
      for (auto it = qm_zone.link.begin(); it != qm_zone.link.end(); ++it, ++i_cap) {
        it->qm_charge = partial_charges[qmb_size + i_cap] * this->param->unit_factor_charge;
        tot_link_charge += it->qm_charge;
        DEBUG(15, "qm_charge from NN, capping atom " << it->qm_index << "-" << it->mm_index
                  << " : " << it->qm_charge);
      }

      double total_predicted_charge = tot_qm_charge + tot_link_charge;

      const double system_charge =
          sim.param().qmmm.qm_zone.charge +
          sim.param().qmmm.buffer_zone.charge;

      DEBUG(10, "NN charge summary: requested=" << system_charge
                << " predicted=" << total_predicted_charge
                << " (QM=" << tot_qm_charge
                << ", LA=" << tot_link_charge << ")");

    // Homogeneous background charge correction (if needed).
    // We distribute the difference equally over all (IR+BR) atoms and all cap atoms.
      double diff =  system_charge - total_predicted_charge;

      if (total_predicted_charge != system_charge) {
          double delta = diff /  (qm_zone.qm.size()+qm_zone.link.size());

          DEBUG(10, "Applying homogeneous correction of " << delta
                      << " to all QM and link atoms (" << qm_zone.qm.size()+qm_zone.link.size()
                      << " atoms).");

        // Apply correction to QM atoms (IR+BR) in the same order.
        for (int i = 0; i < qmb_size; ++i) {
          const interaction::QM_Atom* a = qm_atoms_order[i];
          const double old = a->qm_charge;
          a->qm_charge = old + delta;
          DEBUG(10, "Charge adjusted for atom " << a->index
                    << " from " << old << " to " << a->qm_charge);
        }

        // Apply correction to link/cap atoms.
        i_cap = 0;
        for (auto it = qm_zone.link.begin(); it != qm_zone.link.end(); ++it, ++i_cap) {
          const double old = it->qm_charge;
          it->qm_charge = old + delta;
          DEBUG(10, "LA charge adjusted for atom " << it->mm_index
                    << " from " << old << " to " << it->qm_charge);
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
