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

  // Special-trajectory/adaptive-sampling controls used for phi_static output.
  py::int_ ntwse = sim.param().write.energy_index;
  py::float_ val_thresh = py::cast(
      sim.param().qmmm.nn.val_thresh / this->param->unit_factor_energy);

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
      py::list perturbed_atomic_numbers_A =
          py::cast(sim.param().qmmm.nn.pertqm_atomic_number_A);
      py::list perturbed_atomic_numbers_B =
          py::cast(sim.param().qmmm.nn.pertqm_atomic_number_B);

      // Initialize mlp_calculator Pert_SchNet_V2_Calculator Python object
      mlp_calculator = schnet_v2.attr("Pert_SchNet_V2_Calculator")(
          model_path, val_models_paths, write_val_step, write_energy_step,
          spin_mult, total_charge, lambda, perturbed_qm_states,
          "endpoint_atomic_numbers_A"_a=perturbed_atomic_numbers_A,
          "endpoint_atomic_numbers_B"_a=perturbed_atomic_numbers_B,
          "endpoint_total_charge_B"_a=(sim.param().qmmm.qm_zone.pert_charge +
                                       sim.param().qmmm.buffer_zone.pert_charge),
          "endpoint_spin_multiplicity_B"_a=(sim.param().qmmm.qm_zone.pert_spin_mult +
                                            sim.param().qmmm.buffer_zone.pert_spin_mult - 1),
          "ntwse"_a=ntwse, "val_thresh"_a=val_thresh);
    }

    else {
      // Initialize SchNet_V2_Calculator
      mlp_calculator = schnet_class(
          model_path, val_models_paths, write_val_step, write_energy_step,
          spin_mult, total_charge,
          "ntwse"_a=ntwse, "val_thresh"_a=val_thresh);
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

// --- copy from q_equilibration.cc (keep identical for consistency) ---
static inline double pair_potential_q_over_r(const math::Vec& qm_pos,
                                            const interaction::MM_Atom& mm_atom)
  {
    if (mm_atom.is_polarisable) {
      double pot = (mm_atom.charge - mm_atom.cos_charge) / math::abs(qm_pos - mm_atom.pos);
      pot += mm_atom.cos_charge / math::abs(qm_pos - mm_atom.pos - mm_atom.cosV);
      return pot;
    } else {
      return mm_atom.charge / math::abs(qm_pos - mm_atom.pos);
    }
  }

  static inline void build_excluded_mm_atoms(const topology::Topology& topo,
                                            const simulation::Simulation& sim,
                                            const interaction::QM_Zone& qm_zone,
                                            const interaction::QM_Atom& qm_atom,
                                            std::set<unsigned>& excluded_mm)
  {
    excluded_mm.clear();
    for (auto li_it = qm_zone.link.begin(); li_it != qm_zone.link.end(); ++li_it) {
      if (li_it->qm_index == qm_atom.index) excluded_mm.insert(li_it->mm_index);
    }

    if (sim.param().qmmm.mopac.link_atom_mode == 2) {
      std::set<unsigned> excluded_cgs;
      for (auto it = excluded_mm.begin(); it != excluded_mm.end(); ++it) {
        int cg_idx = -1;
        for (unsigned cg = 0; cg < topo.num_chargegroups(); ++cg) {
          if (*it < unsigned(topo.chargegroup(cg + 1))) { cg_idx = int(cg); break; }
        }
        assert(cg_idx != -1);
        excluded_cgs.insert(unsigned(cg_idx));
      }
      for (auto it = excluded_cgs.begin(); it != excluded_cgs.end(); ++it) {
        for (int a = topo.chargegroup(*it); a < topo.chargegroup(*it + 1); ++a) {
          excluded_mm.insert(unsigned(a));
        }
      }
    }
  }

  static inline double total_potential_q_over_r(const topology::Topology& topo,
                                              const simulation::Simulation& sim,
                                              const interaction::QM_Zone& qm_zone,
                                              const interaction::QM_Atom& qm_atom)
  {
    double pot = 0.0;

    if (!qm_atom.is_linked) {
      for (auto mm_it = qm_zone.mm.begin(); mm_it != qm_zone.mm.end(); ++mm_it) {
        pot += pair_potential_q_over_r(qm_atom.pos, *mm_it);
      }
      return pot;
    }

    if (sim.param().qmmm.mopac.link_atom_mode == 0) return 0.0;

    std::set<unsigned> excluded_mm;
    build_excluded_mm_atoms(topo, sim, qm_zone, qm_atom, excluded_mm);

    for (auto mm_it = qm_zone.mm.begin(); mm_it != qm_zone.mm.end(); ++mm_it) {
      if (excluded_mm.find(mm_it->index) == excluded_mm.end()) {
        pot += pair_potential_q_over_r(qm_atom.pos, *mm_it);
      }
    }
    return pot;
  }

interaction::NN_Worker::~NN_Worker() = default;

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
  std::vector<std::pair<const interaction::MM_Atom*, bool> > or_sites;

  if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    // Dynamic/B2: build OR arrays + cutoff and pass them
    std::vector<std::vector<double>> or_coordinates_nm;
    std::vector<double> or_charges_e;

    or_coordinates_nm.reserve(2 * qm_zone.mm.size());
    or_charges_e.reserve(2 * qm_zone.mm.size());
    or_sites.reserve(2 * qm_zone.mm.size());

    for (auto mm_it = qm_zone.mm.begin(); mm_it != qm_zone.mm.end(); ++mm_it) {
      // Positions are already in nm (GROMOS internal). Polarizable atoms are
      // represented by two external charge sites, as in the other QM workers.
      if (mm_it->is_polarisable) {
        const double atom_charge = mm_it->charge - mm_it->cos_charge;
        or_coordinates_nm.push_back({mm_it->pos[0], mm_it->pos[1], mm_it->pos[2]});
        or_charges_e.push_back(atom_charge);
        or_sites.push_back(std::make_pair(&(*mm_it), false));

        const math::Vec cos_pos = mm_it->pos + mm_it->cosV;
        or_coordinates_nm.push_back({cos_pos[0], cos_pos[1], cos_pos[2]});
        or_charges_e.push_back(mm_it->cos_charge);
        or_sites.push_back(std::make_pair(&(*mm_it), true));
      } else {
        or_coordinates_nm.push_back({mm_it->pos[0], mm_it->pos[1], mm_it->pos[2]});
        or_charges_e.push_back(mm_it->charge); // e
        or_sites.push_back(std::make_pair(&(*mm_it), false));
      }
    }

    const double cutoff_nm = sim.param().qmmm.cutoff;
    DEBUG(10, "NN B2 OR list: mm_atoms=" << qm_zone.mm.size()
              << " or_sites=" << or_sites.size()
              << " link_atoms=" << n_caps
              << " cutoff_nm=" << cutoff_nm);

    mlp_calculator.attr("calculate_next_step")(
        atomic_numbers,
        system_coordinates,
        step,
        "dynamic_charges"_a=true,
        "or_positions_nm"_a=or_coordinates_nm,
        "or_charges_e"_a=or_charges_e,
        "cutoff_nm"_a=cutoff_nm,
        "n_link_atoms"_a=n_caps
    );
  }

  else {
    // Legacy: no OR arrays needed
    mlp_calculator.attr("calculate_next_step")(
        atomic_numbers,
        system_coordinates,
        step,
        "dynamic_charges"_a=false
    );
  }

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
    const double energy_derivative = mlp_calculator.attr("get_derivative")().cast<double>() * this->param->unit_factor_energy;
    qm_zone.QM_energy_derivative() = energy_derivative;
  }

  // Assign dynamic charges for IR+BR
  if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {

    // fetch OR forces due to Qeq
    auto or_forces = mlp_calculator.attr("get_or_forces")()
        .cast<std::vector<std::vector<double>>>();

    if (or_forces.size() != or_sites.size()) {
      throw std::runtime_error("OR force size mismatch");
    }

    for (unsigned site = 0; site < or_sites.size(); ++site) {
      const interaction::MM_Atom* mm = or_sites[site].first;
      const bool is_cos_site = or_sites[site].second;

      math::Vec force(or_forces[site][0], or_forces[site][1], or_forces[site][2]);
      force *= this->param->unit_factor_force;

      if (is_cos_site) {
        mm->cos_force += force;
        DEBUG(15, "force from NN, OR COS atom " << mm->index << " : " << math::v2s(mm->cos_force));
      } else {
        mm->force += force;
        DEBUG(15, "force from NN, OR atom " << mm->index << " : " << math::v2s(mm->force));
      }
    }

    {
      std::vector<double> partial_charges = mlp_calculator.attr("get_charges")().cast<std::vector<double>>();
      // Assign QM atom charges in the same order as sent (IR+BR).
      double tot_qm_charge = 0.0;
      for (int i = 0; i < qmb_size; ++i) {
        const interaction::QM_Atom* a = qm_atoms_order[i];
        a->qm_charge = partial_charges[i] * this->param->unit_factor_charge;
        tot_qm_charge += a->qm_charge;
        DEBUG(10, "qm_charge from NN, atom " << a->index << " : " << a->qm_charge);
      }
      // Link/cap charges are internal B2/QEq quantities only. Their
      // electrostatic force contribution is already included in it->force and
      // redistributed by the link-atom machinery, so do not assign qm_charge
      // back to the fictitious link atoms.
      double tot_link_charge = 0.0;
      int i_cap = 0;
      for (auto it = qm_zone.link.begin(); it != qm_zone.link.end(); ++it, ++i_cap) {
        const double link_charge = partial_charges[qmb_size + i_cap] * this->param->unit_factor_charge;
        tot_link_charge += link_charge;
        DEBUG(15, "qm_charge from NN, capping atom " << it->qm_index << "-" << it->mm_index
                  << " : " << link_charge << " (not assigned)");
      }

      double total_predicted_charge = tot_qm_charge + tot_link_charge;

      double system_charge = sim.param().qmmm.qm_zone.charge +
                             sim.param().qmmm.buffer_zone.charge;
      if (sim.param().perturbation.perturbation) {
        const double charge_B = sim.param().qmmm.qm_zone.pert_charge +
                                sim.param().qmmm.buffer_zone.pert_charge;
        system_charge += sim.param().perturbation.lambda * (charge_B - system_charge);
      }

      DEBUG(10, "NN charge summary: requested=" << system_charge
                << " predicted=" << total_predicted_charge
                << " (QM=" << tot_qm_charge
                << ", LA=" << tot_link_charge << ")");
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
