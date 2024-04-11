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

//#include "../../../math/periodicity.h"

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
                                      ml_calculator(),
                                      val_calculator() {};

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
  std::vector<std::string> modules = {
      "sys"
    , "ase"
    , "os"
    , "torch"
    , "logging"
    , "schnetpack"
  };
  #ifdef OMP
    // Python module modifies omp_num_threads value, we gonna restore it
    unsigned num_threads;
    #pragma omp parallel
      #pragma omp single
        num_threads = omp_get_num_threads();
  #endif
  // Initialize Python modules
  for (size_t i = 0; i < modules.size(); ++i) {
    py_modules.emplace(modules[i], py::module::import(modules[i].c_str()));
  }
  // Hint: use them as py_modules["sys"].attr()

  // Print modules versions
  std::string version;
  version = py_modules["sys"].attr("version").cast<std::string>();
    io::messages.add("Python version: " + version
            , this->name(), io::message::notice);
  version = py_modules["ase"].attr("__version__").cast<std::string>();
    io::messages.add("ASE module version: " + version
            , this->name(), io::message::notice);
  version = py_modules["torch"].attr("__version__").cast<std::string>();
    io::messages.add("Torch module version: " + version
            , this->name(), io::message::notice);
  version = py_modules["schnetpack"].attr("__file__").cast<std::string>();
    io::messages.add("Schnetpack module path: " + version
            , this->name(), io::message::notice);

  std::string device;
  switch(sim.param().qmmm.nn.device) {
    case simulation::nn_device_auto: {
      const bool cuda_available = py::bool_(py_modules["torch"].attr("cuda").attr("is_available")());
      if (cuda_available) {
        device = "cuda";
        DEBUG(1, "auto: using CUDA");
      } else {
        device = "cpu";
        DEBUG(1, "auto: Using CPU");
      }
      break;
    }
    case simulation::nn_device_cuda:
      device = "cuda";
        DEBUG(1, "Using CUDA");
      break;
    case simulation::nn_device_cpu:
      device = "cpu";
        DEBUG(1, "Using CPU");
      break;
    default:
      io::messages.add("Unknown NN device", this->name(), io::message::critical);
      return 1;
  }
  /** Load the model */
  py::str model_path = sim.param().qmmm.nn.model_path;
  DEBUG(11, "model_path: " << model_path.cast<std::string>());
  py::str ml_args_path = py_modules["os"].attr("path").attr("join")(py::cast('/').attr("join")(model_path.attr("split")('/')[py::slice(0,-1,1)]), "args.json");
  DEBUG(11, "ml_args_path: " << ml_args_path.cast<std::string>());
  py::object ml_model_args = py_modules["schnetpack"].attr("utils").attr("read_from_json")(ml_args_path);
  py::object ml_model = py_modules["torch"].attr("load")(model_path,"map_location"_a=py_modules["torch"].attr("device")(device));
  
  if (py::bool_(ml_model_args.attr("parallel"))) {
    ml_model = ml_model.attr("module");
  }

  py::object ml_environment = py_modules["schnetpack"].attr("utils").attr("script_utils").attr("settings").attr("get_environment_provider")(ml_model_args, "device"_a=device);

  std::vector<int> mlmm_indices;

  /** Initialize the ML calculator */
  if(sim.param().qmmm.nn.learning_type == simulation::nn_learning_type_all) {
    ml_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(ml_model, "energy"_a="energy", "forces"_a="forces", "device"_a=device, "environment_provider"_a=ml_environment);
  } else if(sim.param().qmmm.nn.learning_type == simulation::nn_learning_type_qmonly) {
    // this should be filled with all the atoms that are in the qm_zone, but not in the buffer zone.
    // currently it is supported only for one, the first atom
    // for multi-atom QM region, the atoms should be ordered such that the QM atoms are always first in the list
    mlmm_indices.push_back(0);
    ml_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(ml_model, "energy"_a="energy", "forces"_a="forces", "device"_a=device, "environment_provider"_a=ml_environment, "mlmm"_a=mlmm_indices);
  } else {
      std::ostringstream msg;
      msg << "Learning type value unknown " << sim.param().qmmm.nn.learning_type;
      io::messages.add(msg.str(), this->name(), io::message::error);
  }
  if (!sim.param().qmmm.nn.val_model_path.empty()) {
    py::str val_model_path = sim.param().qmmm.nn.val_model_path;
    py::str val_args_path = py_modules["os"].attr("path").attr("join")(py::cast('/').attr("join")(val_model_path.attr("split")('/')[py::slice(0,-1,1)]), "args.json");
    py::object val_model_args = py_modules["schnetpack"].attr("utils").attr("read_from_json")(val_args_path);
    py::object val_model = py_modules["torch"].attr("load")(val_model_path,"map_location"_a=py_modules["torch"].attr("device")(device));
    py::object val_environment = py_modules["schnetpack"].attr("utils").attr("script_utils").attr("settings").attr("get_environment_provider")(val_model_args, "device"_a=device);
    if (py::bool_(val_model_args.attr("parallel"))) {
      val_model = val_model.attr("module");
    }
    if(sim.param().qmmm.nn.learning_type == simulation::nn_learning_type_all) {
      val_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(val_model, "energy"_a="energy", "forces"_a="forces", "device"_a=device, "environment_provider"_a=val_environment);
    } else if(sim.param().qmmm.nn.learning_type == simulation::nn_learning_type_qmonly) {
      val_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(val_model, "energy"_a="energy", "forces"_a="forces", "device"_a=device, "environment_provider"_a=val_environment, "mlmm"_a=mlmm_indices);
    } 
  }   
   
  if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    py::str charge_model_path = sim.param().qmmm.nn.charge_model_path;
    DEBUG(11, "charge_model_path: " << charge_model_path.cast<std::string>());
    py::str charge_args_path = py_modules["os"].attr("path").attr("join")(py::cast('/').attr("join")(charge_model_path.attr("split")('/')[py::slice(0,-1,1)]), "args.json");
    DEBUG(11, "charge_args_path: " << charge_args_path.cast<std::string>());
    py::object charge_model_args = py_modules["schnetpack"].attr("utils").attr("read_from_json")(charge_args_path);
    py::object charge_model = py_modules["torch"].attr("load")(charge_model_path,"map_location"_a=py_modules["torch"].attr("device")(device));
    py::object charge_environment = py_modules["schnetpack"].attr("utils").attr("script_utils").attr("settings").attr("get_environment_provider")(charge_model_args, "device"_a=device);
    if (py::bool_(charge_model_args.attr("parallel"))) {
      charge_model = charge_model.attr("module");
    }
    charge_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(charge_model, "charges"_a="charges", "device"_a=device, "environment_provider"_a=charge_environment);
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
  // create the molecule object
  py::object molecule_1 = py_modules["ase"].attr("Atoms")();
  py::object molecule_2 = py_modules["ase"].attr("Atoms")();
  double length_to_nn = 1 / this->param->unit_factor_length;

  // perturbed states - specify indices of state A and B
  // To state A, we take everything between <stateA_first ; stateA_last>
  // To state B, we take everything between <stateB_first ; stateB_last>
  // both states contain also atoms <both_states_first ; +inf)
  const unsigned stateA_first = 0;
  const unsigned stateA_last = 5;
  const unsigned stateB_first = 6;
  const unsigned stateB_last = 10;
  const unsigned both_states_first = 11;

  std::set<interaction::QM_Atom>::const_iterator it, to;
  it = qm_zone.qm.begin();
  to = qm_zone.qm.end();
  for (;it != to; ++it) {
    //unsigned atomic_number = it->atomic_number;
    uint32_t atomic_number = it->atomic_number;
    py::list py_coordinates;
    math::Vec nn_pos = it->pos * length_to_nn;
    DEBUG(15, "atom to NN: " << it->index << " : " << math::v2s(nn_pos));
    py_coordinates.attr("append")(nn_pos[0]);
    py_coordinates.attr("append")(nn_pos[1]);
    py_coordinates.attr("append")(nn_pos[2]);
    py::object atom = py_modules["ase"].attr("Atom")("symbol"_a = atomic_number, "position"_a = py_coordinates);
    const bool is_state_A = stateA_first <= it->index && it->index <= stateA_last;
    const bool is_state_B = stateB_first <= it->index && it->index <= stateB_last;
    const bool is_both_states = it->index >= both_states_first;
    if (is_state_A || is_both_states) {
      molecule_1.attr("append")(atom); 
    }
    if (is_state_B || is_both_states) {
      molecule_2.attr("append")(atom);
    }
    //molecule.attr("append")(atom); //original
  }

  const double lambda = sim.param().perturbation.lambda;

  // Run the calculator
  molecule_1.attr("set_calculator")(ml_calculator);
  molecule_2.attr("set_calculator")(ml_calculator);
  
  const double energy_1 = molecule_1.attr("get_potential_energy")().cast<double>() * this->param->unit_factor_energy;
  const double energy_2 = molecule_2.attr("get_potential_energy")().cast<double>() * this->param->unit_factor_energy;

  // Write the energy
  const double energy = (1-lambda) * energy_1 + lambda * energy_2;
  qm_zone.QM_energy() = energy;
 

  // Get the forces
  /*const size_t qm_size = qm_zone.qm.size();
  math::VArray forces;
  forces.resize(qm_size);
  std::array<double, 3> tmp;
  for (unsigned i = 0; i < qm_size; ++i) {
    //double * ptr = &(forces(i)[0]);
    forces(i) = molecule.attr("get_forces")(i).cast<std::vector<double> >();
    //forces(i)[0] = tmp[0];
    //forces(i)[1] = tmp[1];
    //forces(i)[2] = tmp[2];
  }*/
  
  // Write the forces
  std::vector<std::vector<double>> forces_1 = molecule_1.attr("get_forces")().cast<std::vector<std::vector<double>>>();
  std::vector<std::vector<double>> forces_2 = molecule_2.attr("get_forces")().cast<std::vector<std::vector<double>>>();
  it = qm_zone.qm.begin();
  for (unsigned i = 0, j = 0; it != to; ++it) {
    const bool is_state_A = stateA_first <= it->index && it->index <= stateA_last;
    const bool is_state_B = stateB_first <= it->index && it->index <= stateB_last;
    const bool is_both_states = it->index >= both_states_first;

    math::Vec force_1, force_2;
    if (is_state_A || is_both_states) {
      force_1[0] = forces_1[i][0];//molecule_1.attr("get_forces")().attr("item")(i,0).cast<double>();
      force_1[1] = forces_1[i][1];//molecule_1.attr("get_forces")().attr("item")(i,1).cast<double>();
      force_1[2] = forces_1[i][2];//molecule_1.attr("get_forces")().attr("item")(i,2).cast<double>();
      ++i;
    }
    if (is_state_B || is_both_states) {
      force_2[0] = forces_2[j][0];//molecule_2.attr("get_forces")().attr("item")(j,0).cast<double>();
      force_2[1] = forces_2[j][1];//molecule_2.attr("get_forces")().attr("item")(j,1).cast<double>();
      force_2[2] = forces_2[j][2];//molecule_2.attr("get_forces")().attr("item")(j,2).cast<double>();
      ++j;
    }
    it->force = (1 - lambda) * force_1 + lambda * force_2;
    it->force *= this->param->unit_factor_force;
    DEBUG(15, "force from NN, atom " << it->index << " : " << math::v2s(it->force));
  }

  // Run validation, if asked for
  if (!sim.param().qmmm.nn.val_model_path.empty()
      && (sim.steps() % sim.param().qmmm.nn.val_steps == 0
       || sim.steps() % sim.param().write.energy == 0)) {
    //py::object val_molecule(molecule); we don't need to create a new (reference to a) molecule 
    molecule_1.attr("set_calculator")(val_calculator);
    molecule_2.attr("set_calculator")(val_calculator);
    // Energy of validation model
    const double val_energy_1 = molecule_1.attr("get_potential_energy")().cast<double>() * this->param->unit_factor_energy;
    const double val_energy_2 = molecule_2.attr("get_potential_energy")().cast<double>() * this->param->unit_factor_energy;
    const double val_energy = (1-lambda) * val_energy_1 + lambda * val_energy_2;
    const double dev_1 = energy_1 - val_energy_1;
    const double dev_2 = energy_2 - val_energy_2;
    const double dev_overall = energy - val_energy;
    const double dev = fmax(fabs(dev_1), fabs(dev_2));
    conf.current().energies.nn_valid = dev;
    DEBUG(7, "Deviation from validation model: " << dev);
    if (fabs(dev) > sim.param().qmmm.nn.val_thresh) {
      std::ostringstream msg;
      msg << "Deviation from validation model above threshold in step " << sim.steps() << " : " << dev_overall << ", molecule 1 deviation: " << dev_1 << ", molecule 2 deviation: " << dev_2;
      io::messages.add(msg.str(), this->name(), io::message::notice); // Changed to notice
      //if(sim.param().qmmm.nn.val_forceconstant != 0.0){
      //  // add a biasing force between the two NN networks
      //  double dev_squared = (dev*dev - sim.param().qmmm.nn.val_thresh * sim.param().qmmm.nn.val_thresh);
      //  qm_zone.QM_energy() += 0.25 * sim.param().qmmm.nn.val_forceconstant * dev_squared * dev_squared;
      //  math::Vec val_force;
      //  it = qm_zone.qm.begin();
      //  for (unsigned i = 0; it != to; ++it, ++i) {
      //    val_force[0] = molecule.attr("get_forces")().attr("item")(i,0).cast<double >();
      //    val_force[1] = molecule.attr("get_forces")().attr("item")(i,1).cast<double >();
      //    val_force[2] = molecule.attr("get_forces")().attr("item")(i,2).cast<double >();
      //    val_force *= this->param->unit_factor_force;
      //    it->force += sim.param().qmmm.nn.val_forceconstant * dev_squared * dev * (it->force - val_force);
      //  }
      //}
    }
  }
  // Assign charges for the MM calculation, if asked for
  //if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic
  //    && sim.steps() % sim.param().qmmm.nn.charge_steps == 0) {
  //  //py::object val_molecule(molecule); we don't need to create a new (reference to a) molecule 
  //  molecule.attr("set_calculator")(charge_calculator);
  //  // collect the charges
  //  it = qm_zone.qm.begin();
  //  double totcharge=0.0;
  //  for (unsigned i = 0; it != to; ++it, ++i) {
  //    it->qm_charge = molecule.attr("get_charges")().attr("item")(i).cast<double >() * this->param->unit_factor_charge;
  //    totcharge+=it->qm_charge;
  //  }
  //  if(totcharge != sim.param().qmmm.qm_zone.charge){
  //    std::ostringstream msg;
  //    msg << "Charges from NN model do not add up " << sim.steps() 
  //        << ": requested " << sim.param().qmmm.qm_zone.charge 
  //        << " predicted " << totcharge << ". Homogeneously adjusting charges";
  //    io::messages.add(msg.str(), this->name(), io::message::warning);
  //    // adjust?
  //    double q_adjust=(sim.param().qmmm.qm_zone.charge - totcharge) / qm_zone.qm.size();
  //    it = qm_zone.qm.begin();
  //    for(; it!=to; ++it){
  //      it->qm_charge += q_adjust;
  //      DEBUG(10, "Charge adjusted for atom " << it->index << " from " << it->qm_charge - q_adjust << " to " << it->qm_charge);
  //    }
  //  }
  //}

#endif
  return 0;
}
