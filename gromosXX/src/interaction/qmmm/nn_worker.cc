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
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

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


namespace py = pybind11;
using namespace py::literals;

interaction::NN_Worker::NN_Worker() : QM_Worker("NN Worker"),
                                      param(nullptr),
                                      guard(),
                                      ml_calculator(),
                                      val_calculator() {};

int interaction::NN_Worker::init(simulation::Simulation& sim) {
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.nn);
  QM_Worker::param = this->param;

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

  /** Load the model */
  py::str model_path = sim.param().qmmm.nn.model_path;
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
  py::object ml_model = py_modules["torch"].attr("load")(model_path,"map_location"_a=py_modules["torch"].attr("device")(device));
  
  
  //py::object conn = ase.attr("db").attr("connect")(db_path);
  // "map_location" parameter is to convert model trained on gpu to run on cpu, but seems not to work for my model

  /** Initialize the ML calculator */
  ml_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(ml_model, "energy"_a="energy", "forces"_a="forces", "device"_a=device);

  if (!sim.param().qmmm.nn.val_model_path.empty()) {
    py::str val_model_path = sim.param().qmmm.nn.val_model_path;
    py::object val_model = py_modules["torch"].attr("load")(val_model_path,"map_location"_a=py_modules["torch"].attr("device")(device));
    val_calculator = py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(val_model, "energy"_a="energy", "forces"_a="forces", "device"_a=device);
  }

  #ifdef OMP
    omp_set_num_threads(num_threads);
  #endif

  return 0;
}

interaction::NN_Worker::~NN_Worker() {}

int interaction::NN_Worker::run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone) {
  // run NN interface
  // create the molecule object
  py::object molecule = py_modules["ase"].attr("Atoms")();
  double length_to_nn = 1 / this->param->unit_factor_length;
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
    molecule.attr("append")(atom);
  }

  // Run the calculator
  molecule.attr("set_calculator")(ml_calculator);
  
  // Write the energy
  double energy = molecule.attr("get_potential_energy")().cast<double>();
  qm_zone.QM_energy() = energy * this->param->unit_factor_energy;
 

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
  it = qm_zone.qm.begin();
  for (unsigned i = 0; it != to; ++it, ++i) {
    it->force[0] = molecule.attr("get_forces")().attr("item")(i,0).cast<double >();
    it->force[1] = molecule.attr("get_forces")().attr("item")(i,1).cast<double >();
    it->force[2] = molecule.attr("get_forces")().attr("item")(i,2).cast<double >();
    it->force *= this->param->unit_factor_force;
    DEBUG(15, "force from NN, atom " << it->index << " : " << math::v2s(it->force));
  }

  // Run validation, if asked for
  if (!sim.param().qmmm.nn.val_model_path.empty()
      && (sim.steps() % sim.param().qmmm.nn.val_steps == 0
       || sim.steps() % sim.param().write.energy == 0)) {
    //py::object val_molecule(molecule); we don't need to create a new (reference to a) molecule 
    molecule.attr("set_calculator")(val_calculator);
    // Energy of validation model
    double val_energy = molecule.attr("get_potential_energy")().cast<double>();
    double dev = (energy - val_energy) * this->param->unit_factor_energy;
    conf.current().energies.nn_valid = dev;
    DEBUG(7, "Deviation from validation model: " << dev);
    if (fabs(dev) > sim.param().qmmm.nn.val_thresh) {
      std::ostringstream msg;
      msg << "Deviation from validation model above threshold in step " << sim.steps() << " : " << dev;
      io::messages.add(msg.str(), this->name(), io::message::warning);
    }
  }

  return 0;
}
