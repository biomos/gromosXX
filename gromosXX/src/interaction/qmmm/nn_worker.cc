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


#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm


namespace py = pybind11;
using namespace py::literals;

interaction::NN_Worker::NN_Worker() : QM_Worker("NN Worker"), param(nullptr), guard() {};

int interaction::NN_Worker::init(const simulation::Simulation& sim) {
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

  for (size_t i = 0; i < modules.size(); ++i) {
    py_modules.emplace(modules[i], py::module::import(modules[i].c_str()));
  }
  // use as py_modules["sys"].attr()



  //py::print("SPK version:    ",spk.attr("__version__")); //has no version
  py::print("ASE version:    ",py_modules["ase"].attr("__version__"));
  py::print("Torch version:  ",py_modules["torch"].attr("__version__"));
  py::print("Logging version:",py_modules["logging"].attr("__version__"));

  /** Load the model */
  py::str model_path = sim.param().qmmm.nn.model_path;
  ml_model = new py::object(py_modules["torch"].attr("load")(model_path,"map_location"_a=py_modules["torch"].attr("device")("cpu")));
  
  //py::object conn = ase.attr("db").attr("connect")(db_path);
  // "map_location" parameter is to convert model trained on gpu to run on cpu, but seems not to work for my model

  /** Initialize the ML calculator */
  ml_calculator = new py::object(py_modules["schnetpack"].attr("interfaces").attr("SpkCalculator")(*ml_model, "energy"_a="energy", "forces"_a="forces"));


  return 0;
}

interaction::NN_Worker::~NN_Worker() {
  delete ml_model;
  delete ml_calculator;
}

int interaction::NN_Worker::run_QM(topology::Topology& topo
                     , configuration::Configuration& conf
                     , simulation::Simulation& sim, interaction::QM_Zone & qm_zone) {
  // run NN interface
  // create the molecule object
  py::object molecule = py_modules["ase"].attr("Atoms")();
  std::set<interaction::QM_Atom>::const_iterator it, to;
  it = qm_zone.qm.begin();
  to = qm_zone.qm.end();
  for (;it != to; ++it) {
    //unsigned atomic_number = it->atomic_number;
    uint32_t atomic_number = it->atomic_number;
    py::list py_coordinates;
    py_coordinates.attr("append")(it->pos[0]);
    py_coordinates.attr("append")(it->pos[1]);
    py_coordinates.attr("append")(it->pos[2]);
    py::object atom = py_modules["ase"].attr("Atom")("symbol"_a = atomic_number, "position"_a = py_coordinates);
    molecule.attr("append")(atom);
  }

  // Run the calculator
  molecule.attr("set_calculator")(ml_calculator);
  
  // Write the energy
  qm_zone.QM_energy() = molecule.attr("get_potential_energy")().cast<double>();

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
  }

  return 0;
}
