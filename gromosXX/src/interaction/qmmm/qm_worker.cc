/**
 * @file qm_worker.cc
 * implements the factory function for the QM_Worker class
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/interaction.h>

#include <simulation/parameter.h>

#include <util/timing.h>
#include <util/system_call.h>

#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "mndo_worker.h"
//#include "turbomole_worker.h"
//#include "dftb_worker.h"
//#include "mopac_worker.h"
#include "gaussian_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::QM_Worker::QM_Worker(std::string name) : m_name(name)/*, m_timer(nullptr)*/
                                                      , m_timer(name)
                                                      , param(nullptr)
                                                      , symlink_err(0)
                                                      , minimisation(false)
                                                      , using_tmp(false) {}
interaction::QM_Worker::~QM_Worker() {}

interaction::QM_Worker * interaction::QM_Worker::get_instance(const simulation::Simulation& sim) {

  switch (sim.param().qmmm.software) {
    case simulation::qm_mndo :
      return new MNDO_Worker;
    case simulation::qm_turbomole :
      //return new Turbomole_Worker; // Temporarily disabled
      return nullptr; // REMOVE
    case simulation::qm_dftb :
      //return new DFTB_Worker; // Temporarily disabled
      return nullptr; // REMOVE
    case simulation::qm_mopac :
      //return new MOPAC_Worker; // Temporarily disabled
      return nullptr; // REMOVE
    case simulation::qm_gaussian :
      return new Gaussian_Worker;
    default:
      io::messages.add("QM worker not implemented", "QM_Worker", io::message::critical);
      break;
  }
  return nullptr;
}

int interaction::QM_Worker::run_QM(topology::Topology& topo
                                 , configuration::Configuration& conf
                                 , simulation::Simulation& sim
                                 , interaction::QM_Zone& qm_zone) {
  if (qm_zone.mm.empty()) {
    io::messages.add("Cannot deal with zero MM atoms yet.", this->name()
                    , io::message::warning);
    //return 1;
  }
  m_timer.start();

  DEBUG(15,"Running QM Worker");
  int ret = 0;
  m_timer.start("writing input");
  if ((ret = this->write_input(topo, conf, sim, qm_zone)) != 0)
    return ret;
  m_timer.stop("writing input");
  
  m_timer.start("QM program call");
  if ((ret = this->system_call()) != 0)
    return ret;
  m_timer.stop("QM program call");

  m_timer.start("reading output");
  if ((ret = this->read_output(topo, conf, sim, qm_zone)) != 0)
    return ret;

  m_timer.stop("reading output");
  m_timer.stop();
  return 0;
}

int interaction::QM_Worker::write_input(const topology::Topology& topo
                                      , const configuration::Configuration& conf
                                      , const simulation::Simulation& sim
                                      , const interaction::QM_Zone& qm_zone) {
  std::ofstream ifs;
  if (int ret = this->open_input(ifs, this->param->input_file) != 0)
    return ret;
  /**
   * Custom input generating function goes here
   */
  ifs.close();
  return 0;
}

int interaction::QM_Worker::system_call() {
  return util::system_call(this->param->binary, this->param->input_file, this->param->output_file);
};

int interaction::QM_Worker::read_output(topology::Topology& topo
                                      , configuration::Configuration& conf
                                      , simulation::Simulation& sim
                                      , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  if (int ret = this->open_output(ofs, this->param->output_file) != 0)
    return ret;
  /**
   * Custom parsing function goes here
   */
  ofs.close();
  return 0;
}

int interaction::QM_Worker::open_input(std::ofstream& inputfile_stream, const std::string& input_file) {
  inputfile_stream.open(input_file.c_str());
  if (!inputfile_stream.is_open()) {
    io::messages.add("Unable to write to input file: "
            + input_file, this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::QM_Worker::open_output(std::ifstream& outputfile_stream, const std::string& output_file) {
  outputfile_stream.open(output_file.c_str());
  if (!outputfile_stream.is_open()) {
    io::messages.add("Unable to read from output file: "
            + output_file, this->name(), io::message::error);
    return 1;
  }
  return 0;
}