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
 * @file qm_worker.cc
 * implements the factory function for the QM_Worker class
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/interaction.h"

#include "../../../simulation/parameter.h"

#include "../../../util/timing.h"
#include "../../../util/system_call.h"

#include "qm_atom.h"
#include "qm_link.h"
#include "mm_atom.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "ghost_worker.h"
#include "mndo_worker.h"
#include "turbomole_worker.h"
#include "dftb_worker.h"
#include "mopac_worker.h"
#include "gaussian_worker.h"
#include "orca_worker.h"
#include "nn_worker.h"

#ifdef XTB
  #include "xtb_worker.h"
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::QM_Worker::QM_Worker(std::string name) : m_timer(name)
                                                    , m_name(name)
                                                    , param(nullptr)
                                                    , minimisation(false) {}

interaction::QM_Worker::~QM_Worker() {
  // Remove temporary files and links
#ifdef HAVE_UNLINK
  while (!this->tmp_files.empty()) {
    std::set<std::string>::const_iterator it = this->tmp_files.begin();
    unlink(it->c_str());
    this->tmp_files.erase(*it);
  }
#endif

#ifdef HAVE_RMDIR
  // Remove temporary directory
  if (!this->tmp_dir.empty()) {
    rmdir(this->tmp_dir.c_str());
  }
#endif
}

interaction::QM_Worker * interaction::QM_Worker::get_instance(const simulation::Simulation& sim) {

  switch (sim.param().qmmm.software) {
    case simulation::qm_ghost :
      return new Ghost_Worker;
    case simulation::qm_mndo :
      return new MNDO_Worker;
    case simulation::qm_turbomole :
      return new Turbomole_Worker;
    case simulation::qm_dftb :
      return new DFTB_Worker;
    case simulation::qm_mopac :
      return new MOPAC_Worker;
    case simulation::qm_gaussian :
      return new Gaussian_Worker;
    case simulation::qm_nn :
      return new NN_Worker;
    case simulation::qm_orca :
      return new Orca_Worker;
#ifdef XTB
    case simulation::qm_xtb :
      return new XTB_Worker;
#endif
    default:
      io::messages.add("QM worker not implemented", "QM_Worker", io::message::critical);
      break;
  }
  return nullptr;
}

int interaction::QM_Worker::init(const topology::Topology& topo
                   , const configuration::Configuration& conf
                   , simulation::Simulation& sim
                   , const interaction::QM_Zone& qm_zone) {
  // open trajectory streams
  if (sim.param().qmmm.write > 0) {
    // coordinates
    int err = open_input(input_coordinate_stream, this->param->trajectory_input_coordinate_file);
    if (err) return err;
    // gradients
    err = open_input(output_gradient_stream, this->param->trajectory_output_gradient_file);
    if (err) return err;

    if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) { // electrostatic embedding
      // point charges
      err = open_input(input_point_charge_stream, this->param->trajectory_input_pointcharges_file);
      if (err) return err;
      // point charge gradients
      err = open_input(output_point_charge_gradient_stream, this->param->trajectory_output_mm_gradient_file);
      if (err) return err;
    }
    else { // mechanical embedding
      // charges
      err = open_input(output_charges_stream, this->param->trajectory_output_charges_file);
      if (err) return err;
    }
  }
  return 0;
}

int interaction::QM_Worker::run_QM(topology::Topology& topo
                                 , configuration::Configuration& conf
                                 , simulation::Simulation& sim
                                 , interaction::QM_Zone& qm_zone) {
  m_timer.start(sim);

  DEBUG(15,"Running QM Worker");
  int ret = 0;
  m_timer.start_subtimer("writing input");
  if ((ret = this->process_input(topo, conf, sim, qm_zone)) != 0)
    return ret;
  m_timer.stop_subtimer("writing input");

  if ((sim.param().qmmm.write > 0) &&
	  ((sim.steps()) % (sim.param().qmmm.write) == 0)) {
    m_timer.start_subtimer("writing QM trajectory input");
    // steps reported in output are steps finished already
    this->save_input(sim.steps(), sim, qm_zone);
    m_timer.stop_subtimer("writing QM trajectory input");
  }
  
  m_timer.start_subtimer("QM program call");
  if ((ret = this->run_calculation()) != 0)
    return ret;
  m_timer.stop_subtimer("QM program call");

  m_timer.start_subtimer("reading output");
  if ((ret = this->process_output(topo, conf, sim, qm_zone)) != 0)
    return ret;
  m_timer.stop_subtimer("reading output");

  if ((sim.param().qmmm.write > 0) &&
	  ((sim.steps()) % (sim.param().qmmm.write) == 0)) {
    m_timer.start_subtimer("writing QM trajectory output");
    // steps reported in output are steps finished already
    this->save_output(sim.steps(), sim, qm_zone);
    m_timer.stop_subtimer("writing QM trajectory output");
  }

  m_timer.stop();
  return 0;
}

int interaction::QM_Worker::process_input(const topology::Topology& topo
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

int interaction::QM_Worker::run_calculation() {
  return util::system_call(this->param->binary + " < " + this->param->input_file
                                + " 1> " + this->param->output_file + " 2>&1 ");
};

void interaction::QM_Worker::save_input(const unsigned int step
                                      , const simulation::Simulation& sim
                                      , const interaction::QM_Zone & qm_zone) const {
  this->save_input_coord(input_coordinate_stream, step, qm_zone);
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) { // electrostatic embedding
    this->save_input_point_charges(input_point_charge_stream, step, this->get_num_charges(sim, qm_zone), qm_zone);
  }
}

void interaction::QM_Worker::save_output(const unsigned int step
                                       , const simulation::Simulation& sim
                                       , const interaction::QM_Zone & qm_zone) const {
  this->save_output_gradients(output_gradient_stream, step, qm_zone);
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) { // electrostatic embedding
    this->save_output_pc_gradients(output_point_charge_gradient_stream, step, qm_zone);
  }
  else { // mechanical embedding
    this->save_output_charges(output_charges_stream, step, qm_zone);
  }
}

void interaction::QM_Worker::save_input_coord(std::ofstream& ifs
                                            , const unsigned int step
                                            , const interaction::QM_Zone & qm_zone) const {
  // Gromos -> QM length unit is inverse of input value from QM/MM specification file
  const double len_to_qm = 1.0 / this->param->unit_factor_length;

  // write step size
  this->write_step_size(ifs, step);

  // write QM coordinates
  this->write_coordinate_header(ifs, qm_zone);
  DEBUG(15, "Writing QM coordinates");
  for (std::set<QM_Atom>::const_iterator 
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15, it->index << " " << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  // write capping atoms 
  DEBUG(15, "Writing capping atoms coordinates");
  for (std::set<QM_Link>::const_iterator it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; it++) {
    DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << " "
      << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  } 
  this->write_coordinate_footer(ifs);
}

void interaction::QM_Worker::save_input_point_charges(std::ofstream& ifs
                                                    , const unsigned int step
                                                    , const unsigned int ncharges
                                                    , const interaction::QM_Zone & qm_zone) const {
  // MM -> QM length unit is inverse of input value
  const double cha_to_qm = 1.0 / this->param->unit_factor_charge;
  const double len_to_qm = 1.0 / this->param->unit_factor_length;

  // write step size
  this->write_step_size(ifs, step);

  // write QM coordinates
  ifs << ncharges << '\n';
  for (std::set<MM_Atom>::const_iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    if (it->is_polarisable) {
      // MM atom minus COS
      DEBUG(15, it->index << " " << it->atomic_number << " " 
        << (it->charge - it->cos_charge) * cha_to_qm << " " << math::v2s(it->pos * len_to_qm));
      this->write_mm_atom(ifs, it->atomic_number, it->pos * len_to_qm, (it->charge - it->cos_charge) * cha_to_qm);
      // COS
      DEBUG(15, it->index << " " << it->atomic_number << " " 
        << it->cos_charge * cha_to_qm << " " << math::v2s((it->pos + it->cosV) * len_to_qm));
      this->write_mm_atom(ifs, it->atomic_number, it->cosV * len_to_qm, it->cos_charge * cha_to_qm);
    }
    else {
      this->write_mm_atom(ifs, it->atomic_number, it->pos * len_to_qm, it->charge * cha_to_qm);
    }
  }
}

void interaction::QM_Worker::save_output_gradients(std::ofstream& ifs
                                                 , const unsigned int step
                                                 , const interaction::QM_Zone & qm_zone) const {
  // MM -> QM length unit is inverse of input value
  const double energy_to_qm = 1.0 / this->param->unit_factor_energy;
  const double force_to_qm = 1.0 / this->param->unit_factor_force;

  // write step size
  this->write_step_size(ifs, step);

  // Write energy
  ifs.setf(std::ios::fixed, std::ios::floatfield);
  ifs << std::setprecision(12);
  ifs << "ENERGY: " << qm_zone.QM_energy() * energy_to_qm << '\n';

  // Write QM atoms
  for (std::set<QM_Atom>::iterator
         it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    // forces = negative gradient (!)
    DEBUG(15, "Writing gradients of QM atom " << it->index);
    this->write_gradient(-1.0 * it->force * force_to_qm, ifs);
    DEBUG(15, "Force: " << math::v2s(it->force));
  }
  // Write capping atoms (index i keeps running...)
  for (std::set<QM_Link>::iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15, "Writing gradient of capping atom " << it->qm_index << "-" << it->mm_index);
    this->write_gradient(-1.0 * it->force * force_to_qm, ifs);
    DEBUG(15, "Force: " << math::v2s(it->force));
  }
}

void interaction::QM_Worker::save_output_pc_gradients(std::ofstream& ifs
                                                    , const unsigned int step
                                                    , const interaction::QM_Zone & qm_zone) const {
  // MM -> QM length unit is inverse of input value
  const double force_to_qm = 1.0 / this->param->unit_factor_force;

  // write step size
  this->write_step_size(ifs, step);

  // Parse MM atoms
  for (std::set<MM_Atom>::iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    // forces = negative gradient (!)
    DEBUG(15,"Writing gradient of MM atom " << it->index);
    this->write_gradient(-1.0 * it->force * force_to_qm, ifs);
    DEBUG(15, "Force: " << math::v2s(it->force));
    if (it->is_polarisable) {
      DEBUG(15, "Writing gradient of COS of MM atom " << it->index);
      this->write_gradient(-1.0 * it->cos_force * force_to_qm, ifs);
      DEBUG(15, "Force " << math::v2s(it->cos_force));
    }
  }
}

void interaction::QM_Worker::save_output_charges(std::ofstream& ifs
                                               , const unsigned int step
                                               , const interaction::QM_Zone & qm_zone) const {
  // MM -> QM length unit is inverse of input value
  const double cha_to_qm = 1.0 / this->param->unit_factor_charge;

  // write step size
  this->write_step_size(ifs, step);

  // QM atoms
  unsigned int i = 0;
  for (std::set<QM_Atom>::const_iterator
             it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; it++) {
    DEBUG(15, "Writing charge of QM atom " << it->index);
    this->write_charge(ifs, it->atomic_number, it->qm_charge * cha_to_qm); 
    DEBUG(15, "Charge: " << it->qm_charge);
    ++i;
  }
  // Capping atoms (iterator i keeps running)
  for (std::set<QM_Link>::const_iterator
             it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; it++) {
    DEBUG(15, "Writing charge of capping atom " << it->qm_index << "-" << it->mm_index);
    this->write_charge(ifs, it->atomic_number, it->qm_charge * cha_to_qm);
    DEBUG(15, "Charge: " << it->qm_charge);
    ++i;
  }
}

int interaction::QM_Worker::process_output(topology::Topology& topo
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

int interaction::QM_Worker::open_input(std::ofstream& inputfile_stream, const std::string& input_file) const {
  inputfile_stream.open(input_file.c_str()); 
  if (!inputfile_stream.is_open()) {
    io::messages.add("Unable to write to file: "
            + input_file, this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::QM_Worker::open_output(std::ifstream& outputfile_stream, const std::string& output_file) const {
  outputfile_stream.open(output_file.c_str());
  if (!outputfile_stream.is_open()) {
    io::messages.add("Unable to read from file: "
            + output_file, this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::QM_Worker::get_num_charges(const simulation::Simulation& sim
                                          , const interaction::QM_Zone& qm_zone) const {
  unsigned num_charges = 0;
  switch (sim.param().qmmm.qmmm) {
    case simulation::qmmm_mechanical: {
      num_charges = 0;
      break;
    }
    case simulation::qmmm_electrostatic: {
      num_charges = qm_zone.mm.size();
      break;
    }
    case simulation::qmmm_polarisable: {
      num_charges = qm_zone.mm.size();
      for (std::set<MM_Atom>::const_iterator
          it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
        num_charges += int(it->is_polarisable);
      }
      break;
    }
    default: {
      io::messages.add("Unknown QMMM option", this->name(), io::message::error);
    }
  }
  return num_charges;
}

void interaction::QM_Worker::write_step_size(std::ofstream& ifs, 
                                             const unsigned int step) const {
  ifs << "TIMESTEP" << '\n';
  ifs << "    " << step << '\n';
  ifs << "END" << '\n';
}

void interaction::QM_Worker::write_coordinate_header(std::ofstream& inputfile_stream
                                                   , const QM_Zone& qm_zone) const {
  // default: pass  
}

void interaction::QM_Worker::write_coordinate_footer(std::ofstream& inputfile_stream) const {
  // default: pass
}

void interaction::QM_Worker::write_gradient(const math::Vec& gradient, 
                                            std::ofstream& inputfile_stream) const {
  inputfile_stream.setf(std::ios::fixed, std::ios::floatfield);
  inputfile_stream << std::setprecision(12)
                   << std::setw(20) << gradient(0)
                   << std::setw(20) << gradient(1)
                   << std::setw(20) << gradient(2)
                   << '\n';
}

void interaction::QM_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                         , const int atomic_number
                                         , const math::Vec& pos) const {
  inputfile_stream.setf(std::ios::fixed, std::ios::floatfield);
  inputfile_stream << std::setprecision(20)
                   << std::setw(28) << pos(0)
                   << std::setw(28) << pos(1)
                   << std::setw(28) << pos(2)
                   << std::setw(8)  << atomic_number
                   << '\n';
}

void interaction::QM_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                         , const int atomic_number
                                         , const math::Vec& pos
                                         , const double charge) const {
  if (charge != 0.0) {
    inputfile_stream.setf(std::ios::fixed, std::ios::floatfield);
    inputfile_stream << std::setprecision(6)
                     << std::setw(10) << charge
                     << std::setprecision(20)
                     << std::setw(28) << pos(0)
                     << std::setw(28) << pos(1)
                     << std::setw(28) << pos(2) 
                     << std::setw(8)  << atomic_number
                     << '\n';
  }
}

void interaction::QM_Worker::write_charge(std::ofstream& inputfile_stream
                                         , const int atomic_number
                                         , const double charge) const {
  inputfile_stream.setf(std::ios::fixed, std::ios::floatfield);
  inputfile_stream << std::setprecision(6)
                   << std::setw(10) << charge
                   << std::setw(8)  << atomic_number
                   << '\n';
}
