/**
 * @file turbomole_worker.cc
 * interface to the Turbomole software package
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

#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "turbomole_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::Turbomole_Worker::Turbomole_Worker() : QM_Worker("Turbomole Worker"), param(nullptr) {};

int interaction::Turbomole_Worker::init(const topology::Topology& topo
                                      , const configuration::Configuration& conf
                                      , simulation::Simulation& sim
                                      , const interaction::QM_Zone& qm_zone) {
  DEBUG(15, "Initializing " << this->name());
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.turbomole);
  QM_Worker::param = this->param;
  this->cwd = this->getcwd();
  if (this->cwd == "") return 1;
  return 0;
}

int interaction::Turbomole_Worker::process_input(const topology::Topology& topo
                                             , const configuration::Configuration& conf
                                             , const simulation::Simulation& sim
                                             , const interaction::QM_Zone& qm_zone) {
  if (this->chdir(this->param->working_directory) != 0) return 1;
  std::ofstream ifs;
  int err = this->open_input(ifs, this->param->input_coordinate_file);
  if (err) return err;
  
  ifs << "$coord" << std::endl;
  // write QM zone
  double len_to_qm = 1.0 / this->param->unit_factor_length;
  for (std::set<QM_Atom>::const_iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15, "Writing QM atom " << it->index << " " << it->atomic_number << " " << math::v2s(it->pos));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  // write capping atoms
  for (std::set<QM_Link>::const_iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15, "Writing capping atom " << it->qm_index << "-" << it->mm_index << " " << it->atomic_number << " " << math::v2s(it->pos));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  ifs << "$end" << std::endl;
  ifs.close();

  err = this->open_input(ifs, this->param->input_mm_coordinate_file);
  if (err) return err;

  ifs << "$point_charges" << std::endl;
  double cha_to_qm = 1.0 / this->param->unit_factor_charge;
  for (std::set<MM_Atom>::const_iterator
        it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    if (it->is_polarisable) {
      this->write_mm_atom(ifs, it->pos * len_to_qm, (it->charge - it->cos_charge) * cha_to_qm);
      this->write_mm_atom(ifs, (it->pos + it->cosV) * len_to_qm, it->cos_charge * cha_to_qm);
    }
    else {
      this->write_mm_atom(ifs, it->pos * len_to_qm, it->charge * cha_to_qm);
    }
  }
  ifs << "$end" << std::endl;
  ifs.close();
  return 0;
}

void interaction::Turbomole_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                                , const int atomic_number
                                                , const math::Vec& pos) const
  {
  inputfile_stream << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << "      " << std::left
                   << this->param->elements[atomic_number]
                   << std::endl;
}

void interaction::Turbomole_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                                , const math::Vec& pos
                                                , const double charge) const
  {
  // Turbomole skips zero charges in output
  if (charge != 0.0) {
    inputfile_stream << std::scientific << std::setprecision(17)
                     << std::setw(25) << std::right << pos(0)
                     << std::setw(25) << std::right << pos(1)
                     << std::setw(25) << std::right << pos(2)
                     << std::setw(25) << std::right << charge << std::endl;
  }
}

int interaction::Turbomole_Worker::run_calculation()
  {
  // First delete output files, since Turbomole appends to them and we dont need old data
#ifdef HAVE_UNLINK
  unlink(param->output_energy_file.c_str());
  unlink(param->output_gradient_file.c_str());
  unlink(param->output_charge_file.c_str());
#else
  io::messages.add("unlink function is not available on this platform.", 
        this->name(), io::message::error);
  return 1;
#endif

  DEBUG(15, "Calling external Turbomole program");
  for(std::vector<std::string>::const_iterator it = this->param->toolchain.begin(),
          to = this->param->toolchain.end(); it != to; ++it) {
  
    DEBUG(15, "Current toolchain command: " << *it);

    std::string output_file(*it + ".out"), input_file("");
    
    if (*it == "define") {
      input_file = "define.inp";
    }
    if (!input_file.empty()) input_file = " < " + input_file;
    std::string comm = this->param->binary_directory + "/" + *it + input_file + " 1> " + output_file + " 2>&1";
    int err = util::system_call(comm);
    if (err != 0) {
      std::ostringstream msg;
      msg << "Turbomole program " << *it << " failed. ";
      if (err == 127)
        msg << *it << " probably not in PATH. ";
      msg << "See " << output_file << " for details.";
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    // open the output file and check whether the program ended normally
    std::ifstream ofs;
    err = this->open_output(ofs, output_file);
    if (err) return err;
    bool success = false;
    std::string line;
    while(std::getline(ofs, line)) {
      if (line.find(" ended normally") != std::string::npos ||
          line.find(" : all done") != std::string::npos) {
        success = true;
        break;
      }
    }
    if (!success) {
      std::ostringstream msg;
      msg << "Turbomole program " << *it << " ended abnormally. See " 
              << output_file << " for details.";
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
  }
  return 0;
}

int interaction::Turbomole_Worker::process_output(topology::Topology& topo
                                             , configuration::Configuration& conf
                                             , simulation::Simulation& sim
                                             , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  int err = this->open_output(ofs, this->param->output_energy_file);
  if (err) return err;

  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;
  ofs.close();
  ofs.clear();

  err = this->open_output(ofs, this->param->output_gradient_file);
  if (err) return err;

  err = this->parse_qm_gradients(ofs, qm_zone);
  if (err) return err;
  ofs.close();
  ofs.clear();
  
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    err = this->open_output(ofs, this->param->output_mm_gradient_file);
    if (err) return err;

    err = this->parse_mm_gradients(ofs, qm_zone);
    if (err) return err;
    ofs.close();
    ofs.clear();
  } else if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
              && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    err = this->open_output(ofs, this->param->output_charge_file);
    if (err) return err;

    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
    ofs.close();
    ofs.clear();
  }
  if (this->chdir(this->cwd) != 0) return 1;
  return 0;
}

int interaction::Turbomole_Worker::parse_charges(std::ifstream& ofs
                                               , interaction::QM_Zone& qm_zone) const {
  std::string& out = this->param->output_charge_file;
  // Find the block
  {
    std::string line;
    bool got_charges = false;
    while(std::getline(ofs, line)) {
      if (line.find("charges resulting from fit:") != std::string::npos) {
        got_charges = true;
        break;
      }
    }
    if (!got_charges) {
      io::messages.add("Unable to find QM charges in output file " + out,
                        this->name(), io::message::error);
      return 1;
    }
    const int n_lines = 3;
    // skip n lines
    for (int i = 0; i < n_lines; ++i) {
      std::getline(ofs, line);
    }
  }
  // Parse QM charges
  {
    std::string dummy;
    // QM atoms
    for(std::set<QM_Atom>::iterator
          it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
      DEBUG(15,"Parsing charge of QM atom " << it->index);
      ofs >> dummy >> dummy >> dummy >> it->qm_charge;
      DEBUG(15,"Charge: " << it->qm_charge);
      if (ofs.fail()) {
        std::ostringstream msg;
        msg << "Failed to parse charge of QM atom " << (it->index + 1)
            << " in " << out;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
      it->qm_charge *= this->param->unit_factor_charge;
    }
    // Capping atoms
    for(std::set<QM_Link>::iterator
          it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
      DEBUG(15,"Parsing charge of capping atom " << it->qm_index << "-" << it->mm_index);
      ofs >> dummy >> dummy >> dummy >> it->qm_charge;
      DEBUG(15,"Charge: " << it->qm_charge);
      if (ofs.fail()) {
        std::ostringstream msg;
        msg << "Failed to parse charge of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << out;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
      it->qm_charge *= this->param->unit_factor_charge;
    }

  }
  return 0;
}

int interaction::Turbomole_Worker::parse_energy(std::ifstream& ofs
                                              , interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_energy = false;
  while(std::getline(ofs, line)) {
    // get energy section
    if (line.find("$energy") != std::string::npos) {
      std::getline(ofs, line);
      std::istringstream iss(line);
      int dummy = 0;
      iss >> dummy >> qm_zone.QM_energy();
      if (iss.fail()) {
        std::ostringstream msg;
        msg << "Failed to parse energy in " + this->param->output_energy_file;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
      qm_zone.QM_energy() *= this->param->unit_factor_energy;  
      got_energy = true;
      break;
    }
    }
  if (!got_energy) {
    io::messages.add("Unable to find energy in output file "
                      + this->param->output_energy_file
                      , this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::Turbomole_Worker::parse_qm_gradients(std::ifstream& ofs
                                                    , interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_qm_gradients = false;
  while(std::getline(ofs, line)) {
    if (line.find("$grad") != std::string::npos) {
      got_qm_gradients = true;
      break;
    }
  }
  if (!got_qm_gradients) {
    io::messages.add("Unable to find QM gradients in output file "
                      + this->param->output_gradient_file,
                        this->name(), io::message::error);
    return 1;
  }

  // skip n+1 lines
  const unsigned qm_size = qm_zone.qm.size() + qm_zone.link.size();
  for(unsigned int i = 0; i < qm_size+1; ++i) {
    std::getline(ofs, line);
  }
  // Parse QM atoms
  for(std::set<QM_Atom>::iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15,"Parsing gradient of QM atom " << it->index);
    int err = this->parse_gradient(ofs, it->force);
    DEBUG(15,"Force: " << math::v2s(it->force));
    if (err) {
      std::ostringstream msg;
      msg << "Failed to parse gradient line of QM atom " << (it->index + 1)
          << " in " << this->param->output_gradient_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
  }
  // Parse capping atoms
  for(std::set<QM_Link>::iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15,"Parsing gradient of capping atom " << it->qm_index << "-" << it->mm_index);
    int err = this->parse_gradient(ofs, it->force);
    DEBUG(15,"Force: " << math::v2s(it->force));
    if (err) {
      std::ostringstream msg;
      msg << "Failed to parse gradient line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << this->param->output_gradient_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
  }
  return 0;
}

int interaction::Turbomole_Worker::parse_mm_gradients(std::ifstream& ofs
                                                    , interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_mm_gradients = false;
  while(std::getline(ofs, line)) {
    if (line.find("$point_charge_gradients") != std::string::npos) {
      got_mm_gradients = true;
      break;
    }
  }
  if (!got_mm_gradients) {
    io::messages.add("Unable to find MM gradients in output file "
                      + this->param->output_mm_gradient_file,
                        this->name(), io::message::error);
    return 1;
  }
  // Parse MM atoms
  for(std::set<MM_Atom>::iterator
        it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    // Turbomole ignores zero charges, so we have to skip them as well
    if (it->charge != 0.0 || (it->is_polarisable && (it->charge - it->cos_charge) != 0.0)) {
      DEBUG(15,"Parsing gradient of MM atom " << it->index);
      int err = this->parse_gradient(ofs, it->force);
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse gradient line of MM atom " << (it->index + 1)
            << " in " << this->param->output_mm_gradient_file;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
    if (it->is_polarisable && it->cos_charge != 0.0) {
      DEBUG(15,"Parsing gradient of MM atom " << it->index);
      int err = this->parse_gradient(ofs, it->force);
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse gradient line of COS of MM atom " << (it->index + 1)
            << " in " << this->param->output_mm_gradient_file;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
  }
  return 0;
}

int interaction::Turbomole_Worker::parse_gradient(std::ifstream& ofs,
                                                  math::Vec& force) const {
  std::string line;
  if (!std::getline(ofs, line)) {
    io::messages.add("Failed to read gradient line"
                    , this->name(), io::message::error);
    return 1;
  }
  this->defortranize(line);
  std::istringstream iss(line);
  iss >> force(0) >> force(1) >> force(2);
  if (iss.fail()) return 1;
  // force = - gradient
  force *= - this->param->unit_factor_force;
  return 0;
}
