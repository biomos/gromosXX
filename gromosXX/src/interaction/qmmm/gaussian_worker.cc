/**
 * @file gaussian_worker.cc
 * interface to the Gaussian software package
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

// special interactions
#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "gaussian_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::Gaussian_Worker::Gaussian_Worker() : QM_Worker("Gaussian Worker"), param(nullptr) {};

int interaction::Gaussian_Worker::init(const simulation::Simulation& sim) {
  /** Make a copy of parameters, so QM worker can safely customize
   * them and also noone can change them
   */
  this->param = new simulation::Parameter::qmmm_struct::gaussian_param_struct(sim.param().qmmm.gaussian);
  // Aliases to shorten the code
  std::string& inp = this->param->input_file;
  std::string& out = this->param->output_file;

  if (inp.empty()) {
    this->using_tmp = true;
    if(util::create_tmpfile(inp) < 1) {
        io::messages.add("Unable to create temporary input file: " + inp,
        "Gaussian_Worker", io::message::critical);
      return 1;
    }
    else {
        io::messages.add("Using temporary input file: " + inp,
        "Gaussian_Worker", io::message::notice);
    }
  }

  if (out.empty()) {
    if (util::create_tmpfile(out) < 1) {
      io::messages.add("Unable to create temporary output file: " + out,
        "Gaussian_Worker", io::message::critical);
      return 1;
    }
    else {
        io::messages.add("Using temporary output file: " + out,
        "Gaussian_Worker", io::message::notice);
    }
  }
  return 0;
}

int interaction::Gaussian_Worker::write_input(const topology::Topology& topo
                                            , const configuration::Configuration& conf
                                            , const simulation::Simulation& sim
                                            , const interaction::QM_Zone& qm_zone)
  {
  std::ofstream ifs;
  int err;
  err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  std::string header(this->param->input_header);

  // Write header
  ifs << this->param->input_header;
  ifs << this->param->route_section;
  ifs << std::endl;
  ifs << "GROMOS generated input file" << std::endl;
  ifs << std::endl;
  ifs << this->param->chm;

  double len_to_qm = 1.0 / this->param->unit_factor_length;
  double cha_to_qm = 1.0 / this->param->unit_factor_charge;
  ifs.setf(std::ios::fixed, std::ios::floatfield);
  ifs.precision(8);

  DEBUG(15,"Writing QM coordinates");
  // Write QM coordinates
  for (std::set<QM_Atom>::const_iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15,it->index << " " << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  DEBUG(15,"Writing capping atoms coordinates");
  // Write capping atoms
  for (std::set<QM_Link>::const_iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15,it->index << " " << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  ifs << std::endl;
  
  // Write MM coordinates and charges
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
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
    ifs << std::endl;
    for (std::set<MM_Atom>::const_iterator
          it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
      if (it->is_polarisable) {
        this->write_mm_pos(ifs, it->pos * len_to_qm);
        this->write_mm_pos(ifs, (it->pos + it->cosV) * len_to_qm);
      }
      else {
        this->write_mm_pos(ifs, it->pos * len_to_qm);
      }
    }
  }
  ifs << 3 * '\n';
  ifs.close();
  return 0;
}

int interaction::Gaussian_Worker::system_call() {
  int err = util::system_call(this->param->binary, this->param->input_file, this->param->output_file);
  if (err) {
    std::ostringstream msg;
    msg << "Gaussian failed with code " << err;
    msg << ". See output file " << this->param->output_file << " for details.";
    io::messages.add(msg.str(), "Gaussian_Worker", io::message::error);
    return 1;
  }
  return 0;
}
int interaction::Gaussian_Worker::read_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  int err;
  err = this->open_output(ofs, this->param->output_file);
  if (err) return err;

  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;
  
  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical) {
    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
  }

  if (minimisation) {
    err = this->parse_coordinates(ofs, qm_zone);
    if (err) return err;
  }

  err = this->parse_gradients(sim, ofs, qm_zone);
  if (err) return err;

  ofs.close();
  return 0;
}

void interaction::Gaussian_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                        , const int atomic_number
                                        , const math::Vec& pos)
  {
  inputfile_stream << std::setw(4) << std::left << atomic_number
                   << std::setw(20) << std::setprecision(15) << std::right << pos(0)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(1)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(2)
                   << std::endl;
}

void interaction::Gaussian_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                        , const math::Vec& pos
                                        , const double charge)
  {
  inputfile_stream << std::setw(20) << std::setprecision(15) << std::right << pos(0)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(1)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(2)
                   << std::setw(20) << std::setprecision(4) << std::right << charge
                   << std::endl;
}

void interaction::Gaussian_Worker::write_mm_pos(std::ofstream& inputfile_stream
                                        , const math::Vec& pos)
  {
  inputfile_stream << std::setw(20) << std::setprecision(15) << std::right << pos(0)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(1)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(2)
                   << std::endl;
}

int interaction::Gaussian_Worker::parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
  std::string& out = this->param->output_file;
  std::string line;
  /**
   * Parse charges
   * They are used in mechanical embedding
   */
  bool got_charges = false;
  while (std::getline(ofs, line)) {
    if (line.find("ESP charges:") != std::string::npos) {
      got_charges = true;
      break;
    }
  }
  if (!got_charges) {
    io::messages.add("Charges were not found in the output file "
                      + out, this->name(), io::message::error);
    return 1;
  }
  const unsigned skip_lines = 1;
  for (unsigned i = 0; i < skip_lines; ++i)
    std::getline(ofs, line);
  for(std::set<QM_Atom>::iterator
      it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    if(!std::getline(ofs, line)) {
      std::ostringstream msg;
      msg << "Failed to read charge line of atom " << (it->index + 1)
          << " in " << out;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    std::istringstream iss(line);
    std::string dummy;
    iss >> dummy >> dummy >> it->qm_charge;
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse charge line of atom " << (it->index + 1)
          << " in " << out;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->qm_charge *= this->param->unit_factor_charge;
  }
  return 0;
}

int interaction::Gaussian_Worker::parse_coordinates(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
  std::string& out = this->param->output_file;
  std::string line;
  bool got_coordinates = false;
  // Find coordinates block
  // Not implemented yet
  if (!got_coordinates) {
    io::messages.add("Reading coordinates from Gaussian output is not supported.",
                      this->name(), io::message::error);
    return 1;
  }
  
  // Parse coordinates lines
  math::Vec pos;
  int dummy;
  for(std::set<QM_Atom>::iterator
      it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    if(!std::getline(ofs, line)) {
      std::ostringstream msg;
      msg << "Failed to read coordinates line of atom " << (it->index + 1)
          << " in " << out;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    std::istringstream iss(line);
    iss >> dummy >> dummy >> dummy >> pos(0) >> pos(1) >> pos(2);
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse coordinate line of atom " << (it->index + 1)
          << " in " + out;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->pos = pos * this->param->unit_factor_length;
  }
  return 0;
}

int interaction::Gaussian_Worker::parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
  std::string line;
  // Find energy block
  bool got_energy = false;
  while (std::getline(ofs, line)) {
    if (line.find("SCF Done") != std::string::npos) {
      got_energy = true;
      break;
    }
  }
  if (!got_energy) {
    io::messages.add("Unable to find energy in output file "
                      + this->param->output_file
                      , this->name(), io::message::error);
    return 1;
  }
  // Parse energy
  std::istringstream iss(line);
  std::string dummy;
  iss >> dummy >> dummy >> dummy >> dummy >> qm_zone.QM_energy();
  if (iss.fail()) {
    io::messages.add("Failed to parse QM energy in output file"
                      + this->param->output_file
                      , this->name(), io::message::error);
    return 1;
  }
  qm_zone.QM_energy() *= this->param->unit_factor_energy;
  return 0;
}

int interaction::Gaussian_Worker::parse_gradients(const simulation::Simulation& sim
                                            , std::ifstream& ofs
                                            , interaction::QM_Zone& qm_zone) {
  std::string& out = this->param->output_file;
  std::string line;
  int err;
  
  // Find MM gradients - Electric field comes first in the output
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    bool got_mm_gradients = false;
    while (std::getline(ofs, line)) {
      if (line.find("-------- Electric Field --------") != std::string::npos) {
        got_mm_gradients = true;
        break;
      }
    }
    if (!got_mm_gradients) {
      io::messages.add("Unable to find MM gradients in output file " + out,
                          this->name(), io::message::error);
      return 1;
    }
    const unsigned skip_lines = 2 + qm_zone.qm.size() + qm_zone.link.size();
    for (unsigned i = 0; i < skip_lines; ++i)
      std::getline(ofs, line);
    // Parse MM atoms
    err = this->_parse_gradients(ofs, qm_zone.mm);
    if (err) return err;
  }
  {
    // Find QM gradients
    bool got_qm_gradients = false;
    while (std::getline(ofs, line)) {
      if (line.find("Forces (Hartrees/Bohr)") != std::string::npos) {
        got_qm_gradients = true;
        break;
      }
    }
    if (!got_qm_gradients) {
      io::messages.add("Unable to find QM gradients in output file " + out,
                          this->name(), io::message::error);
      return 1;
    }

    const unsigned skip_lines = 2;
    for (unsigned i = 0; i < skip_lines; ++i)
      std::getline(ofs, line);
    // Parse QM atoms
    err = this->_parse_gradients(ofs, qm_zone.qm);
    if (err) return err;
    // Parse link QM atoms
    err = this->_parse_gradients(ofs, qm_zone.link);
    if (err) return err;
  }
  return 0;
}

template<class AtomType>
int interaction::Gaussian_Worker::_parse_gradients(std::ifstream& ofs,
                                               std::set<AtomType>& atom_set) {
  int err = 0;
  typename std::set<AtomType>::iterator it, to;
  for(it = atom_set.begin(), to = atom_set.end(); it != to; ++it) {
    DEBUG(15, "Parsing force of QM atom " << it->index);
    err = this->parse_gradient(ofs, it->index, it->force, this->param->unit_factor_force);
    if (err) return err;
    DEBUG(15, "stored force on QM atom " << it->index << ": " << math::v2s(it->force));
  }
  return 0;
}

template<>
int interaction::Gaussian_Worker::_parse_gradients
                                        (std::ifstream& ofs
                                       , std::set<interaction::MM_Atom>& atom_set) {
  int err = 0;
  std::set<interaction::MM_Atom>::iterator it, to;
  for(it = atom_set.begin(), to = atom_set.end(); it != to; ++it) {
    int i = it->index;
    DEBUG(15, "Parsing electric field on MM atom " << i);
    err = this->parse_gradient(ofs, i, it->force, this->param->unit_factor_force);
    if (err) return err;
    // We still need to multiply by charge, since here we parsed electric field only
    if (it->is_polarisable) {
      it->force *= it->charge - it->cos_charge;
      err = this->parse_gradient(ofs, i, it->cos_force, this->param->unit_factor_force);
      if (err) return err;
      it->cos_force *= it->cos_charge;
    }
    else {
      DEBUG(15, "Charge of atom " << i << ": " << it->charge);
      it->force *= it->charge;
      DEBUG(12, "stored force on MM atom " << i << ": " << math::v2s(it->force));
    }
  }
  return 0;
}

int interaction::Gaussian_Worker::parse_gradient(std::ifstream& ofs,
                                             const int index,
                                             math::Vec& force,
                                             const double unit_factor) {
  std::string line;
  math::Vec gradient;
  std::string dummy;
  if(!std::getline(ofs, line)) {
    std::ostringstream msg;
    msg << "Failed to read gradient line of atom " << (index)
        << " in " << this->param->output_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  std::istringstream iss(line);
  iss >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
  DEBUG(12, "Parsed gradient line of atom " << index << ": " << math::v2s(gradient));
  if (iss.fail()) {
    std::ostringstream msg;
    msg << "Failed to parse gradient line of atom " << (index)
        << " in " << this->param->output_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  force = gradient * unit_factor;
  return 0;
}

interaction::Gaussian_Worker::~Gaussian_Worker() {
#ifdef HAVE_UNLINK
  // Delete temporary files
  if (this->using_tmp) {
    unlink(this->param->input_file.c_str());
    unlink(this->param->output_file.c_str());
  }
  if (this->param != nullptr) {
    delete this->param;
    this->param = nullptr;
  }
#endif
}
