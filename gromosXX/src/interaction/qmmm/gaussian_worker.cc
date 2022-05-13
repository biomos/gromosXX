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

int interaction::Gaussian_Worker::init(const topology::Topology& topo
                                     , const configuration::Configuration& conf
                                     , simulation::Simulation& sim
                                     , const interaction::QM_Zone& qm_zone) {
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.gaussian);
  QM_Worker::param = this->param;
  // Aliases to shorten the code
  std::string& inp = this->param->input_file;
  std::string& out = this->param->output_file;

  if (inp.empty()) {
    if(util::create_tmpfile(inp) < 1) {
      io::messages.add("Unable to create temporary input file: " + inp,
        "Gaussian_Worker", io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(inp);
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
      this->tmp_files.insert(out);
      io::messages.add("Using temporary output file: " + out,
        "Gaussian_Worker", io::message::notice);
    }
  }
  
#ifndef HAVE_UNLINK
  {
    io::messages.add("Unlink function not supported on this platform. "
    + "Please delete temporary files manually.",
    this->name(), io::message::warning);
  }
#endif
  return 0;
}

int interaction::Gaussian_Worker::process_input(const topology::Topology& topo
                                            , const configuration::Configuration& conf
                                            , const simulation::Simulation& sim
                                            , const interaction::QM_Zone& qm_zone)
  {
  std::ofstream ifs;
  int err = 0;
  err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  std::string header(this->param->input_header);

  // Write header
  ifs << this->param->input_header;
  std::string guess;
  if (sim.steps() != 0) {
    guess = "guess=read";
  }
  std::string route_section = io::replace_string(this->param->route_section, "@@GUESS@@", guess);
  ifs << route_section;
  ifs << std::endl;
  ifs << "GROMOS generated input file" << std::endl;
  ifs << std::endl;
  std::string chsm = this->param->chsm;
  chsm = io::replace_string(chsm, "@@CHARGE@@", std::to_string(qm_zone.charge()));
  chsm = io::replace_string(chsm, "@@SPINM@@", std::to_string(qm_zone.spin_mult()));
  ifs << chsm;

  double len_to_qm = 1.0 / this->param->unit_factor_length;
  double cha_to_qm = 1.0 / this->param->unit_factor_charge;

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
    DEBUG(15,"Capping atom " << it->qm_index << "-" << it->mm_index << " " 
            << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
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
  ifs << std::string(3, '\n');
  ifs.close();
  return 0;
}

int interaction::Gaussian_Worker::run_calculation() {
  int err = util::system_call(this->param->binary + " < " + this->param->input_file
                                + " 1> " + this->param->output_file + " 2>&1 ");
  if (err) {
    std::ostringstream msg;
    msg << "Gaussian failed with code " << err;
    msg << ". See output file " << this->param->output_file << " for details.";
    io::messages.add(msg.str(), "Gaussian_Worker", io::message::error);
    return 1;
  }
  return 0;
}
int interaction::Gaussian_Worker::process_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  int err = 0;
  err = this->open_output(ofs, this->param->output_file);
  if (err) return err;

  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;
  
  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
      && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
  }

  if (minimisation) {
    err = this->parse_coordinates(ofs, qm_zone);
    if (err) return err;
  }

  err = this->parse_forces(sim, ofs, qm_zone);
  if (err) return err;

  ofs.close();
  return 0;
}

void interaction::Gaussian_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                        , const int atomic_number
                                        , const math::Vec& pos)
  {
  inputfile_stream << std::setw(4) << std::left << atomic_number
                   << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::endl;
}

void interaction::Gaussian_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                        , const math::Vec& pos
                                        , const double charge)
  {
  inputfile_stream << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::setw(25) << std::right << charge
                   << std::endl;
}

void interaction::Gaussian_Worker::write_mm_pos(std::ofstream& inputfile_stream
                                        , const math::Vec& pos)
  {
  /* Format 3F20.12 */
  inputfile_stream << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::endl;
}

int interaction::Gaussian_Worker::parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
  std::string& out = this->param->output_file;
  // Find the block
  {
  std::string line;
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
    // skip line
    std::getline(ofs, line);
  }
  // Parse charges of QM atoms
  {
    std::string dummy;
    for(std::set<QM_Atom>::iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
      ofs >> dummy >> dummy >> it->qm_charge;
      if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse charge of QM atom " << (it->index + 1)
            << " in " << out;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
      it->qm_charge *= this->param->unit_factor_charge;
    }
    // Do the capping atoms as well
    for(std::set<QM_Link>::iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
      ofs >> dummy >> dummy >> it->qm_charge;
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
  int dummy = 0;
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

int interaction::Gaussian_Worker::parse_forces(const simulation::Simulation& sim
                                            , std::ifstream& ofs
                                            , interaction::QM_Zone& qm_zone) {
  std::string& out = this->param->output_file;
  std::string line;
  int err = 0;
  
  // Find MM forces - in the form of electric field
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    bool got_mm_efield = false;
    while (std::getline(ofs, line)) {
      if (line.find("-------- Electric Field --------") != std::string::npos) {
        got_mm_efield = true;
        break;
      }
    }
    if (!got_mm_efield) {
      io::messages.add("Unable to find electric field at MM charges in output file " + out,
                          this->name(), io::message::error);
      return 1;
    }
    const unsigned skip_lines = 2 + qm_zone.qm.size() + qm_zone.link.size();
    for (unsigned i = 0; i < skip_lines; ++i)
      std::getline(ofs, line);
    // Parse MM atoms
    std::set<interaction::MM_Atom>::iterator it, to;
    for(std::set<interaction::MM_Atom>::iterator
        it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
      const int i = it->index;
      DEBUG(15, "Parsing electric field on MM atom " << i);
      err = this->parse_force(ofs, it->force);
      if (err) return err;
      // We still need to multiply by charge, since here we parsed electric field only
      if (it->is_polarisable) {
        it->force *= it->charge - it->cos_charge;
        err = this->parse_force(ofs, it->cos_force);
        if (err) return err;
        it->cos_force *= it->cos_charge;
      }
      else {
        DEBUG(15, "Charge of atom " << i << ": " << it->charge);
        it->force *= it->charge;
        DEBUG(15, "stored force on MM atom " << i << ": " << math::v2s(it->force));
      }
    }
  }
  {
    // Find QM forces
    bool got_qm_forces = false;
    while (std::getline(ofs, line)) {
      if (line.find("Forces (Hartrees/Bohr)") != std::string::npos) {
        got_qm_forces = true;
        break;
      }
    }
    if (!got_qm_forces) {
      io::messages.add("Unable to find QM forces in output file " + out,
                          this->name(), io::message::error);
      return 1;
    }

    const unsigned skip_lines = 2;
    for (unsigned i = 0; i < skip_lines; ++i)
      std::getline(ofs, line);
    // Parse QM atoms
    for(std::set<QM_Atom>::iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
      err = this->parse_force(ofs, it->force);
      if (err) return err;
    }
    // Parse link QM atoms
    for(std::set<QM_Link>::iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
      err = this->parse_force(ofs, it->force);
      if (err) return err;
    }
  }
  return 0;
}

int interaction::Gaussian_Worker::parse_force(std::ifstream& ofs,
                                              math::Vec& force) {
  std::string dummy;
  ofs >> dummy >> dummy >> force(0) >> force(1) >> force(2);
  if (ofs.fail()) {
    std::ostringstream msg;
    msg << "Failed to parse force line of atom "
        << " in " << this->param->output_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  force *= this->param->unit_factor_force;
  return 0;
}
