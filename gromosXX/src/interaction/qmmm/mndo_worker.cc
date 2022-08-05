/**
 * @file mndo_worker.cc
 * interface to the MNDO software package
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
#include "mndo_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::MNDO_Worker::MNDO_Worker() : QM_Worker("MNDO Worker"), param(nullptr) {};

int interaction::MNDO_Worker::init(const topology::Topology& topo
                                 , const configuration::Configuration& conf
                                 , simulation::Simulation& sim
                                 , const interaction::QM_Zone& qm_zone) {
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.mndo);
  QM_Worker::param = this->param;
  // Aliases to shorten the code
  std::string& inp = this->param->input_file;
  std::string& out = this->param->output_file;
  std::string& out_grad = this->param->output_gradient_file;
  std::string& dens = this->param->density_matrix_file;

  /** Check if optimization in QM is requested, then also read new positions
   * Maybe leave this on user input, since parsing of the header is cumbersome.
   * MNDO input is case- and space-insensitive.
   */

  if (inp.empty()) {
    if(util::create_tmpfile(inp) < 1) {
      io::messages.add("Unable to create temporary input file: " + inp,
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(inp);
      io::messages.add("Using temporary input file: " + inp,
        "MNDO_Worker", io::message::notice);
    }
  }

  if (out.empty()) {
    if (util::create_tmpfile(out) < 1) {
      io::messages.add("Unable to create temporary output file: " + out,
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(out);
      io::messages.add("Using temporary output file: " + out,
        "MNDO_Worker", io::message::notice);
    }
  }

  if (out_grad.empty()) {
    if (util::create_tmpfile(out_grad) < 1) {
      io::messages.add("Unable to create temporary output gradient file: " 
        + out_grad, "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(out_grad);
      io::messages.add("Using temporary output gradient file: " + out_grad,
        "MNDO_Worker", io::message::notice);
    }
  }
  else {
    // Try to create a file specified by user - fort.15 will be symlinked to it
    // If it already exist, we dont want to overwrite it, exit
    std::ofstream of(out_grad.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output gradient file: "
      + out_grad, "MNDO_Worker", io::message::critical);
      return 1;
    }
    of.close();
  }

  if (dens.empty()) {
    if (util::create_tmpfile(dens) < 1) {
      io::messages.add("Unable to create temporary density matrix file: " 
        + dens, "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(dens);
      io::messages.add("Using temporary density matrix file: " + dens,
        "MNDO_Worker", io::message::notice);
    }
  } // This file might exist and might be reused so we dont raise an error if it exists
  
#ifdef HAVE_SYMLINK
  // create fort.15 link for gradients
  if ("fort.15" != out_grad) {
    if (symlink(out_grad.c_str(), "fort.15") != 0) {
      io::messages.add("Unable to create symbolic link from fort.15 to "
        + out_grad + " - check permissions.",
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert("fort.15");
  }
  // create fort.11 link for density matrix
  if ("fort.11" != dens) {
    if (symlink(dens.c_str(), "fort.11") != 0) {
      io::messages.add("Unable to create symbolic link from fort.11 to "
        + dens + " - check permissions.",
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert("fort.11");
  }
#else
  {
    out_grad = "fort.15";
    dens = "fort.11";
    io::messages.add("Symbolic links are not supported in this built. "
      + "Output gradient file is now set to fort.15",
      "MNDO_Worker", io::message::warning);
    io::messages.add("Symbolic links are not supported in this built. "
      + "Density matrix file is now set to fort.11",
      "MNDO_Worker", io::message::warning);
  }
#endif
  
#ifndef HAVE_UNLINK
  {
    io::messages.add("Unlink function not supported on this platform. "
    + "Please delete temporary files manually.",
    this->name(), io::message::warning);
  }
#endif
  return 0;
}

int interaction::MNDO_Worker::process_input(const topology::Topology& topo
                                        , const configuration::Configuration& conf
                                        , const simulation::Simulation& sim
                                        , const interaction::QM_Zone& qm_zone)
  {
  std::ofstream ifs;
  int err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  std::string header(this->param->input_header);
  
  // Get number of links and replace
  unsigned num_links = qm_zone.link.size();
  // Include also linked QM atom - count them
  std::stringstream link_atoms;
  {
    unsigned i = 1;
    for (std::set<QM_Atom>::const_iterator
          it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it, ++i) {
      if (it->is_linked) {
        ++num_links;
        link_atoms << i << " ";
      }
    }
  }
  header = io::replace_string(header, "@@NUM_LINKS@@", std::to_string(num_links));
  
  // Get number of MM charges and replace
  unsigned num_charges = this->get_num_charges(sim, qm_zone);
  header = io::replace_string(header, "@@NUM_CHARGES@@", std::to_string(num_charges));
  header = io::replace_string(header, "@@CHARGE@@", std::to_string(qm_zone.charge()));
  header = io::replace_string(header, "@@SPINM@@", std::to_string(qm_zone.spin_mult()));

  // Write header
  ifs << header << std::endl;

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
    DEBUG(15,"capping atom " << it->qm_index << "-" << it->mm_index << " " 
            << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  // Write termination line
  this->write_qm_atom(ifs, 0, math::Vec(0.0));

  // Write link atom indices to be excluded in QM-MM electrostatics
  ifs << link_atoms.rdbuf();
  // Write capping atom indices - they are last in the geometry
  const unsigned qm_size = qm_zone.qm.size();
  const unsigned last_link = qm_zone.link.size() + qm_size;
  for (unsigned i = qm_size + 1; i <= last_link; ++i) {
    ifs << i << " ";
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
  }
  ifs.close();
  return 0;
}

int interaction::MNDO_Worker::run_calculation() {
  int err = util::system_call(this->param->binary + " < " + this->param->input_file
                                + " 1> " + this->param->output_file + " 2>&1 ");
  if (err) {
    std::ostringstream msg;
    msg << "MNDO failed with code " << err;
    if (err == 127)
      msg << ". mndo command probably not in PATH";
    msg << ". See output file " << this->param->output_file << " for details.";
    io::messages.add(msg.str(), "MNDO_Worker", io::message::error);
    return 1;
  }
  return 0;
}
int interaction::MNDO_Worker::process_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  int err = this->open_output(ofs, this->param->output_file);
  if (err) return err;
  
  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
      && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
  }

  ofs.close();
  ofs.clear();

  err = this->open_output(ofs, this->param->output_gradient_file);
  if (err) return err;

  if (minimisation) {
    err = this->parse_coordinates(ofs, qm_zone);
    if (err) return err;
  }

  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;

  err = this->parse_gradients(sim, ofs, qm_zone);
  if (err) return err;

  ofs.close();
  return 0;
}

void interaction::MNDO_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                        , const int atomic_number
                                        , const math::Vec& pos
                                        , const int opt_flag) const
  {
/*
      a(1,i)    11-20  f10.5  First coordinate.
      la        23-24    i2   Optimization of first coordinate.
                              = 0 a(1,i) is not optimized.
                              = 1 a(1,i) is optimized.
      a(2,i)    31-40  f10.5  Second coordinate.
      lb        43-44   i2    Optimization of second coordinate.
                              = 0 a(2,i) is not optimized.
                              = 1 a(2,i) is optimized.
      a(3,i)    51-60  f10.5  Third coordinate.
      lc        63-64   i2    Optimization of third coordinate.
                              = 0 a(3,i) is not optimized.
                              = 1 a(3,i) is optimized.
*/
  inputfile_stream << std::setw(2) << std::left << atomic_number
                   << std::string(2,' ')
                   << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::string(2,' ') << std::setw(2) << std::right << opt_flag
                   << std::setw(25) << std::right << pos(1)
                   << std::string(2,' ') << std::setw(2) << std::right << opt_flag
                   << std::setw(25) << std::right << pos(2)
                   << std::string(2,' ') << std::setw(2) << std::right << opt_flag
                   << std::endl;
}

void interaction::MNDO_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                        , const math::Vec& pos
                                        , const double charge) const
  {
    /*
      cm(1,m)    1-12  f12.7  x coordinate.
      cm(2,m)   13-24  f12.7  y coordinate.
      cm(3,m)   25-36  f12.7  z coordinate.
      qm(m)     37-44   f8.4  External point charge.
    */
  inputfile_stream << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::setw(25) << std::right << charge
                   << std::endl;
}

int interaction::MNDO_Worker::parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string& out = this->param->output_file;
  std::string line;
  /**
   * Parse charges
   * They are used in mechanical embedding
   */
  bool got_charges = false;
  while (std::getline(ofs, line)) {
    if (line.find("NET ATOMIC CHARGES") != std::string::npos) {
      got_charges = true;
      break;
    }
  }
  if (!got_charges) {
    io::messages.add("Charges were not found in the output file "
                      + out, this->name(), io::message::error);
    return 1;
  }
  // Skip three lines
  for (unsigned i = 0; i < 3; ++i)
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
  // Do the same for capping atoms
  for(std::set<QM_Link>::iterator
      it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    if(!std::getline(ofs, line)) {
      std::ostringstream msg;
      msg << "Failed to read charge line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << out;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    std::istringstream iss(line);
    std::string dummy;
    iss >> dummy >> dummy >> it->qm_charge;
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse charge line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << out;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->qm_charge *= this->param->unit_factor_charge;
  }
  return 0;
}

int interaction::MNDO_Worker::parse_coordinates(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string& out_grad = this->param->output_gradient_file;
  std::string line;
  bool got_coordinates = false;
  // Find coordinates block
  while (std::getline(ofs, line)) {
    if (line.find("CARTESIAN COORDINATES:") != std::string::npos) {
      got_coordinates = true;
      break;
    }
  }
  if (!got_coordinates) {
    io::messages.add("Unable to find coordinates in output file " + out_grad,
                      this->name(), io::message::error);
    return 1;
  }
  // Parse coordinates lines
  int dummy = 0;
  for(std::set<QM_Atom>::iterator
      it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    if(!std::getline(ofs, line)) {
      std::ostringstream msg;
      msg << "Failed to read coordinates line of atom " << (it->index + 1)
          << " in " << out_grad;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    std::istringstream iss(line);
    iss >> dummy >> dummy >> it->pos(0) >> it->pos(1) >> it->pos(2);
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse coordinate line of atom " << (it->index + 1)
          << " in " + out_grad;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->pos *= this->param->unit_factor_length;
  }
  // Do it also for capping atoms
  for(std::set<QM_Link>::iterator
      it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    if(!std::getline(ofs, line)) {
      std::ostringstream msg;
      msg << "Failed to read coordinates line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << out_grad;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    std::istringstream iss(line);
    iss >> dummy >> dummy >> it->pos(0) >> it->pos(1) >> it->pos(2);
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse coordinate line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << out_grad;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->pos *= this->param->unit_factor_length;
  }
  return 0;
}

int interaction::MNDO_Worker::parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string line;
  // Find energy block
  bool got_energy = false;
  while (std::getline(ofs, line)) {
    if (line.find("ENERGY") != std::string::npos) {
      got_energy = true;
      break;
    }
  }
  if (!got_energy) {
    io::messages.add("Unable to find energy in output file "
                      + this->param->output_gradient_file
                      , this->name(), io::message::error);
    return 1;
  }
  // Parse energy
  ofs >> qm_zone.QM_energy();
  if (ofs.fail()) {
    io::messages.add("Failed to parse QM energy in output file"
                      + this->param->output_gradient_file
                      , this->name(), io::message::error);
    return 1;
  }
  qm_zone.QM_energy() *= this->param->unit_factor_energy;
  return 0;
}

int interaction::MNDO_Worker::parse_gradients(const simulation::Simulation& sim
                                            , std::ifstream& ofs
                                            , interaction::QM_Zone& qm_zone) const {
  std::string& out_grad = this->param->output_gradient_file;
  std::string line;
  int err = 0;
  
  {
    // Find QM gradients
    bool got_qm_gradients = false;
    while (std::getline(ofs, line)) {
      if (line.find("CARTESIAN GRADIENT") != std::string::npos
          && line.find("OF MM ATOMS") == std::string::npos) {
        got_qm_gradients = true;
        break;
      }
    }
    if (!got_qm_gradients) {
      io::messages.add("Unable to find QM gradients in output file " + out_grad,
                          this->name(), io::message::error);
      return 1;
    }
    
    // Parse QM atoms
    for(std::set<QM_Atom>::iterator
          it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
      DEBUG(15,"Parsing gradient of QM atom " << it->index);
      err = this->parse_gradient(ofs, it->force);
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse gradient line of QM atom " << (it->index + 1)
            << " in " << out_grad;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
    // Parse capping QM atoms
    for(std::set<QM_Link>::iterator
          it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
      DEBUG(15,"Parsing gradient of capping atom " << it->qm_index << "-" << it->mm_index);
      err = this->parse_gradient(ofs, it->force);
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse gradient line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << out_grad;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
  }
  // Find MM gradients
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    bool got_mm_gradients = false;
    while (std::getline(ofs, line)) {
      if (line.find("CARTESIAN GRADIENT OF MM ATOMS") != std::string::npos) {
        got_mm_gradients = true;
        break;
      }
    }
    if (!got_mm_gradients) {
      io::messages.add("Unable to find MM gradients in output file " + out_grad,
                          this->name(), io::message::error);
      return 1;
    }
    // Parse MM atoms
    for(std::set<MM_Atom>::iterator
          it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
      DEBUG(15,"Parsing gradient of MM atom " << it->index);
      err = this->parse_gradient(ofs, it->force);
      if (it->is_polarisable) {
        DEBUG(15,"Parsing gradient on COS of MM atom " << it->index);
        err += this->parse_gradient(ofs, it->cos_force);
      }
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse gradient line of MM atom " << (it->index + 1)
            << " in " << out_grad;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
  }
  return 0;
}

int interaction::MNDO_Worker::parse_gradient(std::ifstream& ofs,
                                             math::Vec& force) const {
  std::string line;
  if (!std::getline(ofs, line)) {
    std::ostringstream msg;
    msg << "Failed to read gradient line "
        << " in " << this->param->output_gradient_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  std::istringstream iss(line);
  int dummy = 0;
  iss >> dummy >> dummy >> force(0) >> force(1) >> force(2);
  if (ofs.fail()) {
    std::ostringstream msg;
    msg << "Failed to parse gradient line in "
        << this->param->output_gradient_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  // force = - gradient
  force *= - this->param->unit_factor_force;
  return 0;
}
