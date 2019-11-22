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

//#include "../../../math/periodicity.h"

// special interactions
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

int interaction::MNDO_Worker::init(const simulation::Simulation& sim) {
  /** Make a copy of parameters, so QM worker can safely customize
   * them and also noone can change them
   */
  this->param = new simulation::Parameter::qmmm_struct::mndo_param_struct(sim.param().qmmm.mndo);
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
    this->using_tmp = true;
    if(util::create_tmpfile(inp) < 1) {
        io::messages.add("Unable to create temporary input file: " + inp,
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
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
        io::messages.add("Using temporary output gradient file: " + out_grad,
        "MNDO_Worker", io::message::notice);
    }
  }
  else {
    // Try to create a file specified by user
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
        io::messages.add("Using temporary density matrix file: " + dens,
        "MNDO_Worker", io::message::notice);
    }
  }
  
#ifdef HAVE_SYMLINK
  // create fort.15 link for gradients
  if ((this->symlink_err = symlink(out_grad.c_str(), "fort.15")) != 0) {
    io::messages.add("Unable to create symbolic link from fort.15 to "
      + out_grad + " - check permissions.",
      "MNDO_Worker", io::message::critical);
    return 1;
  }
  // create fort.11 link for density matrix
  if ((this->symlink_err = symlink(dens.c_str(), "fort.11")) != 0) {
    io::messages.add("Unable to create symbolic link from fort.11 to "
      + dens + " - check permissions.",
      "MNDO_Worker", io::message::critical);
    return 1;
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
    "MNDO_Worker", io::message::notice);
  }
#endif
  return 0;
}

int interaction::MNDO_Worker::write_input(const topology::Topology& topo
                                        , const configuration::Configuration& conf
                                        , const simulation::Simulation& sim
                                        , const interaction::QM_Zone& qm_zone)
  {
  std::ofstream ifs;
  int err;
  err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  std::string header(this->param->input_header);
  
  // Get number of links and replace
  unsigned num_links = qm_zone.link.size();
  header = io::replace_string(header, "@@NUM_LINKS@@", std::to_string(num_links));
  
  // Get number of MM charges and replace
  unsigned num_charges;
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
      io::messages.add("Uknown QMMM option", io::message::error);
      return 1;
    }
  }
  header = io::replace_string(header, "@@NUM_CHARGES@@", std::to_string(num_charges));

  // Write header
  ifs << header << std::endl;

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
  // Write termination line
  this->write_qm_atom(ifs, 0, math::Vec(0.0));

  // Write capping atom numbers - they are last in the geometry
  const unsigned qm_size = qm_zone.qm.size();
  unsigned last_link = num_links + qm_size;
  for (unsigned i = qm_size + 1; i <= last_link; ++i)
    ifs << i << " ";
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

int interaction::MNDO_Worker::system_call() {
  int err = util::system_call(this->param->binary, this->param->input_file, this->param->output_file);
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
int interaction::MNDO_Worker::read_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  int err;
  err = this->open_output(ofs, this->param->output_file);
  if (err) return err;
  
  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical) {
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
                                        , const math::Vec& pos)
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


  bool opt_flag = 0; // optimization flag - if the atom position should be optimized
  inputfile_stream << std::setw(2) << std::left << atomic_number
                   << std::string(2,' ')
                   << std::setw(20) << std::setprecision(15) << std::right << pos(0)
                   << std::string(2,' ') << std::setw(2) << std::right << opt_flag
                   //<< std::string(6,' ')
                   << std::setw(20) << std::setprecision(15) << std::right << pos(1)
                   << std::string(2,' ') << std::setw(2) << std::right << opt_flag
                   //<< std::string(6,' ')
                   << std::setw(20) << std::setprecision(15) << std::right << pos(2)
                   << std::string(2,' ') << std::setw(2) << std::right << opt_flag
                   << std::endl;
}

void interaction::MNDO_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                        , const math::Vec& pos
                                        , const double charge)
  {
    /*
      cm(1,m)    1-12  f12.7  x coordinate.
      cm(2,m)   13-24  f12.7  y coordinate.
      cm(3,m)   25-36  f12.7  z coordinate.
      qm(m)     37-44   f8.4  External point charge.
    */
  inputfile_stream << std::setw(20) << std::setprecision(15) << std::right << pos(0)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(1)
                   << std::setw(20) << std::setprecision(15) << std::right << pos(2)
                   << std::setw(20) << std::setprecision(4) << std::right << charge
                   << std::endl;
}

int interaction::MNDO_Worker::parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
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
  return 0;
}

int interaction::MNDO_Worker::parse_coordinates(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
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
  math::Vec pos;
  int dummy;
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
    iss >> dummy >> dummy >> pos(0) >> pos(1) >> pos(2);
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse coordinate line of atom " << (it->index + 1)
          << " in " + out_grad;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->pos = pos * this->param->unit_factor_length;
  }
  return 0;
}

int interaction::MNDO_Worker::parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) {
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
  std::getline(ofs, line);
  std::istringstream iss(line);
  iss >> qm_zone.QM_energy();
  if (iss.fail()) {
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
                                            , interaction::QM_Zone& qm_zone) {
  std::string& out_grad = this->param->output_gradient_file;
  std::string line;
  int err;
  
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
    err = this->_parse_gradients(ofs, qm_zone.qm);
    if (err) return err;
    // Parse link QM atoms
    err = this->_parse_gradients(ofs, qm_zone.link);
    if (err) return err;
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
    err = this->_parse_gradients(ofs, qm_zone.mm);
    if (err) return err;
  }
  return 0;
}

template<class AtomType>
int interaction::MNDO_Worker::_parse_gradients(std::ifstream& ofs,
                                               std::set<AtomType>& atom_set) {
  int err = 0;
  const double unit_factor =
      this->param->unit_factor_energy / this->param->unit_factor_length;
  typename std::set<AtomType>::iterator it, to;
  for(it = atom_set.begin(), to = atom_set.end(); it != to; ++it) {
    err = this->parse_gradient(ofs, it->index, it->force, unit_factor);
    if (err) return err;
  }
  return 0;
}

template<>
int interaction::MNDO_Worker::_parse_gradients
                                        (std::ifstream& ofs
                                       , std::set<interaction::MM_Atom>& atom_set) {
  int err = 0;
  const double unit_factor =
      this->param->unit_factor_energy / this->param->unit_factor_length;
  std::set<interaction::MM_Atom>::iterator it, to;
  for(it = atom_set.begin(), to = atom_set.end(); it != to; ++it) {
    unsigned i = it->index;
    err = this->parse_gradient(ofs, i, it->force, unit_factor);
    if (err) return err;
    if (it->is_polarisable) {
      err = this->parse_gradient(ofs, i, it->cos_force, unit_factor);
      if (err) return err;
    }
  }
  return 0;
}

int interaction::MNDO_Worker::parse_gradient(std::ifstream& ofs,
                                             const unsigned index,
                                             math::Vec& force,
                                             const double unit_factor) {
  std::string line;
  math::Vec gradient;
  int dummy;
  if(!std::getline(ofs, line)) {
    std::ostringstream msg;
    msg << "Failed to read gradient line of atom " << (index)
        << " in " << this->param->output_gradient_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  std::istringstream iss(line);
  iss >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
  if (iss.fail()) {
    std::ostringstream msg;
    msg << "Failed to parse gradient line of atom " << (index)
        << " in " << this->param->output_gradient_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  force = - gradient * unit_factor;
  return 0;
}

interaction::MNDO_Worker::~MNDO_Worker() {
#ifdef HAVE_UNLINK
  // Remove symbolic links
  if (this->symlink_err == 0) {
    unlink("fort.11");
    unlink("fort.15");
  }
  // Delete temporary files
  if (this->using_tmp) {
    unlink(this->param->input_file.c_str());
    unlink(this->param->output_file.c_str());
    unlink(this->param->output_gradient_file.c_str());
    unlink(this->param->density_matrix_file.c_str());
  }
  if (this->param != nullptr) {
    delete this->param;
    this->param = nullptr;
  }
#endif
}
