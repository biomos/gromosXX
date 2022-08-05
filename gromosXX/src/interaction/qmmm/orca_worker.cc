/**
 * @file orca_worker.cc
 * interface to the Orca software package
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
#include "orca_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::Orca_Worker::Orca_Worker() : QM_Worker("Orca Worker"), param(nullptr) {};

int interaction::Orca_Worker::init(const topology::Topology& topo
                                 , const configuration::Configuration& conf
                                 , simulation::Simulation& sim
                                 , const interaction::QM_Zone& qm_zone) {

  DEBUG(15, "Initializing " << this->name());

  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.orca);
  QM_Worker::param = this->param;

  // Aliases to shorten the code
  std::string& inp = this->param->input_file;
  std::string& out = this->param->output_file;
  std::string& coordinates = this->param->input_coordinate_file;
  std::string& pointcharges = this->param->input_pointcharges_file;
  std::string& engrad = this->param->output_gradient_file;
  std::string& pcgrad = this->param->output_mm_gradient_file;

  // input file
  if (inp.empty()) {
    if (util::create_tmpfile(inp) < 1) {
      io::messages.add("Unable to create temporary input file: " + inp,
          this->name(), io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(inp);
      io::messages.add("Using temporary input file: " + inp,
        this->name(), io::message::notice);
    }
  }
  else {
    // Try to create a file specified by user
    std::ofstream of(inp.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create input file: "
        + inp, this->name(), io::message::critical);
      return 1;
    }
    else {
      io::messages.add("File initialized: "
        + inp, this->name(), io::message::notice);
    }
    of.close();
  }

  // we need to construct the names of the orca files based on input file name
  int ext_pos = inp.size() - 4;
  std::string fname; // this will be the base file name
  if (ext_pos > 0 && (inp.substr(ext_pos) == ".inp")) {
    fname = inp.substr(0, ext_pos);
  } else {
    fname = inp;
  }

  // generate file names based on input file name
  std::string out_fname = fname + ".out";
  std::string coordinates_fname = fname + ".xyz";
  std::string pointcharges_fname = fname + ".pc";
  std::string engrad_fname = fname + ".engrad";
  std::string pcgrad_fname = fname + ".pcgrad";

  // other files that orca writes but we do not need
  // there are potentially more
  std::string densities_fname = fname + ".densities";
  std::string gbw_fname = fname + ".gbw";
  std::string ges_fname = fname + ".ges";
  std::string properties_fname= fname + "_property.txt";


  // initialize files
  int err = this->initialize_file(out, out_fname);
  if (err) return err;

  err = this->initialize_file(coordinates, coordinates_fname);
  if (err) return err; 

  err = this->initialize_file(pointcharges, pointcharges_fname);
  if (err) return err; 

  err = this->initialize_file(engrad, engrad_fname);
  if (err) return err;

  err = this->initialize_file(pcgrad, pcgrad_fname);
  if (err) return err;

  // let Gromos know about temporary files so that they can be cleaned up
  this->tmp_files.insert(densities_fname);
  this->tmp_files.insert(gbw_fname);
  this->tmp_files.insert(ges_fname);
  this->tmp_files.insert(properties_fname);

  DEBUG(15, "Initialized " << this->name());
  return 0;
}

int interaction::Orca_Worker::initialize_file(std::string& file_name, const std::string& new_name) {
  if (file_name.empty()) {
    // generate temporary file
    file_name = new_name; // update file name
    this->tmp_files.insert(file_name);
    io::messages.add("Using temporary output file: " + file_name,
      this->name(), io::message::notice);
  }
  else {
    // Try to create a file specified by user
    std::ofstream of(file_name.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output file: "
        + file_name, this->name(), io::message::critical);
      return 1;
    }
    else {
      io::messages.add("File initialized: "
        + file_name, this->name(), io::message::notice);
    }
    of.close();
  }
  return 0;
}

int interaction::Orca_Worker::process_input(const topology::Topology& topo
                                        , const configuration::Configuration& conf
                                        , const simulation::Simulation& sim
                                        , const interaction::QM_Zone& qm_zone) {
  std::ofstream ifs;
  
  // First create Orca input file
  int err = this->write_input_parameters(ifs, topo, conf, sim, qm_zone);
  if (err) return err;

  // create xyz file with QM zone coordinates
  err = this->write_input_coordinates(ifs, topo, conf, sim, qm_zone);
  if (err) return err;

  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // create external file with information on point charges
    err = this->write_input_pointcharges(ifs, topo, conf, sim, qm_zone);
    if (err) return err;
  }

  return 0;
}

int interaction::Orca_Worker::write_input_parameters(std::ofstream& ifs
                                                   , const topology::Topology& topo
                                                   , const configuration::Configuration& conf
                                                   , const simulation::Simulation& sim
                                                   , const interaction::QM_Zone& qm_zone) {
  int err = this->open_input(ifs, this->param->input_file);
  if (err) return err;

  // initialize header from QM/MM specification file
  std::string header(this->param->input_header);

  // Replace variables in header
  header = io::replace_string(header, "@@CHARGE@@", std::to_string(qm_zone.charge())); 
  header = io::replace_string(header, "@@SPINM@@", std::to_string(qm_zone.spin_mult())); 
  header = io::replace_string(header, "@@POINTCHARGES@@", this->param->input_pointcharges_file);
  header = io::replace_string(header, "@@COORDINATES@@", this->param->input_coordinate_file);

  // Write header and one empty line
  ifs << header << std::endl;
  ifs.close();
  ifs.clear();

  return 0;
}

int interaction::Orca_Worker::write_input_coordinates(std::ofstream& ifs
                                                    , const topology::Topology& topo
                                                    , const configuration::Configuration& conf
                                                    , const simulation::Simulation& sim
                                                    , const interaction::QM_Zone& qm_zone) {
  int err = this->open_input(ifs, this->param->input_coordinate_file);
  if (err) return err;

  const double len_to_qm = 1.0 / this->param->unit_factor_length;

  DEBUG(15, "Writing QM coordinates");
  // number of atoms and comment line
  ifs << qm_zone.qm.size() + qm_zone.link.size() << '\n' << "Coordinates of the QM zone (in Ångström) generated by GROMOS" << std::endl;

  // write QM coordinates
  for (std::set<QM_Atom>::const_iterator it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15, it->index << " " << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }
  DEBUG(15, "Writing capping atoms coordinates");
  // write capping atoms
  for (std::set<QM_Link>::const_iterator it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; it++) {
    DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << " "
      << it->atomic_number << " " << math::v2s(it->pos * len_to_qm));
    this->write_qm_atom(ifs, it->atomic_number, it->pos * len_to_qm);
  }

  ifs.close();
  ifs.clear();

  return 0;
}

int interaction::Orca_Worker::write_input_pointcharges(std::ofstream& ifs
                                                     , const topology::Topology& topo
                                                     , const configuration::Configuration& conf
                                                     , const simulation::Simulation& sim
                                                     , const interaction::QM_Zone& qm_zone) {
  int err = this->open_input(ifs, this->param->input_pointcharges_file);
  if (err) return err;

  // write number of charges
  ifs << this->get_num_charges(sim, qm_zone) << std::endl;

  const double len_to_qm = 1.0 / this->param->unit_factor_length;
  const double cha_to_qm = 1.0 / this->param->unit_factor_charge;
  for (std::set<MM_Atom>::const_iterator it = qm_zone.mm.begin(), to = qm_zone.mm.end();
  it != to; ++it) {
    if (it->is_polarisable) {
      this->write_mm_atom(ifs, it->pos * len_to_qm, (it->charge - it->cos_charge) * cha_to_qm);
      this->write_mm_atom(ifs, (it->pos + it->cosV) * len_to_qm, it->cos_charge * cha_to_qm);
    }
    else {
      this->write_mm_atom(ifs, it->pos * len_to_qm, it->charge * cha_to_qm);
    }
  }

  ifs.close();
  ifs.clear();

  return 0;
}

void interaction::Orca_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                           , const int atomic_number
                                           , const math::Vec& pos) const {
  inputfile_stream.setf(std::ios::fixed, std::ios::floatfield);
  inputfile_stream << std::setw(6) << this->param->elements[atomic_number]
                   << std::setprecision(17)
                   << std::setw(25) << pos(0)
                   << std::setw(25) << pos(1)
                   << std::setw(25) << pos(2)
                   << std::endl;
}

void interaction::Orca_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                                , const math::Vec& pos
                                                , const double charge) const {
  if (charge != 0.0) {
    inputfile_stream.setf(std::ios::fixed, std::ios::floatfield);
    inputfile_stream << std::setprecision(6)
                     << std::setw(10) << charge
                     << std::setprecision(17)
                     << std::setw(25) << pos(0)
                     << std::setw(25) << pos(1)
                     << std::setw(25) << pos(2) 
                     << std::endl;
  }
}

int interaction::Orca_Worker::run_calculation() {
  DEBUG(15, "Calling external Orca program");
  // orca has to be called with full path in case several cores are requested
  int err = util::system_call(this->param->binary + " " + this->param->input_file + " > " + this->param->output_file);
  if (err) {
    std::ostringstream msg;
    msg << "Orca failed with code " << err;
    if (err == 127)
      msg << ". orca probably not in PATH";
    msg << ". See output file " << this->param->output_file << " for details.";
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::Orca_Worker::process_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  // parse energy
  int err = this->open_output(ofs, this->param->output_gradient_file);
  if (err) return err;
  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;
  ofs.close();
  ofs.clear();

  // parse QM gradients
  err = this->open_output(ofs, this->param->output_gradient_file);
  if (err) return err;
  err = this->parse_qm_gradients(ofs, qm_zone);
  if (err) return err;
  ofs.close();
  ofs.clear();

  // parse MM gradients or charges
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // also parse MM gradients
    err = this->open_output(ofs, this->param->output_mm_gradient_file);
    if (err) return err;
    err = this->parse_mm_gradients(ofs, qm_zone);
    if (err) return err;
  }
  else if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
             && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    // also extract charges
    err = this->open_output(ofs, this->param->output_file);
    if (err) return err;
    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
  }

  return 0;
}

int interaction::Orca_Worker::parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_energy = false;
  while(std::getline(ofs, line)) {
    // get energy section
    if (line.find("# The current total energy in Eh") != std::string::npos) {
      std::getline(ofs, line); // comment line
      std::getline(ofs, line); // read energy value
      std::istringstream iss(line);
      iss >> qm_zone.QM_energy();
      if (iss.fail()) {
        std::ostringstream msg;
        msg << "Failed to parse energy in " + this->param->output_gradient_file;
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
                  + this->param->output_gradient_file
                  , this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::Orca_Worker::parse_qm_gradients(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_qm_gradients = false;
  while (std::getline(ofs, line)) {
    if (line.find("# The current gradient in Eh/bohr") != std::string::npos) {
      got_qm_gradients = true;
      break;
    }
  }
  if (!got_qm_gradients) {
    io::messages.add("Unable to find QM gradients in output file "
                    + this->param->output_gradient_file
                    , this->name(), io::message::error);
    return 1;
  }

  // skip one line (comment)
  std::getline(ofs, line);
  // Parse QM atoms
  for (std::set<QM_Atom>::iterator it = qm_zone.qm.begin(), to = qm_zone.qm.end();
       it != to; ++it) {
    DEBUG(15, "Parsing gradients of QM atom " << it->index);
    int err = this->parse_gradient1D(ofs, it->force);
    DEBUG(15, "Force: " << math::v2s(it->force));
    if (err) {
      std::ostringstream msg;
      msg << "Failed to parse gradient line of QM atom " << (it->index + 1)
          << " in " << this->param->output_gradient_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
  }
  // Parse capping atoms
  for (std::set<QM_Link>::iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15, "Parsing gradient of capping atom " << it->qm_index << "-" << it->mm_index);
    int err = this->parse_gradient1D(ofs, it->force);
    DEBUG(15, "Force: " << math::v2s(it->force));
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

int interaction::Orca_Worker::parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_charges = false;
  while(std::getline(ofs, line)) {
    // get Mulliken charge section
    if (line.find("MULLIKEN ATOMIC CHARGES") != std::string::npos) {
      got_charges = true;
      break;
    }
  }
  if (!got_charges) {
    io::messages.add("Unable to find QM charges in output file "
                        + this->param->output_file,
                          this->name(), io::message::error);
    return 1;
  }
  std::getline(ofs, line); // part of the block header
  std::string dummy;
  // QM atoms
  for (std::set<QM_Atom>::const_iterator
             it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15, "Parsing charge of QM atom " << it->index);
    std::getline(ofs, line); // get the actual line
    std::istringstream iss(line);
    // charge is last element in line - but there may be three or four tokens in the line
    while (!iss.eof()) {
      iss >> dummy;
    }
    // last element is charge
    it->qm_charge = std::stod(dummy);
    DEBUG(15, "Charge: " << it->qm_charge);       
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse charge of QM atom " << (it->index + 1)
          << " in " << this->param->output_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->qm_charge *= this->param->unit_factor_charge;
  }
  // capping atoms
  for (std::set<QM_Link>::const_iterator
         it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15, "Parsing charge of capping atom " << it->qm_index << "-" << it->mm_index);
    std::istringstream iss(line);
    // charge is last element in line - but there may be three or four tokens in the line
    while (!iss.eof()) {
      iss >> dummy;
    }
    // last element is charge
    it->qm_charge = std::stod(dummy);
    DEBUG(15, "Charge: " << it->qm_charge);
    if (iss.fail()) {
      std::ostringstream msg;
      msg << "Failed to parse charge of capping atom " << (it->qm_index + 1)
        << "-" << (it->mm_index + 1) << " in " << this->param->output_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    it->qm_charge *= this->param->unit_factor_charge;
  }

  return 0;
}

int interaction::Orca_Worker::parse_mm_gradients(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string line;
  std::getline(ofs, line); // number of point charges
  // Parse MM atoms
  for (std::set<MM_Atom>::iterator
         it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    DEBUG(15,"Parsing gradient of MM atom " << it->index);
    int err = this->parse_gradient3D(ofs, it->force);
    DEBUG(15, "Force: " << math::v2s(it->force));
    if (err) {
      std::ostringstream msg;
      msg << "Failed to parse gradient line of MM atom " << (it->index + 1)
            << " in " << this->param->output_mm_gradient_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
    if (it->is_polarisable) {
      DEBUG(15,"Parsing gradient of COS of MM atom " << it->index);
      int err = this->parse_gradient3D(ofs, it->cos_force);
      DEBUG(15, "Force: " << math::v2s(it->cos_force));
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

int interaction::Orca_Worker::parse_gradient1D(std::ifstream& ofs,
                                                  math::Vec& force) const {
  std::string line;
  for (unsigned int i = 0; i < 3; ++i) {
      if (!std::getline(ofs, line)) {
      io::messages.add("Failed to read gradient line"
                       , this->name(), io::message::error);
      return 1;
      }
      std::istringstream iss(line);
      iss >> force(i);
      if (iss.fail()) return 1;
    }
    // force = -gradient
    force *= -this->param->unit_factor_force;
    return 0;
}

int interaction::Orca_Worker::parse_gradient3D(std::ifstream& ofs,
                                                  math::Vec& force) const {
  std::string line;
  if (!std::getline(ofs, line)) {
    io::messages.add("Failed to read gradient line"
                    , this->name(), io::message::error);
    return 1;
  }
  std::istringstream iss(line);
  iss >> force(0) >> force(1) >> force(2);
  if (iss.fail()) return 1;
  // force = - gradient
  force *= -this->param->unit_factor_force;
  return 0;
}
