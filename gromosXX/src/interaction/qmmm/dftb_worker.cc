/**
 * @file dftb_worker.cc
 * interface to the DFTB+ software package
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction.h"
#include "../../../util/timing.h"
#include "../../../util/system_call.h"

#include "qm_atom.h"
#include "mm_atom.h"
#include "qm_link.h"
#include "qm_zone.h"
#include "qm_worker.h"
#include "dftb_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::DFTB_Worker::DFTB_Worker() : QM_Worker("DFTB Worker"), param(nullptr) {};

int interaction::DFTB_Worker::init(const topology::Topology& topo
                                 , const configuration::Configuration& conf
                                 , simulation::Simulation& sim
                                 , const interaction::QM_Zone& qm_zone) {
  DEBUG(15, "Initializing " << this->name());
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.dftb);
  QM_Worker::param = this->param;
  this->cwd = this->getcwd();
  if (this->cwd == "") return 1;

  // Set input and output file, these are hardcoded in DFTB+
  this->param->input_file = "dftb_in.hsd";
  this->param->output_file = "results.tag";

  // Change to working directory and create the input file
  if (this->chdir(this->param->working_directory) != 0) return 1;
  std::ofstream ifs;
  int err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  ifs << this->param->input_header;
  // Change back to GromosXX directory
  if (this->chdir(this->cwd) != 0) return 1;
  return 0;
}

int interaction::DFTB_Worker::process_input(const topology::Topology& topo
                                        , const configuration::Configuration& conf
                                        , const simulation::Simulation& sim
                                        , const interaction::QM_Zone& qm_zone) {
  if (this->chdir(this->param->working_directory) != 0) return 1;
  std::ofstream ifs;
  int err = this->open_input(ifs, this->param->input_coordinate_file);
  if (err) return err;

  unsigned qm_size = qm_zone.qm.size() + qm_zone.link.size();
  ifs << std::setw(5) << qm_size << " C" << std::endl;
  std::map<unsigned, unsigned> type_indices;
  this->write_atom_types_list(ifs, type_indices, qm_zone);

  // Write QM atoms
  double len_to_qm = 1.0 / this->param->unit_factor_length;
  {
    unsigned i = 1;
    for (std::set<QM_Atom>::const_iterator
          it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it, ++i) {
      this->write_qm_atom(ifs, i, type_indices[it->atomic_number]
                            , it->pos * len_to_qm);
    }
    // write capping atoms
    for (std::set<QM_Link>::const_iterator
          it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it, ++i) {
      this->write_qm_atom(ifs, i, type_indices[it->atomic_number]
                            , it->pos * len_to_qm);
    }
  }
  ifs.close();
  
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    // write MM atoms
    err = this->open_input(ifs, this->param->input_mm_coordinate_file);
    if (err) return err;

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
    ifs.close();
  }
  return 0;
}

void interaction::DFTB_Worker::write_atom_types_list(std::ofstream& input_file_stream
                                                   , std::map<unsigned, unsigned>& type_indices
                                                   , const interaction::QM_Zone& qm_zone) const {
  // Get all atom types
  std::set<unsigned> atom_types;
  for (std::set<QM_Atom>::const_iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    atom_types.insert(it->atomic_number);
  }
  for (std::set<QM_Link>::const_iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    atom_types.insert(it->atomic_number);
  }
  // Write as element names and store type index
  {
    unsigned i = 1;
    for (std::set<unsigned>::const_iterator
          it = atom_types.begin(), to = atom_types.end(); it != to; ++it, ++i) {
      input_file_stream << "  " << this->param->elements[*it];
      type_indices[*it] = i;
    }
    input_file_stream << std::endl;
  }
}

void interaction::DFTB_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                           , const int id
                                           , const int atomic_type_id
                                           , const math::Vec& pos) const
  {
  inputfile_stream << std::scientific << std::setprecision(17)
                   << std::setw(5) << std::right << id << " "
                   << std::setw(3) << std::right << atomic_type_id << " "
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::endl;
}

void interaction::DFTB_Worker::write_mm_atom(std::ofstream& inputfile_stream
                                           , const math::Vec& pos
                                           , const double charge) const
  {
  inputfile_stream << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::setw(25) << std::right << charge
                   << std::endl;
}

int interaction::DFTB_Worker::run_calculation()
  {
  DEBUG(15, "Calling external DFTB+ program");
  // First delete output files
#ifdef HAVE_UNLINK
  unlink(param->output_file.c_str());
#else
  io::messages.add("unlink function is not available on this platform.", 
        this->name(), io::message::error);
  return 1;
#endif
  std::string comm = this->param->binary + " 1> " + this->param->stdout_file + " 2>&1";
  int err = util::system_call(comm);
  if (err != 0) {
    std::ostringstream msg;
    msg << "DFTB+ program failed. ";
    if (err == 127)
      msg << "Probably incorrect path to DFTB+ binary. ";
    msg << "See " << this->param->stdout_file << " for details.";
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::DFTB_Worker::process_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  
  int err = this->open_output(ofs, this->param->output_file);
  if (err) return err;

  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;

  err = this->parse_qm_forces(ofs, qm_zone);
  if (err) return err;

  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    err = this->parse_mm_forces(ofs, qm_zone);
    if (err) return err;
  }

  if (sim.param().qmmm.qmmm == simulation::qmmm_mechanical
      && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
  }

  ofs.close();
  ofs.clear();

  if (this->chdir(this->cwd) != 0) return 1;
  return 0;
}

int interaction::DFTB_Worker::parse_charges(std::ifstream& ofs
                                          , interaction::QM_Zone& qm_zone) const {
  std::string& out = this->param->output_file;
  std::string line;
  /**
   * Parse charges
   * They are used in mechanical embedding
   */
  bool got_charges = false;
  while (std::getline(ofs, line)) {
    if (line.find("gross_atomic_charges") != std::string::npos) {
      got_charges = true;
      break;
    }
  }
  if (!got_charges) {
    io::messages.add("Charges were not found in the output file "
                      + out, this->name(), io::message::error);
    return 1;
  }
  for(std::set<QM_Atom>::iterator
      it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15,"Parsing charge of QM atom " << it->index);
    ofs >> it->qm_charge;
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
  // Also for capping atoms
  for(std::set<QM_Link>::iterator
      it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15,"Parsing charge of capping atom " << it->qm_index << "-" << it->mm_index);
    ofs >> it->qm_charge;
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
  return 0;
}

int interaction::DFTB_Worker::parse_energy(std::ifstream& ofs
                                         , interaction::QM_Zone& qm_zone) const {
  std::string line;
  // Find energy block and parse
  bool got_energy = false;
  while (std::getline(ofs, line)) {
    if (line.find("total_energy") != std::string::npos) {
      ofs >> qm_zone.QM_energy();
      if (ofs.fail()) {
        io::messages.add("Failed to parse energy in output file "
                          + this->param->output_file
                          , this->name(), io::message::error);
        return 1;
      }
      got_energy = true;
      qm_zone.QM_energy() *= this->param->unit_factor_energy;
      DEBUG(15, "Parsed QM_energy: " << qm_zone.QM_energy());
      break;
    }
  }
  if (!got_energy) {
    io::messages.add("Unable to find energy in output file "
                      + this->param->output_file
                      , this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::DFTB_Worker::parse_qm_forces(std::ifstream& ofs
                                            , interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_qm_forces = false;
  while(std::getline(ofs, line)) {
    if (line.find("forces   ") != std::string::npos) {
      got_qm_forces = true;
      break;
    }
  }
  if (!got_qm_forces) {
    io::messages.add("Unable to find QM forces in output file "
                      + this->param->output_file,
                        this->name(), io::message::error);
    return 1;
  }
  // Parse QM atoms
  for(std::set<QM_Atom>::iterator
        it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
    DEBUG(15,"Parsing force of QM atom " << it->index);
    int err = this->parse_force(ofs, it->force);
    DEBUG(15,"Force: " << math::v2s(it->force));
    if (err) {
      std::ostringstream msg;
      msg << "Failed to parse force line of QM atom " << (it->index + 1)
          << " in " << this->param->output_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
  }
  // Parse capping atoms
  for(std::set<QM_Link>::iterator
        it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
    DEBUG(15,"Parsing force of capping atom " << it->qm_index << "-" << it->mm_index);
    int err = this->parse_force(ofs, it->force);
    DEBUG(15,"Force: " << math::v2s(it->force));
    if (err) {
      std::ostringstream msg;
      msg << "Failed to parse force line of capping atom " << (it->qm_index + 1)
          << "-" << (it->mm_index + 1) << " in " << this->param->output_file;
      io::messages.add(msg.str(), this->name(), io::message::error);
      return 1;
    }
  }
  return 0;
}

int interaction::DFTB_Worker::parse_mm_forces(std::ifstream& ofs
                                               , interaction::QM_Zone& qm_zone) const {
  std::string line;
  bool got_mm_forces = false;
  while(std::getline(ofs, line)) {
    if (line.find("forces_ext_charges") != std::string::npos) {
      got_mm_forces = true;
      break;
    }
  }
  if (!got_mm_forces) {
    io::messages.add("Unable to find MM forces in output file "
                      + this->param->output_file,
                        this->name(), io::message::error);
    return 1;
  }
  // Parse MM atoms
  for(std::set<MM_Atom>::iterator
        it = qm_zone.mm.begin(), to = qm_zone.mm.end(); it != to; ++it) {
    if (it->charge != 0.0 || (it->is_polarisable && (it->charge - it->cos_charge) != 0.0)) {
      DEBUG(15,"Parsing force of MM atom " << it->index);
      int err = this->parse_force(ofs, it->force);
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse force line of MM atom " << (it->index + 1)
            << " in " << this->param->output_file;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
    if (it->is_polarisable && it->cos_charge != 0.0) {
      DEBUG(15,"Parsing force of MM atom " << it->index);
      int err = this->parse_force(ofs, it->force);
      if (err) {
        std::ostringstream msg;
        msg << "Failed to parse force line of COS of MM atom " << (it->index + 1)
            << " in " << this->param->output_file;
        io::messages.add(msg.str(), this->name(), io::message::error);
        return 1;
      }
    }
  }
  return 0;
}

int interaction::DFTB_Worker::parse_force(std::ifstream& ofs,
                                          math::Vec& force) const {
  std::string line;
  ofs >> force(0) >> force(1) >> force(2);
  if(ofs.fail()) {
    std::ostringstream msg;
    msg << "Failed to parse force line in " << this->param->output_file;
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  force *= this->param->unit_factor_force;
  return 0;
}
