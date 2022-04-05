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

int interaction::Orca_Worker::init(simulation::Simulation& sim) {
  DEBUG(15, "Initializing " << this->name());
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.orca);
  QM_Worker::param = this->param;
  // Aliases to shorten the code
  std::string& inp = this->param->input_file;
  std::string& out = this->param->output_file;
  std::string& pointcharges = this->param->input_mm_coordinate_file;
  std::string& engrad = this->param->output_gradient_file;
  std::string& pcgrad = this->param->output_mm_gradient_file;

  // Check set filenames and generate temporary files where necessary
  if (inp.empty()) {
      if(util::create_tmpfile(inp) < 1) {
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

  // create links for output files
  // first we need to construct the names of the Orca output files from the input file name
  int ext_pos = inp.size() - 4;
  std::string fname;
  if (ext_pos > 0 && (inp.substr(ext_pos) == ".out")) {
    fname = inp.substr(0, ext_pos);
  }
  else {
    fname = inp;
  }
  std::string out_link = fname + ".out";

  if (out.empty()) {
    out = out_link;
    this->tmp_files.insert(out);
    io::messages.add("Using temporary output file: " + out,
      this->name(), io::message::notice);
  }
  else {
    // Try to create a file specified by user - .out file will be symlinked to it
    // it it already exists, we don't want to overwrite it, exit
    std::ofstream of(out.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output file: "
        + out, this->name(), io::message::critical);
      return 1;
    }
    of.close();
  }

#ifdef HAVE_SYMLINK
  // create .out link for output file
  if (out_link != out) {
    if (symlink(out.c_str(), out_link.c_str()) != 0) {
      io::messages.add("Unable to create symbolic link from " + out_link + " to "
      + out + " - check permissions.",
      this->name(), io::message::critical);
      return 1;
    }
    this->tmp_files.insert(out_link);
  }
#else
  {
    out = out_link;
    io::messages.add("Symbolic links are not supported in this build. "
    + "Output file is now set to " + out,
    this->name(), io::message::warning);
  }
#endif

#ifndef HAVE_UNLINK
  {
    io::messages.add("Unlink function not supported on this platform. "
    + "Please delete temporary files manually.",
    this->name(), io::message::warning);
  }
#endif
  DEBUG(15, "Initialized " << this->name());
  return 0;
}

int interaction::Orca_Worker::write_input(const topology::Topology& topo
                                        , const configuration::Configuration& conf
                                        , const simulation::Simulation& sim
                                        , const interaction::QM_Zone& qm_zone) {
  std::ofstream ifs;
  int err;
  // First create Orca inout file
  err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  std::string header(this->param->input_header);

  // Replace variables in header
  header = io::replace_string(header, "@@CHARGE@@", std::to_string(qm_zone.charge()));  

  double len_to_qm = 1.0 / this->param->unit_factor_length;

  DEBUG(15, "Writing QM coordinates");
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

  return 0;
}

int interaction::Orca_Worker::system_call() {
  DEBUG(15, "Calling external Orca program");
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

int interaction::Orca_Worker::read_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {

}

void interaction::Orca_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                           , const int atomic_number
                                           , const math::Vec& pos
                                           , const int var_flag) const {
  inputfile_stream << std::setw(4) << std::left << atomic_number
                   << std::scientific << std::setprecision(17)
                   << std::setw(25)   << std::right << pos(0)
                   << std::setw(4)    << std::right << var_flag
                   << std::setw(25)   << std::right << pos(1)
                   << std::setw(4)    << std::right << var_flag
                   << std::setw(25)   << std::right << pos(2)
                   << std::setw(4)    << std::right << var_flag
                   << std::endl;
}

void interaction::Orca_Worker::write_mm_potential(std::ofstream& inputfile_stream
                                                , const int atomic_number
                                                , const math::Vec& pos
                                                , double potential) const {

}

double interaction::Orca_Worker::total_potential(const topology::Topology& topology
                                               , const simulation::Simulation& simulation
                                               , const QM_Zone& qm_zone
                                               , const QM_Atom& qm_atom) const {

}

void interaction::Orca_Worker::get_excluded_mm(const topology::Topology& topo
                                             , const simulation::Simulation& sim
                                             , const QM_Zone& qm_zone
                                             , const QM_Atom& qm_atom
                                             , std::set<unsigned> excluded) const {
}

double interaction::Orca_Worker::pair_potential(const math::Vec& pos
                                              , const MM_Atom& mm_atom) const {

}

int interaction::Orca_Worker::parse_charges(std::ifstream& ofs
                                          , interaction::QM_Zone& qm_zone) const {

}

int interaction::Orca_Worker::parse_coordinates(std::ifstream& ofs
                                              , interaction::QM_Zone& qm_zone) const {

}

int interaction::Orca_Worker::parse_energy(std::ifstream& ofs
                                         , interaction::QM_Zone& qm_zone) const {

}

int interaction::Orca_Worker::parse_qm_gradients(const simulation::Simulation& sim
                                               , std::ifstream& ofs
                                               , interaction::QM_Zone& qm_zone) const {

}

void interaction::Orca_Worker::calculate_mm_forces(const topology::Topology& topo
                                                 , const simulation::Simulation& sim
                                                 , interaction::QM_Zone& qm_zone) const {

}

void interaction::Orca_Worker::calculate_pair_force(const math::Vec& qm_pos
                                                  , const math::Vec& mm_pos
                                                  , math::Vec& qm_force
                                                  , math::Vec& mm_force
                                                  , const double qmq_mmq_four_pi_eps_i) const {

}