/**
 * @file MOPAC_worker.cc
 * interface to the MOPAC software package
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
#include "mopac_worker.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE qmmm

interaction::MOPAC_Worker::MOPAC_Worker() : QM_Worker("MOPAC Worker"), param(nullptr) {};

int interaction::MOPAC_Worker::init(const topology::Topology& topo
                                  , const configuration::Configuration& conf
                                  , simulation::Simulation& sim
                                  , const interaction::QM_Zone& qm_zone) {
  DEBUG(15, "Initializing " << this->name());
  // Get a pointer to simulation parameters
  this->param = &(sim.param().qmmm.mopac);
  QM_Worker::param = this->param;
  // Aliases to shorten the code
  std::string& inp = this->param->input_file;
  std::string& out = this->param->output_file;
  std::string& aux = this->param->output_aux_file;
  std::string& arc = this->param->output_arc_file;
  std::string& stdout = this->param->stdout_file;
  std::string& den = this->param->output_dens_file;
  std::string& molin = this->param->molin_file;
  //".mop", ".dat", or ".arc"

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
  // first we need to construct the names of the MOPAC output files from the input file name
  int ext_pos = inp.size() - 4;
  std::string fname;
  if (ext_pos > 0 && (inp.substr(ext_pos) == ".mop" || inp.substr(ext_pos) == ".dat" || inp.substr(ext_pos) == ".arc")) {
    fname = inp.substr(0, ext_pos);
  } else {
    fname = inp;
  }
  std::string out_link = fname + ".out";
  std::string aux_link = fname + ".aux";
  std::string arc_link = fname + ".arc";
  std::string den_link = fname + ".den";
  
  std::string molin_link;
  {
    size_t last_slash_pos = inp.find_last_of("/");
    if (last_slash_pos == std::string::npos) {
      molin_link = "mol.in";
    }
    else {
      std::string inp_path = inp.substr(0, last_slash_pos);
      molin_link = inp_path + "/mol.in";
    }
  }

  if (out.empty()) {
    out = out_link;
    this->tmp_files.insert(out);
    io::messages.add("Using temporary output file: " + out,
      this->name(), io::message::notice);
  }
  else {
    // Try to create a file specified by user - .out file will be symlinked to it
    // If it already exist, we dont want to overwrite it, exit
    std::ofstream of(out.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output file: "
        + out, this->name(), io::message::critical);
      return 1;
    }
    of.close();
  }

  if (aux.empty()) {
    aux = aux_link;
    this->tmp_files.insert(aux);
    io::messages.add("Using temporary output aux file: " + aux,
      this->name(), io::message::notice);
  }
  else {
    // Try to create a file specified by user - .aux file will be symlinked to it
    // If it already exist, we dont want to overwrite it, exit
    std::ofstream of(aux.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output aux file: "
        + aux, this->name(), io::message::critical);
      return 1;
    }
    of.close();
  }

  if (arc.empty()) {
    arc = arc_link;
    this->tmp_files.insert(arc);
    io::messages.add("Using temporary output arc file: " + arc,
      this->name(), io::message::notice);
  }
  else {
    // Try to create a file specified by user - .arc file will be symlinked to it
    // If it already exist, we dont want to overwrite it, exit
    std::ofstream of(arc.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output arc file: "
        + arc, this->name(), io::message::critical);
      return 1;
    }
    of.close();
  }

  if (stdout.empty()) {
    if(util::create_tmpfile(stdout) < 1) {
        io::messages.add("Unable to create temporary stdout file: " + stdout,
        this->name(), io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(stdout);
      io::messages.add("Using temporary stdout file: " + stdout,
        this->name(), io::message::notice);
    }
  }

  if (den.empty()) {
    den = den_link;
    this->tmp_files.insert(den);
    io::messages.add("Using temporary density matrix file: " + den,
      this->name(), io::message::notice);
  } // This file might exist and might be reused thus we dont raise an error if it exists
  

  if (molin.empty()) {
    if (util::create_tmpfile(molin) < 1) {
      io::messages.add("Unable to create temporary mol.in file: " + molin,
        this->name(), io::message::critical);
      return 1;
    }
    else {
      this->tmp_files.insert(molin);
      io::messages.add("Using temporary mol.in file: " + molin,
        this->name(), io::message::notice);
    }
  }


  
#ifdef HAVE_SYMLINK
  // create .out link for output file
  if (out_link != out) {
    if (symlink(out.c_str(), out_link.c_str()) != 0) {
      io::messages.add("Unable to create symbolic link from " + out_link + " to "
        + out + " - check permissions.",
        "MOPAC_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert(out_link);
  }
  // create .aux link for auxilliary output file
  if (aux_link != aux) {
    if (symlink(aux.c_str(), aux_link.c_str()) != 0) {
      io::messages.add("Unable to create symbolic link from " + aux_link + " to "
        + aux + " - check permissions.",
        "MOPAC_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert(aux_link);
  }
  // create .arc link for archive output file
  if (arc_link != arc) {
    if (symlink(arc.c_str(), arc_link.c_str()) != 0) {
      io::messages.add("Unable to create symbolic link from " + arc_link + " to "
        + arc + " - check permissions.",
        "MOPAC_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert(arc_link);
  }
  // create .den link for density file
  if (den_link != den) {
    if (symlink(den.c_str(), den_link.c_str()) != 0) {
      io::messages.add("Unable to create symbolic link from " + den_link + " to "
        + den + " - check permissions.",
        "MOPAC_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert(den_link);
  }
  // create mol.in link for QMMM input (mol.in)
  if (molin_link != molin) {
    if (symlink(molin.c_str(), molin_link.c_str()) != 0) {
      io::messages.add("Unable to create symbolic link from " + molin_link + " to "
        + molin + " - check permissions.",
        "MOPAC_Worker", io::message::critical);
      return 1;
    }
    this->tmp_files.insert(molin_link);
  }
#else
  {
    out = out_link;
    aux = aux_link;
    arc = arc_link;
    den = den_link;
    molin = molin_link;
    io::messages.add("Symbolic links are not supported in this built. "
      + "Output file is now set to " + out,
      "MOPAC_Worker", io::message::warning);
    io::messages.add("Symbolic links are not supported in this built. "
      + "Output aux file is now set to " + aux,
      "MOPAC_Worker", io::message::warning);
    io::messages.add("Symbolic links are not supported in this built. "
      + "Density matrix file is now set to " + den,
      "MOPAC_Worker", io::message::warning);
    io::messages.add("Symbolic links are not supported in this built. "
      + "mol.in file is now set to " + molin,
      "MOPAC_Worker", io::message::warning);
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

int interaction::MOPAC_Worker::process_input(const topology::Topology& topo
                                         , const configuration::Configuration& conf
                                         , const simulation::Simulation& sim
                                         , const interaction::QM_Zone& qm_zone)
  {
  std::ofstream ifs;
  // First do MOPAC input file
  int err = this->open_input(ifs, this->param->input_file);
  if (err) return err;
  std::string header(this->param->input_header);
  
  // Replace variables in header
  header = io::replace_string(header, "@@CHARGE@@", std::to_string(qm_zone.charge()));
  std::string oldens;
  if (sim.steps()) {
    oldens = "OLDENS DENOUT";
  } else {
    // in the first step we dont read the old density, only write it
    oldens = "DENOUT";
  }
  header = io::replace_string(header, "@@OLDENS@@", oldens);

  // Write header and one empty line
  ifs << header << std::endl;

  double len_to_qm = 1.0 / this->param->unit_factor_length;

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
  ifs.close();
  // Now do mol.in file (only for electrostatic and polarizable embedding)
  if (sim.param().qmmm.qmmm > simulation::qmmm_mechanical) {
    err = this->open_input(ifs, this->param->molin_file);
    if (err) return err;
    // Write empty line
    ifs << std::endl;
    // Write number of QM atoms and number of capping atoms (the latter has actually no effect)
    ifs << qm_zone.qm.size() << " " << qm_zone.link.size() << std::endl;
    // Write electrostatic potentials at the QM atoms from MM atoms
    // Calculate potential and convert to MOPAC units
    // four_pi_eps_i with conversion factors included
    double four_pi_eps_i_conv = math::four_pi_eps_i * this->param->unit_factor_charge
                              / this->param->unit_factor_energy;
    DEBUG(15,"four_pi_eps_i_conv: " << four_pi_eps_i_conv);
    // potential = four_psi_eps_i*SUM(qj/rij)
    
    // Do this for QM atoms
    for (std::set<QM_Atom>::const_iterator
          qm_it = qm_zone.qm.begin(), qm_to = qm_zone.qm.end(); qm_it != qm_to; ++qm_it) {
      double potential = four_pi_eps_i_conv * this->total_potential(topo, sim, qm_zone, *qm_it);
      DEBUG(15,"writing potential on QM atom: " << qm_it->index << " " << qm_it->atomic_number << " "
            << math::v2s(qm_it->pos * len_to_qm) << " " << potential);
      this->write_mm_potential(ifs, qm_it->atomic_number, qm_it->pos * len_to_qm, potential);
    }

    // In capping atoms, we should exclude all MM atoms
    for (std::set<QM_Link>::const_iterator
          li_it = qm_zone.link.begin(), li_to = qm_zone.link.end(); li_it != li_to; ++li_it) {
      double potential = 0.0;
      DEBUG(15,"writing potential on capping atom "
                << li_it->qm_index << "-" << li_it->mm_index << ": " << potential);
      this->write_mm_potential(ifs, li_it->atomic_number, li_it->pos * len_to_qm, potential);
    }
    ifs.close();
  }
  return 0;
}

double interaction::MOPAC_Worker::total_potential(const topology::Topology& topo
                                                , const simulation::Simulation& sim
                                                , const QM_Zone& qm_zone
                                                , const QM_Atom& qm_atom) const {
  double potential = 0.0;
  // In linked atoms, we can exclude some or all MM atoms
  /* sim.param().qmmm.mopac.link_atom_mode value:
  *   0: Link atoms see no MM atoms
  *   1: Link atoms see all MM atoms except the linked MM atom
  *   2: Link atoms see all MM atoms except the linked MM charge group
  *   3: Link atoms see all MM atoms
  */
  DEBUG(15,"QM atom: " << qm_atom.index);
  DEBUG(15,"is linked?: " << qm_atom.is_linked);
  if (!qm_atom.is_linked || sim.param().qmmm.mopac.link_atom_mode == 3) {
    for (std::set<MM_Atom>::const_iterator
            mm_it = qm_zone.mm.begin(), mm_to = qm_zone.mm.end(); mm_it != mm_to; ++mm_it) {
      DEBUG(15,"Adding potential from " << mm_it->index << ": " << this->pair_potential(qm_atom.pos, *mm_it));
      potential += this->pair_potential(qm_atom.pos, *mm_it);
    }
  }
  else if (sim.param().qmmm.mopac.link_atom_mode == 0) {
    DEBUG(15,"QM atom is link atom, mode 0, potential 0.0");
    potential = 0.0;
  }
  else {
    // get linked MM atoms
    std::set<unsigned> excluded_mm;
    DEBUG(15,"Generating a set of excluded MM atoms");
    this->get_excluded_mm(topo, sim, qm_zone, qm_atom, excluded_mm);
    // Calculate potential on link atom with exclusions
    for (std::set<MM_Atom>::const_iterator
            mm_it = qm_zone.mm.begin(), mm_to = qm_zone.mm.end(); mm_it != mm_to; ++mm_it) {
      if (excluded_mm.find(mm_it->index) == excluded_mm.end()) {
        DEBUG(15, "Adding potential from " << mm_it->index << ": " << this->pair_potential(qm_atom.pos, *mm_it));
        potential += this->pair_potential(qm_atom.pos, *mm_it);
      } else {
        DEBUG(15, "MM atom " << mm_it->index << " excluded");
      }
    }
  }
  return potential;
}

void interaction::MOPAC_Worker::get_excluded_mm(const topology::Topology& topo
                                              , const simulation::Simulation& sim
                                              , const QM_Zone& qm_zone
                                              , const QM_Atom& qm_atom
                                              , std::set<unsigned> excluded) const {
  for (std::set<QM_Link>::const_iterator
        li_it = qm_zone.link.begin(), li_to = qm_zone.link.end(); li_it != li_to; ++li_it) {
    if (li_it->qm_index == qm_atom.index) {
      excluded.insert(li_it->mm_index);
      DEBUG(15, "Excluding MM atom " << li_it->mm_index);
    }
  }
  // Optionally also collect charge groups
  if (sim.param().qmmm.mopac.link_atom_mode == 2) {
    std::set<unsigned> excluded_cgs;
    for (std::set<unsigned>::const_iterator
          it = excluded.begin(), to = excluded.end(); it != to; ++it) {
      // Find chargegroups of the linked MM atoms
      int it_cg = -1;
      for (unsigned cg = 0; cg < topo.num_chargegroups(); ++cg) {
        if (*it < unsigned(topo.chargegroup(cg + 1))) {
          it_cg = cg;
          break;
        }
      }
      assert(it_cg != -1);
      excluded_cgs.insert(it_cg);
    }
    // Translate chargegroup indices into MM atom indices
    for (std::set<unsigned>::const_iterator
          it = excluded_cgs.begin(), to = excluded_cgs.end(); it != to; ++it) {
      for (int i = topo.chargegroup(*it); i < topo.chargegroup(*it + 1); ++i) {
        excluded.insert(i);
      }
    }
  }
  return;
}
  

double interaction::MOPAC_Worker::pair_potential(const math::Vec& pos
                                               , const MM_Atom& mm_atom) const {
  if (mm_atom.is_polarisable) {
    // MM charge
    double potential = (mm_atom.charge - mm_atom.cos_charge) / math::abs(pos - mm_atom.pos);
    // COS charge
    potential += mm_atom.cos_charge / math::abs(pos - mm_atom.pos - mm_atom.cosV);
    return potential;
  }
  else {
    // MM charge only
    return mm_atom.charge / math::abs(pos - mm_atom.pos);
  }
}

int interaction::MOPAC_Worker::run_calculation() {
  // MOPAC writes automatically to inputfilename.out
  // We redirect stderr to stdout and store in separate file
  DEBUG(15, "Calling external MOPAC program");
  int err = util::system_call(this->param->binary + " " + this->param->input_file
                                + " 1> " + this->param->stdout_file + " 2>&1 ");
  if (err) {
    std::ostringstream msg;
    msg << "MOPAC failed with code " << err;
    if (err == 127)
      msg << ". mopac command probably not in PATH";
    msg << ". See output file " << this->param->output_file << " for details.";
    io::messages.add(msg.str(), this->name(), io::message::error);
    return 1;
  }
  return 0;
}

int interaction::MOPAC_Worker::process_output(topology::Topology& topo
                                        , configuration::Configuration& conf
                                        , simulation::Simulation& sim
                                        , interaction::QM_Zone& qm_zone) {
  std::ifstream ofs;
  int err = this->open_output(ofs, this->param->output_aux_file);
  if (err) return err;

  err = this->parse_energy(ofs, qm_zone);
  if (err) return err;

  if (minimisation) {
    err = this->parse_coordinates(ofs, qm_zone);
    if (err) return err;
  }
  
  // here we always parse charges to calculate QM-MM interactions manually
  // except for mechanical embedding and constant QM charges
  if (!( sim.param().qmmm.qmmm == simulation::qmmm_mechanical
          && sim.param().qmmm.qm_ch == simulation::qm_ch_constant )) {
    err = this->parse_charges(ofs, qm_zone);
    if (err) return err;
  }

  err = this->parse_qm_gradients(sim, ofs, qm_zone);
  if (err) return err;

  this->calculate_mm_forces(topo, sim, qm_zone);

  ofs.close();
  return 0;
}

void interaction::MOPAC_Worker::write_qm_atom(std::ofstream& inputfile_stream
                                        , const int atomic_number
                                        , const math::Vec& pos
                                        , const int var_flag) const
  {
  inputfile_stream << std::setw(4) << std::left << atomic_number
                   << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(4)  << std::right << var_flag
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(4)  << std::right << var_flag
                   << std::setw(25) << std::right << pos(2)
                   << std::setw(4)  << std::right << var_flag
                   << std::endl;
}

void interaction::MOPAC_Worker::write_mm_potential(std::ofstream& inputfile_stream
                                                 , const int atomic_number
                                                 , const math::Vec& pos
                                                 , double potential) const
  {
  inputfile_stream << std::setw(4) << std::left << atomic_number
                   << std::scientific << std::setprecision(17)
                   << std::setw(25) << std::right << pos(0)
                   << std::setw(25) << std::right << pos(1)
                   << std::setw(25) << std::right << pos(2)
                   << std::setw(25) << std::right << potential
                   << std::endl;
}

int interaction::MOPAC_Worker::parse_charges(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string& out_aux = this->param->output_aux_file;
  std::string line;
  /**
   * Parse charges
   * They are used in all embedding schemes
   */
  while (std::getline(ofs, line)) {
    if (line.find("ATOM_CHARGES") != std::string::npos) {
      for(std::set<QM_Atom>::iterator
            it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
        ofs >> it->qm_charge;
        if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse charge of QM atom " << (it->index + 1)
              << " in " << out_aux;
          io::messages.add(msg.str(), this->name(), io::message::error);
          return 1;
        }
        it->qm_charge *= this->param->unit_factor_charge;
        DEBUG(15, "QM atom " << it->index << " new charge: " << it->qm_charge);
      }
      // Also parse capping atoms
      for(std::set<QM_Link>::iterator
            it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
        ofs >> it->qm_charge;
        if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse charge of capping atom " << (it->qm_index + 1)
              << "-" << (it->mm_index + 1) << " in " << out_aux;
          io::messages.add(msg.str(), this->name(), io::message::error);
          return 1;
        }
        it->qm_charge *= this->param->unit_factor_charge;
        DEBUG(15, "Capping atom " << it->qm_index << "-" << it->mm_index << " new charge: " << it->qm_charge);
      }
      return 0;
    }
  }
  io::messages.add("Charges were not found in the output aux file "
                    + out_aux, this->name(), io::message::error);
  return 1;
}

int interaction::MOPAC_Worker::parse_coordinates(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string& out_aux = this->param->output_aux_file;
  std::string line;
  // Find coordinates block
  while (std::getline(ofs, line)) {
    if (line.find("ATOM_X_OPT") != std::string::npos) {
      for(std::set<QM_Atom>::iterator
            it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
        ofs >> it->pos(0) >> it->pos(1) >> it->pos(2);
        if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse coordinates of QM atom " << (it->index + 1)
              << " in " << out_aux;
          io::messages.add(msg.str(), this->name(), io::message::error);
          return 1;
        }
        it->pos *= this->param->unit_factor_length;
      }
      // Also parse capping atoms
      for(std::set<QM_Link>::iterator
            it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
        ofs >> it->pos(0) >> it->pos(1) >> it->pos(2);
        if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse coordinates of capping atom " << (it->qm_index + 1)
              << "-" << (it->mm_index + 1) << " in " << out_aux;
          io::messages.add(msg.str(), this->name(), io::message::error);
          return 1;
        }
        it->pos *= this->param->unit_factor_length;
      }
      return 0;
    }
  }
  io::messages.add("Coordinates were not found in the output aux file "
                    + out_aux, this->name(), io::message::error);
  return 1;
}

int interaction::MOPAC_Worker::parse_energy(std::ifstream& ofs, interaction::QM_Zone& qm_zone) const {
  std::string line;
  // Find energy line
  while (std::getline(ofs, line)) {
    if (line.find("TOTAL_ENERGY") != std::string::npos) {
      line = line.substr(line.find("=") + 1);
      this->defortranize(line);
      qm_zone.QM_energy() = this->param->unit_factor_energy * atof(line.c_str());
      return 0;
    }
  }
  io::messages.add("Unable to find energy in output aux file "
                    + this->param->output_aux_file
                    , this->name(), io::message::error);
  return 1;
}

int interaction::MOPAC_Worker::parse_qm_gradients(const simulation::Simulation& sim
                                            , std::ifstream& ofs
                                            , interaction::QM_Zone& qm_zone) const {
  std::string& out_aux = this->param->output_aux_file;
  std::string line;
  
  // MOPAC gradients dont contain MM contribution - this has to be calculated from
  // the QM and MM charges
  // Find QM gradients
  while (std::getline(ofs, line)) {
    if (line.find("GRADIENTS") != std::string::npos) {
      for(std::set<QM_Atom>::iterator
            it = qm_zone.qm.begin(), to = qm_zone.qm.end(); it != to; ++it) {
        ofs >> it->force(0) >> it->force(1) >> it->force(2);
        if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse gradient of QM atom " << (it->index + 1)
              << " in " << out_aux;
          io::messages.add(msg.str(), this->name(), io::message::error);
          return 1;
        }
        it->force *= - this->param->unit_factor_force;
      }
      // Also parse capping atoms
      for(std::set<QM_Link>::iterator
            it = qm_zone.link.begin(), to = qm_zone.link.end(); it != to; ++it) {
        ofs >> it->force(0) >> it->force(1) >> it->force(2);
        if (ofs.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse gradient of capping atom " << (it->qm_index + 1)
              << "-" << (it->mm_index + 1) << " in " << out_aux;
          io::messages.add(msg.str(), this->name(), io::message::error);
          return 1;
        }
        it->force *= - this->param->unit_factor_force;
      }
      return 0;
    }
  }
  io::messages.add("Gradients were not found in the output aux file "
                    + out_aux, this->name(), io::message::error);
  return 1;
}

void interaction::MOPAC_Worker::calculate_mm_forces(const topology::Topology& topo
                                                  , const simulation::Simulation& sim
                                                  , interaction::QM_Zone& qm_zone) const {
  for (std::set<QM_Atom>::iterator qm_it = qm_zone.qm.begin(), qm_to = qm_zone.qm.end();
        qm_it != qm_to; ++qm_it) {
    const double qmq_four_pi_eps_i = math::four_pi_eps_i * qm_it->qm_charge;
    
    DEBUG(15,"QM atom: " << qm_it->index);
    DEBUG(15,"is linked?: " << qm_it->is_linked);
    if (!qm_it->is_linked || sim.param().qmmm.mopac.link_atom_mode == 3) {
      for (std::set<MM_Atom>::iterator
            mm_it = qm_zone.mm.begin(), mm_to = qm_zone.mm.end(); mm_it != mm_to; ++mm_it) {
        DEBUG(15, "Calculating force with MM atom " << mm_it->index);
        this->calculate_pair_force(qm_it->pos, mm_it->pos, qm_it->force, mm_it->force
                                  , qmq_four_pi_eps_i * mm_it->charge);        
      }
    }
    else if (sim.param().qmmm.mopac.link_atom_mode == 0) {
      DEBUG(15,"QM atom is link atom, mode 0, QM-MM force 0.0");
    }
    else {
      // get linked MM atoms
      std::set<unsigned> excluded_mm;
      DEBUG(15,"Generating a set of excluded MM atoms");
      this->get_excluded_mm(topo, sim, qm_zone, *qm_it, excluded_mm);
      // Calculate forces on link atom with exclusions
      for (std::set<MM_Atom>::iterator
              mm_it = qm_zone.mm.begin(), mm_to = qm_zone.mm.end(); mm_it != mm_to; ++mm_it) {
        if (excluded_mm.find(mm_it->index) == excluded_mm.end()) {
          DEBUG(15, "Calculating force with MM atom " << mm_it->index);
          this->calculate_pair_force(qm_it->pos, mm_it->pos, qm_it->force, mm_it->force
                                    , qmq_four_pi_eps_i * mm_it->charge);
        } else {
          DEBUG(15, "MM atom " << mm_it->index << " excluded");
        }
      }
    }
  }
}

void interaction::MOPAC_Worker::calculate_pair_force(const math::Vec& qm_pos
                                                   , const math::Vec& mm_pos
                                                   , math::Vec& qm_force
                                                   , math::Vec& mm_force
                                                   , const double qmq_mmq_four_pi_eps_i) const {
  const math::Vec r = qm_pos - mm_pos;
  const double d2 = math::abs2(r);
  const double d = sqrt(d2);
  const double d3 = d2 * d;
  const double qq_four_pi_eps_i_d3 = qmq_mmq_four_pi_eps_i / d3;
  const math::Vec force = qq_four_pi_eps_i_d3 * r;
  DEBUG(15, "Force:" << math::v2s(force));
  qm_force += force;
  mm_force -= force;
  return;
}
