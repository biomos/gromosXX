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
#include "../../../util/timing.h"

//#include "../../../math/periodicity.h"

// special interactions
#include "qm_storage.h"
#include "mm_atom.h"
#include "qm_worker.h"
#include "mndo_worker.h"
#include "../../../util/system_call.h"
#include "../../../io/blockinput.h"
#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special


int interaction::MNDO_Worker::run_QM(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        const math::VArray & qm_pos,
        const std::vector<interaction::MM_Atom> & mm_atoms,
        interaction::QM_Storage& storage,
        interaction::QM_Storage& LA_storage,
        const configuration::Configuration& qmmm_conf){
     

/*
  if (mm_atoms.empty()) {
    io::messages.add("Cannot deal with zero MM atoms yet.", "MNDO_Worker", 
            io::message::error);
    return 1;
  }
  */
  bool verbose=false;

  std::ofstream inp(this->input_file.c_str());
  if (!inp.is_open()) {
    io::messages.add("Unable to write to input file: "
            + this->input_file, "MNDO_Worker", io::message::critical);
    return 1;
  }
  std::string header(this->input_header);
  // get the number of point charges
  unsigned num_charge = mm_atoms.size();
  for(unsigned i = 0; i < mm_atoms.size(); ++i) {
    if (topo.is_polarisable(mm_atoms[i].index)) {
      ++num_charge;
    }
  }
  std::ostringstream oss;
  oss << num_charge;
  header = io::replace_string(header, "@@NUM_ATOMS@@", oss.str());

  // get the number of links
  std::vector<unsigned int> link_atoms;
  unsigned int li = 1;
  for (std::set<topology::qm_atom_struct>::const_iterator
    it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++li) {
    if (it->link) {
      link_atoms.push_back(li);
    }
  }
  //std::ostringstream oss; oss << num_charge ;
  header = io::replace_string(header, "@@NUM_ATOMS@@", oss.str());


  oss.str("");
  oss.clear();
  oss << link_atoms.size();
  header = io::replace_string(header, "@@NUM_LINKS@@", oss.str());

  // write header
  inp << header << std::endl;

  // write QM zone
  double len_to_qm = 1.0 / sim.param().qmmm.unit_factor_length;
  unsigned int pi = 0;
  for (std::set<topology::qm_atom_struct>::const_iterator
    it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++pi) {
    

    inp.setf(std::ios::fixed, std::ios::floatfield);
    inp.precision(8);

    inp << std::setw(2) << std::left << it->atomic_number;
    for (unsigned int i = 0; i < 3; ++i) {
      inp << std::setw(14) << std::right << qm_pos(pi)(i) * len_to_qm << " 0";
    }
    inp << std::endl;
  }

 pi=0;
  
  // append list of link-atoms to the QM coordinates
  //   the order of the capping H atoms follows the same
  //  order as the QM atoms in the "qmmm" file.
  //
  math::Vec posCap,dR;
  unsigned int m1,q1;
  const double rch = sim.param().qmmm.cap_dist; //Carbon-Hydrogen bond length
  //create a conf for capped system
  //configuration::Configuration qmmm_conf = conf;
  //qmmm_conf.current().force=0.0;

  if (verbose) {
     std::cout << "number of link atoms: " << topo.qm_mm_pair().size()  << std::endl;
  }
  for (std::vector< std::pair<unsigned int,unsigned int> >::const_iterator
    it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it,++pi  )
  {  
     q1=it->first;
     m1=it->second;  
     //dR =conf.current(q1).pos - conf.current(m1).pos;
     dR = conf.current().pos(m1) - conf.current().pos(q1);
     posCap = (rch/abs(dR)    )    * dR + conf.current().pos(q1) ;
    // std::cout << "modifying position of " << m1 << " from " << qmmm_conf.current().pos(m1)[0] << "to  "  <<
    //           posCap(0) << std::endl;
   //  qmmm_conf.current().pos(m1) = posCap;
     inp << std::setw(5) << std::left << "1 " ;  //link_atoms[i];
    for (unsigned int i = 0; i < 3; ++i) {
      inp << std::setw(14) << std::right << posCap(i) * len_to_qm << " 0"  ;
    }
    if (++pi % 16 == 0 || pi == link_atoms.size()) {
      inp << std::endl;
    }
    inp << std::endl ;
  }  

  // the termination line. just a couple of zeros
  inp << std::setw(2) << std::left << 0;
  for (unsigned int i = 0; i < 3; ++i) {
    inp << std::setw(14) << std::right << 0.0 << " 0";
  }
  inp << std::endl;
  //write list of link atoms


  pi=0;
  for (std::set<topology::qm_atom_struct>::const_iterator
       it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++pi) {
  {
    if (it->link == true ) {
      inp << it->index - topo.qm_zone().begin()->index + 1 << " " ;
    }
  }
  }

  inp << std::endl;
  // write point charges
  //  comment out for simulations in vacuum
  double chg_to_qm = 1.0 / sim.param().qmmm.unit_factor_charge;
  
  for (unsigned int i = 0; i < mm_atoms.size(); ++i) {
    if (topo.is_polarisable(mm_atoms[i].index)) {
      inp.setf(std::ios::fixed, std::ios::floatfield);
      inp.precision(8);
      for (unsigned int j = 0; j < 3; ++j) {
        inp << std::setw(14) << std::right << mm_atoms[i].pos(j) * len_to_qm;
      }
      inp << std::setw(14) << std::right <<
              (mm_atoms[i].charge - topo.coscharge(mm_atoms[i].index)) * chg_to_qm << std::endl;
      const math::Vec pos_cos(mm_atoms[i].pos + conf.current().posV(mm_atoms[i].index));
      for (unsigned int j = 0; j < 3; ++j) {
        inp << std::setw(14) << std::right << pos_cos(j) * len_to_qm;
      }
      inp << std::setw(14) << std::right << topo.coscharge(mm_atoms[i].index) * chg_to_qm << std::endl;
    } else {
      inp.setf(std::ios::fixed, std::ios::floatfield);
      inp.precision(8);
      for (unsigned int j = 0; j < 3; ++j) {
        inp << std::setw(14) << std::right << mm_atoms[i].pos(j) * len_to_qm;
      }
      inp << std::setw(14) << std::right << mm_atoms[i].charge * chg_to_qm << std::endl;
    }
  }
  inp.close();
  
  // starting MNDO
  int result = util::system_call(sim.param().qmmm.mndo.binary,
          this->input_file, this->output_file);
  if (result != 0) {
    std::ostringstream msg;
    msg << "MNDO failed with code " << result;
    if (result == 127)
      msg << ". mndo command probably not in PATH";
    msg << ". See output file " << this->output_file << " for details.";
    io::messages.add(msg.str(), "MNDO_Worker", io::message::error);
    return 1;
  }
  
  // read output
  std::ifstream output(this->output_file.c_str());
  if (!output.is_open()) {
    io::messages.add("Cannot open MNDO output file", "MNDO_Worker", 
            io::message::error);
    return 1;
  }
  while(true) {
    std::string line;
    std::getline(output, line);
    if (output.eof() || output.fail()) {
      break;
    }
    else if (line.find("FATAL INPUT ERROR") != std::string::npos) {
      // MNDO doesnt return 1 if it fails on input, so this is a workaround
      io::messages.add("MNDO input error. See output file "
        + this->output_file + " for details.", "MNDO_Worker", io::message::error);
      return 1;
    }
    // get charges
    if (line.find("NET ATOMIC CHARGES") != std::string::npos) {
      // skip 3 lines
      for (unsigned i = 0; i < 3; ++i) std::getline(output, line);
      // read atoms
      for (unsigned i = 0; i < topo.qm_zone().size(); ++i) {
        std::getline(output, line);
        if (output.fail()) {
          io::messages.add("Failed to read charge line of " + std::to_string(i+1) + "th atom.",
            "MNDO_Worker", io::message::error);
          return 1;
        }
        
        std::istringstream is(line);
        std::string dummy; int idummy;
        double charge;
        is >> dummy >> idummy >> charge;
        if (is.fail()) {
          io::messages.add("Failed to parse charge line of " + std::to_string(i+1) + "th atom.",
            "MNDO_Worker", io::message::error);
          return 1;
        }
        storage.charge(i) = charge * sim.param().qmmm.unit_factor_charge;
      } // for atoms
    }
  } // while output
  output.close();
  
  // read output in fort.15
  output.open(this->output_gradient_file.c_str());
  if (!output.is_open()) {
    io::messages.add("Cannot open MNDO energy, gradient file (fort.15)",
            "MNDO_Worker", io::message::error);
    return 1;
  }
  while(1) {
    std::string line;
    std::getline(output, line);
    if (output.eof() || output.fail()) {
      break;
    }
    if /* energy */(line.find("ENERGY") != std::string::npos) {
      // next line
      std::getline(output, line);
      double energy;
      std::istringstream is(line);
      is >> energy;
      if (is.fail()) {
        io::messages.add("Failed to read energy result", "MNDO_Worker", 
                io::message::error);
        return 1;
      }
      storage.energy = energy * sim.param().qmmm.unit_factor_energy;
    } /* force QM atoms */ else if (line.find("CARTESIAN GRADIENT") != std::string::npos &&
                           line.find("OF MM ATOMS") == std::string::npos) {
      // read atoms
      for (std::set<topology::qm_atom_struct>::const_iterator
        it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
        std::getline(output, line);
        if (output.fail()) {
          io::messages.add("Failed to read gradient line of QM atom " + std::to_string(it->index+1),
            "MNDO_Worker", io::message::error);
          return 1;
        }
        
        std::istringstream is(line);
        int dummy;
        math::Vec gradient;
        // skip atom number, Z 
        is >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
        if (is.fail()) {
          io::messages.add("Failed to parse gradient line of QM atom " + std::to_string(it->index+1),
            "MNDO_Worker", io::message::error);
          return 1;
        }
        storage.force(it->index) = gradient * (-(sim.param().qmmm.unit_factor_energy) / 
                sim.param().qmmm.unit_factor_length);
      } // for QM atoms
        //read now the forces on the capping atom
        int dummy;
        math::Vec gradient;
        for (unsigned int i = 0; i < topo.qm_mm_pair().size(); i++) {
            std::istringstream is(line);
            std::getline(output, line);
            is >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
            if (output.fail()) {
                std::ostringstream msg;
                msg << "Failed to read charge line of " << (i + 1) << "th linkatom.";
                io::messages.add(msg.str(), "MNDO_Worker", io::message::error);
                return 1;
            }
            LA_storage.force(i) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
                                              sim.param().qmmm.unit_factor_length);
        } // for atoms
        

    } /* force MM atoms */ else if (line.find("CARTESIAN GRADIENT") != std::string::npos &&
                           line.find("OF MM ATOMS") != std::string::npos) {
      // read atoms
      for (unsigned int i = 0; i < mm_atoms.size(); ++i) {        
        std::getline(output, line);
        if (output.fail()) {
          io::messages.add("Failed to read gradient line of MM atom " + std::to_string(i+1),
            "MNDO_Worker", io::message::error);
          return 1;
        }
        
        std::istringstream is(line);
        int dummy;
        math::Vec gradient;
        // skip atom number, Z 
        is >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
        if (is.fail()) {
          io::messages.add("Failed to parse gradient line of MM atom " + std::to_string(i+1),
            "MNDO_Worker", io::message::error);
          return 1;
        }
        storage.force(mm_atoms[i].index) += gradient * (-(sim.param().qmmm.unit_factor_energy) /
                sim.param().qmmm.unit_factor_length);

        // get force on the charge-on-spring
        if (topo.is_polarisable(mm_atoms[i].index)) {
          std::getline(output, line);
          if (output.fail()) {
            io::messages.add("Failed to read gradient line of COS on MM atom " + std::to_string(i+1),
              "MNDO_Worker", io::message::error);
            return 1;
          }

          std::istringstream is(line);
          // skip atom number, Z 
          is >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
          if (is.fail()) {
            io::messages.add("Failed to parse gradient line of COS MM atom " + std::to_string(i+1),
              "MNDO_Worker", io::message::error);
            return 1;
          }
          storage.cos_force(mm_atoms[i].index) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
                  sim.param().qmmm.unit_factor_length);
        }
      } // for MM atoms 
    }
  }
  output.close();


/*
    for(std::vector<topology::three_body_term_struct>::const_iterator
            it=sim.param().qmmm.QM_CAP_topo->solute().angles().begin(),
            to=sim.param().qmmm.QM_CAP_topo->solute().angles().end();
            it!=to;++it){
        std::cout << "angle:" << it->i << "  " << it->j << "  "  << it->k <<
                  "  "  << it->type << std::endl;

    }
    for(std::vector<topology::two_body_term_struct>::const_iterator
                it=sim.param().qmmm.QM_CAP_topo->solute().bonds().begin(),
                to=sim.param().qmmm.QM_CAP_topo->solute().bonds().end();
        it!=to;++it){
        std::cout << "bond:" << it->i << "  " << it->j <<
                  "  "  << it->type << std::endl;

    }


    if (topo.qm_mm_pair().size()!=0) {

        sim.param().qmmm.qmmm_newdihed_cap->calculate_interactions(*sim.param().qmmm.QM_CAP_topo,
                                                                   qmmm_conf, sim);
        sim.param().qmmm.qmmm_qbond_cap->calculate_interactions(*sim.param().qmmm.QM_CAP_topo,
                                                                qmmm_conf,
                                                                sim);//this should be zero, since we set the C-H distance to the equilibrium value
        sim.param().qmmm.qmmm_ang_cap->calculate_interactions(*sim.param().qmmm.QM_CAP_topo,
                                                              qmmm_conf, sim);
        //for (std::set<topology::qm_atom_struct>::const_iterator
        //             it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {

        if (verbose) {
            for (unsigned int i = 0; i < 26; i++) {
                std::cout << "index :   " << i << "  " << qmmm_conf.current().force(i)[0] <<
                          "   " << qmmm_conf.current().force(i)[1] << "  "
                          << qmmm_conf.current().force(i)[2] << std::endl;
            }
        }

        for (std::set<topology::qm_atom_struct>::const_iterator
                     it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
            storage.force(it->index) -= qmmm_conf.current().force(it->index);
        }
    }
    /
*/
  if (verbose) {
      std::cout << "QM/MM Forces " << std::endl;
      //print out the gradient
      for (std::set<topology::qm_atom_struct>::const_iterator
                   it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
          std::cout << "index :   " << it->index << "  " << storage.force(it->index)[0] << "  "
                    << storage.force(it->index)[1] << "   " <<
                    storage.force(it->index)[2] << std::endl;
      }
  }
  return 0;
  
}

int interaction::MNDO_Worker::init(topology::Topology& topo, 
        configuration::Configuration& conf, simulation::Simulation& sim) {
  this->input_file = sim.param().qmmm.mndo.input_file;
  this->output_file = sim.param().qmmm.mndo.output_file;
  this->output_gradient_file = sim.param().qmmm.mndo.output_gradient_file;
  this->density_matrix_file = sim.param().qmmm.mndo.density_matrix_file;
  this->input_header=sim.param().qmmm.mndo.input_header;
  if (this->input_file.empty()) {
    using_tmp = true;
    if(util::create_tmpfile(this->input_file) < 1) {
        io::messages.add("Unable to create temporary input file: " + this->input_file,
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
        io::messages.add("Using temporary input file: " + this->input_file,
        "MNDO_Worker", io::message::notice);
    }
  }

  if (this->output_file.empty()) {
    if (util::create_tmpfile(this->output_file) < 1) {
      io::messages.add("Unable to create temporary output file: " + this->output_file,
        "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
        io::messages.add("Using temporary output file: " + this->output_file,
        "MNDO_Worker", io::message::notice);
    }
  }

  if (this->output_gradient_file.empty()) {
    if (util::create_tmpfile(this->output_gradient_file) < 1) {
      io::messages.add("Unable to create temporary output gradient file: " 
        + this->output_gradient_file, "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
        io::messages.add("Using temporary output gradient file: " + this->output_gradient_file,
        "MNDO_Worker", io::message::notice);
    }
  }
  else {
    // Create file specified by user
    std::ofstream of(this->output_gradient_file.c_str());
    if (!of.is_open()) {
      io::messages.add("Unable to create output gradient file: "
      + this->output_gradient_file, "MNDO_Worker", io::message::critical);
      return 1;
    }
    of.close();
  }

  if (this->density_matrix_file.empty()) {
    if (util::create_tmpfile(this->density_matrix_file) < 1) {
      io::messages.add("Unable to create temporary density matrix file: " 
        + this->density_matrix_file, "MNDO_Worker", io::message::critical);
      return 1;
    }
    else {
        io::messages.add("Using temporary density matrix file: " + this->density_matrix_file,
        "MNDO_Worker", io::message::notice);
    }
  }
  
#ifdef HAVE_SYMLINK
  // create fort.15 link for gradients
  if (symlink_err += symlink(this->output_gradient_file.c_str(), "fort.15") != 0) {
    io::messages.add("Unable to create symbolic link from fort.15 to "
      + this->output_gradient_file + " - check permissions.",
      "MNDO_Worker", io::message::critical);
    return 1;
  }
  // create fort.11 link for density matrix
  if (symlink_err += symlink(this->density_matrix_file.c_str(), "fort.11") != 0) {
    io::messages.add("Unable to create symbolic link from fort.11 to "
      + this->density_matrix_file + " - check permissions.",
      "MNDO_Worker", io::message::critical);
    return 1;
  }
#else
  {
    this->output_gradient_file = "fort.15";
    this->density_matrix_file = "fort.11";
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

interaction::MNDO_Worker::~MNDO_Worker() {
#ifdef HAVE_UNLINK
  // Remove symbolic links
  if (symlink_err == 0) {
    unlink("fort.11");
    unlink("fort.15");
  }
  // Delete temporary files
  if (using_tmp) {
    unlink(this->input_file.c_str());
    unlink(this->output_file.c_str());
    unlink(this->output_gradient_file.c_str());
    unlink(this->density_matrix_file.c_str());
  }
#endif
}
