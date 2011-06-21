/**
 * @file turbomole_worker.cc
 * interface to the Trubomole software package
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
#include "turbomole_worker.h"
#include "../../../util/system_call.h"
#include "../../../io/blockinput.h"
#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

#define MAXPATH 10240

interaction::Turbomole_Worker::~Turbomole_Worker() {
}

int interaction::Turbomole_Worker::init(topology::Topology& topo, 
        configuration::Configuration& conf, simulation::Simulation& sim) {
  return 0;
}

int impl(std::string func) {
  std::ostringstream msg;
  msg << "Function " << func << " is missing on this platform but required.";
  io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
  return 1;
}

void defortranize(std::string & str) {
  for(std::string::iterator it = str.begin(), to = str.end(); it != to; ++it) {
    if (*it == 'D') {
      *it = 'E';
    }
  }
}

int interaction::Turbomole_Worker::run_QM(topology::Topology& topo, 
        configuration::Configuration& conf, simulation::Simulation& sim, 
        const math::VArray& qm_pos, const std::vector<MM_Atom>& mm_atoms, 
        interaction::QM_Storage& storage) {
  char current_working_directory[MAXPATH];
#ifdef HAVE_GETCWD
  if (getcwd(current_working_directory, MAXPATH) == NULL) {
    io::messages.add("Cannot get current working directory. "
            "Increase MAXPATH in turbomole_worker.cc", 
            "Turbomole_Worker", io::message::error);
    return 1;
  }
#else
  return impl("getcwd");
#endif
  
#ifdef HAVE_CHDIR
  if (chdir(sim.param().qmmm.turbomole.working_directory.c_str()) != 0) {
    io::messages.add("Cannot change into Turbomole working directory",
            "Turbomole_Worker", io::message::error);
    return 1;
  }
#else
  return impl("chdir");
#endif

  std::ofstream input_coord(sim.param().qmmm.turbomole.input_coordinate_file.c_str());
  if (!input_coord.is_open()) {
    io::messages.add("Could not create input coordinate file for Turbomole at "
            " the location " + sim.param().qmmm.turbomole.input_coordinate_file,
            "Turbomole_Worker", io::message::critical);
    return 1;
  }
  
  input_coord << "$coord" << std::endl;
  // write QM zone
  double len_to_qm = 1.0 / sim.param().qmmm.unit_factor_length;
  unsigned int pi = 0;
  for (std::set<topology::qm_atom_struct>::const_iterator
    it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++pi) {
    
    input_coord.setf(std::ios::fixed, std::ios::floatfield);
    input_coord.precision(14);
    for (unsigned int i = 0; i < 3; ++i) {
      input_coord << std::setw(20) << std::right << qm_pos(pi)(i) * len_to_qm;
    }
    input_coord << "      " << std::left << 
            sim.param().qmmm.turbomole.elements[it->atomic_number];
    input_coord << std::endl;
  }
  input_coord << "$end" << std::endl;
  input_coord.close();
  
  
  std::ofstream input_mm_coord(sim.param().qmmm.turbomole.input_mm_coordinate_file.c_str());
  if (!input_mm_coord.is_open()) {
    io::messages.add("Could not create input MM coordinate file for Turbomole at "
            " the location " + sim.param().qmmm.turbomole.input_mm_coordinate_file,
            "Turbomole_Worker", io::message::critical);
    return 1;
  }
  
  input_mm_coord << "$point_charges" << std::endl;
  double chg_to_qm = 1.0 / sim.param().qmmm.unit_factor_charge;
  for (unsigned int i = 0; i < mm_atoms.size(); ++i) {
    if (topo.is_polarisable(mm_atoms[i].index)) {
      input_mm_coord.setf(std::ios::fixed, std::ios::floatfield);
      input_mm_coord.precision(14);
      for (unsigned int j = 0; j < 3; ++j) {
        input_mm_coord << std::setw(20) << std::right << mm_atoms[i].pos(j) * len_to_qm;
      }
      input_mm_coord.precision(8);
      input_mm_coord << std::setw(15) << std::right <<
              (mm_atoms[i].charge - topo.coscharge(mm_atoms[i].index)) * chg_to_qm;
      input_mm_coord << std::endl;
      input_mm_coord.setf(std::ios::fixed, std::ios::floatfield);
      input_mm_coord.precision(14);
      const math::Vec pos_cos(mm_atoms[i].pos + conf.current().posV(mm_atoms[i].index));
      for (unsigned int j = 0; j < 3; ++j) {
        input_mm_coord << std::setw(20) << std::right << pos_cos(j) * len_to_qm;
      }
      input_mm_coord.precision(8);
      input_mm_coord << std::setw(15) << std::right <<
              topo.coscharge(mm_atoms[i].index) * chg_to_qm;
      input_mm_coord << std::endl;
    } else {
      input_mm_coord.setf(std::ios::fixed, std::ios::floatfield);
      input_mm_coord.precision(14);
      for (unsigned int j = 0; j < 3; ++j) {
        input_mm_coord << std::setw(20) << std::right << mm_atoms[i].pos(j) * len_to_qm;
      }
      input_mm_coord.precision(8);
      input_mm_coord << std::setw(15) << std::right << mm_atoms[i].charge;
      input_mm_coord << std::endl;
    }
  }
  
  input_mm_coord << "$end" << std::endl;
  input_mm_coord.close();
  
  for(std::vector<std::string>::const_iterator it = sim.param().qmmm.turbomole.toolchain.begin(),
          to = sim.param().qmmm.turbomole.toolchain.end(); it != to; ++it) {
    std::string output_file(*it + ".out"), input_file("");
    
    if (*it == "define") {
      input_file = "define.inp";
    }
    
    if (util::system_call(sim.param().qmmm.turbomole.binary_directory + "/" + *it,
            input_file, output_file) != 0) {
      std::ostringstream msg;
      msg << "Turbomole program " << *it << " failed. See " 
              << sim.param().qmmm.turbomole.working_directory << "/" << output_file
              << " for details.";
      io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
      return 1;
    }
    // open the output file and check whether the program ended normally
    std::ifstream output(output_file.c_str());
    if (!output.is_open()) {
      std::ostringstream msg;
      msg << "Could not open output file " << output_file;
      io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
      return 1;
    }
    bool success = false;
    while(true) {
      std::string line;
      std::getline(output, line);
      if (output.fail()) {
        break;
      }
      if (line.find(" ended normally") != std::string::npos ||
          line.find(" : all done") != std::string::npos) {
        success = true;
        break;
      }
    }
    if (!success) {
      std::ostringstream msg;
      msg << "Turbomole program " << *it << " failed. See " 
              << sim.param().qmmm.turbomole.working_directory << "/" << output_file
              << " for details.";
      io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
      return 1;
    } 
  }
  // get energy
  std::ifstream output(sim.param().qmmm.turbomole.output_energy_file.c_str());
  if (!output.is_open()) {
    io::messages.add("Could not read energy output file at "
            " the location " + sim.param().qmmm.turbomole.output_energy_file,
            "Turbomole_Worker", io::message::critical);
    return 1;
  }
  
  while(1) {
    std::string line;
    std::getline(output, line);
    if (output.eof() || output.fail()) {
      break;
    }
    // get energy section
    if (line.find("$energy") != std::string::npos) {
      std::getline(output, line);
      std::istringstream is(line);
      int dummy;
      double energy;
      is >> dummy >> energy;
      if (output.fail() || is.fail()) {
        io::messages.add("Could not parse energy output file.",
            "Turbomole_Worker", io::message::error);
        return 1;
      }
      storage.energy = energy * sim.param().qmmm.unit_factor_energy;  
    }
  }
  output.close();
  output.clear();
  // get gradient
  output.open(sim.param().qmmm.turbomole.output_gradient_file.c_str());
  if (!output.is_open()) {
    io::messages.add("Could not read gradient output file at "
            " the location " + sim.param().qmmm.turbomole.output_gradient_file,
            "Turbomole_Worker", io::message::critical);
    return 1;
  }
  
  while(1) {
    std::string line;
    std::getline(output, line);
    if (output.eof() || output.fail()) {
      break;
    }
    // get gradient
    if (line.find("$grad") != std::string::npos) {
      // skip n+1 lines
      for(unsigned int i = 0; i < topo.qm_zone().size()+1; ++i) {
        std::getline(output, line);
      }
      
      for (std::set<topology::qm_atom_struct>::const_iterator
        it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
        std::getline(output, line);
        if (output.fail()) {
          std::ostringstream msg;
          msg << "Failed to read gradient line of QM atom " << it->index+1 << " atom.";
          io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
          return 1;
        }
        
        defortranize(line);
        std::istringstream is(line);  
        math::Vec gradient;
        // skip atom number, Z 
        is >> gradient(0) >> gradient(1) >> gradient(2);
        if (is.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse gradient line of QM atom " << it->index+1 << " atom.";
          io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
          return 1;
        }
        storage.force(it->index) = gradient * (-(sim.param().qmmm.unit_factor_energy) / 
                sim.param().qmmm.unit_factor_length);
      }
    }
  }
  output.close();
  output.clear();
  // get gradient
  output.open(sim.param().qmmm.turbomole.output_mm_gradient_file.c_str());
  if (!output.is_open()) {
    io::messages.add("Could not read MM gradient output file at "
            " the location " + sim.param().qmmm.turbomole.output_mm_gradient_file,
            "Turbomole_Worker", io::message::critical);
    return 1;
  }
  
  while(1) {
    std::string line;
    std::getline(output, line);
    if (output.eof() || output.fail()) {
      break;
    }
    // get gradient
    if (line.find("$point_charge_gradients") != std::string::npos) { 
      for (std::vector<interaction::MM_Atom>::const_iterator
        it = mm_atoms.begin(), to = mm_atoms.end(); it != to; ++it) {
        std::getline(output, line);
        if (output.fail()) {
          std::ostringstream msg;
          msg << "Failed to read gradient line of MM atom " << it->index+1 << ".";
          io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
          return 1;
        }
        
        defortranize(line);
        std::istringstream is(line);  
        math::Vec gradient;
        // skip atom number, Z 
        is >> gradient(0) >> gradient(1) >> gradient(2);
        if (is.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse gradient line of MM atom " << it->index+1 << ".";
          io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
          return 1;
        }
        storage.force(it->index) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
                sim.param().qmmm.unit_factor_length);
        
        // get the gradient on the charge-on-spring
        if (topo.is_polarisable(it->index)) {
          std::getline(output, line);
          if (output.fail()) {
            std::ostringstream msg;
            msg << "Failed to read gradient line of COS of MM atom " << it->index + 1 << ".";
            io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
            return 1;
          }

          defortranize(line);
          std::istringstream is(line);
          // skip atom number, Z 
          is >> gradient(0) >> gradient(1) >> gradient(2);
          if (is.fail()) {
            std::ostringstream msg;
            msg << "Failed to parse gradient line of COS of MM atom " << it->index + 1 << ".";
            io::messages.add(msg.str(), "Turbomole_Worker", io::message::error);
            return 1;
          }
          storage.cos_force(it->index) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
                  sim.param().qmmm.unit_factor_length);
        } 
      } // for mm atoms
    }
  }
  output.close();
  output.clear();
  
  // remove files again
#ifdef HAVE_REMOVE
  remove(sim.param().qmmm.turbomole.output_energy_file.c_str());
  remove(sim.param().qmmm.turbomole.output_gradient_file.c_str());
  remove(sim.param().qmmm.turbomole.output_mm_gradient_file.c_str());
#else
  return impl("remove");
#endif
  
#ifdef HAVE_CHDIR
  if (chdir(current_working_directory) != 0) {
    io::messages.add("Cannot change back from into Turbomole working directory.",
            "Turbomole_Worker", io::message::error);
    return 1;
  }
#else
  return impl("chdir");
#endif
          
  return 0;
}
