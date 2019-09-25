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
//#include "../../../math/periodicity.h"

// special interactions
#include "qm_storage.h"
#include "mm_atom.h"
#include "qm_worker.h"
#include "dftb_worker.h"
#include "../../../util/system_call.h"
#include "../../../io/blockinput.h"
#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special
#define MAXPATH 10240

int interaction::DFTB_Worker::run_QM(topology::Topology& topo,
        configuration::Configuration& conf,
        simulation::Simulation& sim,
        const math::VArray & qm_pos,
        const std::vector<interaction::MM_Atom> & mm_atoms,
        interaction::QM_Storage& storage,
        interaction::QM_Storage& LA_storage,
        const configuration::Configuration& qmmm_conf){

char current_working_directory[MAXPATH];
#ifdef HAVE_GETCWD
    if (getcwd(current_working_directory, MAXPATH) == NULL) {
        io::messages.add("Cannot get current working directory. "
                "Increase MAXPATH in turbomole_worker.cc",
                "dftb_Worker", io::message::error);
        return 1;
    }
#else
    return impl("getcwd");
#endif
    std::string tout="bla";
   
#ifdef HAVE_CHDIR
    if (chdir(sim.param().qmmm.dftb.working_directory.c_str()) != 0) {

        io::messages.add("Cannot change into DFTB working directory",
                "dftb_Worker", io::message::error);
        return 1;
    }
#else
    return impl("chdir");
#endif


  
  if (mm_atoms.empty()) {
    io::messages.add("Cannot deal with zero MM atoms yet.", "DFTB_Worker",
            io::message::error);
    return 1;
  }

  std::ofstream inp(input_file.c_str());
  
  if (!inp.is_open()) {
    io::messages.add("Could not create input file for DFTB at the location "
            + input_file, "DFTB_Worker", io::message::critical);
    return 1;
  }
  
  std::string header(sim.param().qmmm.dftb.input_header);
  // get the number of point charges
  unsigned int num_charge = mm_atoms.size();
  for(unsigned int i = 0; i < mm_atoms.size(); ++i) {
    if (topo.is_polarisable(mm_atoms[i].index)) {
      ++num_charge;
    }
  }
  std::ostringstream oss; oss << num_charge;
  header = io::replace_string(header, "@@NUM_ATOMS@@", oss.str());

  // write header
  inp << header << std::endl;
  inp.close();
  std::ofstream input_coord(sim.param().qmmm.dftb.geom_file.c_str());
  // write QM zone
  double len_to_qm = 1.0 / sim.param().qmmm.unit_factor_length3;
  double len_to_qm2 = 1.0 / sim.param().qmmm.unit_factor_mmlen ; //
  unsigned int pi = 0;
   //map Z to spec_number of dftb
  std::map<unsigned int,unsigned int >  map_Z;
  input_coord << topo.qm_zone().size() << " C" << std::endl;
   // input_coord << topo.qm_zone().size()+topo.qm_mm_pair().size() << " C" << std::endl;
  int spec=0;
  for (unsigned int i =0 ; i !=sim.param().qmmm.dftb.elements.size()  ;i++)
  {    
      spec=sim.param().qmmm.dftb.elements[i];
      input_coord << sim.param().qmmm.qmmm_at_to_num[spec]<< " ";
      map_Z[spec ] = i+1;
  } 
  input_coord << std::endl;

  for (std::set<topology::qm_atom_struct>::const_iterator
    it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++pi) {
    input_coord.setf(std::ios::fixed, std::ios::floatfield);
    input_coord.precision(8);
    input_coord << std::setw(3) << std::left << pi+1 << 
                      map_Z[it->atomic_number]<< "  ";
    for (unsigned int i = 0; i < 3; ++i) {
      input_coord << std::setw(9) << std::right << qm_pos(pi)(i) * len_to_qm2 << " ";
    }
    input_coord << std::endl;
  }
//write link atoms
    math::Vec posCap, dR;
    unsigned int m1, q1;
    const double rch = sim.param().qmmm.cap_dist; //Carbon-Hydrogen bond length
    for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
                 it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it, ++pi) {
        q1 = it->first;
        m1 = it->second;
        dR = conf.current().pos(m1) - conf.current().pos(q1);
        posCap = (rch / abs(dR)) * dR + conf.current().pos(q1);
        input_coord << std::setw(3) << std::left << pi+1 <<
                    map_Z[1]<< "  ";
        for (unsigned int i = 0; i < 3; ++i) {
            input_coord << std::setw(9) << std::right << posCap(i) * len_to_qm2 << " ";
        }
        input_coord << std::endl;
    }
    input_coord.close();
    // write point charges

  // write point charges
  std::ofstream input_mm(sim.param().qmmm.dftb.output_charg_file.c_str());
  double chg_to_qm = 1.0 / sim.param().qmmm.unit_factor_charge3;
 // std::cout << "number of mm atoms: " << mm_atoms.size() << std::endl; 
  for (unsigned int i = 0; i < mm_atoms.size(); ++i) {
      input_mm.setf(std::ios::fixed, std::ios::floatfield);
      input_mm.precision(8);
      for (unsigned int j = 0; j < 3; ++j) {
        input_mm << std::setw(14) << std::right << mm_atoms[i].pos(j) * len_to_qm;
      }
      input_mm << std::setw(14) << std::right << mm_atoms[i].charge * chg_to_qm << std::endl;
    }
  //}
  input_mm.close();
  
  // starting DFTB

  tout = "dftb.out";
//  std::ostringstream command_to_launch;
//  command_to_launch << sim.param().qmmm.dftb.binary;
// 
//   command_to_launch << " < " << input_file;
// command_to_launch << " 1> " << output_file << " 2>&1 ";
// std::cout << "starting dftb  "<<  command_to_launch.str().c_str() << std::endl;
//  system(command_to_launch.str().c_str());
//  
  if (util::system_call(sim.param().qmmm.dftb.binary ,
            input_file, output_file) != 0) {
      std::ostringstream msg;
      msg << "dftb program " << " failed. See " 
              << sim.param().qmmm.dftb.binary   << " for details.";
      io::messages.add(msg.str(), "dftb_Worker", io::message::error);
      return 1;
  }
  
  // read output
   std::ifstream output(sim.param().qmmm.dftb.output_file.c_str());
  if (!output.is_open()) {
    io::messages.add("Cannot open DFTB output file", "DFTB_Worker", 
            io::message::error);
    return 1;
  }
//  while(1) {
//    std::string line;
//    std::getline(output, line);
//    if (output.eof() || output.fail()) {
//      break;
//    }
//    // get charges
////    if (line.find("NET ATOMIC CHARGES") != std::string::npos) {
////      // skip 3 lines
////      for (unsigned int i = 0; i < 3; ++i) std::getline(output, line);
////      // read atoms
////      for (unsigned int i = 0; i < topo.qm_zone().size(); ++i) {
////        std::getline(output, line);
////        if (output.fail()) {
////          std::ostringstream msg;
////          msg << "Failed to read charge line of " << (i+1) << "th atom.";
////          io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
////          return 1;
////        }
////        
////        std::istringstream is(line);
////        std::string dummy; int idummy;
////        double charge;
////        is >> dummy >> idummy >> charge;
////        if (is.fail()) {
////          std::ostringstream msg;
////          msg << "Failed to parse charge line of " << (i+1) << "th atom.";
////          io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
////          return 1;
////        }
////        storage.charge(i) = charge * sim.param().qmmm.unit_factor_charge;
////      } // for atoms
////    }
//  } // while output
//  output.close();
  // get gradient
//  std::ifstream of(sim.param().qmmm.dftb.output_file.c_str());
//  of.open(sim.param().qmmm.dftb.output_gradient_file.c_str());
//  
//  of.close();
  
  // read output in fort.15

  //output.open(output_file.c_str());
  //output.clear();

  while(1) {
    std::string line;
    std::string dummy;
    std::getline(output, line);
    if (output.eof() || output.fail()) {
      break;
    }
    if /* energy */(line.find("Total energy:") != std::string::npos) {
      // next line
      double energy;
      int  niter;
      std::istringstream is(line);
      is >> dummy >> dummy >> energy;
      if (is.fail()) {
        io::messages.add("Failed to read energy result", "DFTB_Worker", 
                io::message::error);
        return 1;
      }
      storage.energy = energy * sim.param().qmmm.unit_factor_energy3;
    } /* force QM atoms */
    else if (line.find("Total Forces") != std::string::npos ) 
    {
       // std::cout << " foundQM forces" << std::endl;
      // read atoms
      for (std::set<topology::qm_atom_struct>::const_iterator
        it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
        std::getline(output, line);
        if (output.fail()) {
          std::ostringstream msg;
          msg << "Failed to read gradient line of QM atom " << it->index+1 << " atom.";
          io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
          return 1;
        }
        
        std::istringstream is(line);
        int dummy;
        math::Vec gradient;
        is >> gradient(0) >> gradient(1) >> gradient(2);
      //  std::cout << gradient(0) <<"   " << gradient(1) << " "  << gradient(2) << std::endl;
        if (is.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse gradient line of QM atom " << it->index+1 << " atom.";
          io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
          return 1;
        }
        storage.force(it->index) += gradient * ((sim.param().qmmm.unit_factor_energy3) /
                sim.param().qmmm.unit_factor_length3);
        
      } // for QM atoms 
        //read now the forces on the capping atom
        math::Vec gradient;
        for (unsigned int i = 0; i < topo.qm_mm_pair().size(); i++) {
            std::istringstream is(line);
            std::getline(output, line);
            is >> gradient(0) >> gradient(1) >> gradient(2);
            if (output.fail()) {
                std::ostringstream msg;
                msg << "Failed to read charge line of " << (i + 1) << "th linkatom.";
                io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
                return 1;
            }
            LA_storage.force(i) = gradient * ((sim.param().qmmm.unit_factor_energy3) /
                                              sim.param().qmmm.unit_factor_length3);
        }
        
    } /* force MM atoms */ else if (line.find("Forces on external charges") != std::string::npos ) {
      // read atoms
    //    std::cout << " found MM forces" << std::endl;
      for (unsigned int i = 0; i < mm_atoms.size(); ++i) {        
        std::getline(output, line);
        if (output.fail()) {
          std::ostringstream msg;
          msg << "Failed to read gradient line of MM atom " << i+1 << ".";
          io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
          return 1;
        }
        
        std::istringstream is(line);
        int dummy;
        math::Vec gradient;
        // skip atom number, Z 
        is >>  gradient(0) >> gradient(1) >> gradient(2);
        //std::cout << gradient(0) <<"   " << gradient(1) << " "  << gradient(2) << std::endl;
        if (is.fail()) {
          std::ostringstream msg;
          msg << "Failed to parse gradient line of MM atom " << i+1 << ".";
          io::messages.add(msg.str(), "DFTB_Worker", io::message::error);
          return 1;
        }
        storage.force(mm_atoms[i].index) = gradient * ((sim.param().qmmm.unit_factor_energy3) / 
                sim.param().qmmm.unit_factor_length3);

        // get force on the charge-on-spring
//        if (topo.is_polarisable(mm_atoms[i].index)) {
//          std::getline(output, line);
//          if (output.fail()) {
//            std::ostringstream msg;
//            msg << "Failed to read gradient line of COS on MM atom " << i + 1 << ".";
//            io::messages.add(msg.str(), "MNDO_Worker", io::message::error);
//            return 1;
//          }
//
//          std::istringstream is(line);
//          // skip atom number, Z 
//          is >> dummy >> dummy >> gradient(0) >> gradient(1) >> gradient(2);
//          if (is.fail()) {
//            std::ostringstream msg;
//            msg << "Failed to parse gradient line of COS MM atom " << i + 1 << ".";
//            io::messages.add(msg.str(), "MNDO_Worker", io::message::error);
//            return 1;
//          }
//          storage.cos_force(mm_atoms[i].index) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
//                  sim.param().qmmm.unit_factor_length);
//        }
      } // for MM atoms 
    }
    
  }
  output.close();
  #ifdef HAVE_REMOVE
  //remove(sim.param().qmmm.dftb.output_file.c_str());
#else
  return impl("remove");
#endif
#ifdef HAVE_CHDIR
  if (chdir(current_working_directory) != 0) {
    io::messages.add("Cannot change back from into DFTB working directory.",
            "dftb worker", io::message::error);
    return 1;
  }
#else
  return impl("chdir");
#endif
          
  return 0;
  
}

int interaction::DFTB_Worker::init(topology::Topology& topo, 
        configuration::Configuration& conf, simulation::Simulation& sim) { 
      input_file = sim.param().qmmm.dftb.input_file;
      output_file = sim.param().qmmm.dftb.output_file;
      output_charg_file = sim.param().qmmm.dftb.output_charg_file;
      geom_file = sim.param().qmmm.dftb.geom_file;
      // Make a global periodic table with atomic number, type, mass ...
      sim.param().qmmm.qmmm_at_to_num[1]="H";
      sim.param().qmmm.qmmm_at_to_num[6]="C";
      sim.param().qmmm.qmmm_at_to_num[7]="N";
      sim.param().qmmm.qmmm_at_to_num[8]="O";
      sim.param().qmmm.qmmm_at_to_num[9]="F";
      sim.param().qmmm.qmmm_at_to_num[15]="P";
      sim.param().qmmm.qmmm_at_to_num[16]="S";
      sim.param().qmmm.qmmm_at_to_num[17]="Cl";
      sim.param().qmmm.qmmm_at_to_num[35]="Br";
      sim.param().qmmm.qmmm_at_to_num[53]="I";
      if (output_file.empty()) {
          std::cout << "OUTPUTFILE EMPTY?!?!?!?" << std::endl;
#ifdef HAVE_TMPNAM
    char tmp[TMP_MAX];
    if (tmpnam(tmp) == NULL) {
      io::messages.add("Could not get temporary file",
              "dftb_Worker", io::message::error);
      return 1;
    }
    output_file = std::string(tmp);
#else
    io::messages.add("Temporary files are not supported on this platform. "
            "Please provide the file names manually using the MNDOFILES "
            "block in the QM/MM specification file",
            "MNDO_Worker", io::message::critical);
    return 1;
#endif
  }

    return 0;
}

interaction::DFTB_Worker::~DFTB_Worker() {
this->del_qmID();
}
