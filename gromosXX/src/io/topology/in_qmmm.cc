
/**
 * @file in_posres.cc
 * implements methods of In_Posresspec and In_Posres.
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../interaction/interaction_types.h"
#include "../../configuration/configuration.h"

#include "../../io/instream.h"
#include "../../io/blockinput.h"
#include "../../io/configuration/in_configuration.h"
#include "../../io/configuration/out_configuration.h"

#include "in_qmmm.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

/**
 * @section qmzone QMZONE block
 * The QMZONE block specifies the atoms which are treated quantum mechanically
 *
 * The block is read from the QM/MM specification file
 * (\@qmmm).
 *
 * @verbatim
QMZONE
# QMI:  index of the QM atom
# QMZ:  atomic number of the QM atom
# QMLI: 0,1 atom is a link atom
#
# Warning: the first 17 characters are ignored!
# RESIDUE   ATOM     QMI   QMZ   QMLI
    1 H2O   OW         1     8      0
    1 H2O   HW1        2     1      0
    1 H2O   HW2        3     1      0
END
@endverbatim
 *
 * @section qmunbit QMUNIT block
 * The QMUNIT block specifies conversion factors for units
 *
 * The block is read from the QM/MM specification file
 * (\@qmmm).
 *
 * @verbatim
QMUNIT
# QMULEN: Conversion factor to convert the QM length unit to the GROMOS one
# QMUENE: Conversion factor to convert the QM energy unit to the GROMOS one
# QMUCHR: Conversion factor to convert the QM charge unit to the GROMOS one
#
# QMULEN    QMUENE    QMUCHR
     0.1     4.184       1.0
END
@endverbatim
 *
 * @section MNDO blocks for the MNDO worker
 * The MNDOFILES blocks specifies where MNDO writes the input and output files
 *
 * Temporary files are used if this block is omitted.
 *
 * @verbatim
MNDOFILES
/path/to/mndo/binary
/path/to/mndo.in
/path/to/mndo.out
/path/to/mndo_gradient.out
END
@endverbatim
 *
 * The MNDOHEADER block specifies the header part of the MNDO input file. Variables
 * are allowed. Implemented are
 * - NUM_ATOMS: the number of atoms
 * - NUM_LINK: Number of link atoms
 * 
@verbatim
MNDOHEADER
kharge=0 iop=-8 +
kitscf=200 +
ktrial=11 +
igeom=1 iform=1 nsav15=4 ipubo=1 jop=-2 +
mminp=2 mmcoup=2 mmlink=2 nlink=@@NUM_LINK@@ numatm=@@NUM_ATOMS@@
Title line
END
@endverbatim
 * 
 * @section Turbomole blocks for the Turbomole worker
 * 
 * The TURBOMOLEFILES blocks specifies where Turbomole writes the input and output files.
 * The first line contains the directory which contains the Turbomole binaries.
 * The second line is the turbomole working directory containing the control file.
 * In this control file the relative paths for the coordinate, point charges coordinates,
 * energy, gradients and point charges gradients are defined. The next lines contains
 * these file names for GROMOS.
 *
 * @verbatim
TURBOMOLEFILES
/path/to/turbomole/binary/directory
/path/to/working/directory/containing/control/file
coordinate.in
mm_coordinate.in
energy.out
gradient.out
mm_gradient.out
END
@endverbatim
 * 
 * The TURBOMOLETOOLCHAIN block specifies the Turbomole programs that are executed.
 * Each line contains one program that is called. By default, it is assumed that
 * the control file is static, i.e. that the number of QM atoms cannot change. To 
 * modify the control file during the simulation the TURBOMOLE program define is 
 * to be used. The input for this program has to be given as a file named "define.inp" 
 * which is in the same directory as the "control" file. 
 * Note: You still have to provide an initial control file as at the time the 
 * program define is only cannot include the $point_charges directives. 
 * 
 * @verbatim
 TURBOMOLETOOLCHAIN
 ridft
 rdgrad
 END
 @endverbatim
 *
 * The TURBOMOLEELEMENTS block specifies the element name used in Turbomole.
 * It is determined by the atomic number given in the QMZONE block.
 * 
@verbatim
TURBOMOLEELEMENTS
1 h
6 c
7 n
8 o
END
@endverbatim
 */
void
io::In_QMMM::read(topology::Topology& topo,
        simulation::Simulation & sim,
        std::ostream & os) {

  DEBUG(7, "reading in a QM/MM specification file");
  std::vector<std::string> buffer;
  { // QMZONE
    buffer = m_block["QMZONE"];
    DEBUG(10, "QMZONE block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no QMZONE block in QM/MM specification file",
              "In_QMMM", io::message::error);
      return;
    }

    std::vector<std::string>::const_iterator it = buffer.begin() + 1,
            to = buffer.end() - 1;

    DEBUG(10, "reading in QMZONE data");

    unsigned int i, nr, z;
    int link;
    for (i = 0; it != to; ++i, ++it) {

      DEBUG(11, "\tnr " << i);

      std::string line(*it);
      if (line.length() < 17) {
        io::messages.add("line too short in QMZONE block", "In_QMMM",
                io::message::error);
      }

      // the first 17 chars are ignored
      line.erase(line.begin(), line.begin() + 17);

      _lineStream.clear();
      _lineStream.str(line);

      _lineStream >> nr >> z >> link;

      DEBUG(11, "\t" << nr << "\t" << z);

      if (_lineStream.fail()) {
        io::messages.add("bad line in QMZONE block",
                "In_QMMM", io::message::error);
        return;
      }

      if (nr < 1 || nr > topo.num_atoms()) {
        io::messages.add("QMZONE block: atom out of range",
                "In_QMMM", io::message::error);
        return;
      }

      if (link < 0 || link > 1) {
        io::messages.add("QMZONE block: QMLI has to be 0 or 1",
                "In_QMMM", io::message::error);
        return;
      }

      --nr;
      topo.qm_zone().insert(topology::qm_atom_struct(nr, z, link, topo.charge(nr)));
    }
  } // QMZONE
  { // QMUNIT
    buffer = m_block["QMUNIT"];
    DEBUG(10, "QMUNIT block : " << buffer.size());

    if (!buffer.size()) {
      io::messages.add("no QMUNIT block in QM/MM specification file",
              "In_QMMM", io::message::error);
      return;
    }
    std::string s;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    _lineStream >> sim.param().qmmm.unit_factor_length
            >> sim.param().qmmm.unit_factor_energy
            >> sim.param().qmmm.unit_factor_charge;

    if (_lineStream.fail()) {
      io::messages.add("bad line in QMUNIT block.",
              "In_QMMM", io::message::error);
      return;
    }
  } // QMUNIT

  // check for MNDO specific data
  if (sim.param().qmmm.software == simulation::qmmm_software_mndo) {
    { // MNDOFILES
      buffer = m_block["MNDOFILES"];
      DEBUG(10, "MNDOFILES block : " << buffer.size());

      if (!buffer.size()) {
        io::messages.add("Using temporary files for MNDO input/output and assuming that the binary is in the PATH",
                "In_QMMM", io::message::notice);
        sim.param().qmmm.mndo.binary = "mndo";
      } else {
        if (buffer.size() != 6) {
          io::messages.add("MNDOFILES block corrupt. Provide 4 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.mndo.binary = buffer[1];
        sim.param().qmmm.mndo.input_file = buffer[2];
        sim.param().qmmm.mndo.output_file = buffer[3];
        sim.param().qmmm.mndo.output_gradient_file = buffer[4];
      }
    } // MNDOFILES
    { // MNDOHEADER
      buffer = m_block["MNDOHEADER"];
      DEBUG(10, "MNDOHEADER block : " << buffer.size());

      if (!buffer.size()) {
        io::messages.add("no MNDOHEADER block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }

      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.mndo.input_header);
    } // MNDOHEADER
  } else if (sim.param().qmmm.software == simulation::qmmm_software_turbomole) {
    { // TURBOMOLEELEMENTS
      buffer = m_block["TURBOMOLEELEMENTS"];
      DEBUG(10, "TURBOMOLEELEMENTS block : " << buffer.size());

      if (!buffer.size()) {
        io::messages.add("no TURBOMOLEELEMENTS block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      std::string s;
      _lineStream.clear();
      _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
      unsigned int Z;
      std::string element;
      while(_lineStream >> Z >> element) 
        sim.param().qmmm.turbomole.elements[Z] = element;
      
      // check whether all elements have been provided
      bool error = false;
      for (std::set<topology::qm_atom_struct>::const_iterator
        it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
        if (sim.param().qmmm.turbomole.elements.find(it->atomic_number) == 
                sim.param().qmmm.turbomole.elements.end()) {
          std::ostringstream msg;
          msg << "Please provide element name for atomic number " << it->atomic_number
                  << " in TURBOMOLEELEMENTS block.";
          io::messages.add(msg.str(), "In_QMMM", io::message::error);
        }
      }
      if (error) 
        return;
    } // TURBOMOLEELEMENTS
    { // TURBOMOLEFILES
      buffer = m_block["TURBOMOLEFILES"];
      DEBUG(10, "TURBOMOLEFILES block : " << buffer.size());

      if (!buffer.size()) {
        io::messages.add("TURBOMOLEFILES block missing",
                "In_QMMM", io::message::error);
        return;
      } else {
        if (buffer.size() != 9) {
          io::messages.add("TURBOMOLEFILES block corrupt. Provide 7 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.turbomole.binary_directory = buffer[1];
        sim.param().qmmm.turbomole.working_directory = buffer[2];
        sim.param().qmmm.turbomole.input_coordinate_file = buffer[3];
        sim.param().qmmm.turbomole.input_mm_coordinate_file = buffer[4];
        sim.param().qmmm.turbomole.output_energy_file = buffer[5];
        sim.param().qmmm.turbomole.output_gradient_file = buffer[6];
        sim.param().qmmm.turbomole.output_mm_gradient_file = buffer[7];
      }
    } // TURBOMOLEFILES
    { // TURBOMOLETOOLCHAIN
      buffer = m_block["TURBOMOLETOOLCHAIN"];
      DEBUG(10, "TURBOMOLETOOLCHAIN block : " << buffer.size());

      if (!buffer.size()) {
        io::messages.add("TURBOMOLETOOLCHAIN block missing",
                "In_QMMM", io::message::error);
        return;
      } 
      for(unsigned int i = 1; i < buffer.size()-1; ++i) {
        _lineStream.clear();
        _lineStream.str(buffer[i]);
        std::string tool;
        _lineStream >> tool;
        if (_lineStream.fail()) {
          io::messages.add("bad line in TURBOMOLETOOLCHAIN block",
                "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.turbomole.toolchain.push_back(tool);
      }
    } // TURBOMOLETOOLCHAIN  
  }
}

