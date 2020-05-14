
/**
 * @file in_qmmm.cc
 * implements qm/mm methods.
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

#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

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
 * @section bufferzone BUFFERZONE block
 * The BUFFERZONE block specifies the atoms which are treated in both quantum and
 * classical way. They are added to QMZONE atoms to form the full QM zone. Energies
 * and interactions are then calculated as a contribution to forcefield deltaE and
 * deltaF
 * deltaE = E(fullQM) - E(buffer)
 * deltaF = F(fullQM) - F(buffer)
 * 
 * 
 * 
 *
 * The block is read from the QM/MM specification file
 * (\@qmmm).
 *
 * @verbatim
BUFFERZONE
# BUFCUT: cutoff of the adaptive buffer zone (default = 0.0, no adaptive buffer zone)
          If the specified atoms occur within BUFCUT of QM atom, they are considered
          as buffer atoms. Otherwise they are considered solely as MM atoms. BUFCUT = 0.0
          means that they are always considered as buffer atoms.
# BUFCUT
    1.4
# NETCH:
# SPINM:
# QMI:   index of the QM atom
# QMZ:   atomic number of the QM atom
# QMLI:  0,1 atom is a link atom
#
# Warning: the first 17 characters are ignored!
# RESIDUE   ATOM     QMI   QMZ   QMLI
    1 H2O   OW         1     8      0
    1 H2O   HW1        2     1      0
    1 H2O   HW2        3     1      0
END
@endverbatim
 *
 * @section qmunit QMUNIT block
 * The QMUNIT block specifies conversion factors for units
 *
 * The block is read from the QM/MM specification file
 * (\@qmmm).
 *
 * @verbatim
QMUNIT
# QMULEN: Conversion factor to convert the QM length unit to the GROMOS one
# QMUENE: Conversion factor to convert the QM energy unit to the GROMOS one
# QMUFOR: Conversion factor to convert the QM force unit to the GROMOS one
# QMUCHR: Conversion factor to convert the QM charge unit to the GROMOS one
#
# QMULEN    QMUENE    QMUFOR    QMUCHR
     0.1     4.184     41.84       1.0
END
@endverbatim
 * @section cap_length CAPLEN block
 * The CAPLEN block specifies the distance between QM linking atom and its capping atom
 *
 * The block is read from the QM/MM specification file
 * (\@qmmm).
 *
 * @verbatim
CAPLEN
# CAPLEN:  distance between QM atom and capping atom in nm (default is 0.109 nm)
#
# CAPLEN
    0.109
END
@endverbatim
 *
 * @section MNDOBINARY block for the MNDO worker
 * The MNDOBINARY block specifies path to MNDO binary
 *
 * This block is optional. If unspecified, mndo command from PATH environment variable
 * is used.
 *
 * @verbatim
MNDOBINARY
/path/to/mndo/binary
END
@endverbatim
 * @section MNDOFILES block for the MNDO worker
 * The MNDOFILES block specifies input and output files to exchange data with MNDO
 *
 * This block is optional. If unspecified, temporary files are created using TMPDIR
 * environment variable. User-specified files are not deleted after use.
 *
 * @verbatim
MNDOFILES
/path/to/mndo.in
/path/to/mndo.out
/path/to/mndo_gradient.out   ##fort.15 file
/path/to/mndo_density.bin    ##fort.11 file
END
@endverbatim
 *
 * The MNDOHEADER block specifies the header part of the MNDO input file. Variables
 * are allowed. Implemented are
 * - CHARGE: net charge of the QM zone
 * - SPINM: spin multiplicity of the QM zone (be aware, that MNDO uses 0 for ground-state)
 * - NUM_CHARGES: the number of MM atoms
 * - NUM_LINK: Number of link atoms
 * 
@verbatim
MNDOHEADER
kharge=@@CHARGE@@ imult=@@SPINM@@ +
iop=-8 +
kitscf=2000 +
idiis=1 +
ktrial=11 +
igeom=1 iform=1 nsav15=4 ipubo=1 jop=-2 +
mminp=2 mmcoup=2 mmlink=1 nlink=@@NUM_LINK@@ numatm=@@NUM_CHARGES@@
title line
END
@endverbatim
 *
 * @section GAUSSIANBINARY block for the Gaussian worker
 * The GAUSSIANBINARY block specifies path to GAUSSIAN binary
 *
 * This block is optional. If unspecified, g16 command from PATH environment variable
 * is used.
 *
 * @verbatim
GAUSSIANBINARY
/path/to/gaussian/binary
END
@endverbatim
 * @section GAUSSIANFILES block for the Gaussian worker
 * The GAUSSIANFILES block specifies input and output files to exchange data with Gaussian
 *
 * This block is optional. If unspecified, temporary files are created using TMPDIR
 * environment variable. User-specified files are not deleted after use.
 *
 * @verbatim
GAUSSIANFILES
/path/to/gaussian.gjf
/path/to/gaussian.out
END
@endverbatim
 *
 * The GAUSSIANHEADER block specifies the header part of the Gaussian input file.
 * 
@verbatim
GAUSSIANHEADER
%nproc=8
%mem=2GB
%NoSave
%chk=tmp
END
@endverbatim
 *
 * The GAUSSIANROUTE block specifies the route section of the Gaussian input file.
 * hashsign (#) should be omitted. It is beneficial to generate an initial checkpoint file
 * and reuse it in subsequent steps with guess=read option.
 * 
@verbatim
GAUSSIANROUTE
N hf/STO-3G nosymm pop(mk) charge(angstroms) force charge(angstroms) prop=(field,read)
END
@endverbatim
 *
 * The GAUSSIANCHSM block specifies the net charge and the spin multiplicity of the system.
 * Variables are allowed. Implemented are
 * - CHARGE: net charge of the QM zone
 * - SPINM: spin multiplicity of the QM zone
@verbatim
GAUSSIANCHSM
@@CHARGE@@ @@SPINM@@
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
 * @verbatim
TURBOMOLEELEMENTS
1 h
6 c
7 n
8 o
END
@endverbatim
 *
 * The NNMODEL block specifies the PyTorch NN model to use.
 * 
 * @verbatim
NNMODEL
/path/to/schnetpack/model
END
@endverbatim
 *
 * The NNDEVICE block specifies the device to use.
 * Allowed are
 * - cpu: Use CPU, a model trained on CUDA will be mapped to CPU
 * - cuda: Use CUDA
 * - auto: Use CUDA, if available, otherwise use CPU (this is default)
 * 
@verbatim
NNDEVICE
auto
END
@endverbatim
 */
void
io::In_QMMM::read(topology::Topology& topo,
        simulation::Simulation& sim,
        std::ostream & os) {
  io::messages.add("Reading QM/MM specification file",
                   "In_QMMM", io::message::notice);
  //std::string blockname = "QMZONE";
  this->read_zone(topo, sim, "QMZONE");
  //blockname = "BUFFERZONE";
  if (!sim.param().qmmm.qm_constraint) {
    DEBUG(15, "Removing constraints");
    this->remove_constraints(topo);
    io::messages.add("QM-QM distance constraints removed",
        "In_QMMM", io::message::notice);
  }
  this->read_zone(topo, sim, "BUFFERZONE");
  
  std::vector<std::string> buffer;

  const simulation::qm_software_enum sw = sim.param().qmmm.software;
  /**
   * MNDO
   */
  if (sw == simulation::qm_mndo) {
    this->read_units(sim, &sim.param().qmmm.mndo);
    { // MNDOBINARY

      DEBUG(15, "Reading MNDOBINARY");
      buffer = m_block["MNDOBINARY"];

      if (!buffer.size()) {
        io::messages.add("Assuming that the mndo binary is in the PATH",
                "In_QMMM", io::message::notice);
        sim.param().qmmm.mndo.binary = "mndo";
      } else {
        if (buffer.size() != 3) {
          io::messages.add("MNDOBINARY block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.mndo.binary = buffer[1];
      }
    } // MNDOBINARY
    { // MNDOFILES

      DEBUG(15, "Reading MNDOFILES");
      buffer = m_block["MNDOFILES"];

      if (!buffer.size()) {
        io::messages.add("Using temporary files for MNDO input/output",
                "In_QMMM", io::message::notice);
      } else {
        if (buffer.size() != 6) {
          io::messages.add("MNDOFILES block corrupt. Provide 4 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.mndo.input_file = buffer[1];
        sim.param().qmmm.mndo.output_file = buffer[2];
        sim.param().qmmm.mndo.output_gradient_file = buffer[3];
        sim.param().qmmm.mndo.density_matrix_file = buffer[4];
      }
    } // MNDOFILES
    { // MNDOHEADER
      buffer = m_block["MNDOHEADER"];

      if (!buffer.size()) {
        io::messages.add("no MNDOHEADER block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.mndo.input_header);
    } // MNDOHEADER
  }

  /**
   * Turbomole
   */
  else if (sw == simulation::qm_turbomole) {
    this->read_units(sim, &sim.param().qmmm.turbomole);
    this->read_elements(topo, &sim.param().qmmm.turbomole);

    { // TURBOMOLEFILES
      buffer = m_block["TURBOMOLEFILES"];

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

  /**
   * DFTB
   */
  else if (sw == simulation::qm_dftb) {
    this->read_units(sim, &sim.param().qmmm.dftb);
    this->read_elements(topo, &sim.param().qmmm.turbomole);
    { // DFTBFILES
      buffer = m_block["DFTBFILES"];
      if (!buffer.size()) {
        io::messages.add("Using temporary files for DFTB input/output and assuming that the binary is in the PATH",
                         "In_QMMM", io::message::notice);
        sim.param().qmmm.dftb.binary = "dftb";
      } else {
        if (buffer.size() != 8) {
          io::messages.add("DFTB block corrupt. Provide 6 lines.",
                           "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.dftb.binary = buffer[1];
        sim.param().qmmm.dftb.working_directory = buffer[2];
        sim.param().qmmm.dftb.input_file = buffer[3];
        sim.param().qmmm.dftb.output_file = buffer[4];
        sim.param().qmmm.dftb.output_charge_file = buffer[5];
        sim.param().qmmm.dftb.geom_file = buffer[6];
      }
    } // DFTBFILES

    { // DFTBHEADER
      buffer = m_block["DFTBHEADER"];
      if (!buffer.size()) {
        io::messages.add("no DFTBHEADER block in QM/MM specification file",
                         "In_QMMM", io::message::error);
        return;
      }

      concatenate(buffer.begin() + 1, buffer.end() - 1,
                  sim.param().qmmm.dftb.input_header);
    } // DFTBHEADER
  }

  /**
   * MOPAC
   */
  else if (sw == simulation::qm_mopac) {
    this->read_units(sim, &sim.param().qmmm.mopac);
    {
      //MOPACFILES
      buffer = m_block["MOPACFILES"];

      if (!buffer.size()) {
        io::messages.add("Using temporary files for MOPAC input/output and assuming that the binary is in the PATH",
                         "In_QMMM", io::message::notice);
        sim.param().qmmm.mopac.binary = "mopac";
      } else {
        if (buffer.size() != 7) {
          io::messages.add("MOPAC block corrupt. Provide 4 lines.",
                           "In_QMMM", io::message::error);
          return;
        }
      sim.param().qmmm.mopac.binary = buffer[1];
      sim.param().qmmm.mopac.input_file = buffer[2];
      sim.param().qmmm.mopac.output_file = buffer[3];
      sim.param().qmmm.mopac.output_gradient_file = buffer[4];
      sim.param().qmmm.mopac.molin_file = buffer[5];
      }
    } // MOPACFILES
    { // MOPACHEADER
      buffer = m_block["MOPACHEADER"];
      if (!buffer.size()) {
        io::messages.add("no MOPACHEADER block in QM/MM specification file",
                         "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
                  sim.param().qmmm.mopac.input_header);
    } // MOPACHEADER
  }

  /**
   * Gaussian
   */
  else if (sw == simulation::qm_gaussian) {
    this->read_units(sim, &sim.param().qmmm.gaussian);
    { // GAUSSIANBINARY

      DEBUG(15, "Reading GAUSSIANBINARY");
      buffer = m_block["GAUSSIANBINARY"];

      if (!buffer.size()) {
        io::messages.add("Assuming that the g16 binary is in the PATH",
                "In_QMMM", io::message::notice);
        sim.param().qmmm.gaussian.binary = "g16";
      } else {
        if (buffer.size() != 3) {
          io::messages.add("GAUSSIANBINARY block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.gaussian.binary = buffer[1];
      }
    } // GAUSSIANBINARY
    { // GAUSSIANFILES

      DEBUG(15, "Reading GAUSSIANFILES");
      buffer = m_block["GAUSSIANFILES"];

      if (!buffer.size()) {
        io::messages.add("Using temporary files for Gaussian input/output",
                "In_QMMM", io::message::notice);
      } else {
        if (buffer.size() != 4) {
          io::messages.add("GAUSSIANFILES block corrupt. Provide 2 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.gaussian.input_file = buffer[1];
        sim.param().qmmm.gaussian.output_file = buffer[2];
      }
    } // GAUSSIANFILES
    { // GAUSSIANHEADER
      buffer = m_block["GAUSSIANHEADER"];

      if (!buffer.size()) {
        io::messages.add("no GAUSSIANHEADER block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.gaussian.input_header);
      DEBUG(1, "sim.param().qmmm.gaussian.input_header:");
      DEBUG(1, sim.param().qmmm.gaussian.input_header);
    } // GAUSSIANHEADER
    { // GAUSSIANROUTE
      buffer = m_block["GAUSSIANROUTE"];

      if (!buffer.size()) {
        io::messages.add("no GAUSSIANROUTE block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.gaussian.route_section);
      sim.param().qmmm.gaussian.route_section = "#" + sim.param().qmmm.gaussian.route_section;
      DEBUG(1, "sim.param().qmmm.gaussian.route_section:");
      DEBUG(1, sim.param().qmmm.gaussian.route_section);
    } // GAUSSIANROUTE
    { // GAUSSIANCHSM
      buffer = m_block["GAUSSIANCHSM"];

      if (!buffer.size()) {
        io::messages.add("no GAUSSIANCHSM block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      if (buffer.size() != 3) {
        io::messages.add("GAUSSIANCHSM block corrupt. Provide 1 line.",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.gaussian.chsm);
      DEBUG(1, "sim.param().qmmm.gaussian.chsm:");
      DEBUG(1, sim.param().qmmm.gaussian.chsm);
    } // GAUSSIANCHSM
  }

  /**
   * Schnetpack NN
   */
  else if (sw == simulation::qm_nn) {
    this->read_units(sim, &sim.param().qmmm.nn);
    //this->read_elements(topo, &sim.param().qmmm.nn);

    { // NNMODEL
      buffer = m_block["NNMODEL"];

      if (!buffer.size()) {
        io::messages.add("NNMODEL block missing",
                "In_QMMM", io::message::error);
        return;
      } else {
        if (buffer.size() != 3) {
          io::messages.add("NNMODEL block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.nn.model_path = buffer[1];
      }
    } // NNMODEL
    { // NNDEVICE
      buffer = m_block["NNDEVICE"];
      const std::set<std::string> allowed = {"auto", "cpu", "cuda"};

      if (!buffer.size()) {
        io::messages.add("NNDEVICE not specified.",
            "In_QMMM", io::message::notice);
      } else {
        if (buffer.size() != 3) {
          io::messages.add("NNDEVICE block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        const std::string device = buffer[1];
        if (device == "auto") sim.param().qmmm.nn.device = simulation::nn_device_auto;
        else if (device == "cuda") sim.param().qmmm.nn.device = simulation::nn_device_cuda;
        else if (device == "cpu") sim.param().qmmm.nn.device = simulation::nn_device_cpu;
        else {
          std::ostringstream msg;
          msg << "NNDEVICE block corrupt. Unknown option - " << device;
          io::messages.add(msg.str(), "In_QMMM", io::message::error);
          return;
        }
      }
    } // NNDEVICE
  }

  // Cap length definition
  if(topo.qmmm_link().size() > 0 ) {
    _lineStream.clear();
    buffer = m_block["CAPLEN"];
    if (!buffer.size()) {
      std::ostringstream msg;
      msg << "Using default capping atom bond length - "
          << sim.param().qmmm.cap_length;
      io::messages.add(msg.str(),"In_QMMM", io::message::notice);
    }
    else {
      double caplen;
      std::string line(*(buffer.begin() + 1));
      _lineStream.clear();
      _lineStream.str(line);
      _lineStream >> caplen;
      if (_lineStream.fail()) {
        io::messages.add("bad line in CAPLEN block.",
                          "In_QMMM", io::message::error);
        return;
        sim.param().qmmm.cap_length = caplen;
      }
    }
  }
}

void io::In_QMMM::read_elements(const topology::Topology& topo
    , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param)
  {
  std::vector<std::string> buffer = m_block["ELEMENTS"];

  if (!buffer.size()) {
    io::messages.add("No ELEMENTS block in QM/MM specification file",
            "In_QMMM", io::message::error);
    return;
  }
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1));
  unsigned Z;
  std::string element;
  while(_lineStream >> Z >> element) 
    qm_param->elements[Z] = element;
  if (_lineStream.fail())
    io::messages.add("Cannot read ELEMENTS block", "In_QMMM", io::message::error);

  // check whether all elements were provided
  const std::vector<unsigned>& atomic_number = topo.qm_atomic_number();
  for (std::vector<unsigned>::const_iterator
        it = atomic_number.begin(), to = atomic_number.end(); it != to; ++it)
    {
    if ((*it != 0) && (qm_param->elements.find(*it) == qm_param->elements.end())) {
      std::ostringstream msg;
      msg << "ELEMENTS block: No element name provided for atomic number " << *it;
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }
  }
}

void io::In_QMMM::read_units(const simulation::Simulation& sim
                  , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param)
  {
  std::map<simulation::qm_software_enum, std::array<double, 4> >
  unit_factor_defaults = {
    {simulation::qm_mndo,
                    { math::angstrom /* A */
                    , math::kcal /* kcal */
                    , math::kcal / math::angstrom /* kcal/A*/
                    , math::echarge /* e */}},
    {simulation::qm_turbomole,
                    { math::bohr /* a.u. */
                    , math::hartree * math::avogadro /* a.u. */
                    , math::hartree * math::avogadro / math::bohr /* a.u. */
                    , math::echarge /* e */}},
    {simulation::qm_dftb,
                    { math::bohr /* a.u. */
                    , math::hartree * math::avogadro /* a.u. */
                    , math::hartree * math::avogadro / math::bohr /* a.u. */
                    , math::echarge /* e */}},
    {simulation::qm_mopac,
                    { math::angstrom /* A */
                    , math::kcal /* kcal */
                    , math::kcal / math::angstrom /* kcal/A */
                    , math::echarge /* e */}},
    {simulation::qm_gaussian,
                    { math::angstrom /* A */
                    , math::hartree * math::avogadro /* a.u. */
                    , math::hartree * math::avogadro / math::bohr /* a.u. */
                    , math::echarge /* e */}},
    {simulation::qm_nn,
                    { math::bohr /* a.u. */
                    , math::hartree * math::avogadro /* a.u. */
                    , math::hartree * math::avogadro / math::bohr /* a.u. */
                    , math::echarge /* e */}}
  };

  std::vector<std::string> buffer = m_block["QMUNIT"];
  if (!buffer.size()) {
    std::array<double, 4> defaults = unit_factor_defaults[sim.param().qmmm.software];
    qm_param->unit_factor_length = defaults[0];
    qm_param->unit_factor_energy = defaults[1];
    qm_param->unit_factor_force  = defaults[2];
    qm_param->unit_factor_charge = defaults[3];
    std::ostringstream msg;
    msg << "Using default QMUNIT: "
        << qm_param->unit_factor_length << ", "
        << qm_param->unit_factor_energy << ", "
        << qm_param->unit_factor_force << ", "
        << qm_param->unit_factor_charge;
    io::messages.add(msg.str(),"In_QMMM", io::message::notice);
    return;
  }
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1));

  _lineStream >> qm_param->unit_factor_length
              >> qm_param->unit_factor_energy
              >> qm_param->unit_factor_force
              >> qm_param->unit_factor_charge;

  if (_lineStream.fail()) {
    io::messages.add("Bad line in QMUNIT block.",
            "In_QMMM", io::message::error);
    return;
  }
  DEBUG(15, "QM units read done");
}

void io::In_QMMM::read_zone(topology::Topology& topo
                           , simulation::Simulation& sim
                           , const std::string& blockname)
  {
  std::vector<std::string> buffer;
  buffer = m_block[blockname];

  if (!buffer.size()) {
    if (blockname == "QMZONE") {
      io::messages.add("No QMZONE block in QM/MM specification file",
              "In_QMMM", io::message::error);
      return;
    }
    else if (blockname == "BUFFERZONE") {
      return;
    }
  }

  std::string line(buffer[1]);
  _lineStream.clear();
  _lineStream.str(line);

  int charge, spin_mult;
  
  _lineStream >> charge
              >> spin_mult;
  if (blockname == "BUFFERZONE") {
    _lineStream >> sim.param().qmmm.buffer_zone.cutoff;
    sim.param().qmmm.buffer_zone.charge = charge;
    sim.param().qmmm.buffer_zone.spin_mult = spin_mult;
  }
  else if (blockname == "QMZONE") {
    sim.param().qmmm.qm_zone.charge = charge;
    sim.param().qmmm.qm_zone.spin_mult = spin_mult;
  }

  if (_lineStream.fail()) {
    std::ostringstream msg;
    msg << "Bad first line in " << blockname << " block";
    io::messages.add(msg.str(), "In_QMMM", io::message::error);
    return;
  }

  unsigned qmi, qmz, qmli;
  for (std::vector<std::string>::const_iterator it = buffer.begin() + 2
                                              , to = buffer.end() - 1
                                              ; it != to; ++it) {
    std::string line(*it);
    if (line.length() < 17) {
      std::ostringstream msg;
      msg << "Line too short in " << blockname << " block";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
    }

    // the first 17 chars are ignored
    line.erase(line.begin(), line.begin() + 17);

    _lineStream.clear();
    _lineStream.str(line);

    _lineStream >> qmi >> qmz >> qmli;

    if (_lineStream.fail()) {
      std::ostringstream msg;
      msg << "Bad line in " << blockname << " block";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }

    if (qmi < 1 || qmi > topo.num_atoms()) {
      std::ostringstream msg;
      msg << blockname << " block: atom out of range";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }

    if ((blockname == "QMZONE")
        && (qmi > topo.num_solute_atoms())) {
      std::ostringstream msg;
      msg << blockname << " block: QM atom should be in solute";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }

    if (qmz < 1) {
      std::ostringstream msg;
      msg << blockname << " block: wrong atomic number (QMZ)";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }

    if (qmli < 0 ) {
      std::ostringstream msg;
      msg << blockname << " block: QMLI has to be 0 or index of linked MM atom";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }
    topo.is_qm(qmi - 1) = topo.is_qm(qmi - 1) || (blockname == "QMZONE");
    const bool is_qm_buffer = (blockname == "BUFFERZONE");
    topo.is_qm_buffer(qmi - 1) = topo.is_qm_buffer(qmi - 1) || is_qm_buffer;
    sim.param().qmmm.use_qm_buffer = sim.param().qmmm.use_qm_buffer
                                      || is_qm_buffer;
    if (topo.is_qm(qmi - 1) && topo.is_qm_buffer(qmi - 1)) {
      std::ostringstream msg;
      msg << blockname << " block: atom " << qmi
          << " cannot be both in QM and buffer zone";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }

    topo.qm_atomic_number(qmi - 1) = qmz;
    if (qmli > 0 ) {
      DEBUG(15, "Linking " << qmi << " to " << qmli);
      topo.qmmm_link().insert(std::make_pair(qmi - 1, qmli - 1));
    }
  }

  for (std::set< std::pair<unsigned,unsigned> >::const_iterator
      it = topo.qmmm_link().begin(), to = topo.qmmm_link().end();
      it != to; ++it)
    {
    if (topo.is_qm(it->second)) {
      std::ostringstream msg;
      msg << blockname << " block: Invalid link - QMLI should be MM atom";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
    }
  }
}

void io::In_QMMM::remove_constraints(topology::Topology& topo)
  {
  // Remove distance constraints between QM atoms
  std::vector<topology::two_body_term_struct> & dc = topo.solute().distance_constraints();
  std::vector<topology::two_body_term_struct>::const_iterator it = dc.begin();
  while (it != dc.end()) {
    if (topo.is_qm(it->i) && topo.is_qm(it->j)) {
      DEBUG(15, "Removing distance constraint: " << it->i << "-" << it->j);
      it = dc.erase(it);
    }
    else ++it;
  }
}
