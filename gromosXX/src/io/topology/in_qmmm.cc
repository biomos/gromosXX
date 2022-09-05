
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
# NETCH: net charge of the QM zone
# SPINM: net spin multiplicity of the QM zone
# QMI:  index of the QM atom
# QMZ:  atomic number of the QM atom
# QMLI: 0,1 atom is a link atom
#
# NETCH SPINM
      0     1
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
 * The block is read from the QM/MM specification file
 * (\@qmmm).
 *
 * @verbatim
BUFFERZONE
# NETCH: net charge of the buffer zone
# SPINM: net spin multiplicity of the buffer zone
# BUFCUT: cutoff of the adaptive buffer zone (default = 0.0, no adaptive buffer zone)
          If the specified atoms occur within BUFCUT of QM atom, they are considered
          as buffer atoms. Otherwise they are considered solely as MM atoms. BUFCUT = 0.0
          means that they are always considered as buffer atoms.
# QMI:   index of the QM atom
# QMZ:   atomic number of the QM atom
# QMLI:  0,1 atom is a link atom
# NETCH SPINM BUFCUT
      0     1    1.4
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
 * @section The ELEMENTS block specifies the element names used in QM packages.
 * It is determined by the atomic number given in the QMZONE block. This block
 * is required only for Turbomole, DFTB, and ORCA
 * 
@verbatim
ELEMENTS
1 h
6 c
7 n
8 o
END
@endverbatim
 *
 * @section The IAC block specifies the atomic number that should be used for each 
 * IAC atom type used. This block is only requiered for XTB to model point charges.
 * Default parameters are generated for the 54a7 ff atom types if the block is omitted.
 * IAC numbering start from "1" to match the documentation.
 * 
@verbatim
IAC
# IAC   ATMS   ATMN
  1     O      8
  2     O      8
  3     O      8
  4     O      8
  5     O      8
  6     N      7
  7     N      7
  8     N      7
  9     N      7
  10    N      7
  11    N      7
  12    C      6
  13    C      6
  14    C      6
  15    C      6
  16    C      6
  17    C      6
  18    C      6
  19    C      6
  20    H      1
  21    H      1
  22    Z      0
  23    S      16
  24    Cu     29
  25    Cu     29
  26    Fe     26
  27    Zn     30
  28    Mg     12
  29    Ca     20
  30    P      15
  31    Ar     18
  32    F      9
  33    Cl     17
  34    Br     35
  35    C      6
  36    O      8
  37    Na     11
  38    Cl     17
  39    C      6
  40    Cl     17
  41    H      1
  42    S      16
  43    C      6
  44    O      8
  45    C      6
  46    Cl     17
  47    F      9
  48    C      6
  49    C      6
  50    O      8
  51    C      6
  52    O      8
  53    N      7
  54    C      6
  55    I      53
  56    Cl     17
  57    B      5
  58    Se     34
  59    H      1
  60    Cl     17
  61    Br     35
  62    O      8
  63    N      7
  64    C      6
  65    C      6
  66    N      7
  67    N      7
  68    O      8
  69    Si     14
  70    P      15
END
@endverbatim
 *
 * @section MNDO blocks for the MNDO worker
 * 
 * MNDOBINARY block for the MNDO worker
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
 *
 * MNDOFILES block for the MNDO worker
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
 * - NUM_LINK: Number of link and capping atoms atoms
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
 * @section Turbomole blocks for the Turbomole worker
 * 
 * The TMOLEFILES block specifies where Turbomole writes the input and output files.
 * The first line contains the directory which contains the Turbomole binaries.
 * The second line is the turbomole working directory containing the control file.
 * In this control file the relative paths for the coordinate, point charges coordinates,
 * energy, gradients and point charges gradients are defined. Next lines contains
 * these filenames. Last line is a filename containing point charges of QM atoms
 * calculated using ESP. This is necessary only if GromosXX should obtain QM charges
 * from the QM calculation (NTQMMM = 1). The naming format of the standard output from
 * the Turbomole program is [program].out. Depending on your TMOLETOOLCHAIN definition
 * you have to decide which program outputs the charges and provide its output filename.
 * If (NTQMMM != 1), this line can be omitted.
 *
 * @verbatim
TMOLEFILES
/path/to/turbomole/binary/directory
/path/to/working/directory/containing/control/file
coordinate.in
mm_coordinate.in
energy.out
gradient.out
mm_gradient.out
ridft.out
END
@endverbatim
 * 
 * The TMOLETOOLCHAIN block specifies the Turbomole programs that are executed.
 * Each line contains one program that is called. By default, it is assumed that
 * the control file is static, i.e. that the number of QM atoms cannot change. To 
 * modify the control file during the simulation the TURBOMOLE program define is 
 * to be used. The input for this program has to be given as a filename "define.inp" 
 * which is in the same directory as the "control" file. Standard output of every
 * program is written to separate file called [program].out
 * Note: You still have to provide an initial control file as at the time the 
 * program define it only cannot include the $point_charges directives. 
 * 
 * @verbatim
 TMOLETOOLCHAIN
 ridft
 rdgrad
 END
 @endverbatim
 *
 * @section DFTB blocks for the DFTB worker
 * 
 * The DFTBFILES block specifies paths necessary to run the DFTB+ program. First line is
 * the path to DFTB binary. Second line defines path to the working directory where DFTB+
 * will run and write all the output files.
 * In this directory GromosXX creates dftb_in.hsd input file and fills it with
 * the content of the DFTBINPUT block. Also hardcoded DFTB output files are written here
 * including the details.out file used to transfer data back to GromosXX.
 * The third line defines the coordinate file. The fourth line is the MM coordinate
 * file containing the MM charges and is used only with electrostatic and polarisable
 * embedding. The input coordinate file has to be specified in the DFTBINPUT block.
 * This applies also to the MM coordinate file if used.
 * The fifth line specifies the path, where DFTB standard output is written. This file is
 * not used for data exchange and is usually needed for debugging purposes only.
 * 
 * @verbatim
DFTBFILES
/path/to/dftb/binary/directory
/path/to/working/directory
coordinate.in
mm_coordinate.in
stdout.out
END
@endverbatim 
 *
 * The DFTBINPUT block contains the complete content of the dftb_in.hsd file created 
 * by GromosXX in the working directory specified in the  DFTBFILES block. The example
 * block contains minimal settings that must be specified for proper data exchange
 * between GromosXX and DFTB+. User may provide additional DFTB+ blocks, e.g. 
 * Slater-Koster files path and/or other settings. Please refer to DFTB+ manual
 * for additional information on DFTB+ input.
 * 
 * @verbatim
DFTBINPUT
Geometry = genFormat {
 <<< "coordinate.in"   ## Always equired
}
Hamiltonian = DFTB {

  ## For electrostatic and polarisable embedding also use this block:
  ElectricField = {
    PointCharges = {
      CoordsAndCharges = {
        <<< "mm_coordinate.in"
      }
    }
  }
}
Analysis = {
  ## ... optional custom settings ...

  CalculateForces = Yes   ## Always required

  ## For mechanical embedding also specify:
  MullikenAnalysis = Yes
}
END
@endverbatim
 *
 * @section MOPAC blocks for the MOPAC worker
 * 
 * MOPACBINARY block for the MOPAC worker
 * The MOPACBINARY block specifies path to MOPAC binary
 *
 * This block is optional. If unspecified, mopac command from PATH environment variable
 * is used.
 *
 * @verbatim
MOPACBINARY
/path/to/mopac/binary
END
@endverbatim
 *
 * MOPACFILES block for the MOPAC worker
 * The MOPACFILES block specifies input and output files to exchange data with MOPAC
 *
 * This block is optional. If any line is not specified, temporary file is created
 * using TMPDIR environment variable and deleted after use.
 * User-specified files are not deleted.
 *
 * @verbatim
MOPACFILES
/path/to/mopac.mop    ##input file
/path/to/mopac.out    ##output file
/path/to/mopac.aux    ##auxiliary output file - MOPAC output is read from this file
/path/to/mopac.arc    ##archive output file
/path/to/stdout.out   ##standard output file - STDOUT and STDERR are redirected here
/path/to/mopac.den    ##density matrix file to be reused in subsequent step
/path/to/mol.in       ##mol.in file with electric potentials from external charges
END
@endverbatim
 *
 * The MOPACHEADER block specifies the header part of the MOPAC input file. Variables
 * are allowed. Implemented are
 * - CHARGE: net charge of the QM zone
 * - OLDENS: will generate density file in every step using DENOUT and then reuses
 *           it using OLDENS in subsequent step
 * 
@verbatim
MOPACHEADER
PM7 1SCF CHARGE=@@CHARGE@@ GRAD QMMM AUX(PRECISION=9) PRECISE @@OLDENS@@
title line
END
@endverbatim
 *
 * The MOPACLINKATOM block specifies the way the link atoms are treated. Link atoms are atoms
 * that are bonded to MM atoms.
@verbatim
MOPACLINKATOM
#  mode of treatment:
#    0(none) : link atoms see no MM atom (default)
#    1(exclude_atom) : Link atoms see all MM atoms except the MM atoms involved in the link
#    2(exclude_chargegroup) : Link atoms see all MM atoms except the MM chargegroup involved in the link
#    3(all) : Link atoms see all MM atoms
# 
   0
END
@endverbatim
 *
 * @section Gaussian blocks for the Gaussian worker
 * 
 * GAUBINARY block for the Gaussian worker
 * The GAUBINARY block specifies path to GAUSSIAN binary
 *
 * This block is optional. If unspecified, g16 command from PATH environment variable
 * is used.
 *
 * @verbatim
GAUBINARY
/path/to/gaussian/binary
END
@endverbatim
 * GAUFILES block for the Gaussian worker
 * The GAUFILES block specifies input and output files to exchange data with Gaussian
 *
 * This block is optional. If unspecified, temporary files are created using TMPDIR
 * environment variable. User-specified files are not deleted after use.
 *
 * @verbatim
GAUFILES
/path/to/gaussian.gjf
/path/to/gaussian.out
END
@endverbatim
 *
 * The GAUHEADER block specifies the header part of the Gaussian input file.
 * 
@verbatim
GAUHEADER
%nproc=8
%mem=2GB
%NoSave
%chk=tmp
END
@endverbatim
 *
 * The GAUROUTE block specifies the route section of the Gaussian input file.
 * hashsign (#) should be omitted. It is beneficial to generate an initial checkpoint file
 * and reuse it in subsequent steps with guess=read option. Allowed variables
 * - GUESS: replaced by 'guess=read' after first step
 * 
@verbatim
GAUROUTE
N hf/STO-3G @@GUESS@@ nosymm pop(mk) force charge(angstroms) prop=(field,read)
END
@endverbatim
 *
 * The GAUCHSM block specifies the net charge and the spin multiplicity of the system.
 * Variables are allowed. Implemented are
 * - CHARGE: net charge of the QM zone
 * - SPINM: spin multiplicity of the QM zone
@verbatim
GAUCHSM
@@CHARGE@@ @@SPINM@@
END
@endverbatim
 *
 * @section XTB blocks for the XTB worker
 * 
 * XTBOPTIONS block for the XTB worker
 * The XTBOPTIONS block specifies options for the XTB worker
 * Supported are selection of the Hamiltonian (GFNHAM): 1 or 2 for 
 * GFN1-xTB or GFN2-xTB, respectively and verbosity level: 
 * 0 (none), 1 (minimimal), and 2 (full) (XVERBO)
 *
 @verbatim
XTBOPTIONS
# GFNHAM    XVERBO
       1         1
END
 @endverbatim
 */
void io::In_QMMM::read(topology::Topology& topo,
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
    this->read_iac_elements(topo, &sim.param().qmmm.mndo);
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
    this->read_iac_elements(topo, &sim.param().qmmm.turbomole);

    { // TMOLEFILES
      buffer = m_block["TMOLEFILES"];

      if (!buffer.size()) {
        io::messages.add("TMOLEFILES block missing",
                "In_QMMM", io::message::error);
        return;
      } else {
        if (buffer.size() > 10
          || buffer.size() < 9) {
          io::messages.add("TMOLEFILES block corrupt. Provide 7 or 8 lines.",
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
        if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
          if (buffer.size() != 10) {
            io::messages.add("output charge filename missing in TMOLEFILES block",
                  "In_QMMM", io::message::error);
            return;
          }
          sim.param().qmmm.turbomole.output_charge_file = buffer[8];
        }
      }
    } // TMOLEFILES
    { // TMOLETOOLCHAIN
      buffer = m_block["TMOLETOOLCHAIN"];

      if (!buffer.size()) {
        io::messages.add("TMOLETOOLCHAIN block missing",
                "In_QMMM", io::message::error);
        return;
      } 
      for(unsigned int i = 1; i < buffer.size()-1; ++i) {
        _lineStream.clear();
        _lineStream.str(buffer[i]);
        std::string tool;
        _lineStream >> tool;
        if (_lineStream.fail()) {
          io::messages.add("bad line in TMOLETOOLCHAIN block",
                "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.turbomole.toolchain.push_back(tool);
      }
    } // TMOLETOOLCHAIN  
  }

  /**
   * DFTB
   */
  else if (sw == simulation::qm_dftb) {
    this->read_units(sim, &sim.param().qmmm.dftb);
    this->read_elements(topo, &sim.param().qmmm.dftb);
    this->read_iac_elements(topo, &sim.param().qmmm.dftb);
    { // DFTBFILES
      buffer = m_block["DFTBFILES"];
      if (!buffer.size()) {
        io::messages.add("DFTBFILES block missing",
                "In_QMMM", io::message::error);
        return;
      }
      if (buffer.size() != 7) {
        io::messages.add("DFTBFILES block corrupt. Provide 5 lines.",
                          "In_QMMM", io::message::error);
        return;
      }
      sim.param().qmmm.dftb.binary = buffer[1];
      sim.param().qmmm.dftb.working_directory = buffer[2];
      sim.param().qmmm.dftb.input_coordinate_file = buffer[3];
      sim.param().qmmm.dftb.input_mm_coordinate_file = buffer[4];
      sim.param().qmmm.dftb.stdout_file = buffer[5];
    } // DFTBFILES

    { // DFTBINPUT
      buffer = m_block["DFTBINPUT"];
      if (!buffer.size()) {
        io::messages.add("no DFTBINPUT block in QM/MM specification file",
                         "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
                  sim.param().qmmm.dftb.input_header);
    } // DFTBINPUT
  }

  /**
   * MOPAC
   */
  else if (sw == simulation::qm_mopac) {
    this->read_units(sim, &sim.param().qmmm.mopac);
    this->read_iac_elements(topo, &sim.param().qmmm.mopac);
    { // MOPACBINARY

      DEBUG(15, "Reading MOPACBINARY");
      buffer = m_block["MOPACBINARY"];

      if (!buffer.size()) {
        io::messages.add("Assuming that the mopac binary is in the PATH",
                "In_QMMM", io::message::notice);
        sim.param().qmmm.mopac.binary = "mopac";
      } else {
        if (buffer.size() != 3) {
          io::messages.add("MOPACBINARY block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.mopac.binary = buffer[1];
      }
    } // MOPACBINARY
    {
      //MOPACFILES
      buffer = m_block["MOPACFILES"];

      if (!buffer.size()) {
        io::messages.add("Using temporary files for MOPAC input/output",
                         "In_QMMM", io::message::notice);
      } else {
        if (buffer.size() != 9) {
          io::messages.add("MOPACFILES block corrupt. Provide 7 lines.",
                           "In_QMMM", io::message::error);
          return;
        }
      sim.param().qmmm.mopac.input_file = buffer[1];
      sim.param().qmmm.mopac.output_file = buffer[2];
      sim.param().qmmm.mopac.output_aux_file = buffer[3];
      sim.param().qmmm.mopac.output_arc_file = buffer[4];
      sim.param().qmmm.mopac.stdout_file = buffer[5];
      sim.param().qmmm.mopac.output_dens_file = buffer[6];
      sim.param().qmmm.mopac.molin_file = buffer[7];
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
    { // MOPACLINKATOM
      buffer = m_block["MOPACLINKATOM"];
      if (!buffer.size()) {
        io::messages.add("MOPAC link atoms see no MM atoms (default)",
                         "In_QMMM", io::message::notice);
        sim.param().qmmm.mopac.link_atom_mode = 0;
      } else {
        if (buffer.size() != 3) {
          io::messages.add("MOPACLINKATOM block corrupt. Provide 1 line.",
                           "In_QMMM", io::message::error);
          return;
        }
        std::string& s = buffer[1];
        if (s == "default" || s == "none" || s == "0") sim.param().qmmm.mopac.link_atom_mode = 0;
        else if (s == "exclude_atom" || s == "1") sim.param().qmmm.mopac.link_atom_mode = 1;
        else if (s == "exclude_chargegroup" || s == "2") sim.param().qmmm.mopac.link_atom_mode = 2;
        else if (s == "all" || s == "3") sim.param().qmmm.mopac.link_atom_mode = 3;
        else sim.param().qmmm.mopac.link_atom_mode = 0;
      }
    } // MOPACLINKATOM
  }

  /**
   * Gaussian
   */
  else if (sw == simulation::qm_gaussian) {
    this->read_units(sim, &sim.param().qmmm.gaussian);
    this->read_iac_elements(topo, &sim.param().qmmm.gaussian);
    { // GAUBINARY

      DEBUG(15, "Reading GAUBINARY");
      buffer = m_block["GAUBINARY"];

      if (!buffer.size()) {
        io::messages.add("Assuming that the g16 binary is in the PATH",
                "In_QMMM", io::message::notice);
        sim.param().qmmm.gaussian.binary = "g16";
      } else {
        if (buffer.size() != 3) {
          io::messages.add("GAUBINARY block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.gaussian.binary = buffer[1];
      }
    } // GAUBINARY
    { // GAUFILES

      DEBUG(15, "Reading GAUFILES");
      buffer = m_block["GAUFILES"];

      if (!buffer.size()) {
        io::messages.add("Using temporary files for Gaussian input/output",
                "In_QMMM", io::message::notice);
      } else {
        if (buffer.size() != 4) {
          io::messages.add("GAUFILES block corrupt. Provide 2 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.gaussian.input_file = buffer[1];
        sim.param().qmmm.gaussian.output_file = buffer[2];
      }
    } // GAUFILES
    { // GAUHEADER
      buffer = m_block["GAUHEADER"];

      if (!buffer.size()) {
        io::messages.add("no GAUHEADER block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.gaussian.input_header);
      DEBUG(1, "sim.param().qmmm.gaussian.input_header:");
      DEBUG(1, sim.param().qmmm.gaussian.input_header);
    } // GAUHEADER
    { // GAUROUTE
      buffer = m_block["GAUROUTE"];

      if (!buffer.size()) {
        io::messages.add("no GAUROUTE block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.gaussian.route_section);
      sim.param().qmmm.gaussian.route_section = "#" + sim.param().qmmm.gaussian.route_section;
      DEBUG(1, "sim.param().qmmm.gaussian.route_section:");
      DEBUG(1, sim.param().qmmm.gaussian.route_section);
    } // GAUROUTE
    { // GAUCHSM
      buffer = m_block["GAUCHSM"];

      if (!buffer.size()) {
        io::messages.add("no GAUCHSM block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      if (buffer.size() != 3) {
        io::messages.add("GAUCHSM block corrupt. Provide 1 line.",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.gaussian.chsm);
      DEBUG(1, "sim.param().qmmm.gaussian.chsm:");
      DEBUG(1, sim.param().qmmm.gaussian.chsm);
    } // GAUCHSM
  }

  /*
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
        if (buffer.size() != 4) {
          io::messages.add("NNMODEL block corrupt. Provide 2 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.nn.model_path = buffer[1];
        std::string line(buffer[2]);
        _lineStream.clear();
        _lineStream.str(line);
        int model_type;
        int learning_type;
        _lineStream >> model_type >> learning_type;
        if (model_type == 0) sim.param().qmmm.nn.model_type = simulation::nn_model_type_burnn;
        else if (model_type == 1) sim.param().qmmm.nn.model_type = simulation::nn_model_type_standard;
        else {
          std::ostringstream msg;
          msg << "NNMODEL block corrupt. Unknown option for model type - " << model_type;
          io::messages.add(msg.str(), "In_QMMM", io::message::error);
          return;
        }
        if (learning_type == 1) sim.param().qmmm.nn.learning_type = simulation::nn_learning_type_all;
        else if (learning_type == 2) sim.param().qmmm.nn.learning_type = simulation::nn_learning_type_qmonly;
        if (_lineStream.fail()) {
          io::messages.add("bad line in NNMODEL block",
                "In_QMMM", io::message::error);
          return;
        }
      }
    } // NNMODEL
    { // NNVALID
      buffer = m_block["NNVALID"];
      if (buffer.size()) {
        if (buffer.size() != 4) {
          io::messages.add("NNVALID block corrupt. Provide 2 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.nn.val_model_path = buffer[1];
        std::string line(buffer[2]);
        _lineStream.clear();
        _lineStream.str(line);
        unsigned val_steps;
        double val_thresh, val_forceconstant;
        _lineStream >> val_steps >> val_thresh >> val_forceconstant;
        sim.param().qmmm.nn.val_steps = val_steps;
        sim.param().qmmm.nn.val_thresh = val_thresh;
        sim.param().qmmm.nn.val_forceconstant = val_forceconstant;
        if (_lineStream.fail()) {
          io::messages.add("bad line in NNVALID block",
                "In_QMMM", io::message::error);
          return;
        }
      }
    } // NNVALID
    { // NNCHARGE
      buffer = m_block["NNCHARGE"];
      if (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic) {
        if (buffer.size()) {
          if (buffer.size() != 4) {
            io::messages.add("NNCHARGE block corrupt. Provide 2 lines.",
                    "In_QMMM", io::message::error);
            return;
          }
          sim.param().qmmm.nn.charge_model_path = buffer[1];
          std::string line(buffer[2]);
          _lineStream.clear();
          _lineStream.str(line);
          unsigned charge_steps;
          _lineStream >> charge_steps;
          sim.param().qmmm.nn.charge_steps = charge_steps;
          if (_lineStream.fail()) {
            io::messages.add("bad line in NNCHARGE block",
                  "In_QMMM", io::message::error);
            return;
          }
        } else {
          io::messages.add("Dynamic QM charges requested but NNCHARGE is not specified.",
              "In_QMMM", io::message::error);
        }
      }
    } // NNCHARGE
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

  /**
   * Orca
   */
  else if (sw == simulation::qm_orca) {
    this->read_units(sim, &sim.param().qmmm.orca);
    this->read_elements(topo, &sim.param().qmmm.orca);
    this->read_iac_elements(topo, &sim.param().qmmm.orca);
    { // ORCABINARY

      DEBUG(15, "Reading ORCABINARY");
      buffer = m_block["ORCABINARY"];

      if (!buffer.size()) {
        io::messages.add("Assuming that the orca binary is in the PATH",
                "In_QMMM", io::message::notice);
        sim.param().qmmm.orca.binary = "orca";
      } else {
        if (buffer.size() != 3) {
          io::messages.add("ORCABINARY block corrupt. Provide 1 line.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.orca.binary = buffer[1];
      }
    } // ORCABINARY
    { // ORCAFILES

      DEBUG(15, "Reading ORCAFILES");
      buffer = m_block["ORCAFILES"];

      if (!buffer.size()) {
        io::messages.add("Using temporary files for Orca input/output",
                "In_QMMM", io::message::notice);
      } else {
        if (buffer.size() != 8) {
          io::messages.add("ORCAFILES block corrupt. Provide 6 lines.",
                  "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.orca.input_file = buffer[1];
        sim.param().qmmm.orca.output_file = buffer[2];
        sim.param().qmmm.orca.input_coordinate_file = buffer[3];
        sim.param().qmmm.orca.input_pointcharges_file = buffer[4];
        sim.param().qmmm.orca.output_gradient_file = buffer[5];
        sim.param().qmmm.orca.output_mm_gradient_file = buffer[6];
      }
    } // ORCAFILES
    { // ORCAHEADER
      buffer = m_block["ORCAHEADER"];

      if (!buffer.size()) {
        io::messages.add("no ORCAHEADER block in QM/MM specification file",
                "In_QMMM", io::message::error);
        return;
      }
      concatenate(buffer.begin() + 1, buffer.end() - 1,
              sim.param().qmmm.orca.input_header);
      DEBUG(1, "sim.param().qmmm.orca.input_header:");
      DEBUG(1, sim.param().qmmm.orca.input_header);
    } // ORCAHEADER
  }

  /**
   * XTB
   */
  else if (sw == simulation::qm_xtb) {
    this->read_units(sim, &sim.param().qmmm.xtb);
    this->read_iac_elements(topo, &sim.param().qmmm.xtb);
    { // XTBOPTIONS
      buffer = m_block["XTBOPTIONS"];

      if (!buffer.size()) {
        io::messages.add("no XTBOPTIONS block in QM/MM specification file",
                "In_QMMM", io::message::error);
      }
      else {
        unsigned int hamiltonian, verbosity;
        std::string line(buffer[1]);
        _lineStream.clear();
        _lineStream.str(line);
        _lineStream >> hamiltonian >> verbosity;
        if (_lineStream.fail()) {
          io::messages.add("bad line in XTBOPTIONS block.",
                            "In_QMMM", io::message::error);
          return;
        }
        sim.param().qmmm.xtb.hamiltonian = hamiltonian;
        sim.param().qmmm.xtb.verbosity = verbosity;
      }
      DEBUG(1, "sim.param().qmmm.xtb.hamiltonian:");
      DEBUG(1, sim.param().qmmm.xtb.hamiltonian);
      DEBUG(1, "sim.param().qmmm.xtb.verbosity:");
      DEBUG(1, sim.param().qmmm.xtb.verbosity);
    } // XTBOPTIONS
  }

  // Cap length definition
  { // CAPLEN
    if(topo.qmmm_link().size() > 0 ) {
      _lineStream.clear();
      buffer = m_block["CAPLEN"];
      if (!buffer.size()) {
        io::messages.add("no CAPLEN block in QM/MM specification file",
                "In_QMMM", io::message::error);
      }
      else {
        double caplen = 0.0;
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
  } // CAPLEN

  // check if NN charge model is defined
  if (sw == simulation::qm_nn
        && sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic
        && sim.param().qmmm.nn.charge_model_path.empty()) {
    io::messages.add("dynamic QM charge requested but no NN charge model specified",
                "In_QMMM", io::message::error);
  }

  // allow learning_type == nn_learning_type_qmonly only for single-atom QM region
  if (sw == simulation::qm_nn
        && sim.param().qmmm.nn.learning_type == simulation::nn_learning_type_qmonly) {
    size_t qm_size = 0;
    for (unsigned i = 0; i < topo.num_atoms(); ++i) {
      if (topo.is_qm(i)) ++qm_size;
      if (qm_size > 1) {
        io::messages.add("multi-atom QM region with learning type other than 1 is not implemented",
                "In_QMMM", io::message::error);
        return;
      }
    }
  }

  // if dynamic charges for buffer will be used
  sim.param().qmmm.dynamic_buffer_charges =
      sim.param().qmmm.use_qm_buffer && (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic);
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
  std::string bstr = concatenate(buffer.begin() + 1, buffer.end() - 1);
  // Strip away the last newline character
  bstr.pop_back();
  _lineStream.str(bstr);
  unsigned Z = 0;
  std::string element;
  while(!_lineStream.eof()) {
    _lineStream >> Z >> element;
    if (_lineStream.fail()) {
      io::messages.add("Cannot read ELEMENTS block", "In_QMMM", io::message::error);
      return;
    }
    qm_param->elements[Z] = element;
  }

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

void io::In_QMMM::read_iac_elements(topology::Topology& topo
    , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param)
  {
  std::vector<std::string> buffer = m_block["IAC"];

  if (!buffer.size()) { // provide default matching
    io::messages.add("No IAC block in QM/MM specification file. Using Gromos 54a7_ff definitions.",
            "In_QMMM", io::message::notice);
    // indices start with "0" to match the Gromos C++ coding style
    qm_param->iac_elements[0]  = 8;
    qm_param->iac_elements[1]  = 8;
    qm_param->iac_elements[2]  = 8;
    qm_param->iac_elements[3]  = 8;
    qm_param->iac_elements[4]  = 8;
    qm_param->iac_elements[5]  = 7;
    qm_param->iac_elements[6]  = 7;
    qm_param->iac_elements[7]  = 7;
    qm_param->iac_elements[8]  = 7;
    qm_param->iac_elements[9] = 7;
    qm_param->iac_elements[10] = 7;
    qm_param->iac_elements[11] = 6;
    qm_param->iac_elements[12] = 6;
    qm_param->iac_elements[13] = 6;
    qm_param->iac_elements[14] = 6;
    qm_param->iac_elements[15] = 6;
    qm_param->iac_elements[16] = 6;
    qm_param->iac_elements[17] = 6;
    qm_param->iac_elements[18] = 6;
    qm_param->iac_elements[19] = 1;
    qm_param->iac_elements[20] = 1;
    qm_param->iac_elements[21] = 0; // dummy atom
    qm_param->iac_elements[22] = 16;
    qm_param->iac_elements[23] = 29;
    qm_param->iac_elements[24] = 29;
    qm_param->iac_elements[25] = 26;
    qm_param->iac_elements[26] = 30;
    qm_param->iac_elements[27] = 12;
    qm_param->iac_elements[28] = 20;
    qm_param->iac_elements[29] = 14; // P oder Si
    qm_param->iac_elements[30] = 18;
    qm_param->iac_elements[31] = 9;
    qm_param->iac_elements[32] = 17;
    qm_param->iac_elements[33] = 35;
    qm_param->iac_elements[34] = 6;
    qm_param->iac_elements[35] = 8;
    qm_param->iac_elements[36] = 11;
    qm_param->iac_elements[37] = 17;
    qm_param->iac_elements[38] = 6;
    qm_param->iac_elements[39] = 17;
    qm_param->iac_elements[40] = 1;
    qm_param->iac_elements[41] = 16;
    qm_param->iac_elements[42] = 6;
    qm_param->iac_elements[43] = 8;
    qm_param->iac_elements[44] = 6;
    qm_param->iac_elements[45] = 17;
    qm_param->iac_elements[46] = 9;
    qm_param->iac_elements[47] = 6;
    qm_param->iac_elements[48] = 6;
    qm_param->iac_elements[49] = 8;
    qm_param->iac_elements[50] = 6;
    qm_param->iac_elements[51] = 8;
    qm_param->iac_elements[52] = 7;
    qm_param->iac_elements[53] = 6;
  }
  else {
    io::messages.add("Reading IAC block in QMMM specification file.",
            "In_QMMM", io::message::notice);
    _lineStream.clear();
    std::string bstr = concatenate(buffer.begin() + 1, buffer.end() - 1);
    // Strip away the last newline character
    bstr.pop_back();
    _lineStream.str(bstr);
    unsigned int iac;
    std::string atomic_symbol; // unused
    unsigned int atomic_number;
    while(!_lineStream.eof()) {
      _lineStream >> iac >> atomic_symbol >> atomic_number;
      if (_lineStream.fail()) {
        io::messages.add("Cannot read IAC block", "In_QMMM", io::message::error);
        return;
      }
      iac--;
      // IAC indices in QM/MM specification file start from "1" to match documentation
      // IAC indices in Gromos C++ standard start from "0" to match C++ data structure layouts
      qm_param->iac_elements[iac] = atomic_number;
    }
  }
  
  // check whether all atom types of the force field used have a matching atom number
  const std::map<std::string, int>& atomic_names = topo.atom_names();
  for (std::map<std::string, int>::const_iterator
        it = atomic_names.begin(), to = atomic_names.end(); it != to; ++it) {
    if (qm_param->iac_elements.find(it->second) == qm_param->iac_elements.end()) {
      std::ostringstream msg;
      msg << "IAC block: No atomic number provided for IAC " << (it->second + 1); // +1: match documentation for IAC atoms
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }
    else {
      DEBUG(15, "IAC atom \"" << it->first << "\" (#" << (it->second + 1) << ") is mapped to element #" 
        << qm_param->iac_elements.at(it->second)); // +1: match documentation for IAC atoms
    }
  }
  // assign atoms in topology - note that this is only relevant for MM atoms
  // this procedure allows usage of QM atom types for which no force field parametrers
  // are available (e.g. Pd can be modelled with Si) 
  for (unsigned int i = 0; i < topo.iac().size(); ++i) {
    if (!topo.is_qm(i)) { // QM atoms are assigned in QMZONE block
      DEBUG(15, "Atom #" << (i+1) << " with IAC #" << topo.iac(i) + 1 << " is assigned element symbol #" // +1: match documentation
        << qm_param->iac_elements.at(topo.iac(i))); 
      topo.qm_atomic_number(i) = qm_param->iac_elements.at(topo.iac(i));
    }
  }
}

void io::In_QMMM::read_units(const simulation::Simulation& sim
                  , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param)
  {
  std::vector<std::string> buffer = m_block["QMUNIT"];
  if (!buffer.size()) {
    io::messages.add("no QMUNIT block in QM/MM specification file",
            "In_QMMM", io::message::error);
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

  int charge = 0, spin_mult = 0;
  
  _lineStream >> charge
              >> spin_mult;
  if (blockname == "BUFFERZONE") {
    _lineStream >> sim.param().qmmm.buffer_zone.cutoff;
    if (sim.param().pairlist.skip_step > 1) {
    // Adaptive QM buffer and skipping pairlist may suddenly appear/disappear particles
      io::messages.add("BUFFERZONE block: With adaptive QM buffer, the pairlist should be generated every step"
                      , "In_QMMM", io::message::warning);
    }
    sim.param().qmmm.buffer_zone.charge = charge;
    sim.param().qmmm.buffer_zone.spin_mult = spin_mult;
  }
  else if (blockname == "QMZONE") {
    sim.param().qmmm.qm_zone.charge = charge;
    sim.param().qmmm.qm_zone.spin_mult = spin_mult;
  }
  else {
    std::ostringstream msg;
    msg << blockname << " block not implemented";
    io::messages.add(msg.str(), "In_QMMM", io::message::critical);
    return;
  }

  if (_lineStream.fail()) {
    std::ostringstream msg;
    msg << "Bad first line in " << blockname << " block";
    io::messages.add(msg.str(), "In_QMMM", io::message::error);
    return;
  }
  // Count number of electrons and check consistency with charge and multiplicity
  int num_elec = -charge;
  unsigned qmi = 0, qmz = 0, qmli = 0;
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
      msg << blockname << " block: solvent in QMZONE is not supported (atom " << qmi << ")";
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

    if (is_qm_buffer
        && (sim.param().qmmm.qm_ch == simulation::qm_ch_dynamic)
        && (sim.param().innerloop.method != simulation::sla_off)
        && (qmi > topo.num_solute_atoms())) {
      std::ostringstream msg;
      msg << blockname 
          << " block: dynamic QM charge for solvent atoms cannot be used"
          << " with non-default solvent loops";
      io::messages.add(msg.str(), "In_QMMM", io::message::error);
      return;
    }

    topo.qm_atomic_number(qmi - 1) = qmz;
    if (qmli > 0 ) {
      DEBUG(15, "Linking " << qmi << " to " << qmli);
      topo.qmmm_link().insert(std::make_pair(qmi - 1, qmli - 1));
    }

    num_elec += qmz + (bool)qmli;
  }

  // Check charge and multiplicity consistency
  bool open_shell = num_elec % 2;
  // sm=0 is closed-shell in MNDO program and is in principle sm=1
  if (!spin_mult) spin_mult = 1;
  if (open_shell == spin_mult % 2) {
      std::ostringstream msg;
      msg << blockname << " block: " << num_elec << " electrons but multiplicity " << spin_mult;
      io::messages.add(msg.str(), "In_QMMM", io::message::warning);
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
  std::vector<topology::two_body_term_struct>::iterator it = dc.begin();
  while (it != dc.end()) {
    if (topo.is_qm(it->i) && topo.is_qm(it->j)) {
      DEBUG(15, "Removing distance constraint: " << it->i << "-" << it->j);
      it = dc.erase(it);
    }
    else ++it;
  }
}
