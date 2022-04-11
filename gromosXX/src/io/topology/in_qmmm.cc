
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
 * @section ORCA blocks for the ORCA worker
 * 
 * ORCABINARY block for the ORCA worker
 * The ORCABINARY block specifies path to ORCA binary
 *
 * This block is optional. If unspecified, orca command from PATH environment variable
 * is used.
 *
 * @verbatim
ORCABINARY
/path/to/orca/binary
END
@endverbatim
 *
 * ORCAFILES block for the ORCA worker
 * The ORCAFILES block specifies input and output files to exchange data with ORCA
 *
 * This block is optional. If unspecified, temporary files are created using TMPDIR
 * environment variable. User-specified files are not deleted after use.
 *
 * @verbatim
ORCAFILES
/path/to/orca.inp
/path/to/orca.out
/path/to/orca.xyz
/path/to/pointcharges.pc
/path/to/orca.engrad
/path/to/orca.pcgrad
END
@endverbatim
 *
 * The ORCAHEADER block specifies the header part of the ORCA input file. Variables
 * are allowed. Implemented are
 * - CHARGE: net charge of the QM zone
 * - SPINM: spin multiplicity of the QM zone 
 * - POINTCHARGES: file with information on pointcharges
 * - COORDINATES: coordinates of the QM zone
 * 
@verbatim
ORCAHEADER
! BP86 def2-SVP defgrid3 EnGrad TightSCF
%pal nprocs 4 end
%scf MaxIter 1500 end
%pointcharges "@@POINTCHARGES@@"
* xyzfile @@CHARGE@@ @@SPINM@@ @@COORDINATES@@
END
@endverbatim
 * Alternatively, for sempi-empiral methods:
 @verbatim
ORCAHEADER
! ZINDO/1 EnGrad TightSCF
%scf MaxIter 1500 end
%pointcharges "@@POINTCHARGES@@"
* xyzfile @@CHARGE@@ @@SPINM@@ @@COORDINATES@@
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

  /**
   * Orca
   */
  else if (sw == simulation::qm_orca) {
    this->read_units(sim, &sim.param().qmmm.orca);
    this->read_elements(topo, &sim.param().qmmm.orca);
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

  // Cap length definition
  if(topo.qmmm_link().size() > 0 ) {
    _lineStream.clear();
    buffer = m_block["CAPLEN"];
    if (!buffer.size()) {
      io::messages.add("no CAPLEN block in QM/MM specification file",
              "In_QMMM", io::message::error);
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
  std::string bstr = concatenate(buffer.begin() + 1, buffer.end() - 1);
  // Strip away the last newline character
  bstr.pop_back();
  _lineStream.str(bstr);
  unsigned Z;
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

void io::In_QMMM::read_units(const simulation::Simulation& sim
                  , simulation::Parameter::qmmm_struct::qm_param_struct* qm_param)
  {
  // std::map<simulation::qm_software_enum, std::array<double, 4> >
  // unit_factor_defaults = {
  //   {simulation::qm_mndo,
  //                   { math::angstrom /* A */
  //                   , math::kcal /* kcal */
  //                   , math::kcal / math::angstrom /* kcal/A*/
  //                   , math::echarge /* e */}},
  //   {simulation::qm_turbomole,
  //                   { math::bohr /* a.u. */
  //                   , math::hartree * math::avogadro /* a.u. */
  //                   , math::hartree * math::avogadro / math::bohr /* a.u. */
  //                   , math::echarge /* e */}},
  //   {simulation::qm_dftb,
  //                   { math::angstrom /* A */
  //                   , math::hartree * math::avogadro /* a.u. */
  //                   , math::hartree * math::avogadro / math::bohr /* a.u. */
  //                   , math::echarge /* e */}},
  //   {simulation::qm_mopac,
  //                   { math::angstrom /* A */
  //                   , math::kcal /* kcal */
  //                   , math::kcal / math::angstrom /* kcal/A */
  //                   , math::echarge /* e */}},
  //   {simulation::qm_gaussian,
  //                   { math::angstrom /* A */
  //                   , math::hartree * math::avogadro /* a.u. */
  //                   , math::hartree * math::avogadro / math::bohr /* a.u. */
  //                   , math::echarge /* e */}}
  //   {simulation::qm_orca,
  //                   { math::angstrom /* A */
  //                   , math::hartree * math::avogadro /* a.u. */
  //                   , math::hartree * math::avogadro / math::bohr /* a.u. */
  //                   , math::echarge /* e */}}
  // };

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
  // Count number of electrons and check consistency with charge and multiplicity
  int num_elec = -charge;
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
