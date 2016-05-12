/**
 * @file in_parameter.cc
 * implements methods of In_Parameter
 */
#include <stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>

#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include <io/instream.h>
#include <io/blockinput.h>
#include <io/parameter/in_parameter.h>

#include <math/random.h>

#include <configuration/energy.h>

#include <string>
#include <iterator>

#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#define MODULE io
#undef SUBMODULE
#define SUBMODULE parameter

static std::set<std::string> block_read;

/**
 * Store standard parameters in the Parameter
 */
void io::In_Parameter::read(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(7, "reading input");

  if (!quiet)
    os << "INPUT\n"
          << title;

  // store the title...
  param.title = title;

  read_ENERGYMIN(param);
  read_SYSTEM(param);
  read_COMTRANSROT(param); // has to be read before INITIALISE
  read_INITIALISE(param);
  read_STEP(param);
  read_BOUNDCOND(param);
  read_REPLICA(param); // has to be read in before MULTIBATH
  read_MULTIBATH(param);
  read_PRESSURESCALE(param);
  read_PRINTOUT(param);
  read_WRITETRAJ(param);
  read_CONSTRAINT(param);
  read_FORCE(param);
  read_COVALENTFORM(param);
  read_CGRAIN(param);
  read_HOOMD(param);
  read_PAIRLIST(param);
  read_NONBONDED(param);
  read_POSITIONRES(param);
  read_DISTANCERES(param);
  read_DISTANCEFIELD(param); 
  read_DIHEDRALRES(param); // needs to be called after CONSTRAINT!
  read_PERTURBATION(param);
  read_JVALUERES(param);
  read_ORDERPARAMRES(param);
  read_RDCRES(param);
  read_XRAYRES(param);
  read_PERSCALE(param);
  read_ROTTRANS(param);
  read_INNERLOOP(param);
  read_MULTICELL(param);
  read_READTRAJ(param);
  read_INTEGRATE(param);
  read_STOCHDYN(param);
  read_EWARN(param);
  read_MULTISTEP(param);
  read_MONTECARLO(param);
  read_POLARISE(param);
  read_RANDOMNUMBERS(param);
  read_EDS(param);
  read_LAMBDAS(param); // needs to be called after FORCE
  read_LOCALELEV(param);
  read_BSLEUS(param);
  read_ELECTRIC(param);
  read_SASA(param);
  read_ADDECOUPLE(param); // needs to be called after MULTIBATH and FORCE
  read_NEMD(param);
  read_MULTIGRADIENT(param);
  read_QMMM(param);
  read_SYMRES(param);

  read_known_unsupported_blocks();

  DEBUG(7, "input read...");

  for (std::map<std::string, std::vector<std::string> >::const_iterator
    it = m_block.begin(),
          to = m_block.end();
          it != to;
          ++it) {

    if (block_read.count(it->first) == 0 && it->second.size()) {
      io::messages.add("block " + it->first + " not supported!",
              "In_Parameter",
              io::message::error);
    }
  }

  if (!quiet)
    os << "END\n";

}

/**
 * @section system SYSTEM block
 * @verbatim
SYSTEM
# NPM : protein molecules (0, 1)
# NSM : solvent molecules (>= 0)
#      NPM      NSM 
         1        0
END
@endverbatim
 */
void io::In_Parameter::read_SYSTEM(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "reading SYSTEM");

  std::vector<std::string> buffer;
  buffer = m_block["SYSTEM"];
  std::string s;

  if (!buffer.size()) {
    io::messages.add("no SYSTEM block in input", "In_Parameter", io::message::error);
    param.system.nsm = 0;
    param.system.npm = 0;
    return;
  }

  block_read.insert("SYSTEM");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  _lineStream >> param.system.npm
          >> param.system.nsm;

  if (_lineStream.fail())
    io::messages.add("bad line in SYSTEM block",
          "In_Parameter", io::message::error);

  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in SYSTEM block, but should have been: \n" + _lineStream.str() +  "\n",
               "In_Parameter", io::message::warning);
   */

  // we might need to also allow for 0...
  if (param.system.npm != 1 && param.system.npm != 0)
    io::messages.add("SYSTEM: currently only NPM=1 allowed (NPM=0 experimental)",
          "io::In_Parameter::read_SYSTEM",
          io::message::error);
  if (param.system.nsm < 0)
    io::messages.add("SYSTEM: NSM should be >0",
          "io::In_Parameter::read_SYSTEM",
          io::message::error);

}

/**
 * @section energymin ENERGYMIN block
 * @verbatim
ENERGYMIN
# NTEM: 0..2 controls energy minimisation mode
#       0: do not do energy minimisation (default)
#       1: steepest-descent minimisation
#       2: conjugate-gradient minimisation (promd only)
# NCYC: >0 number of steps before resetting of conjugate-gradient search direction
# DELE: >0.0 energy threshold for convergence
# DX0: > 0.0 initial step size
# DXM: > 0.0 maximum step size
# NMIN > 0 minimum number of minimisation steps
# FLIM >= 0.0 limit force to maximum value (FLIM > 0.0 is not recommended).
#     NTEM    NCYC    DELE    DX0     DXM    NMIN    FLIM
         1       0     0.1   0.01    0.05       1       0
END
@endverbatim
 */
void io::In_Parameter::read_ENERGYMIN(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "reading ENERGYMIN");

  std::vector<std::string> buffer;
  buffer = m_block["ENERGYMIN"];
  std::string s;

  if (!buffer.size()) {
    // no minimisation
    return;
  }

  buffer = m_block["ENERGYMIN"];

  block_read.insert("ENERGYMIN");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  _lineStream >> param.minimise.ntem
          >> param.minimise.ncyc
          >> param.minimise.dele
          >> param.minimise.dx0
          >> param.minimise.dxm
          >> param.minimise.nmin
          >> param.minimise.flim;

  if (_lineStream.fail())
    io::messages.add("bad line in ENERGYMIN block",
          "In_Parameter", io::message::error);

  // allow 0 to disable feature...
  if (param.minimise.nmin == 0)
    param.minimise.nmin = 1;

  if (param.minimise.ntem != 0 && param.minimise.ntem != 1)
    io::messages.add("ENERGYMIN block: currently only steepest descent implemented",
          "io::In_Parameter",
          io::message::error);
  if (param.minimise.ntem == 1 && param.minimise.ncyc > 0)
    io::messages.add("ENERGYMIN block: NCYC > 0 has no effect for steepest descent",
          "io::In_Parameter",
          io::message::warning);

  if (param.minimise.ncyc < 0)
    io::messages.add("ENERGYMIN block: NCYC should be >0",
          "io::In_Parameter",
          io::message::error);

  if (param.minimise.dele < 0)
    io::messages.add("ENERGYMIN block: DELE should be >0",
          "io::In_Parameter",
          io::message::error);

  if (param.minimise.dx0 < 0)
    io::messages.add("ENERGYMIN block: DX0 should be >0",
          "io::In_Parameter", io::message::error);

  if (param.minimise.dxm < param.minimise.dx0)
    io::messages.add("ENERGYMIN block: DXM should be > DX0",
          "io::In_Parameter", io::message::error);

  if (param.minimise.nmin <= 0)
    io::messages.add("ENERGYMIN block: NMIN should be >= 0",
          "io::In_Parameter", io::message::error);

  if (param.minimise.flim < 0)
    io::messages.add("ENERGYMIN: FLIM should be >= 0",
          "io::In_Parameter",
          io::message::error);
  if (param.minimise.flim > 0)
    io::messages.add("ENERGYMIN: FLIM > 0 may result in "
          "failure of the minimisation procedure."
          " Only to be used in special cases.",
          "io::In_Parameter",
          io::message::warning);

}

/**
 * @section step STEP block
 * @verbatim
STEP
#   NSTLIM  number of steps
#   T       initial time
#   DT      time step
#
#   NSTLIM         T        DT
       100       0.0     0.005
END
@endverbatim
 */
void io::In_Parameter::read_STEP(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read STEP");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["STEP"];

  if (!buffer.size()) {
    io::messages.add("no STEP block in input", "In_Parameter", io::message::error);
    return;
  }

  block_read.insert("STEP");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  _lineStream >> param.step.number_of_steps
          >> param.step.t0
          >> param.step.dt;

  if (_lineStream.fail())
    io::messages.add("bad line in STEP block", "In_Parameter", io::message::error);

  if (param.step.t0 < 0 && param.step.t0 != -1.0)
    io::messages.add("STEP block: Negative time is not supported",
          "In_Parameter", io::message::error);
  if (param.step.number_of_steps <= 0)
    io::messages.add("STEP block: We want to do at least one step...",
          "In_Parameter", io::message::error);
}

/**
 * @section constraint CONSTRAINT block
 * @verbatim
CONSTRAINT
#	NTC
#		1	"solvent"	solvent only
#		2	"hydrogen"	solvent and solute bonds containing hydrogens
#		3	"all"		solvent and solute all bonds
#		4	"specified"	solvent and constraints in the CONSTRAINTS block
    3
#       NTCP: solute algorithm
#               1        shake
#               2        lincs
#               3        flexshake
#       NTCP
        1
#       NTCP0(1)..(3): algorithm options
#         - shake: tolerance
#         - lincs: order
#         - flexshake: tolerance, readin, order
#       NTCP0(1)   NTCP0(2)   NTCP0(3)
        0.0001
#	NTCS: solvent algorithm
#               1        shake
#               2        lincs
#               3        flexshake
#               4        settle
#               5        m_shake (only implemented for water and methanol!)
#               6        gpu_shake
#       NTCS
        1
#       NTCS0(1):  algorithm options
#         - shake: tolerance
#         - lincs: order
#         - flexshake: tolerance, readin, order
#         - settle: no arguments
#         - m_shake: tolerance
#       NTCS0(1)
        0.0001
#       NTCG: Number of GPUs
#       NTCD: Device number of the GPU
#       NTCG          NTCD
        1             0
END
@endverbatim
 */
void io::In_Parameter::read_CONSTRAINT(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read CONSTRAINT");

  std::vector<std::string> buffer;
  std::string s, salg;

  buffer = m_block["CONSTRAINT"];

  if (!buffer.size()) {
    param.constraint.ntc = 1;
    param.constraint.solute.algorithm = simulation::constr_off;
    param.constraint.solvent.algorithm = simulation::constr_shake;

    io::messages.add("no CONSTRAINT block", "In_Parameter",
            io::message::warning);
    return;
  }

  block_read.insert("CONSTRAINT");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  std::string sntc;
  _lineStream >> sntc;
  if (sntc == "off" || sntc == "0") param.constraint.ntc = 0;
  if (sntc == "solvent" || sntc == "1") param.constraint.ntc = 1;
  else if (sntc == "hydrogen" || sntc == "2") param.constraint.ntc = 2;
  else if (sntc == "all" || sntc == "3") param.constraint.ntc = 3;
  else if (sntc == "specified" || sntc == "4") param.constraint.ntc = 4;
  else {
    std::stringstream ss(sntc);
    if (!(ss >> param.constraint.ntc))
      io::messages.add("CONSTRAINT block: NTC not understood",
            "In_Parameter", io::message::error);
  }

  if (param.constraint.ntc < 0 || param.constraint.ntc > 4) {
    io::messages.add("CONSTRAINT block: NTC out of range",
            "In_Parameter", io::message::error);

    param.constraint.solute.algorithm = simulation::constr_off;
    param.constraint.solvent.algorithm = simulation::constr_off;

    return;
  }

  // SOLUTE
  _lineStream >> salg;

  DEBUG(8, "Constraints (solute): " << salg);

  if (_lineStream.fail())
    io::messages.add("bad line in CONSTRAINT block",
          "In_Parameter", io::message::error);

  std::transform(salg.begin(), salg.end(), salg.begin(), tolower);

  if (salg == "shake" || salg == "1") {
    DEBUG(9, "constraints solute shake");

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_shake;
    else param.constraint.solute.algorithm = simulation::constr_off;

    _lineStream >> param.constraint.solute.shake_tolerance;

    if (param.constraint.solute.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0",
            "In_Parameter", io::message::error);
  } else if (salg == "flexshake" || salg == "3") {
    DEBUG(9, "constraints solute flexshake");

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_flexshake;
    else param.constraint.solute.algorithm = simulation::constr_off;

    _lineStream >> param.constraint.solute.shake_tolerance
            >> param.constraint.solute.flexshake_readin
            >> param.constraint.solute.flexshake_mode;

    if (param.constraint.solute.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0.0",
            "In_Parameter", io::message::error);

    if (param.constraint.solute.flexshake_mode < 0 ||
            param.constraint.solute.flexshake_mode > 3)
      io::messages.add("CONSTRAINT block: flexshake mode should be >= 0 and <= 3",
            "In_Parameter", io::message::error);

  } else if (salg == "lincs" || salg == "2") {
    DEBUG(9, "constraints solute lincs");

    if (param.constraint.ntc > 1)
      param.constraint.solute.algorithm = simulation::constr_lincs;
    else
      param.constraint.solute.algorithm = simulation::constr_off;

    _lineStream >> param.constraint.solute.lincs_order;

    if (param.constraint.solute.lincs_order < 1)
      io::messages.add("CONSTRAINT block: lincs order should be > 1",
            "In_Parameter", io::message::error);

  } else if (salg == "off" || salg == "0") {
    DEBUG(9, "constraints solute off");
    param.constraint.solute.algorithm = simulation::constr_off;
  } else {
    DEBUG(9, "constraints solute error");

    io::messages.add("CONSTRAINT block: unknown algorithm (solute)",
            "In_Parameter", io::message::error);
    param.constraint.solute.algorithm = simulation::constr_off;
  }

  // SOLVENT
  _lineStream >> salg;

  DEBUG(8, "constraints solvent: " << salg);

  if (_lineStream.fail())
    io::messages.add("bad line in CONSTRAINT block",
          "In_Parameter", io::message::error);

  std::transform(salg.begin(), salg.end(), salg.begin(), tolower);

  if (salg == "shake" || salg == "1") {
    DEBUG(9, "constraints solvent shake");

    param.constraint.solvent.algorithm = simulation::constr_shake;
    _lineStream >> param.constraint.solvent.shake_tolerance;

    if (param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0.0",
            "In_Parameter", io::message::error);
  } else if (salg == "flexshake" || salg == "3") {
    DEBUG(9, "constraints solvent flexshake");

    param.constraint.solvent.algorithm = simulation::constr_flexshake;
    _lineStream >> param.constraint.solvent.shake_tolerance;

    if (param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: shake tolerance should be > 0.0",
            "In_Parameter", io::message::error);
  } else if (salg == "lincs" || salg == "2") {
    DEBUG(9, "constraints solvent lincs");

    param.constraint.solvent.algorithm = simulation::constr_lincs;
    _lineStream >> param.constraint.solvent.lincs_order;

    if (param.constraint.solvent.lincs_order < 1)
      io::messages.add("CONSTRAINT block: lincs order should be >1",
            "In_Parameter", io::message::error);

  } else if (salg == "settle" || salg == "4") {
    DEBUG(9, "constraints solvent settle");

    param.constraint.solvent.algorithm = simulation::constr_settle;
  } else if (salg == "m_shake" || salg == "5") {
    DEBUG(9, "constraints solvent m_shake");

    param.constraint.solvent.algorithm = simulation::constr_m_shake;
    _lineStream >> param.constraint.solvent.shake_tolerance;

    if (param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: m_shake tolerance should be > 0.0",
            "In_Parameter", io::message::error);
  }/*else if (salg == "gpu_settle" || salg == "6") {
    DEBUG(9, "constraints solvent settle on GPU");

    param.constraint.solvent.algorithm = simulation::constr_gpu_settle;

  } */ else if (salg == "gpu_shake" || salg == "6") {
    DEBUG(9, "constraints solvent gpu_shake");

    param.constraint.solvent.algorithm = simulation::constr_gpu_shake;
    _lineStream >> param.constraint.solvent.shake_tolerance;

    if (param.constraint.solvent.shake_tolerance <= 0.0)
      io::messages.add("CONSTRAINT block: gpu_shake tolerance should be > 0.0",
            "In_Parameter", io::message::error);

    // Get number of GPUs and their IDs
    _lineStream >> param.constraint.solvent.number_gpus;
    if (_lineStream.fail() || param.constraint.solvent.number_gpus == (unsigned int) 0)
      io::messages.add("CUDA enabled, but number of GPUs is zero",
            "In_Parameter", io::message::error);
    else {
      unsigned int temp;
      bool fail = false;
      for (unsigned int i = 0; i < param.constraint.solvent.number_gpus; i++) {
        _lineStream >> temp;
        if (_lineStream.fail()) {
          fail = true;
          break;
        }
        param.constraint.solvent.gpu_device_number.push_back(temp);
      }
      if (fail) {
        param.constraint.solvent.gpu_device_number.clear();
        param.constraint.solvent.gpu_device_number.resize(param.constraint.solvent.number_gpus, -1);
        io::messages.add("CUDA driver will determine devices for M-SHAKE",
                "In_Parameter", io::message::notice);
      }
    }
  } else if (salg == "off" || salg == "0") {
    DEBUG(9, "constraints solvent off");

    param.constraint.solvent.algorithm = simulation::constr_off;
    io::messages.add("CONSTRAINT block: no constraints for SOLVENT",
            "In_Parameter", io::message::warning);
  } else {

    DEBUG(9, "constraints solvent error");
    io::messages.add("CONSTRAINT block: unknown algorithm (solvent)",
            "In_Parameter", io::message::error);

    param.constraint.solvent.algorithm = simulation::constr_off;
  }
}

/**
 * @section printout PRINTOUT block
 * @verbatim
PRINTOUT
#  NTPR: print out energies, etc. every NTPR steps
#  NTPP: =1 perform dihedral angle transition monitoring
#     NTPR      NTPP
         0         0
END
@endverbatim
 */
void io::In_Parameter::read_PRINTOUT(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read PRINTOUT");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["PRINTOUT"];

  if (!buffer.size()) {
    io::messages.add("no PRINTOUT block", "In_Parameter", io::message::notice);
    return;
  }

  block_read.insert("PRINTOUT");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  _lineStream >> param.print.stepblock
          >> param.print.monitor_dihedrals;

  if (_lineStream.fail())
    io::messages.add("bad line in PRINTOUT block",
          "In_Parameter", io::message::error);

  if (param.print.stepblock < 0)
    io::messages.add("PRINTOUT block: NTPR should be >=0.",
          "In_Parameter", io::message::error);

  if (param.print.monitor_dihedrals < 0 || param.print.monitor_dihedrals > 1)
    io::messages.add("PRINTOUT block: NTPP should be 0 or 1.",
          "In_Parameter", io::message::error);

}

/**
 * @section writetraj WRITETRAJ block
 * @verbatim
WRITETRAJ
# NTWX       controls writing of coordinate trajectory
#       0: no coordinate trajectory is written (default)
#      >0: write solute and solvent coordinates every NTWX steps
#      <0: write solute coordinates every |NTWX| steps
# NTWSE >= 0 selection criteria for coordinate trajectory writing
#       0: write normal coordinate trajectory
#      >0: write minimum-energy coordinate and energy trajectory (based on the
#          energy entry selected by NTWSE and as blocks of length NTWX)
#          (see configuration/energy.cc or ene_ana library for indices)
# NTWV       controls writing of velocity trajectory
#       0: no velocity trajectory is written (default)
#      >0: write solute and solvent velocities every NTWV steps
#      <0: write solute velocities every |NTWV| steps
# NTWF       controls writing of force trajectory
#       0: no force trajectory is written (default)
#      >0: write solute and solvent forces every NTWF steps
#      <0: write solute forces every |NTWF| steps
# NTWE >= 0 controls writing of energy trajectory
#       0: no energy trajectory is written (default)
#      >0: write energy trajectory every NTWE steps
# NTWG >= 0 controls writing of free energy trajectory
#       0: no free energy trajectory is written (default)
#      >0: write free energy trajectory every NTWG steps
# NTWB >= 0 controls writing of block-averaged energy trajectory
#       0: no block averaged energy trajectory is written (default)
#      >0: write block-averaged energy variables every |NTWB| steps
#          (and free energies if NTWG > 0) trajectory
#
#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
      100          0         0         0       100         0       100
END
@endverbatim
 */
void io::In_Parameter::read_WRITETRAJ(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read WRITETRAJ");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["WRITETRAJ"];

  if (!buffer.size()) {
    io::messages.add("no WRITETRAJ block", "In_Parameter", io::message::notice);
    return;
  }

  block_read.insert("WRITETRAJ");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  _lineStream >> param.write.position
          >> param.write.energy_index
          >> param.write.velocity
          >> param.write.force
          >> param.write.energy
          >> param.write.free_energy
          >> param.write.block_average;

  if (_lineStream.fail())
    io::messages.add("bad line in WRITETRAJ block",
          "In_Parameter", io::message::error);

  if (param.write.energy_index < 0 ||
          param.write.energy_index > int(configuration::Energy::MAX_ENERGY_INDEX)) {
    std::ostringstream msg;
    msg << "WRITETRAJ block: NTWSE must be 0 to "
            << configuration::Energy::MAX_ENERGY_INDEX;
    io::messages.add(msg.str(), "In_Parameter", io::message::error);
  }

  if (param.write.energy_index > 0) {
    if (param.write.position == 0) {
      io::messages.add("WRITETRAJ block: NTWX must be a block size >= 0 for "
              "minimum-energy trajectory.",
              "In_Parameter", io::message::error);
    }
    if (param.write.energy != 0 && param.write.energy != abs(param.write.position)) {
      io::messages.add("WRITETRAJ block: NTWE must be 0 or abs(NTWX) for "
              "minimum-energy trajectory.",
              "In_Parameter", io::message::error);
    }
    // from the documentation all others needs to be zero
    if (param.write.velocity != 0 || param.write.force != 0 ||
            param.write.free_energy != 0 || param.write.block_average != 0) {
      io::messages.add("WRITETRAJ block: NTWV, NTWF, NTWG and NTWB must be 0 for "
              "minimum-energy trajectory.",
              "In_Parameter", io::message::error);
    }
  }

  if (param.write.position < 0) {
    param.write.position_solute_only = true;
    param.write.position = -param.write.position;
    io::messages.add("writing solute only position trajectory",
            "In_Parameter", io::message::notice);
  }

  if (param.write.velocity < 0) {
    param.write.velocity_solute_only = true;
    param.write.velocity = -param.write.velocity;
    io::messages.add("writing solute only velocity trajectory",
            "In_Parameter", io::message::notice);
  }

  if (param.write.force < 0) {
    param.write.force_solute_only = true;
    param.write.force = -param.write.force;
    io::messages.add("writing solute only force trajectory",
            "In_Parameter", io::message::notice);
  }

  if (param.write.energy < 0)
    io::messages.add("WRITETRAJ block: NTWE should be >= 0",
          "In_Parameter", io::message::error);
  if (param.write.free_energy < 0)
    io::messages.add("WRITETRAJ block: NTWG should be >= 0",
          "In_Parameter", io::message::error);
  if (param.write.block_average < 0)
    io::messages.add("WRITETRAJ block: NTWB should be >= 0",
          "In_Parameter", io::message::error);
}

/**
 * @section pressurescale PRESSURESCALE block
 * @verbatim
PRESSURESCALE
#	COUPLE:	off(0), calc(1), scale(2)
#	SCALE:  off(0), iso(1), aniso(2), full(3), semianiso(4)
#	VIRIAL: none(0), atomic(1), group(2)
#
#   COUPLE  SCALE   COMP        TAUP    VIRIAL
    calc    iso     4.575E-4    0.5     atomic
#   SEMI (semianisotropic couplings: X, Y, Z)
#       e.g. 1 1 2: x and y jointly coupled and z separately coupled 
#       e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath
    1 1 2
#   reference pressure
    0.06102     0.00000     0.00000
    0.00000     0.06102     0.00000
    0.00000     0.00000     0.06102
END
@endverbatim
 */
void io::In_Parameter::read_PRESSURESCALE(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read PRESSURESCALE");

  std::vector<std::string> buffer;
  std::string s;


  // first try for a PRESSURESCALE block
  buffer = m_block["PRESSURESCALE"];

  if (buffer.size()) {

    block_read.insert("PRESSURESCALE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    std::string s1, s2, s3;
    int xs, ys, zs;

    _lineStream >> s1 >> s2
            >> param.pcouple.compressibility
            >> param.pcouple.tau
            >> s3;

    _lineStream >> xs >> ys >> zs;

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        _lineStream >> param.pcouple.pres0(i, j);
      }
    }

    if (_lineStream.fail())
      io::messages.add("bad line in PRESSURESCALE block",
            "In_Parameter", io::message::error);

    //std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    //std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    //std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    if (s1 == "off" || s1 == "0") {
      param.pcouple.calculate = false;
      param.pcouple.scale = math::pcouple_off;
    } else if (s1 == "calc" || s1 == "1") {
      param.pcouple.calculate = true;
      param.pcouple.scale = math::pcouple_off;
    } else if (s1 == "scale" || s1 == "2") {
      param.pcouple.calculate = true;

      if (s2 == "off" || s2 == "0") {
        io::messages.add("PRESSURESCALE block: requesting scaling but SCALE set to OFF",
                "In_Parameter", io::message::error);
        param.pcouple.scale = math::pcouple_off;
      } else if (s2 == "iso" || s2 == "1")
        param.pcouple.scale = math::pcouple_isotropic;
      else if (s2 == "aniso" || s2 == "2")
        param.pcouple.scale = math::pcouple_anisotropic;
      else if (s2 == "full" || s2 == "3")
        param.pcouple.scale = math::pcouple_full_anisotropic;
      else if (s2 == "semianiso" || s2 == "4") {
        param.pcouple.scale = math::pcouple_semi_anisotropic;
        if (xs < 0 || xs > 2)
          io::messages.add("PRESSURESCALE block: bad value for x component of COUPLINGS FOR SEMIANISOTROPIC"
                "(0,1,2)",
                "In_Parameter", io::message::error);
        if (ys < 0 || ys > 2)
          io::messages.add("PRESSURESCALE block: bad value for y component of COUPLINGS FOR SEMIANISOTROPIC"
                "(0,1,2)",
                "In_Parameter", io::message::error);
        if (zs < 0 || zs > 2)
          io::messages.add("PRESSURESCALE block: bad value for z component of COUPLINGS FOR SEMIANISOTROPIC"
                "(0,1,2)",
                "In_Parameter", io::message::error);
        if (xs == ys && ys == zs)
          io::messages.add("PRESSURESCALE block: x_semi = y_semi = z_semi in SEMI, maybe you want isotropic pressure scaling?"
                "(0,1,2)",
                "In_Parameter", io::message::error);
        param.pcouple.x_semi = xs;
        param.pcouple.y_semi = ys;
        param.pcouple.z_semi = zs;
      } else {
        io::messages.add("PRESSURESCALE block: bad value for SCALE switch "
                "(off,iso,aniso,full)",
                "In_Parameter", io::message::error);
        param.pcouple.scale = math::pcouple_off;
      }

    } else {
      io::messages.add("bad value for calc switch in PRESSURESCALE block\n"
              "(off,calc,scale)",
              "In_Parameter", io::message::error);
      param.pcouple.calculate = false;
    }

    if (param.pcouple.calculate) {
      if (s3 == "none" || s3 == "0") {
        io::messages.add("requesting pressure calculation but "
                "no virial specified",
                "In_Parameter", io::message::error);
        param.pcouple.virial = math::no_virial;
      } else if (s3 == "atomic" || s3 == "1")
        param.pcouple.virial = math::atomic_virial;
      else if (s3 == "molecular" || s3 == "group" || s3 == "2")
        param.pcouple.virial = math::molecular_virial;
      else {
        io::messages.add("bad value for virial switch in PRESSURESCALE block\n"
                "(none,atomic,molecular)",
                "In_Parameter", io::message::error);
        param.pcouple.virial = math::no_virial;
      }
    } else
      param.pcouple.virial = math::no_virial;

  } // PRESSURESCALE block

  if (param.pcouple.calculate == false && param.pcouple.scale != math::pcouple_off)
    io::messages.add("PRESSURESCALE block: pressure coupling activated but "
          "not calculating pressure",
          "In_Parameter",
          io::message::error);
  if (param.pcouple.calculate == true && param.pcouple.virial == math::no_virial)
    io::messages.add("PRESSURESCALE block: pressure calculation requested but"
          " no virial specified!", "In_Parameter",
          io::message::error);
  if (param.pcouple.compressibility <= 0)
    io::messages.add("PRESSURESCALE block: compressibility should be >0 ",
          "In_Parameter", io::message::error);
  if (param.pcouple.tau <= 0)
    io::messages.add("PRESSURESCALE block: tau should be >0 ",
          "In_Parameter", io::message::error);
}

/**
 * @section boundcond BOUNDCOND block
 * @verbatim
BOUNDCOND
#  NTB: boundary conditions
#       -1 : truncoct
#        0 : vacuum
#        1 : rectangular
#        2 : triclinic
#  NDFMIN: number of degrees of freedin subtracted for temperature
#
#         NTB    NDFMIN
            1         0
END
@endverbatim
 */
void io::In_Parameter::read_BOUNDCOND(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read BOUNDCOND");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["BOUNDCOND"];

  if (!buffer.size()) {
    io::messages.add("no BOUNDCOND block", "In_Parameter", io::message::error);
    param.boundary.boundary = math::vacuum;
    return;
  }

  block_read.insert("BOUNDCOND");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int ntb;
  _lineStream >> ntb >> param.boundary.dof_to_subtract;

  if (_lineStream.fail())
    io::messages.add("bad line in BOUNDCOND block",
          "In_Parameter", io::message::error);

  if (ntb == 0) param.boundary.boundary = math::vacuum;
  else if (ntb == 1) param.boundary.boundary = math::rectangular;
  else if (ntb == 2) param.boundary.boundary = math::triclinic;
  else if (ntb == -1) {
    param.boundary.boundary = math::truncoct;
    io::messages.add("The box is converted to triclinic boundary coniditons internally", "In_Parameter", io::message::notice);
  } else {
    std::ostringstream msg;
    msg << "BOUNDCOND block: wrong value for NTB "
            << ntb << "\nvacuum (0), rectangular (1), triclinic (2), truncoct (-1)";
    io::messages.add(msg.str(), "In_Parameter", io::message::error);
    param.boundary.boundary = math::vacuum;
  }

  if (param.boundary.dof_to_subtract < 0) {
    io::messages.add("BOUNDCOND block: NDFMIN must be >= 0.",
            "In_Parameter", io::message::error);
    param.boundary.dof_to_subtract = 0;
  }
}

/**
 * @section perturbation PERTURBATION block
 * @verbatim
PERTURBATION
#    NTG: 0..1 controls use of free-energy calculation.
#         0: no free-energy calculation (default)
#         1: calculate dH/dRLAM
#  NRDGL: 0,1 controls reading of initial value for RLAM.
#         0: use initial RLAM parameter from PERTURBATION block
#         1: read from configuration
#   RLAM: 0.0..1.0 initial value for lambda
#  DLAMT: >= 0.0 rate of lambda increase in time.
# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter
#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter
#   NLAM: > 0 power dependence of lambda coupling
# NSCALE: 0..2 controls use of interaction scaling
#         0: no interaction scaling
#         1: interaction scaling
#         2: perturbation for all atom pairs with scaled
#            interactions. No perturbation for others.
#
#     NTG   NRDGL    RLAM   DLAMT
        0       0     0.0     0.0
#  ALPHLJ   ALPHC    NLAM  NSCALE
      0.0     0.0       1       0
END
@endverbatim
 */
void io::In_Parameter::read_PERTURBATION(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read PERTURBATION");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["PERTURBATION"];
  if (!buffer.size())
    return;

  block_read.insert("PERTURBATION");


  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  std::string b, s1, s2;
  int ntg, scale;
  int nrdgl;

  _lineStream >> ntg
          >> nrdgl
          >> param.perturbation.lambda
          >> param.perturbation.dlamt
          >> param.perturbation.soft_vdw
          >> param.perturbation.soft_crf
          >> param.perturbation.lambda_exponent
          >> scale;

  if (_lineStream.fail())
    io::messages.add("bad line in PERTURBATION block",
          "In_Parameter", io::message::error);

  switch (ntg) {
    case 0:
      param.perturbation.perturbation = false;
      break;
    case 1:
      param.perturbation.perturbation = true;
      break;
    default:
      io::messages.add("PERTURBATION block: NTG must be 0 or 1.",
              "In_Parameter", io::message::error);
  }

  switch (nrdgl) {
    case 0: // use from input file
      param.perturbation.read_initial = false;
      break;
    case 1: // use from configuration
      param.perturbation.read_initial = true;
      break;
    default:
      io::messages.add("PERTURBATION block: NRDGL must be 0 or 1.",
              "In_Parameter", io::message::error);
  }

  if (param.perturbation.lambda < 0.0 ||
          param.perturbation.lambda > 1.0) {
    io::messages.add("PERTURBATION block: RLAM must be 0.0 to 1.0.",
            "In_Parameter", io::message::error);
  }

  if (param.perturbation.dlamt < 0.0) {
    io::messages.add("PERTURBATION block: DLAMT must be >= 0.",
            "In_Parameter", io::message::error);
  }

  if (param.perturbation.lambda_exponent <= 0) {
    io::messages.add("PERTURBATION block: NLAM must be > 0.",
            "In_Parameter", io::message::error);
  }

  if (param.perturbation.soft_vdw < 0.0) {
    io::messages.add("PERTURBATION block: ALPHLJ must be >= 0.",
            "In_Parameter", io::message::error);
  }

  if (param.perturbation.soft_crf < 0.0) {
    io::messages.add("PERTURBATION block: ALPHC must be >= 0.",
            "In_Parameter", io::message::error);
  }

  switch (scale) {
    case 0: // no scaling
      param.perturbation.scaling = false;
      param.perturbation.scaled_only = false;
      break;
    case 1: // scaling on
      param.perturbation.scaling = true;
      param.perturbation.scaled_only = false;
      break;
    case 2: // scaled only
      param.perturbation.scaling = true;
      param.perturbation.scaled_only = true;
      break;
    default:
      io::messages.add("PERTURBATION block: NSCALE must be 0 to 2.",
              "In_Parameter", io::message::error);
  }
}

/**
 * @section force FORCE block
 * @verbatim
FORCE
# NTF(1..6): 0,1 determines terms used in force calculation
#             0: do not include terms
#             1: include terms
# NEGR: ABS(NEGR): number of energy groups
#             > 0: use energy groups
#             < 0: use energy and force groups
# NRE(1..NEGR): >= 1.0 last atom in each energy group
# NTF(1)    NTF(2)    NTF(3)    NTF(4)    NTF(5)        NTF(6)
# bonds     angles    improper  dihedral  electrostatic vdW
  0         1         1         1         1             1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     1        60           
END
@endverbatim
 */
void io::In_Parameter::read_FORCE(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read FORCE");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["FORCE"];

  if (!buffer.size()) {
    DEBUG(8, "no force block found???");
    io::messages.add("no FORCE block", "In_Parameter", io::message::error);
    return;
  }

  block_read.insert("FORCE");

  // NTS
  _lineStream.clear();
  _lineStream.str(buffer[1]);
  
  int a;
  std::vector<int> ntf;
  while (_lineStream >> a) {
    ntf.push_back(a);
  }
  
  int bondH, angleH, impH, dihedralH;
  if (ntf.size() == 10) { // old FORCE block
    io::messages.add("Old FORCE block used.", "In_Parameter",
            io::message::warning);
    bondH = ntf[0];
    param.force.bond = ntf[1];
    angleH = ntf[2];
    param.force.angle = ntf[3];
    impH = ntf[4];
    param.force.improper = ntf[5];
    dihedralH = ntf[6];
    param.force.dihedral = ntf[7];
    param.force.nonbonded_crf = ntf[8];
    param.force.nonbonded_vdw = ntf[9];
    
    // checks
    if (bondH ^ param.force.bond)
    io::messages.add("FORCE block: switch for bond and bond H has to be equal",
          "In_Parameter", io::message::error);

    if (angleH ^ param.force.angle)
      io::messages.add("FORCE block: switch for angle and angle H has to be equal",
            "In_Parameter", io::message::error);

    if (impH ^ param.force.improper)
      io::messages.add("FORCE block: switch for improper and improper H has to be equal",
            "In_Parameter", io::message::error);

    if (dihedralH ^ param.force.dihedral)
      io::messages.add("FORCE block: switch for dihedral and dihedral H has to be equal",
            "In_Parameter", io::message::error);
  } else if (ntf.size() == 6) { // new FORCE block
    param.force.bond = ntf[0];
    param.force.angle = ntf[1];
    param.force.improper = ntf[2];
    param.force.dihedral = ntf[3];
    param.force.nonbonded_crf = ntf[4];
    param.force.nonbonded_vdw = ntf[5];
  } else {
    io::messages.add("FORCE block: bad line in NTS",
            "In_Parameter", io::message::error);
  }
  
  // energy groups
  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 2, buffer.end() - 1, s));

  int snum;
  unsigned int num, e, old_e = 0;

  _lineStream >> snum;
  
  if (snum < 0) {
#ifdef XXFORCEGROUPS
    param.force.force_groups = true;
    num = abs(snum);
#else
    io::messages.add("Force groups requested but not compiled with support for them."
    "Use --enable-forcegroups for configure.", "In_Parameter",
            io::message::error);
    return;
#endif
  }
  num = abs(snum);
  if (!num) {
    io::messages.add("number of energy groups must not be zero.", "In_Parameter",
            io::message::error);
    return;
  }

  for (unsigned int i = 0; i < num; ++i) {
    _lineStream >> e;
    DEBUG(10, "\tadding energy group " << e - 1);
    param.force.energy_group.push_back(e - 1);
    if (e <= old_e) {
      DEBUG(10, "energy groups not in order...");
      io::messages.add("FORCE block: energy groups are not in order",
              "In_Parameter", io::message::error);
      return;
    }
    old_e = e;
  }
  // Now that we have the energy groups, we initialize the  	 	 
  // LAMBDAS parameters that depend on them.  	 	 
  // NOTE: lambdas vectors may be resized again in in_topology.cc
  // in case an extra energy group is added. This will be the case
  // if the last atom of the last energy group in the force array 
  // is not the last atom of the system.

  int maxnilg = param.force.energy_group.size();
  std::vector< double > one(maxnilg, 1.0);
  std::vector< double > zero(maxnilg, 0.0);
  for (unsigned int i = 0; i < param.lambdas.a.size(); i++) {
    param.lambdas.a[i].resize(maxnilg, zero);
    param.lambdas.b[i].resize(maxnilg, zero);
    param.lambdas.c[i].resize(maxnilg, zero);
    param.lambdas.d[i].resize(maxnilg, one);
    param.lambdas.e[i].resize(maxnilg, zero);
  }

  DEBUG(10, "number of energy groups: " << param.force.energy_group.size());

  if (_lineStream.fail())
    io::messages.add("FORCE block: bad line in ENERGYGROUP",
          "In_Parameter", io::message::error);

  if ((!param.force.nonbonded_crf) && param.force.nonbonded_vdw)
    io::messages.add("FORCE block: setting charges to zero",
          "In_Parameter", io::message::notice);

  if (param.force.nonbonded_crf && (!param.force.nonbonded_vdw))
    io::messages.add("FORCE block: setting atom types to dummy",
          "In_Parameter", io::message::notice);

  if (_lineStream.fail())
    io::messages.add("bad line in FORCE block", "In_Parameter", io::message::error);

  /*
  if (!_lineStream.eof())
    io::messages.add("End of line not reached in FORCE block, but should have been: \n" + s +  "\n",
             "In_Parameter", io::message::warning);
   */

  if (param.force.bond < 0 || param.force.bond > 1)
    io::messages.add("FORCE block: Illegal value for force switch for bond",
          "In_Parameter", io::message::error);
  if (param.force.angle < 0 || param.force.angle > 1)
    io::messages.add("FORCE block: Illegal value for force switch for angle",
          "In_Parameter", io::message::error);
  if (param.force.improper < 0 || param.force.improper > 1)
    io::messages.add("FORCE block: Illegal value for force switch for improper dihedral",
          "In_Parameter", io::message::error);
  if (param.force.dihedral < 0 || param.force.dihedral > 1)
    io::messages.add("FORCE block: Illegal value for force switch for dihedral",
          "In_Parameter", io::message::error);
  if (param.force.nonbonded_vdw < 0 || param.force.nonbonded_vdw > 1)
    io::messages.add("FORCE block: Illegal value for force switch for nonbonded (vdw)",
          "In_Parameter", io::message::error);
  if (param.force.nonbonded_crf < 0 || param.force.nonbonded_crf > 1)
    io::messages.add("FORCE block: Illegal value for force switch for nonbonded (crf)",
          "In_Parameter", io::message::error);
}

/**
 * @section covalentform COVALENTFORM block
 * @verbatim
COVALENTFORM
# NTBBH: 0,1 controls bond-stretching potential
#        0: quartic form (default)
#        1: harmonic form
# NTBAH: 0,1 controls bond-angle bending potential
#        0: cosine-harmonic (default)
#        1: harmonic
# NTBDN: 0,1 controls torsional dihedral potential
#        0: arbitrary phase shifts (default)
#        1: phase shifts limited to 0 and 180 degrees.
#   NTBBH    NTBAH    NTBDN
        0        0        0
END
@endverbatim
 */
void io::In_Parameter::read_COVALENTFORM(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read COVALENTFORM");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["COVALENTFORM"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("COVALENTFORM");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int bond, angle, dihedral;

  _lineStream >> bond
          >> angle
          >> dihedral;

  if (bond != 0 && bond != 1) {
    io::messages.add("COVALENTFORM block: NTBBH must be 0 (quartic) "
            "or 1 (harmonic).",
            "In_Parameter", io::message::error);
  } else {
    if (param.force.bond != 0) {
      switch (bond) {
        case 1: param.force.bond = 2;
          break;
        case 0:
        default: param.force.bond = 1;
      }
    }
  }

  if (angle != 0 && angle != 1) {
    io::messages.add("COVALENTFORM block: NTBAH must be 0 (quartic) "
            "or 1 (harmonic).",
            "In_Parameter", io::message::error);
  } else {
    if (param.force.angle != 0) {
      switch (angle) {
        case 1: param.force.angle = 2;
          break;
        case 0:
        default: param.force.angle = 1;
      }
    }
  }

  if (dihedral != 0 && dihedral != 1) {
    io::messages.add("COVALENTFORM block: NTBDN must be 0 (arbitrary "
            "phase shifts) or 1 (phase shifts limited).",
            "In_Parameter", io::message::error);
  } else {
    if (param.force.dihedral != 0) {
      switch (dihedral) {
        case 1: param.force.dihedral = 2;
          break;
        case 0:
        default: param.force.dihedral = 1;
      }
    }
  }
}

/**
 * @section initialise INITIALISE block
 * @verbatim
INITIALISE
# NTIVEL: 0,1 controls generation of initial velocities.
#         0: read from configuration (default)
#         1: generate from Maxell distribution at temperature TEMPI
# NTISHK: 0..3 controls shaking of initial configuration
#         0: no intial SHAKE (default)
#         1: initial SHAKE on coordinates only
#         2: initial SHAKE on velocities only
#         3: initial SHAKE on coordinates and velocities
# NTINHT: 0,1 controls generation of initial Nose-Hoover chain variables
#         0: read from configuration (default)
#         1: reset variables to zero.
# NTINHB: 0,1 controls generation of initial Nose-Hoover (chain) barostat
#             variables
#         0: read from strartup file (if applicable) (default)
#         1: reset variables to zero
# NTISHI: 0,1 controls initial setting for lattice shift vectors
#         0: read from configuration (default)
#         1: reset shifts to zero.
# NTIRTC: 0,1 controls initial setting of positions and orientations for
#             roto-translational constraints
#         0: read from configuration (default)
#         1: reset based on initial configuraion of startup file
# NTICOM: 0,1 controls initial removal of COM motion
#         0: no initial system COM motion removal (default)
#         1: initial COM translation is removed
#         2: initial COM rotation is removed
# NTISTI: 0,1 controls generation of stochastic integrals
#         0: read stochastic integrals and IG from configuration (default)
#         1: set stochastic integrals to zero and use IG from here.
# IG:     random number generator seed
# TEMPI:  initial temperature
#
#  NTIVEL  NTISHK  NTINHT  NTINHB
        0       0       0       0
#  NTISHI  NTIRTC  NTICOM
        0       0       0
#  NTISTI
        0
#      IG   TEMPI
        0     0.0
END
@endverbatim
 */
void io::In_Parameter::read_INITIALISE(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read INITIALISE");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["INITIALISE"];

  if (!buffer.size()) {
    io::messages.add("no INITIALISE block", "In_Parameter", io::message::error);
    return;
  }

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  block_read.insert("INITIALISE");

  int ntivel, ntishk, ntinht, ntinhb, ntishi, ntirtc, nticom, ntisti;
  _lineStream >> ntivel >> ntishk >> ntinht >> ntinhb
          >> ntishi >> ntirtc >> nticom
          >> ntisti
          >> param.start.ig >> param.start.tempi;

  if (_lineStream.fail())
    io::messages.add("bad line in INITIALISE block",
          "In_Parameter", io::message::error);

  // generation of initial velocities
  switch (ntivel) {
    case 0: param.start.generate_velocities = false;
      break;
    case 1: param.start.generate_velocities = true;
      break;
    default: io::messages.add("INITIALISE block: NTIVEL must be 0 or 1",
              "In_Parameter", io::message::error);
  }

  // controls initial SHAKE
  switch (ntishk) {
    case 0: // no initial SHAKE
      param.start.shake_pos = false;
      param.start.shake_vel = false;
      break;
    case 1: // SHAKE coordinates
      param.start.shake_pos = true;
      param.start.shake_vel = false;
      break;
    case 2: // SHAKE velocities
      param.start.shake_pos = false;
      param.start.shake_vel = true;
      break;
    case 3: // SHAKE coordinates & velocities
      param.start.shake_pos = true;
      param.start.shake_vel = true;
      break;
    default: io::messages.add("INITIALISE block: NTISHK must be 0 to 3",
              "In_Parameter", io::message::error);
  }

  // controls reading of Nose-Hoover chain variables.
  switch (ntinht) {
    case 0: param.start.read_nosehoover_chains = true;
      break;
    case 1: param.start.read_nosehoover_chains = false;
      break; // reset them
    default: io::messages.add("INITIALISE block: NTINHT must be 0 or 1",
              "In_Parameter", io::message::error);
  }

  // controls reading of Nose-Hoover chain barostat variables: not implemented.
  switch (ntinhb) {
    case 0: param.start.read_nosehoover_barostat = true;
      break;
    case 1: param.start.read_nosehoover_barostat = false;
      break; // reset them
    default: io::messages.add("INITIALISE block: NTINHB must be 0 or 1",
              "In_Parameter", io::message::error);
  }

  switch (ntishi) {
    case 0: param.start.read_lattice_shifts = true;
      break;
    case 1: param.start.read_lattice_shifts = false;
      break;
    default: io::messages.add("INITIALISE block: NTISHI must be 0 or 1",
              "In_Parameter", io::message::error);
  }

  // controls reading of restart data for roto-translational constraints
  switch (ntirtc) {
    case 0: param.start.read_rottrans = true;
      break;
    case 1: param.start.read_rottrans = false;
      break;
    default: io::messages.add("INITIALISE block: NTIRTC must be 0 or 1",
              "In_Parameter", io::message::error);
  }

  // controls removal of COM translation and rotation.
  switch (nticom) {
    case 0:
      param.start.remove_com_rotation = false;
      param.start.remove_com_translation = false;
      break;
    case 1:
      param.start.remove_com_rotation = false;
      param.start.remove_com_translation = true;
      break;
    case 2:
      param.start.remove_com_rotation = true;
      param.start.remove_com_translation = true;
      break;
    default: io::messages.add("INITIALISE block: NTICOM must be 0 to 2",
              "In_Parameter", io::message::error);
  }

  // controls reading of stochastic integrals
  switch (ntisti) {
    case 0:
      param.stochastic.generate_integral = false;
      break;
    case 1:
      param.stochastic.generate_integral = true;
      break;
    default: io::messages.add("INITIALISE block: NTISTI must be 0 or 1",
              "In_Parameter", io::message::error);
  }

  if (param.start.tempi < 0)
    io::messages.add("Illegal value for TEMPI in INITIALISE block (>=0)",
          "In_Parameter", io::message::error);
}

/**
 * @section comtransrot COMTRANSROT block
 * @verbatim
COMTRANSROT
#    NSCM : controls system centre-of-mass (com) motion removal
#           0: no com motion removal (default)
#         < 0: com translation and rotation are removed every abs(NSCM)
#              steps.
#         > 0: com tranlsation is removed every NSCM steps.
#     NSCM
         0
END
@endverbatim
 */
void io::In_Parameter::read_COMTRANSROT(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read COMTRANSROT");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "reading COMTRANSROT block");
  buffer = m_block["COMTRANSROT"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("COMTRANSROT");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
  int nscm;

  _lineStream >> nscm;

  if (_lineStream.fail())
    io::messages.add("bad line in COMTRANSROT block",
          "In_Parameter", io::message::error);

  if (nscm > 0) {
    param.centreofmass.skip_step = nscm;
    param.centreofmass.remove_rot = false;
    param.centreofmass.remove_trans = true;
  } else if (nscm < 0) {
    param.centreofmass.skip_step = -nscm;
    param.centreofmass.remove_rot = true;
    param.centreofmass.remove_trans = true;
  } else { // nscm == 0;
    param.centreofmass.skip_step = 0;
    param.centreofmass.remove_rot = false;
    param.centreofmass.remove_trans = false;
  }
}

/**
 * @section hoomd HOOMD block (optional)
 * @verbatim
HOOMD
#       PROCESSOR: cpu gpus
#
#       PROCESSOR  
    gpus
END
@endverbatim
 */
void io::In_Parameter::read_HOOMD(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read HOOMD");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "hoomd block");

  // try a HOOMD
  buffer = m_block["HOOMD"];
  if (buffer.size()) {
    block_read.insert("HOOMD");

    std::string s1;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    _lineStream >> s1;

    if (_lineStream.fail()) {
      io::messages.add("bad line in HOOMD block",
              "In_Parameter", io::message::error);
    }

#ifndef HAVE_HOOMD
    io::messages.add("HOOMD block not supported. Recompile with --with-hoomd=DIR", "In_Parameter", io::message::error);
#else
    if (s1 == "cpu") {
      param.hoomd.processor = simulation::cpu;
    } else if (s1 == "gpus") {
      param.hoomd.processor = simulation::gpus;
    } else {
      io::messages.add("HOOMD block: wrong processor chosen",
              "In_Parameter", io::message::error);
    }
#endif
  }
}

/**
 * @section pairlist PAIRLIST block
 * @verbatim
PAIRLIST
#       ALGORITHM: standard(0) (gromos96 like pairlist)
#                  grid(1) (md++ grid pairlist)
#                  grid_cell(2) (creates a mask)
#       SIZE:      grid cell size (or auto = 0.5 * RCUTP)
#       TYPE:      chargegoup(0) (chargegroup based cutoff)
#                  atomic(1) (atom based cutoff)
#
#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE
        grid            5       0.8     1.4     auto    chargegroup
#
END
@endverbatim
 */
void io::In_Parameter::read_PAIRLIST(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read PAIRLIST");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "pairlist block");

  // try a PAIRLIST
  buffer = m_block["PAIRLIST"];
  if (buffer.size()) {
    block_read.insert("PAIRLIST");

    std::string s1, s2, s3;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    _lineStream >> s1
            >> param.pairlist.skip_step
            >> param.pairlist.cutoff_short
            >> param.pairlist.cutoff_long
            >> s2
            >> s3;

    if (_lineStream.fail()) {
      io::messages.add("bad line in PAIRLIST block",
              "In_Parameter", io::message::error);
    }

    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);
    std::transform(s2.begin(), s2.end(), s2.begin(), tolower);
    std::transform(s3.begin(), s3.end(), s3.begin(), tolower);

    if (s1 == "grid" || s1 == "1") param.pairlist.grid = 1;
    else if (s1 == "standard" || s1 == "0") param.pairlist.grid = 0;
    else if (s1 == "grid_cell" || s1 == "2") {
      param.pairlist.grid = 2;
    }
    else {
      io::messages.add("PAIRLIST block: wrong pairlist algorithm chosen",
              "In_Parameter", io::message::error);
      param.pairlist.grid = false;
    }

    if (param.pairlist.grid) {
      if (s2 == "auto")
        param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
      else {
        std::istringstream css;
        css.str(s2);
        css >> param.pairlist.grid_size;
        if (!param.pairlist.grid_size)
          param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
        // param.pairlist.grid_size = atof(s2.c_str());
        if (css.fail()) {
          io::messages.add("PAIRLIST block: wrong pairlist grid size chosen",
                  "In_Parameter", io::message::error);
          param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
        }
      }
    } else param.pairlist.grid_size = 0;

    if (s3 == "atomic" || s3 == "1") param.pairlist.atomic_cutoff = true;
    else if (s3 == "chargegroup" || s3 == "0") param.pairlist.atomic_cutoff = false;
    else {
      io::messages.add("PAIRLIST block: wrong cutoff type chosen "
              "(allowed: atomic(1), chargegroup(0)",
              "In_Parameter", io::message::error);
      param.pairlist.atomic_cutoff = false;
    }
  }

  if (param.pairlist.grid && param.pairlist.grid_size <= 0)
    io::messages.add("PAIRLIST block: Illegal value for grid size (>0)",
          "In_Parameter", io::message::error);
  if (param.pairlist.cutoff_short < 0) {
    io::messages.add("PAIRLIST block: Illegal value for short range cutoff (>0)",
            "In_Parameter", io::message::error);
  }
  if (param.pairlist.cutoff_long < param.pairlist.cutoff_short) {
    io::messages.add("PAIRLIST block: Illegal value for long range cutoff (>=RCUTP)",
            "In_Parameter", io::message::error);
  }
}

/**
 * @section cgrain CGRAIN block
 * @verbatim
CGRAIN
# NTCGRAN 0..3 coarse grain selection
#         0: atomistic (off)
#         1: coarse-grained using MARTINI model (on)
#         2: coarse-grained using GROMOS model (on)
#         3: mixed-grained using GROMOS model (on)
#     EPS >= 0.0 dielectric constant for coarse grained coulombic interaction
#    EPSM >= 0.0 dielectric constant for mixed CG-FG coulombic interaction
# NTCGRAN     EPS     EPSM
        1      20        1 
END
@endverbatim
 */
void io::In_Parameter::read_CGRAIN(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read CGRAIN");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["CGRAIN"];

  if (buffer.size()) {

    block_read.insert("CGRAIN");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int ntcgran;
    _lineStream >> ntcgran >> param.cgrain.EPS >> param.cgrain.EPSM;

    if (_lineStream.fail())
      io::messages.add("bad line in CGRAIN block",
            "In_Parameter", io::message::error);

    switch (ntcgran) {
      case 0:
        param.cgrain.level = 0;
        break;
      case 1:
        param.cgrain.level = 1;
        param.force.interaction_function = simulation::cgrain_func;
        break;
      case 2:
        param.cgrain.level = 2;
        param.force.interaction_function = simulation::cggromos_func;
        break;
      case 3:
        param.cgrain.level = 3;
        param.force.interaction_function = simulation::cggromos_func;
        break;
      default:
        param.cgrain.level = 0;
        io::messages.add("CGRAIN block: NTCGRAN must be 0 to 3.",
                "In_Parameter", io::message::error);
    }
    DEBUG(6, "coarse graining level = " << param.cgrain.level);

    if (param.cgrain.EPS < 0)
      io::messages.add("CGRAIN block: EPS must be >= 0.0.",
            "In_Parameter", io::message::error);
    if (param.cgrain.EPSM < 0)
      io::messages.add("CGRAIN block: EPSM must be >= 0.0.",
            "In_Parameter", io::message::error);
  }
}

/**
 * @section multibath MULTIBATH block
 * @verbatim
MULTIBATH
# ALGORITHM: temperature coupling algorithm
#		weak-coupling(0)
#		nose-hoover(1)
#		nose-hoover-chains(2)	num
#		(where num is the number of chains to use)
#   ALGORITHM           NUM
    nose-hoover-chains	3
#   NBATHS: number of baths
    2
#   TEMP0  TAU
    300    0.10
    300    0.10
#   DOFSET: number of different couplings
    1
#   LAST   COM-BATH  IR-BATH
    60     1         2
#   (this couples the first 60 atoms com motion to bath 1 and
#    the internal / rotational motion to bath 2)
END
@endverbatim
 */
void io::In_Parameter::read_MULTIBATH(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read MULTIBATH");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "TEMPERATURE COUPLING block");

  param.multibath.couple = false;

  // is there a MULTIBATH block
  buffer = m_block["MULTIBATH"];

  if (buffer.size()) {

    block_read.insert("MULTIBATH");

    param.multibath.found_multibath = true;
    param.multibath.found_tcouple = false;

    DEBUG(11, "MULTIBATH present");
    /*io::messages.add("using MULTIBATH block",
             "In_Parameter", io::message::notice);*/

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    std::string alg;

    // the algorithm
    _lineStream >> alg;
    std::transform(alg.begin(), alg.end(), alg.begin(), tolower);

    if (alg == "weak-coupling")
      param.multibath.nosehoover = 0;
    else if (alg == "nose-hoover")
      param.multibath.nosehoover = 1;
    else if (alg == "nose-hoover-chains")
      param.multibath.nosehoover = 2;
    else {
      std::stringstream s(alg);
      if (!(s >> param.multibath.nosehoover) ||
              param.multibath.nosehoover < 0 || param.multibath.nosehoover > 2) {
        io::messages.add("MUTLIBATH block: algorithm not understood",
                "In_Parameter", io::message::error);

        param.multibath.nosehoover = 0;
        return;
      }
    }

    if (param.multibath.nosehoover == 2) {
      // read in the number of chains
      int num;
      _lineStream >> num;

      if (num < 2) {
        io::messages.add("MULTIBATH block: wrong number of Nose-Hoover chains",
                "In_Parameter", io::message::error);
        param.multibath.nosehoover = 0;
        return;
      }
      param.multibath.nosehoover = num;
    }

    int num;
    unsigned int last;
    unsigned int com_bath, ir_bath;
    double temp, tau;

    // the baths
    _lineStream >> num;

    for (int i = 0; i < num; ++i) {
      _lineStream >> temp >> tau;

      if (temp < 0.0 || (tau <= 0.0 && tau != -1)) {
        io::messages.add("MULTIBATH block: illegal value for temp or tau",
                "In_Parameter", io::message::error);
      }

      param.multibath.multibath.add_bath(temp, tau);
      if (tau != -1) param.multibath.couple = true;
    }

    if (_lineStream.fail()) {
      io::messages.add("bad line in MULTIBATH block",
              "In_Parameter", io::message::error);
    }

    // now the ranges
    _lineStream >> num;

    if (param.multibath.multibath.size() == 0 &&
            num > 0) {

      io::messages.add("MULTIBATH block: no baths but coupling groups specified",
              "In_Parameter", io::message::error);
      num = 0;
    }

    for (int i = 0; i < num; ++i) {
      _lineStream >> last >> com_bath >> ir_bath;
      // let it figure out the last molecule on its own

      if (last < 1 || com_bath < 1 || ir_bath < 1) {
        io::messages.add("MULTIBATH block: range parameter < 1",
                "In_Parameter", io::message::error);
        if (last < 1) last = 1;
        if (com_bath < 1) com_bath = 1;
        if (ir_bath < 1) ir_bath = 1;
      }

      if (com_bath > param.multibath.multibath.size() ||
              ir_bath > param.multibath.multibath.size()) {
        io::messages.add("MULTIBATH block: ir bath or com bath index too large",
                "In_Parameter", io::message::error);
        if (com_bath > param.multibath.multibath.size()) com_bath = param.multibath.multibath.size();
        if (ir_bath > param.multibath.multibath.size()) ir_bath = param.multibath.multibath.size();
      }

      param.multibath.multibath.add_bath_index(last - 1, 0, com_bath - 1, ir_bath - 1);
    }

    if (_lineStream.fail()) {
      io::messages.add("bad line in MULTIBATH block",
              "In_Parameter", io::message::error);
    }

  }
} // MULTIBATH

/**
 * @section positionres POSITIONRES block
 * @verbatim
POSITIONRES
#    NTPOR 0..3 controls atom positions re(con)straining.
#          0: no position re(con)straints.
#          1: restraining with force constant CPOR
#          2: restraining with force constant CPOR wieghted by
#             atomic B-factors
#          3: constraining
#   NTPORB 0,1 controls reading of reference positions and
#              B-factors
#          0: read reference positions from startup file.
#          1: read reference positions and B-factors from
#             special file
#   NTPORS 0,1 controls scaling of reference positions upon 
#              pressure scaling
#          0: do not scale reference positions
#          1: scale reference positions
#     CPOR >= 0 position restraining force constant
#
#   NTPOR  NTPORB  NTPORS    CPOR
        0       0       0   2.5E4
END
@endverbatim
 */
void io::In_Parameter::read_POSITIONRES(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read POSITIONRES");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "positionres block");
  buffer = m_block["POSITIONRES"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("POSITIONRES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int ntpor, ntpors, read;
  _lineStream >> ntpor
          >> read
          >> ntpors
          >> param.posrest.force_constant;

  if (_lineStream.fail())
    io::messages.add("bad line in POSITIONRES block",
          "In_Parameter", io::message::error);

  switch (ntpor) {
    case 0:
      param.posrest.posrest = simulation::posrest_off;
      break;
    case 1:
      param.posrest.posrest = simulation::posrest_on;
      break;
    case 2:
      param.posrest.posrest = simulation::posrest_bfactor;
      break;
    case 3:
      param.posrest.posrest = simulation::posrest_const;
      break;
    default:
      io::messages.add("POSITIONRES block: NTPOR must be 0 to 3.",
              "In_Parameter", io::message::error);
  }

  switch (read) {
    case 0:
      param.posrest.read = false;
      break;
    case 1:
      param.posrest.read = true;
      break;
    default:
      param.posrest.read = false;
      io::messages.add("POSITIONRES block: NTPORS must be 0 or 1.",
              "In_Parameter", io::message::error);
  }

  switch (ntpors) {
    case 0:
      param.posrest.scale_reference_positions = false;
      break;
    case 1:
      param.posrest.scale_reference_positions = true;
      break;
    default:
      param.posrest.scale_reference_positions = false;
      io::messages.add("POSITIONRES block: NTPORS must be 0 or 1.",
              "In_Parameter", io::message::error);
  }

  if (param.posrest.force_constant < 0.0)
    io::messages.add("POSITIONRES block: Illegal value for CPOR.",
          "In_Parameter", io::message::error);

} // POSITIONRES

/**
 * @section xrayres XRAYRES block
 * @verbatim
XRAYRES
#    NTXR   -2: time-averaged electron density restraints
#           -1: instantaneous electron density restraints
#            0: no xray restraints.
#            1: instantaneous structure factor restraints
#            2: time-averaged structure factor restraints
#            3: biquadratic/timeaveraged structure factor restraints
#    NTXLE   0: do not perform local elevation
#            1: do perform local elevation
#    CXR     >= 0 xray restraining force constant
#    NTWXR   >= 0 write xray data to output file
#            0: don't write xray data
#            > 0 write every NTPXRth step
#    NTWDE   0..3 write density-maps
#            0: write nothing
#            1: write electron densitiy map
#            2: write asymmetric-unit-only electron densitiy map
#            3: write both
#    NTWXM   >= 0 write every NTWXMth step electron density map(s) to external file
#    CXTAU   >=0 xray time-average restraining memory-constant
#    RDAVG   0/1 read sf-timeaverages (from job to job)
#
#   NTXR   NTXLE     CXR   NTWXR   NTWDE   NTWXM   CXTAU   RDAVG
       0       0     0.0       0       0      0      0.0       0
END
@endverbatim
 */
void io::In_Parameter::read_XRAYRES(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read XRAYRES");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "xrayres block");
  buffer = m_block["XRAYRES"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("XRAYRES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int ntxr, ntxle;
  _lineStream >> ntxr >> ntxle
          >> param.xrayrest.force_constant
          >> param.xrayrest.write
          >> param.xrayrest.writedensity
          >> param.xrayrest.writexmap
          >> param.xrayrest.tau
          >> param.xrayrest.readavg;

  if (_lineStream.fail())
    io::messages.add("bad line in XRAYRES block",
          "In_Parameter", io::message::error);

  if (ntxr < 0)
    param.xrayrest.mode = simulation::xrayrest_mode_electron_density;
  else
    param.xrayrest.mode = simulation::xrayrest_mode_structure_factor;

  switch (abs(ntxr)) {
    case 0:
      param.xrayrest.xrayrest = simulation::xrayrest_off;
      break;
    case 1:
      param.xrayrest.xrayrest = simulation::xrayrest_inst;
      break;
    case 2:
      param.xrayrest.xrayrest = simulation::xrayrest_avg;
      break;
    case 3:
      param.xrayrest.xrayrest = simulation::xrayrest_biq;
      break;
    default:
      io::messages.add("XRAYRES block: NTXR must be -2 to 3.",
              "In_Parameter", io::message::error);
  }

  if (ntxr == -3) {
    io::messages.add("XRAYRES block: NTXR must be -2 to 3.",
            "In_Parameter", io::message::error);
  }

  if (param.xrayrest.xrayrest == simulation::xrayrest_off) {
    // abort reading of rest
    return;
  }

  switch (ntxle) {
    case 0:
      param.xrayrest.local_elevation = false;
      break;
    case 1:
      param.xrayrest.local_elevation = true;
      break;
    default:
      io::messages.add("XRAYRES block: NTXLE must be 0 or 1.",
              "In_Parameter", io::message::error);

  }

  if (param.xrayrest.local_elevation) {
    if (param.xrayrest.mode != simulation::xrayrest_mode_electron_density ||
            param.xrayrest.xrayrest == simulation::xrayrest_biq) {
      io::messages.add("XRAYRES block: NTXLE 1 requires NTXR -2 or -1.",
              "In_Parameter", io::message::error);
    }
  }

  if (param.xrayrest.force_constant < 0.0) {
    io::messages.add("XRAYRES block: Illegal value for CXR.",
            "In_Parameter", io::message::error);
  }
  if ((ntxr == 2 || ntxr == 3) && param.xrayrest.force_constant <= 0.0) {
    io::messages.add("XRAYRES block: Illegal value for CXTAU",
            "In_Parameter", io::message::error);
  }
  if (param.xrayrest.writedensity > 3) {
    io::messages.add("XRAYRES block: Illegal value for NTWDE (0..3)",
            "In_Parameter", io::message::error);
  }
} // XRAYRES

/**
 * @section distanceres DISTANCERES block
 * @verbatim
DISTANCERES
#   NTDIR -2..2 controls distance restraining
#         0: no distrance restraining (default)
#         1: instantaneous, using force constant CDIR
#         2: instantaneous, using force constant CDIR x W0
#        -1: time-averaged, using force constant CDIR
#        -2: time-averaged, using force constant CDIR x W0
#  NTDIRA 0,1 controls values for initial distance averages
#         0: generate initial averages
#         1: read from configuration
#    CDIR >= 0.0 force constant for distance restraining
#    DIR0 > 0.0 distance offset in restraining function
#  TAUDIR >= 0.0 coupling time for time averaging
# FORCESCALE 0..2 controls approximation of force scaling
#         0: approximate d<r>/dr = 1
#         1: approximate d<r>/dr = (1.0 - exp(-Dt/tau))
#         2: use d<r>/dr = (1.0 - exp(-Dt/tau))*(<r>/r)^4
#    VDIR 0,1 controls contribution to virial
#         0: no contribution
#         1: distance restraints contribute to virial
#  NTWDIR >= 0 write every NTWDIRth step dist. restr. information to external file
#   NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  FORCESCALE  VDIR  NTWDIR
        0       0     0.0     1.0     0.0           0     0       0
END
@endverbatim
 */
void io::In_Parameter::read_DISTANCERES(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read DISTANCERES");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "distanceres block");
  buffer = m_block["DISTANCERES"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("DISTANCERES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int ntdira;
  _lineStream >> param.distanceres.distanceres
	      >> ntdira
	      >> param.distanceres.K
	      >> param.distanceres.r_linear
	      >> param.distanceres.tau
	      >> param.distanceres.forcescale
	      >> param.distanceres.virial
	      >> param.distanceres.write;

  if (_lineStream.fail())
    io::messages.add("bad line in DISTANCERES block",
          "In_Parameter", io::message::error);


  if (param.distanceres.distanceres < -2 || param.distanceres.distanceres > 2) {
    io::messages.add("DISTANCERES block: NTDIR must be -2..2.",
            "In_Parameter", io::message::error);
  }

  switch (ntdira) {
    case 0: param.distanceres.read = false;
      break;
    case 1: param.distanceres.read = true;
      break;
    default: param.distanceres.read = false;
      io::messages.add("DISTANCERES block: NTDIRA must be 0 or 1.",
              "In_Parameter", io::message::error);
  }

  if (param.distanceres.r_linear <=  0.0) {
    io::messages.add("DISTANCERES block: DIR0 must be > 0.0.",
            "In_Parameter", io::message::error);
  }

  if (param.distanceres.tau < 0.0) {
    io::messages.add("DISTANCERES block: TAUDIR must be >= 0.0.",
            "In_Parameter", io::message::error);
  }

  if (param.distanceres.K < 0) {
    io::messages.add("DISTANCERES block: CDIR must be >= 0.0.",
            "In_Parameter", io::message::error);
  }

  if (param.distanceres.read > 0 && param.distanceres.distanceres > 0) {
    io::messages.add("DISTANCERES block: NTDIRA > 0 but NTDIR > 0 - DISRESEXPAVE ignored",
            "In_Parameter", io::message::warning);
  }

  if(param.distanceres.virial !=0 && param.distanceres.virial != 1){
    io::messages.add("DISTANCERES block: VDIR must be 0 or 1",
		     "In_Parameter", io::message::error);
  }
  if(param.distanceres.forcescale < 0 || param.distanceres.forcescale > 2){
    io::messages.add("DISTANCERES block: FORCESCALE must be 0, 1 or 2",
		     "In_Parameter", io::message::error);
  }
} // DISTANCERES

/**
 * @section distancefield DISTANCEFIELD block
 * @verbatim
DISTANCEFIELD
#   NTDFR 0,1         controls distance field restraining
#         0: no distance field restraining
#         1: apply distance field restraining
#   GRID  > 0.0       grid size for distance field
#   PROTEINOFFSET > 0 penalty for distances through the host
#   PROTEINCUTOFF > 0 distance to protein atoms to be considered inside
#   UPDATE > 0        update frequency for grid
#   RL >= 0           linearize forces for distances larger than RL
#   SMOOTH >= 0       smoothen the protein boundary after grid construction
#                     by SMOOTH layers
#   NTWDF >= 0        write every NTWDF step disfield information to external file 
#   PRINTGRID = 0,1   write grid to final configuration file      
#
#   NTDFR  
        1
#    GRID   PROTEINOFFSET  PROTEINCUTOFF
      0.2   15             0.2         
#  UPDATE   SMOOTH   RL    NTWDF   PRINTGRID
      100   1        1.0      50           0
END  
@endverbatim
 */
void io::In_Parameter::read_DISTANCEFIELD(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read DISTANCEFIELD");

  std::vector<std::string> buffer;
  std::string s;

  DEBUG(10, "distancefield block");
  buffer = m_block["DISTANCEFIELD"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("DISTANCEFIELD");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int printgrid = -1;
  _lineStream >> param.distancefield.distancefield
	      >> param.distancefield.grid
	      >> param.distancefield.proteinoffset
	      >> param.distancefield.proteincutoff
	      >> param.distancefield.update
	      >> param.distancefield.smooth
              >> param.distancefield.r_l
	      >> param.distancefield.write
	      >> printgrid;

  if (_lineStream.fail())
    io::messages.add("bad line in DISTANCEFIELD block",
          "In_Parameter", io::message::error);


  if (param.distancefield.distancefield != 0 && 
      param.distancefield.distancefield != 1) {
    io::messages.add("DISTANCEFIELD block: NTDRF must be 0,1.",
		     "In_Parameter", io::message::error);
  }
  if (param.distancefield.grid <= 0.0) {
    io::messages.add("DISTANCEFIELD block: GRID must be >= 0.0.",
		     "In_Parameter", io::message::error);
  }
  if (param.distancefield.proteinoffset < 0.0) {
    io::messages.add("DISTANCEFIELD block: PROTEINOFFSET must be >= 0.0.",
		     "In_Parameter", io::message::error);
  }
  if (param.distancefield.proteincutoff < 0.0) {
    io::messages.add("DISTANCEFIELD block: PROTEINCUTOFFmust be >= 0.0.",
		     "In_Parameter", io::message::error);
  }
  if (param.distancefield.update < 0) {
    io::messages.add("DISTANCEFIELD block: UPDATE must be > 0.",
		     "In_Parameter", io::message::error);
  }
  if (param.distancefield.smooth < 0) {
    io::messages.add("DISTANCEFIELD block: SMOOTH must be >= 0.",
		     "In_Parameter", io::message::error);
  }
  if (printgrid!=0 && printgrid!=1) {
    io::messages.add("DISTANCEFIELD block: PRINTGRID must be 0 or 1.",
		     "In_Parameter", io::message::error);
  }
  if(printgrid==1) param.distancefield.printgrid = true;
  else param.distancefield.printgrid = false;
  
} // DISTANCEFIELD

/**
 * @section dihedralres DIHEDRALRES block
 * @verbatim
DIHEDRALRES
# NTDLR   0...3 controls dihedral-angle restraining and constraining 
#         0:    off [default]
#         1:    dihedral restraining using CDLR
#         2:    dihedral restraining using CDLR * WDLR
#         3:    dihedral constraining
#    
# CDLR    >=0.0 force constant for dihedral restraining
# PHILIN  >0.0  deviation after which the potential energy function is linearized 
#
# NTDLR  CDLR      PHILIN
  1      100.0     180.0
END
@endverbatim
 */
void io::In_Parameter::read_DIHEDRALRES(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read DIHEDRALRES");

  std::vector<std::string> buffer;
  std::string s;
  double phi_lin;

  DEBUG(10, "DIHEDRALRES block");
  buffer = m_block["DIHEDRALRES"];

  if (!buffer.size()) {
    return;
  }

  block_read.insert("DIHEDRALRES");

  _lineStream.clear();
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  /*
   # NTDLR   0...3 controls dihedral-angle restraining and constraining 
   #         0:    off [default]
   #         1:    dihedral restraining using CDLR
   #         2:    dihedral restraining using CDLR * WDLR
   #         3:    dihedral constraining
   */

  int dihrest;
  _lineStream >> dihrest
          >> param.dihrest.K
          >> phi_lin;

  if (_lineStream.fail())
    io::messages.add("bad line in DIHEDRALRES block",
          "In_Parameter", io::message::error);

  switch (dihrest) {
    case 0:
      param.dihrest.dihrest = simulation::dihedral_restr_off;
      break;
    case 1:
      param.dihrest.dihrest = simulation::dihedral_restr_inst;
      break;
    case 2:
      param.dihrest.dihrest = simulation::dihedral_restr_inst_weighted;
      break;
    case 3:
      param.dihrest.dihrest = simulation::dihedral_constr;
      param.setDevelop("Dihedral constraining is under development."); 
      break;
    default:
      io::messages.add("DIHEDRALRES block: NTDLR must be 0...3.",
              "In_Parameter", io::message::error);
  }

  if (phi_lin <= 0.0)
    io::messages.add("DIHEDRALRES block: Illegal value for force constant (>=0)",
    "In_Parameter", io::message::error);

  param.dihrest.phi_lin = phi_lin * 2 * math::Pi / 360;

  if (param.dihrest.K < 0)
    io::messages.add("DIHEDRALRES block: Illegal value for force constant (>=0)",
          "In_Parameter", io::message::error);

  if (param.dihrest.dihrest == simulation::dihedral_constr) {
    if (param.constraint.ntc == 1 && param.constraint.solute.algorithm == simulation::constr_off)
      param.constraint.solute.algorithm = simulation::constr_shake;

    if (param.constraint.solute.algorithm != simulation::constr_shake) {
      io::messages.add("DIHEDRALRES block: needs SHAKE as (solute) constraints algorithm",
              "In_Parameter",
              io::message::error);
    }
  }

} // DIHEDRALRES

/**
 * @section jval JVALUERES block
 * @verbatim
JVALUERES
# NTJVR    -3..2
#          -3                biquadratic using CJVR * WJVR
#          -2                time-averaged using CJVR * WJVR
#          -1                time-avaraged using CJVR
#           0                no J-value restraints [default]
#           1                instantaneous using CJVR
#           2                instantaneous using CJVR * WJVR
# NTJVRA    0                controls reading of averages from startup file
#           0                start from initial values of J0 [default]
#           1                read time averages from startup file (for continuation time-averaged run)
# CJVR   >= 0                J-value restraining force constant 
#                            (weighted by individual WJVR)
# TAUJVR >= 0                coupling time for time-averaging
# NJVRTARS  0,1              omits or includes force scaling by memory decay factor in case of time-averaging
#           0                omit factor (set (1 - exp(-Dt/tau)) = 1)
#           1                scale force by (1 - exp(-Dt/tau))                                      
# NJVRBIQW  0..2             controls weights (X in Eq. 19 of MD98.17) of the two terms in biquadratic restraining
#           0                X = 1
#           1                X = (1 - exp(-Dt/tau)) 
#           2                X = 0
# LE        0,1              local elevation restraining
#           0                local elevation off [default]
#           1                local elevation on
# NGRID   > 0                number of grid points in local elevation restraining
# DELTA  >= 0.0              no elevation of potential if J is within DELTA of J0
# NTWJV  >= 0                write J-value averages and LE grid to special trajectory
#           0                don't write [default]
#         > 0                write every NTWJVth step
#
#       NTJVR  NTJVRA  CJVR   TAUJVR  NJVRTARS   NJVRBIQW   LE    NGRID   DELTA  NTWJV
           -3  0       10.0      5.0     0          0       1       16     0.5      0
END
@endverbatim
 */
void io::In_Parameter::read_JVALUERES(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read JVALUERES");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["JVALUERES"];
  if (buffer.size()) {

    block_read.insert("JVALUERES");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int ntjvr;
    _lineStream >> ntjvr // NTJVR
            >> param.jvalue.read_av // NTJVRA
            >> param.jvalue.K // CJVR
            >> param.jvalue.tau // TAUJVR
            >> param.jvalue.tarfscale // NJVRTARS
            >> param.jvalue.biqweight // NJVRBIQW
            >> param.jvalue.le // LE
            >> param.jvalue.ngrid // NGRID
            >> param.jvalue.delta // DELTA
            >> param.jvalue.write; // NTJVW

    if (_lineStream.fail())
      io::messages.add("bad line in JVALUERES block",
            "In_Parameter", io::message::error);

    switch (ntjvr) {
      case -3:
        param.jvalue.mode = simulation::jvalue_restr_biq_weighted;
        break;
      case -2:
        param.jvalue.mode = simulation::jvalue_restr_av_weighted;
        break;
      case -1:
        param.jvalue.mode = simulation::jvalue_restr_av;
        break;
      case 0:
        param.jvalue.mode = simulation::jvalue_restr_off;
        break;
      case 1:
        param.jvalue.mode = simulation::jvalue_restr_inst;
        break;
      case 2:
        param.jvalue.mode = simulation::jvalue_restr_inst_weighted;
        break;
      default:
        io::messages.add("JVALUERES block: NTJVR must be -3..2.",
                "In_Parameter", io::message::error);
    }

    if (param.jvalue.tau < 0 ||
            (param.jvalue.tau == 0 && (param.jvalue.mode != simulation::jvalue_restr_off ||
            param.jvalue.mode != simulation::jvalue_restr_inst ||
            param.jvalue.mode != simulation::jvalue_restr_inst_weighted))) {
      io::messages.add("JVALUERES block: bad value for TAU, should be > 0.0",
              "In_Parameter", io::message::error);
    }
    if (param.jvalue.tarfscale < 0 || param.jvalue.tarfscale > 1){
        io::messages.add("JVALUERES block: NJVRTARS must be 0 or 1",
                        "In_Parameter", io::message::error);
    }
    if (param.jvalue.biqweight < 0 || param.jvalue.biqweight > 2){
        io::messages.add("JVALUERES block: NJVRBIQW must be 0, 1 or 2",
                        "In_Parameter", io::message::error);
    }
    if (param.jvalue.mode != simulation::jvalue_restr_off && param.jvalue.K < 0.0) {
      io::messages.add("JVALUERES block: bad value for K in JVALUERES block,"
              "should be > 0.0",
              "In_Parameter", io::message::error);
    }
    if (param.jvalue.le > 0) {
      if (param.jvalue.ngrid < 1) {
        io::messages.add("JVALUERES block: bad value for NGRID in JVALUERES block, "
                "should be > 1",
                "In_Parameter", io::message::error);
      }
    }
    if (param.jvalue.read_av && (param.jvalue.mode != simulation::jvalue_restr_av &&
            param.jvalue.mode != simulation::jvalue_restr_av_weighted
            && param.jvalue.mode != simulation::jvalue_restr_biq_weighted
            && !param.jvalue.le)) {
      io::messages.add("JVALUERES block: Continuation only needed "
              "with time-averaging or LE.",
              "In_Parameter", io::message::error);
    }
  } // JVALUERES
} // JVALUE

/**
 * @section orderparamres ORDERPARAMRES block
 * @verbatim
ORDERPARAMRES
# NTOPR    -2..0
#          -2                time-averaged using COPR * WOPR
#          -1                time-avaraged using COPR
#           0                no order-parameter restraints [default]
# NTOPRA    0                controls reading of averages from startup file
#           0                start from initial values of S0 [default]
#           1                read time averages from startup file (for continuation time-averaged run)
# COPR   >= 0.0              order-parameter restraining force constant
#                            (weighted by individual WOPR)
# TAUOPR >= 0.0              coupling time for time-averaging
# UPDOPR  > 0                update average every UPDOPRth step
# NTWOP  >= 0                write order-parameter to special trajectory
#           0                don't write [default]
#         > 0                write every NTWOP step
#
#       NTOPR  NTOPRA  COPR   TAUOPR   UPDOPR  NTWOP
           -2  0       10.0      5.0        1      0
END
@endverbatim
 */
void io::In_Parameter::read_ORDERPARAMRES(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read ORDERPARAMRES");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["ORDERPARAMRES"];
  if (buffer.size()) {
   
    block_read.insert("ORDERPARAMRES");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int ntopr;
    _lineStream >> ntopr // NTJVR
            >> param.orderparamrest.read // NTOPRA
            >> param.orderparamrest.K // COPR
            >> param.orderparamrest.tau // TAUOPR
            >> param.orderparamrest.update_step // UPDOPR
            >> param.orderparamrest.write; // NTOPR

    if (_lineStream.fail())
      io::messages.add("bad line in ORDERPARAMRES block",
            "In_Parameter", io::message::error);

    switch (ntopr) {
      case -2:
        param.orderparamrest.orderparamrest = simulation::oparam_restr_av_weighted;
        break;
      case -1:
        param.orderparamrest.orderparamrest = simulation::oparam_restr_av;
        break;
      case 0:
        param.orderparamrest.orderparamrest = simulation::oparam_restr_off;
        break;
      case 1:
        param.orderparamrest.orderparamrest = simulation::oparam_restr_winav;
        break;
      case 2:
        param.orderparamrest.orderparamrest = simulation::oparam_restr_winav_weighted;
        break;
      default:
        io::messages.add("ORDERPARAMRES block: NTOPR must be -2..0.",
                "In_Parameter", io::message::error);
    }

    if (param.orderparamrest.orderparamrest != simulation::oparam_restr_off && param.orderparamrest.tau < 0.0) {
      io::messages.add("ORDERPARAMRES block: bad value for TAUOPR, should be >= 0.0",
              "In_Parameter", io::message::error);
    }
    if (param.orderparamrest.orderparamrest != simulation::oparam_restr_off && param.orderparamrest.K < 0.0) {
      io::messages.add("ORDERPARAMRES block: bad value for COPR, should be > 0.0",
              "In_Parameter", io::message::error);
    }
    
    if (param.orderparamrest.orderparamrest != simulation::oparam_restr_off && 
        param.orderparamrest.update_step == 0) {
      io::messages.add("ORDERPARAMRES block: bad value for UPDOPR, should be > 0",
              "In_Parameter", io::message::error);
    }

    if ((param.orderparamrest.orderparamrest == simulation::oparam_restr_av ||
        param.orderparamrest.orderparamrest == simulation::oparam_restr_av_weighted) &&
        param.orderparamrest.update_step != 1) {
      io::messages.add("ORDERPARAMRES block: bad value for UPDOPR, should be 1 for exponential averaging",
              "In_Parameter", io::message::error);
    }
    
    if (param.orderparamrest.orderparamrest == simulation::oparam_restr_winav ||
        param.orderparamrest.orderparamrest == simulation::oparam_restr_winav_weighted) {
      unsigned int window = int(param.orderparamrest.tau / param.step.dt);
      if (window / param.orderparamrest.update_step == 0) {
        io::messages.add("ORDERPARAMRES block: bad value for UPDOPR smaller than TAUOPR / DT.",
                "In_Parameter", io::message::error);
      }
    }
    
    
  } // ORDERPARAMRES
} // ORDERPARAMRES


/**
 * @section rdcres RDCRES block
 * @verbatim
RDCRES
# NTRDCR -4..2                 RDC restraining
#        -4:                   biquadratic using CRDCR * WRDCR
#        -3:                   biquadratic using CRDCR
#        -2:                   time averaged using CRDCR * WRDCR
#        -1:                   time averaged using CRDCR
#         0:                   no RDC restraints [default]
#         1:                   instantaneous using CRDCR
#         2:                   instantaneous using CRDCR * WRDCR
# NTRDCRA 0,1                  controls reading of average RDCs
#         0:                   take initial values RDC0 from the RDC restraint file
#         1:                   read time averages from initial coordinate file 
#                              (for continuation run)
#
# NTRDCT  0..2                 Type of alignment representation
#         0:                   cartesian magnetic field vectors
#         1:                   alignment tensor
#         2:                   spherical harmonics
# NTALR   0,1                  controls reading of values in the chosen representation
#         0:                   start from values given in RDC restraint file
#         1:                   read values from initial coordinate file (for continuation run)
#
# METHOD  0..3                 Method of updating the magnetic field vectors
#         0:                   Energy minimisation
#         1:                   Stochastic dynamics
#         2:                   Molecular dynamics
# EMGRAD  > 0.0                (METHOD = 0, EM) stop minimisation if gradient is below EMGRAD 
# EMDX0   > 0.0                (METHOD = 0, EM) initial step size
# EMNMAX  > 0                  (METHOD = 0, EM) maximum number of minimisation steps
# SDCFRIC >= 0.0               (METHOD = 1, SD) global friction coefficient gamma
# TEMP  >= 0.0                 temperature of stochastic bath (SD) and temperature used for initial velocities (MD, SD)
# DELTA   >= 0                 the flatbottom potential is 2 DELTA wide
# CRDCR   >= 0                 RDC restraining force constant
#                              (weighted by individual WRDCR)
# TAU     >= 0                 coupling time for time averaging
# NRDCRTARS 0,1                omits or includes force scaling by memory decay factor in case of time-averaging
#           0                  omit factor (set (1-exp(-dt/tau))=1 )
#           1                  scale force by (1-exp(-dt/tau))
# NRDCRBIQW 0..2               controls weights of the two terms in biquadratic restraining
#           0                  X = 1
#           1                  X = (1 - exp(-dt/tau))
#           2                  X = 0
# NTWRDC   >= 0                write output to special trajectory
#           0:                 don't write
#          >0:                 write every NTWRDCth step.
#
#      NTRDCR  NTRDCRA  NTRDCT  NTALR  METHOD  
            2        0       0      0       0
#      EMGRAD  EMDX0  EMNMAX  SDCFRIC    TEMP    DELTA  CRDCR  TAU   NRDCRBIQW   NTWRDC   NTWRDC 
        0.001   0.01    1000       20     300      0      1     1     0          0        10000
END
@endverbatim
 */

void io::In_Parameter::read_RDCRES(simulation::Parameter &param,
                                   std::ostream & os)
{
  DEBUG(8, "read RDCRES");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["RDCRES"];
  if (buffer.size()){

    //set at develop flag
    param.setDevelop("RDC restraining is under development!");

    block_read.insert("RDCRES");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin()+1, buffer.end()-1, s));

    int ntrdcr, ntrdct, method;
    _lineStream >> ntrdcr // NTRDCR
                >> param.rdc.read_av // NTRDCRA
                >> ntrdct    // NTRDCT
                >> param.rdc.read_align // NTALR
                >> method // METHOD
                >> param.rdc.emgradient // EMGRAD
                >> param.rdc.emstepsize // EMDX0
                >> param.rdc.emmaxiter // EMNMAX
                >> param.rdc.sdfric // SDCFRIC
                >> param.rdc.temp // TEMP
                >> param.rdc.delta // DELTA [ps^-1]
                >> param.rdc.K // CRDCR // [kJ*ps^2]
                >> param.rdc.tau // TAU // [ps]
                >> param.rdc.tAVfactor // NRDCRTARS
                >> param.rdc.biqfactor // NRDCRBIQW
                >> param.rdc.write; // NTWRDC

    if (_lineStream.fail()){
      io::messages.add("bad line in RDCRES block", "In_Parameter", io::message::error);
    }

    switch (ntrdct) {
      case 0:
        param.rdc.type = simulation::rdc_mf;
        break;
      case 1:
        param.rdc.type = simulation::rdc_t;
        break;
      case 2:
        param.rdc.type = simulation::rdc_sh;
        break;
      default:
        io::messages.add("RDCRES block: NTRDCT must be 0..2", "In_Parameter", io::message::error);
    }
    
    switch (method) {
      case 0:
        param.rdc.method = simulation::rdc_em;
        break;
      case 1:
        param.rdc.method = simulation::rdc_sd;
        break;
      case 2:
        param.rdc.method = simulation::rdc_md;
        break;
      default:
        io::messages.add("RDCRES block: METHOD must be 0..2", "In_Parameter", io::message::error);
    }

    if(method != simulation::rdc_em){
      io::messages.add("RDCRES block: Your choice of METHOD is legal. However, be aware, that rdc_em is faster and always a better choice except when you are experimenting.",
          "In_Parameter", io::message::warning);
    }
    if(ntrdct != simulation::rdc_t){
      io::messages.add("RDCRES block: Your choice of NTRDCT is legal. However, be aware, that rdc_t is faster and always a better choice except when you are experimenting.",
          "In_Parameter", io::message::warning);
    }

    switch (ntrdcr) {
      case -4:
        param.rdc.mode = simulation::rdc_restr_biq_weighted;
        param.rdc.K *= 10e20; // multiply with 10e20 ps^2  to roughly compensate for second quadratic term
        io::messages.add("Multyplying RDC force constant K with 10e20 ps^2 to compensate for second quadratic term", "In_Parameter", io::message::notice);
        break;
      case -3:
        param.rdc.mode = simulation::rdc_restr_biq;
        param.rdc.K *= 10e20; // multiply with 10e20 ps^2  to roughly compensate for second quadratic term
        io::messages.add("Multyplying RDC force constant K with 10e20 ps^2 to compensate for second quadratic term", "In_Parameter", io::message::notice);
        break;
      case -2:
        param.rdc.mode = simulation::rdc_restr_av_weighted;
        break;
      case -1:
        param.rdc.mode = simulation::rdc_restr_av;
        break;
      case 0:
        param.rdc.mode = simulation::rdc_restr_off;
        break;
      case 1:
        param.rdc.mode = simulation::rdc_restr_inst;
        break;
      case 2:
        param.rdc.mode = simulation::rdc_restr_inst_weighted;
        break;
      default:
        io::messages.add("RDCRES block: NTRDCR must be -4..2.", "In_Parameter", io::message::error);
    }

    if (param.rdc.mode != simulation::rdc_restr_off){
      if (param.rdc.tau <= 0.0 && ( param.rdc.mode != simulation::rdc_restr_inst && param.rdc.mode != simulation::rdc_restr_inst_weighted)){
        io::messages.add("RDCRES block: bad value for TAU, should be > 0.0",
                 "In_Parameter", io::message::error);
      }
      if (param.rdc.K < 0.0){
        io::messages.add("RDCRES block: bad value for CRDCR(K) in RDCRES block, should be > 0.0",
                 "In_Parameter", io::message::error);
      }
      if (param.rdc.read_av && (param.rdc.mode != simulation::rdc_restr_av
          && param.rdc.mode != simulation::rdc_restr_av_weighted
          && param.rdc.mode != simulation::rdc_restr_biq
          && param.rdc.mode != simulation::rdc_restr_biq_weighted )){
        io::messages.add("RDCRES block: Continuation only needed with averaging.",
                 "In_Parameter", io::message::error);
      }
    } // rdcs not switched off

    if(!(param.rdc.tAVfactor == 0 || param.rdc.tAVfactor == 1 )) 
        io::messages.add("RDCRES block: The only possible values for NRDCRTARS are 0 and 1.", "In_Parameter", io::message::error);

    if(!(param.rdc.biqfactor == 0 || param.rdc.biqfactor == 1 || param.rdc.biqfactor == 2)) 
        io::messages.add("RDCRES block: The only possible values for NRDCRBIQW are 0, 1 and 2.", "In_Parameter", io::message::error);

  } // RDCRES block
} // read RDCs



/**
 * @section perscale PERSCALE block
 * @verbatim
PERSCALE
# RESTYPE		special energy term to which periodic scaling should
#         		be applied
#	  off		don't apply periodic scaling
#	jrest		apply periodic scaling to J-value restraints
#
# parameters for periodic scaling of J-value restraints
# KDIH	>= 0		maximum scaling factor for dihedral angle potential
# KJ	>= 0		maximum scaling factor for J-Value restraint potential
# T	>= 0		period of cosine scaling function
# DIFF	>= 0		minimum deviation from target value at which to start
#            		scaling period
# RATIO	>= 0		minimum fraction of T that needs to be passed before
#            		starting a new scaling period
# READ	   0,1		controls reading of scaling parameters
#          0		reset scaling parameters
#          1		read from configuration
#
# RESTYPE
      off
#    KDIH      KJ       T   DIFF    RATIO    READ
      0.1     0.1     0.2    0.7      1.0       0
END
@endverbatim
 */
void io::In_Parameter::read_PERSCALE(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read PERSCALE");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["PERSCALE"];
  if (buffer.size()) {

    block_read.insert("PERSCALE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    std::string s1;
    int read;
    _lineStream >> s1
            >> param.pscale.KDIH >> param.pscale.KJ
            >> param.pscale.T >> param.pscale.difference
            >> param.pscale.ratio >> read;
    if (_lineStream.fail()) {
      io::messages.add("bad line in PERSCALE block",
              "In_Parameter", io::message::error);
      return;
    }

    std::transform(s1.begin(), s1.end(), s1.begin(), tolower);

    if (s1 == "jrest" || s1 == "1") {
      param.pscale.jrest = true;
    } else if (s1 == "off" || s1 == "0") {
      param.pscale.jrest = false;
    } else {
      io::messages.add("PERSCALE block: RESTYPE must be jrest or off.",
              "In_Parameter", io::message::error);
    }

    if (param.pscale.KDIH < 0.0)
      io::messages.add("PERSCALE block: KDIH must be >= 0.0.",
            "In_Parameter", io::message::error);
    if (param.pscale.KJ < 0.0)
      io::messages.add("PERSCALE block: KJ must be >= 0.0.",
            "In_Parameter", io::message::error);
    if (param.pscale.T < 0.0)
      io::messages.add("PERSCALE block: T must be >= 0.0.",
            "In_Parameter", io::message::error);
    if (param.pscale.difference < 0.0)
      io::messages.add("PERSCALE block: DIFF must be >= 0.0.",
            "In_Parameter", io::message::error);

    switch (read) {
      case 0:
        param.pscale.read_data = false;
        break;
      case 1:
        param.pscale.read_data = true;
        break;
      default:
        io::messages.add("PERSCALE block: READ must be 0 or 1.",
                "In_Parameter", io::message::error);
    }
  }
} // PERSCALE

/**
 * @section rottrans ROTTRANS block
 * @verbatim
ROTTRANS
# roto-translational constraints
# use either centre of mass removal or roto-translational constraints
# not both!
#
#     RTC: 0,1 controls roto-translational constraints
#          0 don't use roto-translational constraints
#          1 use roto-translational constraints
# RTCLAST: last atom to be roto-translationally constrained
#     RTC RTCLAST
        0       0
END
@endverbatim
 */
void io::In_Parameter::read_ROTTRANS(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read ROTTRANS");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["ROTTRANS"];

  if (buffer.size()) {
    block_read.insert("ROTTRANS");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int rtc;
    _lineStream >> rtc >> param.rottrans.last;

    if (_lineStream.fail())
      io::messages.add("bad line in ROTTRANS block",
            "In_Parameter", io::message::error);

    switch (rtc) {
      case 0:
        param.rottrans.rottrans = false;
        break;
      case 1:
        param.rottrans.rottrans = true;
        break;
      default:
        io::messages.add("ROTTRANS block: RTC must be 0 or 1",
                "In_Parameter", io::message::error);
    }

    if (param.rottrans.last <= 0)
      io::messages.add("ROTTRANS block: RTCLAST must be > 0.",
            "In_Parameter", io::message::error);
  }
}

/**
 * @section spc_loops INNERLOOP block
 * @verbatim
INNERLOOP
# NTILM: 0..3, acceleration method used
#        0: use standard solvent loops [default]
#        1: use fast generic solvent loops
#        2: use solvent loops with hardcoded parameters
#        3: use solvent loops with tabulated forces and energies
#        4: use solvent loops with CUDA library
# NTILS: 0..1, solvent used
#        0: use topology [default]
#        1: use SPC
# NGPUS: number of GPUs to use
# NDEVG: Which GPU device number to use. If not given driver will determine.

# NTILM NTILS NGPUS NDEVG
      4     0     2   0 1
END
@endverbatim
 */
void io::In_Parameter::read_INNERLOOP(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read INNERLOOP");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["INNERLOOP"];

  if (buffer.size()) {
    block_read.insert("INNERLOOP");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int method, solvent;
    _lineStream >> method >> solvent;

    if (_lineStream.fail())
      io::messages.add("bad line in INNERLOOP block",
            "In_Parameter", io::message::error);

    switch (method) {
      case 0:
      {
        // standard solvent loops
        param.innerloop.method = simulation::sla_off;
        break;
      }
      case 1:
      {
        // fast solvent loops
        param.innerloop.method = simulation::sla_generic;
        break;
      }
      case 2:
      {
        param.innerloop.method = simulation::sla_hardcode;
        break;
      }
      case 3:
      {
        // tabulated forces and energies
        param.innerloop.method = simulation::sla_table;
        break;
      }
      case 4:
      {
#ifdef HAVE_LIBCUKERNEL
        // cuda library
        param.innerloop.method = simulation::sla_cuda;
#else
        param.innerloop.method = simulation::sla_off;
        io::messages.add("INNERLOOP block: CUDA solvent loops are not available "
                "in your compilation. Use --with-cukernel for compiling.",
                "In_Parameter", io::message::error);
#endif
        break;
      }
      default:
      {
        param.innerloop.method = simulation::sla_off;
        io::messages.add("INNERLOOP block: bad value for NTILM (0..4)",
                "In_Parameter", io::message::error);
      }
    }

    switch (solvent) {
      case 0:
      {
        // use topology
        param.innerloop.solvent = simulation::sls_topo;
        break;
      }
      case 1:
      {
        // use SPC
        param.innerloop.solvent = simulation::sls_spc;
        break;
      }
      default:
      {
        param.innerloop.solvent = simulation::sls_topo;
        io::messages.add("INNERLOOP block: bad value for NTILS (0,1)",
                "In_Parameter", io::message::error);
      }
    }

    // Get the number of gpus and their device number
    if (param.innerloop.method == simulation::sla_cuda) {
      _lineStream >> param.innerloop.number_gpus;

      if (_lineStream.fail() || param.innerloop.number_gpus == (unsigned int) 0)
        io::messages.add("CUDA enabled, but number of GPUs is zero",
              "In_Parameter", io::message::error);
      else {
        unsigned int temp;
        bool fail = false;
        for (unsigned int i = 0; i < param.innerloop.number_gpus; i++) {
          _lineStream >> temp;
          if (_lineStream.fail()) {
            fail = true;
            break;
          }
          param.innerloop.gpu_device_number.push_back(temp);
        }
        if (fail) {
          param.innerloop.gpu_device_number.clear();
          param.innerloop.gpu_device_number.resize(param.innerloop.number_gpus, -1);
          io::messages.add("CUDA driver will determine devices for nonbonded interaction evaluation.",
                  "In_Parameter", io::message::notice);
        }
      }
    }
  }
}

/**
 * @section replica REPLICA block
 * @verbatim
REPLICA
#     NRET >= 0 number of replica exchange temperatures
#    RET() >= 0.0 temperature for each replica
# LRESCALE 0,1 controls temperature scaling
#          0 don't scale temperatures after exchange trial
#          1 scale temperatures after exchange trial
#   NRELAM >= 0 number of replica exchange lambda values
#  RELAM() >= 0.0 lambda value for each lambda-replica
#   RETS() >= 0.0 timestep of each lambda-replica
# NRETRIAL >= 0 number of overall exchange trials
#  NREQUIL >= 0 number of exchange periods to equilibrate
#               (disallow switches)
#     CONT >= 0 continuation run
#             0 start from one configuration file
#             1 start from multiple configuration files
#   
# NRET
  10
# RET(1..NRET)
  300.0  320.0  340.0 360.0 380.0
  400.0  420.0  440.0 460.0 480.0
# LRESCALE
  1
# NRELAM
  10
# RELAM(1..NRELAM)
  0.0    0.1    0.2   0.3   0.4
  0.5    0.6    0.7   0.8   0.9
# RETS(1..NRELAM)
  0.002  0.001  0.001 0.001 0.002
  0.003  0.002  0.001 0.001 0.002
# NERTRIAL
  100
# NREQUIL
  10
# CONT
  0
END
@endverbatim
 */
void io::In_Parameter::read_REPLICA(simulation::Parameter &param,
        std::ostream & os) {
  DEBUG(8, "read REPLICA");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["REPLICA"];
  if (buffer.size()) {
    block_read.insert("REPLICA");

    bool error = false;

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    _lineStream >> param.replica.num_T;

    if (_lineStream.fail() || param.replica.num_T < 0) {
      io::messages.add("REPLICA block: NRET must be >= 0.",
              "In_Parameter", io::message::error);
      error = true;
    }

    if (!error)
      param.replica.temperature.resize(param.replica.num_T, 0.0);

    for (int i = 0; i < param.replica.num_T; ++i) {
      _lineStream >> param.replica.temperature[i];
      if (_lineStream.fail() || param.replica.temperature[i] < 0.0) {
        std::ostringstream msg;
        msg << "REPLICA block: RET(" << i + 1 << ") must be >= 0.0";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        error = true;
      }
    }

    int scale;
    _lineStream >> scale;
    switch (scale) {
      case 0:
        param.replica.scale = false;
        break;
      case 1:
        param.replica.scale = true;
        break;
      default:
        io::messages.add("REPLICA block: LRSCALE must be 0 or 1",
                "In_Parameter", io::message::error);
        error = true;
    }

    if (_lineStream.fail() || error) {
      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
      param.replica.dt.clear();
      error = false;
    }

    _lineStream >> param.replica.num_l;

    if (_lineStream.fail() || param.replica.num_l < 0) {
      io::messages.add("REPLICA block: NRELAM must be >= 0.",
              "In_Parameter", io::message::error);
      error = true;
    }

    if (!error) {
      param.replica.lambda.resize(param.replica.num_l, 0.0);
      param.replica.dt.resize(param.replica.num_l, 0.0);
    }

    for (int i = 0; i < param.replica.num_l; ++i) {
      _lineStream >> param.replica.lambda[i];
      if (_lineStream.fail() || param.replica.lambda[i] < 0.0) {
        std::ostringstream msg;
        msg << "REPLICA block: RELAM(" << i + 1 << ") must be >= 0.0";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        error = true;
      }
    }
    for (int i = 0; i < param.replica.num_l; ++i) {
      _lineStream >> param.replica.dt[i];
      if (_lineStream.fail() || param.replica.dt[i] < 0.0) {
        std::ostringstream msg;
        msg << "REPLICA block: RETS(" << i + 1 << ") must be >= 0.0";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        error = true;
      }
    }

    if (_lineStream.fail() || error) {
      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
      param.replica.dt.clear();
      error = false;
    }

    _lineStream >> param.replica.trials;
    if (_lineStream.fail() || param.replica.trials < 0) {
      io::messages.add("REPLICA block: NRETRIAL must be >= 0.",
              "In_Parameter", io::message::error);
      error = true;
    }
    _lineStream >> param.replica.equilibrate;
    if (_lineStream.fail() || param.replica.equilibrate < 0) {
      io::messages.add("REPLICA block: NREQUIL must be >= 0.",
              "In_Parameter", io::message::error);
      error = true;
    }
    // do continuation run
    _lineStream >> param.replica.cont;
    if (_lineStream.fail() || param.replica.cont < 0 || param.replica.cont > 1 ) {
      io::messages.add("REPLICA block: CONT must be 0 or 1",
              "In_Parameter", io::message::error);
      error = true;
    }

    if (_lineStream.fail() || error) {
      io::messages.add("bad line in REPLICA block",
              "In_Parameter", io::message::error);

      param.replica.num_T = 0;
      param.replica.num_l = 0;

      param.replica.temperature.clear();
      param.replica.lambda.clear();
      param.replica.dt.clear();
    }
  }
}

/**
 * @section multicell MULTICELL block
 * @verbatim
MULTICELL
#  NTM: 0,1 switch for multiple-unit-cell simulation.
#       0 : single-unit-cell simulation [default]
#       1 : multiple-unit-cell simulation
#         NTM
            0
#  number of subdivisions along axis
#   NCELLA    NCELLB    NCELLC  
         1         1         1
#  periodicity checks (relative tolerance)
#  not available in md++ -> 0.0
#    TOLPX     TOLPV     TOLPF    TOLPFW
       0.0       0.0       0.0       0.0
END
@endverbatim
 */
void io::In_Parameter::read_MULTICELL(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read MULTICELL");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["MULTICELL"];

  if (buffer.size()) {

    block_read.insert("MULTICELL");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    int ntm;
    double tolpx, tolpv, tolpf, tolpfw;
    _lineStream >> ntm
            >> param.multicell.x >> param.multicell.y >> param.multicell.z
            >> tolpx >> tolpv >> tolpf >> tolpfw;

    if (_lineStream.fail()) {
      io::messages.add("bad line in MULTICELL block",
              "In_Parameter", io::message::error);

      param.multicell.multicell = false;
    }

    switch (ntm) {
      case 1:
        param.multicell.multicell = true;
        break;
      case 0:
        param.multicell.multicell = false;
        break;
      default:
        param.multicell.multicell = false;
        io::messages.add("MULTICELL block: NTM must be 0 or 1.",
                "In_Parameter", io::message::error);
    }

    if (param.multicell.multicell) {
      /*
      // disable because broken
      param.multicell.multicell = false;
      io::messages.add("MULTICELL simulations are broken in MD++",
                         "In_Parameter", io::message::error);      */


      // do these checks only if mutlicell is really used.
      if (param.multicell.x < 1 || param.multicell.y < 1 ||
              param.multicell.z < 1) {
        io::messages.add("MULTICELL block: NCELLA, NCELLB and NCELLC "
                "must be >= 1.", "In_Parameter", io::message::error);
      }

      if (param.multicell.x == 1 && param.multicell.y == 1 &&
              param.multicell.z == 1) {
        io::messages.add("MULTICELL block: NCELLA, NCELLB and NCELLC are all 1.\n"
                "disabling MULTICELL simulation.", "In_Parameter",
                io::message::warning);
        param.multicell.multicell = false;
      }

      if (tolpx || tolpv || tolpf || tolpfw) {
        io::messages.add("MULTICELL block: Periodicity checks are not needed "
                "in the MD++ implementation of MULTICELL.",
                "In_Parameter", io::message::warning);
      }
    } else {
      // make sure all values are set to one in a normal, non-multicell, simulation!
      param.multicell.x = param.multicell.y = param.multicell.z = 1;
    }
  }
}

/**
 * @section readtraj READTRAJ block
 * @verbatim
READTRAJ
# NTRD  0,1 controls trajectory-reevaluation mode
#       0: do not use trajectory-reevaluation mode
#       1: use trajectory-reevaluation mode
# NTRN  number of files (ignored)
# NTRB  read box (must be 1)
# NTSHK 0,1 controls SHAKE on old coordinates
#       0 perform SHAKE with respect to previous coordinates
#       1 perform SHAKE with respect to current coordinates
#
#   NTRD    NTRN    NTRB   NTSHK
       0       0       1       0    
END
@endverbatim
 */
void io::In_Parameter::read_READTRAJ(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read READTRAJ");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["READTRAJ"];

  if (buffer.size()) {

    block_read.insert("READTRAJ");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int ntrd, ntrn, ntrb, ntshk;
    _lineStream >> ntrd >> ntrn >> ntrb >> ntshk;

    if (_lineStream.fail()) {
      io::messages.add("bad line in READTRAJ block",
              "In_Parameter", io::message::error);

      param.analyze.analyze = false;
      return;
    }

    switch (ntrd) {
      case 1:
        param.analyze.analyze = true;
        break;
      case 0:
        param.analyze.analyze = false;
        break;
      default:
        param.analyze.analyze = false;
        io::messages.add("READTRAJ block: NTRD must be 0 or 1",
                "In_Parameter", io::message::error);
    }

    if (ntrn)
      io::messages.add("READTRAJ block: NTRN was ignored", "In_Parameter",
            io::message::warning);

    if (ntrb != 1)
      io::messages.add("READTRAJ block: NTRB must be 1.", "In_Parameter",
            io::message::error);

    switch (ntshk) {
      case 1:
        param.analyze.copy_pos = true;
        break;
      case 0:
        param.analyze.copy_pos = false;
        break;
      default:
        param.analyze.copy_pos = false;
        io::messages.add("READTRAJ block: NTSHK must be 0 or 1",
                "In_Parameter", io::message::error);
    }

  }
}

/**
 * @section integrate INTEGRATE block
 * @verbatim
INTEGRATE
#  NINT 0..1 selects integration method
#	0: no integration performed
#	1: leap-frog integration scheme performed
#
#    NINT
        1
END
@endverbatim
 */
void io::In_Parameter::read_INTEGRATE(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read INTEGRATE");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["INTEGRATE"];

  if (buffer.size()) {

    block_read.insert("INTEGRATE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int nint;
    _lineStream >> nint;

    if (_lineStream.fail()) {
      io::messages.add("bad line in INTEGRATE block",
              "In_Parameter", io::message::error);

    }

    switch (nint) {
      case 0:
        param.integrate.method = simulation::integrate_off;
        break;
      case 1:
        param.integrate.method = simulation::integrate_leap_frog;
        break;
      default:
        io::messages.add("INTEGRATE block: NINT must be 0 or 1",
                "In_Parameter", io::message::error);
    }
  }
}

/**
 * @section stochdyn STOCHDYN block
 * @verbatim
STOCHDYN
# NTSD    0,1 controls stochastic dynamics mode
#         0: do not do stochastic dynamics (default)
#         1: do stochastic dynamics
# NTFR    0..3 defines atomic friction coefficients gamma
#         0: set gamma to 0.0 (default)
#         1: set gamma to CFRIC
#         2: set gamma to CFRIC*GAM0
#         3: calculate gamma using subroutine FRIC (based on CFRIC)
# NSFR    > 0 recalculate gamma every NSFR steps
# NBREF   > 0 threshold number of neighbour atoms for a buried atom
# RCUTF   >= 0.0 interatomic distance considered when calculating gamma
# CFRIC   >= 0.0 global weighting for gamma
# TEMPSD  >= 0.0 temperature of stochastic bath
#
#     NTSD     NTFR     NSFR   NBREF  RCUTF    CFRIC    TEMPSD
         0        1        0       6    0.3     91.0     300.0
END
@endverbatim
 */
void io::In_Parameter::read_STOCHDYN(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read STOCHDYN");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["STOCHDYN"];

  if (buffer.size()) {

    block_read.insert("STOCHDYN");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    _lineStream >> param.stochastic.sd >> param.stochastic.ntfr
            >> param.stochastic.nsfr >> param.stochastic.nbref
            >> param.stochastic.rcutf >> param.stochastic.cfric
            >> param.stochastic.temp;

    if (_lineStream.fail()) {
      io::messages.add("bad line in STOCHDYN block",
              "In_Parameter", io::message::error);
      return;
    }

    if (param.stochastic.sd < 0 || param.stochastic.sd > 1)
      io::messages.add("STOCHDYN block: NTSD must be 0 or 1",
            "In_Parameter", io::message::error);

    if (param.stochastic.ntfr < 0 || param.stochastic.ntfr > 3)
      io::messages.add("STOCHDYN block: NTFR must be 0 to 3",
            "In_Parameter", io::message::error);

    if (param.stochastic.nsfr <= 0)
      io::messages.add("STOCHDYN block: NSFR must be > 0",
            "In_Parameter", io::message::error);

    if (param.stochastic.nbref <= 0)
      io::messages.add("STOCHDYN block: NBREF must be > 0",
            "In_Parameter", io::message::error);

    if (param.stochastic.rcutf < 0)
      io::messages.add("STOCHDYN block: RCUTF must be >= 0",
            "In_Parameter", io::message::error);

    if (param.stochastic.cfric < 0)
      io::messages.add("STOCHDYN block: CFRIC must be >= 0",
            "In_Parameter", io::message::error);

    if (param.stochastic.temp < 0)
      io::messages.add("STOCHDYN block: TEMPSD must be >= 0",
            "In_Parameter", io::message::error);
  }
}

/**
 * @section ewarn EWARN block
 * @verbatim
EWARN
# MAXENER issue a warning if total energy is larger than this value
#
# MAXENER
   100000
END
@endverbatim
 */
void io::In_Parameter::read_EWARN(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read EWARN");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["EWARN"];

  if (buffer.size()) {

    block_read.insert("EWARN");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    _lineStream >> param.ewarn.limit;

    if (_lineStream.fail()) {
      io::messages.add("bad line in EWARN block",
              "In_Parameter", io::message::error);
      return;
    }
  } else {
#ifndef HAVE_ISNAN
    io::messages.add("std::isnan() is not available for your compilation. "
            "Consider using the EWARN block.",
            "In_Parameter", io::message::warning);
#endif
  }
}

/**
 * @section multistep MULTISTEP block
 * @verbatim
MULTISTEP
#   STEPS calculate non-bonded every STEPSth step.
#   BOOST 0,1 
#         0: stored forces of STEPSth step are added every step
#         1: stored forces of STEPSth setp are multiplied by STEPS 
#            and added every STEPSth step.
#
#   STEPS   BOOST
        0       0
END
@endverbatim
 */
void io::In_Parameter::read_MULTISTEP(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read MULTISTEP");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["MULTISTEP"];

  if (buffer.size()) {
    block_read.insert("MULTISTEP");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    int boost;
    _lineStream >> param.multistep.steps >> boost;

    if (_lineStream.fail()) {
      io::messages.add("bad line in MULTISTEP block",
              "In_Parameter", io::message::error);
      return;
    }

    if (param.multistep.steps < 0) {
      io::messages.add("MULTISTEP block: STEPS must be >= 0.",
              "In_Parameter", io::message::error);
    }
    switch (boost) {
      case 0:
        param.multistep.boost = false;
        break;
      case 1:
        param.multistep.boost = true;
        break;
      default:
        io::messages.add("MULTISTEP block: BOOST must be 0 or 1",
                "In_Parameter", io::message::error);
    }
  }
}

void io::In_Parameter::read_MONTECARLO(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read MONTECARLO");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["MONTECARLO"];

  if (buffer.size()) {
    block_read.insert("MONTECARLO");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    int mc;
    _lineStream >> mc >> param.montecarlo.steps >> param.montecarlo.dlambda;

    if (_lineStream.fail()) {
      io::messages.add("bad line in MONTECARLO block",
              "In_Parameter", io::message::error);
      return;
    }

    switch (mc) {
      case 0: param.montecarlo.mc = mc;
        break;
      case 1: param.montecarlo.mc = mc;
        break;
      default:
      {
        io::messages.add("MONTECARLO block: MC must be 0 or 1",
                "In_Parameter", io::message::error);
        break;
      }
    }

    // parameters are all positive
    if (param.montecarlo.steps < 0 || param.montecarlo.dlambda < 0) {
      io::messages.add("MONTECARLO block: Negative parameter",
              "In_Parameter", io::message::error);
      return;
    }
    // perturbation
    if (param.montecarlo.mc && !param.perturbation.perturbation) {
      io::messages.add("Chemical MONTECARLO only possible if perturbation is on."
              " Set NTG to 1.",
              "In_Parameter", io::message::error);
      return;
    }
  }
}

/**
 * @section polarise POLARISE block
 * @verbatim
POLARISE
# COS      0,1 use polarisation
#          0: don't use polarisation
#          1: use charge-on-spring model for dipolar polarisation
#          2: use charge-on-spring model for dipolar polarisation with off atom site
# EFIELD   0,1 controls evaluation site for electric field
#          0: evaluate at atom position
#          1: evaluate at cos position
# MINFIELD >0.0 convergence criterium
# DAMP     0,1 controls polarisability damping
#          0: don't damp polarisability
#          1: damp polarisability (with paramaters from topology)
# WRITE    > 0 write COS positions to special trajectory
#          0: don't write
#         >0: write COS positions every WRITEth step
#
#     COS  EFIELD MINFIELD    DAMP  WRITE
        0       0      2.5       0      0
END
@endverbatim
 */
void io::In_Parameter::read_POLARISE(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read POLARISE");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["POLARISE"];

  if (buffer.size()) {

    //param.setDevelop("Polarisation is under development.");

    block_read.insert("POLARISE");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int cos, damp, efield;
    _lineStream >> cos >> efield >> param.polarise.minfield >> damp
            >> param.polarise.write;

    if (_lineStream.fail()) {
      io::messages.add("bad line in POLARISE block",
              "In_Parameter", io::message::error);
      return;
    }

    switch (cos) {
      case 0:
      {
        param.polarise.cos = 0;
        param.polarise.write = 0;
        break;
      }
      case 1:
      {
        param.polarise.cos = 1;
        param.force.interaction_function = simulation::pol_lj_crf_func;
        break;
      }
      case 2:
      {
        param.polarise.cos = 2;
        param.force.interaction_function = simulation::pol_off_lj_crf_func;
        break;
      }

      default:
        io::messages.add("POLARISE block: COS must be 0, 1 or 2.",
                "In_Parameter", io::message::error);
    }

    switch (efield) {
      case 0: param.polarise.efield_site = simulation::ef_atom;
        break;
      case 1: param.polarise.efield_site = simulation::ef_cos;
        break;
      default:
        io::messages.add("POLARISE block: EFIELD must be 0 or 1.",
                "In_Parameter", io::message::error);
    }

    if (param.polarise.minfield <= 0.0) {
      io::messages.add("POLARISE block: MINFIELD must be > 0.0",
              "In_Parameter", io::message::error);
    }

    switch (damp) {
      case 0: param.polarise.damp = false;
        break;
      case 1: param.polarise.damp = true;
        break;
      default:
        io::messages.add("POLARISE block: DAMP must be 0 or 1.",
                "In_Parameter", io::message::error);
    }

    if (!param.polarise.cos)
      param.polarise.write = 0;

    if (param.polarise.write < 0) {
      io::messages.add("POLARISE block: WRITE must be >= 0",
              "In_Parameter", io::message::error);
    }

    if (param.polarise.damp && !param.polarise.cos) {
      io::messages.add("POLARISE block: DAMP is ignored if no polarisation is used",
              "In_Parameter", io::message::warning);
    }

    if (param.polarise.cos && param.multicell.multicell) {
      io::messages.add("multiple unit cell simulation using polarisation may "
              "converge differently. Use a small MINFIELD parameter.",
              "In_Parameter", io::message::warning);
    }
  }
}

/**
 * @section randomnumbers RANDOMNUMBERS block
 * @verbatim
RANDOMNUMBERS
# NTRNG 0,1 random number generator
#         0 use G96 algorithm (default)
#         1 use GSL library
# NTGSL -1.. GSL random number generation algorithm
#         -1: use default algorithm (mt19937)
#       >=0 : run contrib/rng_gsl for a list of possible arguments
#
#   NTRNG   NTGSL
        1      -1
END
@endverbatim
 */
void io::In_Parameter::read_RANDOMNUMBERS(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read RANDOMNUMBERS");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["RANDOMNUMBERS"];

  if (buffer.size()) {
    block_read.insert("RANDOMNUMBERS");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int rng;
    _lineStream >> rng >> param.rng.gsl_rng;

    if (_lineStream.fail()) {
      io::messages.add("bad line in RANDOMNUMBERS block",
              "In_Parameter", io::message::error);
      return;
    }

    switch (rng) {
      case 0:
        param.rng.rng = simulation::random_g96;
        break;
      case 1:
        param.rng.rng = simulation::random_gsl;
        break;
      default:
        io::messages.add("RANDOMNUMBERS block: NTRNG must be 0 (G96) "
                "or 1 (GSL)", "In_Parameter", io::message::error);
    }

    math::RandomGenerator::check(param);
  }
}

/**
 * @section EDS EDS block
 * @verbatim
EDS
# EDS        0,1  
#              0: no enveloping distribution sampling (EDS) [default]
#              1: enveloping distribution sampling
# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter
#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter
# FORM       1-3
#              1: Single s Hamiltonian
#              2: Hamiltonian with NUMSTATES*(NUMSTATES-1)/2 (pairwise) s parameters
#              3: Hamiltonian with (NUMSTATES-1) s parameters
# NUMSTATES >1  : number of states
# if NUMSTATES != 3:
# S         >0.0: smoothness parameter(s)
# if NUMSTATES == 3:
# i   j   S     : state pair i j and associated s parameter
# EIR           : energy offsets for states
#
# EDS          
  1
# SOFT CORE
# ALPHLJ  ALPHC
  0.0     0.0
# FUNCTIONAL FORM
  3
# NUMSTATES
  3
# S
  0.2  0.01 0.1
# EIR 
  0   20   40
END
# example for FORM = 3
EDS
  1    
# SOFT CORE
# ALPHLJ  ALPHC
  0.0     0.0
# FUNCTIONAL FORM
  3
# NUMSTATES
  3
# i  j  S
  1  2  0.1
  2  3  0.5
# EIR 
  0   20   40
END
@endverbatim
 */
void io::In_Parameter::read_EDS(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read EDS");

  std::vector<std::string> buffer;
  std::string s;
  int form;
  double soft_lj, soft_crf;

  buffer = m_block["EDS"];

  if (buffer.size()) {
    block_read.insert("EDS");

    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int eds;
    _lineStream >> eds >> soft_lj >> soft_crf >> form >> param.eds.numstates;
    if (_lineStream.fail()) {
      io::messages.add("bad line in EDS block",
              "In_Parameter", io::message::error);
      return;
    }


    switch (eds) {
      case 0:
        param.eds.eds = false;
        param.eds.numstates = 0;
        break;
      case 1:
        param.eds.eds = true;
        break;
      default:
        io::messages.add("Error in EDS block: EDS must be 0 (no EDS) "
                "or 1 (EDS)", "In_Parameter", io::message::error);
    }
    if (param.eds.eds) {
      if (param.eds.soft_vdw < 0.0) {
        io::messages.add("EDS block: ALPHLJ must be >= 0.",
                "In_Parameter", io::message::error);
      }

      if (param.eds.soft_crf < 0.0) {
        io::messages.add("EDS block: ALPHC must be >= 0.",
                "In_Parameter", io::message::error);
      }
      param.eds.soft_vdw = soft_lj;
      param.eds.soft_crf = soft_crf;
      if (soft_lj > 0.0 || soft_crf > 0.0)
        param.setDevelop("Soft-core EDS is under development.");
      
      switch (form) {
        case 1: {
          param.eds.form = simulation::single_s;
          // read in 1 s value
          param.eds.s.resize(1, 1.0);
          _lineStream >> param.eds.s[0];
          //std::cerr << " s[0] = " << param.eds.s[0] <<  std::endl;
          if (param.eds.s[0] <= 0) {
            io::messages.add("Error in EDS block: S must be >0",
                    "In_Parameter", io::message::error);
          }
          break;
        }
        case 2: {
          param.eds.form = simulation::multi_s;
          const unsigned int n = param.eds.numstates;
          param.eds.s.resize((n * (n - 1)) / 2, 1.0);
          for (unsigned int pair = 0; pair < param.eds.s.size(); pair++) {
            _lineStream >> param.eds.s[pair];
            if (param.eds.s[pair] <= 0) {
              io::messages.add("Error in EDS block: S must be >0",
                      "In_Parameter", io::message::error);
            }
          }
          break;
        }
        case 3: {
          param.eds.form = simulation::pair_s;
          const unsigned int n = param.eds.numstates;
          param.eds.s.resize(n - 1, 1.0);
          param.eds.pairs.resize(n - 1);
          for (unsigned int pair = 0; pair < param.eds.s.size(); pair++) {
            _lineStream >> param.eds.pairs[pair].i
                    >> param.eds.pairs[pair].j
                    >> param.eds.s[pair];
            if (param.eds.s[pair] <= 0) {
              io::messages.add("Error in EDS block: S must be >0",
                      "In_Parameter", io::message::error);
            }
          }
          break;
        }
        default:
          io::messages.add("Error in EDS block: functional form must be 1 (single s) "
                  "or 2 (multiple s)", "In_Parameter", io::message::error);
      }
      if (_lineStream.fail()) {
        io::messages.add("Error when reading S from EDS block",
                "In_Parameter", io::message::error);
        return;
      }
   
      param.eds.eir.resize(param.eds.numstates, 0.0);
      for (unsigned int i = 0; i < param.eds.numstates; i++) {
        _lineStream >> param.eds.eir[i];
        //std::cerr << "eir = " << param.eds.eir[i] << std::endl;
      }
      if (_lineStream.fail()) {
        io::messages.add("Error when reading EIR from EDS block",
                "In_Parameter", io::message::error);
        return;
      }

      if (param.eds.numstates < 2) {
        io::messages.add("Error in EDS block: NUMSTATES must be >=1.",
                "In_Parameter", io::message::error);
      }
      // make sure we simulate at a given temperature (unambiguous kT)
      if (!param.multibath.couple) {
        io::messages.add("Error in EDS block: EDS requires temperature coupling.",
                "In_Parameter", io::message::error);
      }
      // check whether all baths have the same temperature (unambiguous kT)
      for (unsigned int i = 1; i < param.multibath.multibath.size(); i++) {
        if (param.multibath.multibath.bath(i).temperature !=
                param.multibath.multibath.bath(0).temperature) {
          io::messages.add("Error in EDS block: all baths must have the same temperature.",
                  "In_Parameter", io::message::error);
        }
      }
    }
    //std::cerr << "eds (at end) = " <<  eds << ", form = " << form << ", numstates=" << param.eds.numstates;
  }
}

/**
 * @section LAMBDAS LAMBDAS block
 * @verbatim
LAMBDAS
# NTIL    off(0), on(1)
#         0: no special treatment of interactions with individual lambda-values
#         1: interactions are treated with special individual lambda-values
# NTLI(1..)  interaction type to treat with individual lambda: 
#            bond(1), angle(2), dihedral(3), improper(4), vdw(5), vdw_soft(6),
#            crf(7), crf_soft(8), distanceres(9), distancefield(10), 
#            dihedralres(11), mass(12)
# NILG1, NILG2 energy groups of interactions that are treated with individual 
#              lambda values
# ALI, BLI, CLI, DLI, ELI polynomial coefficients linking the individual lambda-
#                         values to the overall lambda-value
# NTIL
   1
# NTLI
  7
# NILG1  NILG2
  1      3
END
@endverbatim
 */
void io::In_Parameter::read_LAMBDAS(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read LAMBDAS");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["LAMBDAS"];

  if (buffer.size()) {

    block_read.insert("LAMBDAS");

    int num = buffer.size() - 3;
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    std::string nm;
    simulation::interaction_lambda_enum j;
    int n1, n2;
    double a, b, c, d, e;
    _lineStream >> nm;
    DEBUG(10, "read NTIL " << nm);

    if (nm == "on" || nm == "1")
      param.lambdas.individual_lambdas = true;
    else if (nm == "off" || nm == "0") {
      param.lambdas.individual_lambdas = false;
      return;
    } else {
      io::messages.add("illegal value for NTIL in LAMBDAS block (on,off,1,0)",
              "In_Parameter", io::message::error);
      return;
    }

    if (param.perturbation.perturbation == false)
      io::messages.add("LAMBDAS block without perturbation is ignored",
            "In_Parameter", io::message::warning);

    // the full matrices for the energy groups were created and
    // filled with a, b, c, e = 0 and d=1 when the FORCE block was read in
    // that way, they are also defined if you do not give the LAMBDAS block
    // we cannot initialize them earlier, because they depend on the 
    // energy groups

    int maxnilg = param.force.energy_group.size();
    for (int i = 0; i < num; ++i) {
      _lineStream >> nm >> n1 >> n2 >> a >> b >> c >> d >> e;
      DEBUG(10, "read : " << nm << n1 << n2 << a << b << c << d << e);

      if (_lineStream.fail()) {
        io::messages.add("bad line in LAMBDAS block" + s,
                "In_Parameter", io::message::error);
        return;
      }
      if (n2 < n1) {
        io::messages.add("only give NILG2 >= NILG1 in LAMBDAS BLOCK",
                "In_Parameter", io::message::error);
        return;
      }
      if (n1 > maxnilg) {
        io::messages.add("NILG1 larger than number of energy groups in FORCE block",
                "In_Parameter", io::message::error);
        return;
      }
      if (n2 > maxnilg) {
        io::messages.add("NILG2 larger than number of energy groups in FORCE block",
                "In_Parameter", io::message::error);
        return;
      }
      n1--;
      n2--;

      if (nm == "bond" || nm == "1")
        j = simulation::bond_lambda;
      else if (nm == "angle" || nm == "2")
        j = simulation::angle_lambda;
      else if (nm == "dihedral" || nm == "3")
        j = simulation::dihedral_lambda;
      else if (nm == "improper" || nm == "4")
        j = simulation::improper_lambda;
      else if (nm == "vdw" || nm == "5")
        j = simulation::lj_lambda;
      else if (nm == "vdw_soft" || nm == "6")
        j = simulation::lj_softness_lambda;
      else if (nm == "crf" || nm == "7")
        j = simulation::crf_lambda;
      else if (nm == "crf_soft" || nm == "8")
        j = simulation::crf_softness_lambda;
      else if (nm == "distanceres" || nm == "9")
        j = simulation::disres_lambda;
      else if (nm == "distancefield" || nm == "10")
	j = simulation::disfield_lambda;
      else if (nm == "dihedralres" || nm == "11")
        j = simulation::dihres_lambda;
      else if (nm == "mass" || nm == "12")
        j = simulation::mass_lambda;
      else {
        io::messages.add("unknown lambda type in LAMBDAS block: " + nm,
                "In_Parameter", io::message::error);
        return;
      }

      // and now replace the matrix with the numbers we just read in
      if (j != simulation::lj_lambda &&
              j != simulation::lj_softness_lambda &&
              j != simulation::crf_lambda &&
              j != simulation::crf_softness_lambda &&
              n1 != n2)
        io::messages.add("NILG1 != NILG2 in LAMBDAS block only allowed for nonbonded interactions",
              "In_Parameter", io::message::warning);

      param.lambdas.a[j][n1][n2] = a;
      param.lambdas.a[j][n2][n1] = a;
      param.lambdas.b[j][n1][n2] = b;
      param.lambdas.b[j][n2][n1] = b;
      param.lambdas.c[j][n1][n2] = c;
      param.lambdas.c[j][n2][n1] = c;
      param.lambdas.d[j][n1][n2] = d;
      param.lambdas.d[j][n2][n1] = d;
      param.lambdas.e[j][n1][n2] = e;
      param.lambdas.e[j][n2][n1] = e;
    }
  }
}

/**
 * @section nonbonded NONBONDED block
 * @verbatim
NONBONDED
# NLRELE    1-3 method to handle electrostatic interactions
#     1 : reaction-field
#     2 : Ewald method
#     3 : P3M method
# APPAK     >= 0.0 reaction-field inverse Debye screening length
# RCRF      >= 0.0 reaction-field radius 
#   0.0 : set to infinity
# EPSRF     = 0.0 || > 1.0 reaction-field permittivity
#   0.0 : set to infinity
# NSLFEXCL  0,1 contribution of excluded atoms to reaction field
      0 : contribution turned off
      1 : contribution considered (default)
# NSHAPE    -1..10 lattice sum charge-shaping function
#    -1 : gaussian
# 0..10 : polynomial 
# ASHAPE    > 0.0 width of the lattice sum charge-shaping function
# NA2CALC   0,2 controls evaluation of lattice sum A2 term
#     0 : set to zero
#     2 : numerical
# TOLA2     > 0.0 tolerance for numerical A2 evaluation 
# EPSLS      = 0.0 || > 1.0 lattice sum permittivity (0.0 = tinfoil)
# NKX, NKY, NKZ > 0 maximum absolute Ewald k-vector components
# KCUT       > 0.0 Ewald k-space cutoff
# NGX, NGY, NGZ > 0 P3M number of grid points
# NASORD    1..5 order of mesh charge assignment function
# NFDORD    0..5 order of the mesh finite difference operator
#     0 : ik - differentiation
#  1..5 : finite differentiation
# NALIAS    > 0 number of mesh alias vectors considered
# NSPORD        order of SPME B-spline functions (ignored! promd specific)
# NQEVAL    >= 0 controls accuracy reevaluation
#     0 : do not reevaluate
#   > 0 : evaluate every NQEVAL steps
# FACCUR    >= 0.0 rms force error threshold to recompute influence function
# NRDGRD    0,1 read influence function
#     0 : calculate influence function at simulation start up
#     1 : read influence function from file (not yet implemented)
# NWRGRD    0,1 write influence function
#     0 : do not write
#     1 : write at the end of the simulation (not yet implemented)
# NLRLJ     0,1 controls long-range Lennard-Jones corrections
#     0 : no corrections
#     1 : do corrections (not yet implemented)
# SLVDNS    > 0.0 average solvent density for long-range LJ correction (ignored)      
# 
#   NLRELE
         1
#    APPAK      RCRF     EPSRF   NSLFEXCL
       0.0       1.4      61.0          1
#   NSHAPE    ASHAPE    NA2CLC      TOLA2    EPSLS
        -1       1.4         2     0.1E-9      0.0
#      NKX       NKY       NKZ       KCUT
        10        10        10      100.0
#      NGX       NGY       NGZ    NASORD    NFDORD    NALIAS    NSPORD
        32        32        32         3         2         3         4
#   NQEVAL    FACCUR    NRDGRD    NWRGDR
    100000       1.6         0         0
#    NLRLJ    SLVDNS
         0      33.3
END
@endverbatim
 */
void io::In_Parameter::read_NONBONDED(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read NONBONDED");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["NONBONDED"];

  if (buffer.size()) {
    block_read.insert("NONBONDED");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    int method;
    int ls_calculate_a2;
    _lineStream >> method >>
            param.nonbonded.rf_kappa >> param.nonbonded.rf_cutoff >> param.nonbonded.rf_epsilon >>
            param.nonbonded.selfterm_excluded_atoms >>
            param.nonbonded.ls_charge_shape >> param.nonbonded.ls_charge_shape_width >>
            ls_calculate_a2 >> param.nonbonded.ls_a2_tolerance >>
            param.nonbonded.ls_epsilon >>
            param.nonbonded.ewald_max_k_x >> param.nonbonded.ewald_max_k_y >>
            param.nonbonded.ewald_max_k_z >> param.nonbonded.ewald_kspace_cutoff >>
            param.nonbonded.p3m_grid_points_x >> param.nonbonded.p3m_grid_points_y >>
            param.nonbonded.p3m_grid_points_z >> param.nonbonded.p3m_charge_assignment >>
            param.nonbonded.p3m_finite_differences_operator >> param.nonbonded.p3m_mesh_alias >>
            param.nonbonded.spme_bspline >>
            param.nonbonded.accuracy_evaluation >> param.nonbonded.influence_function_rms_force_error >>
            param.nonbonded.influence_function_read >> param.nonbonded.influence_function_write >>
            param.nonbonded.lj_correction >> param.nonbonded.lj_solvent_density;

    if (_lineStream.fail()) {
      io::messages.add("bad line in NONBONDED block",
              "In_Parameter", io::message::error);
    }
    bool do_ls = false;
    switch (method) {
      case 0:
        param.force.nonbonded_crf = 0;
        break;
      case 1:
        param.nonbonded.method = simulation::el_reaction_field;
        param.nonbonded.lserf = false;
        break;
      case -1:
        param.nonbonded.method = simulation::el_reaction_field;
        param.nonbonded.lserf = false;
        break;
      case 2:
        param.nonbonded.method = simulation::el_ewald;
        param.setDevelop("Ewald sum is under development");
        param.nonbonded.lserf = false;
        param.nonbonded.rf_excluded = false;
        do_ls = true;
        break;
      case 3:
        param.nonbonded.method = simulation::el_p3m;
        param.nonbonded.lserf = false;
        param.nonbonded.rf_excluded = false;
        do_ls = true;
        break;
      default:
        io::messages.add("NONBONDED block: electrostatic method not implemented",
                "In_Parameter", io::message::error);
    }
    
    if (param.nonbonded.selfterm_excluded_atoms < 0 || param.nonbonded.selfterm_excluded_atoms > 1)
      io::messages.add("NONBONDED block: Illegal value for NSLFEXCL (1,0)",
            "In_Parameter", io::message::error);
    
    // switch to turn of the contribution of the excluded atoms to the reaction field
    if (param.nonbonded.method == simulation::el_reaction_field 
            && param.nonbonded.selfterm_excluded_atoms == 0) {
      param.nonbonded.rf_excluded = false;
      io::messages.add("NONBONDED block: contribution of excluded atoms to the reaction field turned off!",
                "In_Parameter", io::message::warning);
    }

    if (param.nonbonded.method != simulation::el_reaction_field)
      param.force.interaction_function = simulation::lj_ls_func;

    if (param.nonbonded.rf_kappa < 0)
      io::messages.add("NONBONDED block: Illegal value for APPAK (>=0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.rf_cutoff < 0)
      io::messages.add("NONBONDED block: Illegal value for RCRF (>=0.0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.rf_epsilon != 0.0 && param.nonbonded.rf_epsilon < 1.0)
      io::messages.add("NONBONDED block: Illegal value for EPSRF (0.0 / >=1.0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.ls_charge_shape < -1 ||
            param.nonbonded.ls_charge_shape > 10)
      io::messages.add("NONBONDED block: Illegal value for NSHAPE (-1..10)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.ls_charge_shape_width <= 0.0)
      io::messages.add("NONBONDED block: Illegal value for ASHAPE (>0.0)",
            "In_Parameter", io::message::error);

    if (do_ls && param.nonbonded.ls_charge_shape_width > param.pairlist.cutoff_short && 
        fabs(param.nonbonded.ls_charge_shape_width - param.pairlist.cutoff_short) > math::epsilon)
      io::messages.add("NONBONDED block: charge width greater than cutoff! (ASHAPE > RCUTP)",
            "In_Parameter", io::message::warning);


    switch (ls_calculate_a2) {
      case 0:
        param.nonbonded.ls_calculate_a2 = simulation::ls_a2_zero;
        break;
      case 1:
        param.nonbonded.ls_calculate_a2 = simulation::ls_a2t_exact;
        break;
      case 2:
        param.nonbonded.ls_calculate_a2 = simulation::ls_a2_numerical;
        break;
      case 3:
        param.nonbonded.ls_calculate_a2 = simulation::ls_a2t_exact_a2_numerical;
        break;
      case 4:
        param.nonbonded.ls_calculate_a2 = simulation::ls_a2t_exact_a2_numerical;
        if (param.nonbonded.method != simulation::el_p3m) {
          io::messages.add("NONBONDED block: averaged A2~ calculation needs P3M.",
                  "In_Parameter", io::message::error);
        }
        break;
      default:
        io::messages.add("NONBONDED block: A2 calculation method not implemented",
                "In_Parameter", io::message::error);
    }

    if (param.nonbonded.ls_epsilon != 0.0 && param.nonbonded.ls_epsilon < 1.0)
      io::messages.add("NONBONDED block: Illegal value for EPSLS (0.0 / >=1.0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.ewald_max_k_x <= 0 ||
            param.nonbonded.ewald_max_k_y <= 0 ||
            param.nonbonded.ewald_max_k_z <= 0)
      io::messages.add("NONBONDED block: Illegal value for NKX, NKY or NKZ (>0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.ewald_kspace_cutoff <= 0.0)
      io::messages.add("NONBONDED block: Illegal value for NK2 (>0.0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.p3m_grid_points_x <= 0 ||
            param.nonbonded.p3m_grid_points_y <= 0 ||
            param.nonbonded.p3m_grid_points_z <= 0)
      io::messages.add("NONBONDED block: Illegal value for NGA, NGB or NGC (>0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.p3m_grid_points_x % 2 != 0 ||
            param.nonbonded.p3m_grid_points_y % 2 != 0 ||
            param.nonbonded.p3m_grid_points_z % 2 != 0)
      io::messages.add("NONBONDED block: Illegal value for NGA, NGB or NGC (even)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.p3m_charge_assignment < 1 ||
            param.nonbonded.p3m_charge_assignment > 5)
      io::messages.add("NONBONDED block: Illegal value for NASORD (1..5)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.p3m_finite_differences_operator < 0 ||
            param.nonbonded.p3m_finite_differences_operator > 5)
      io::messages.add("NONBONDED block: Illegal value for NFDORD (0..5)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.p3m_mesh_alias <= 0)
      io::messages.add("NONBONDED block: Illegal value for NALIAS (>0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.accuracy_evaluation < 0)
      io::messages.add("NONBONDED block: Illegal value for NQEVAL (>=0)",
            "In_Parameter", io::message::error);

    if (param.pcouple.scale != math::pcouple_off
            && param.nonbonded.method == simulation::el_p3m
            && param.nonbonded.accuracy_evaluation == 0)
      io::messages.add("NONBONDED block: Pressure scaling but no quality evaluation of influence function."
            " Set NQEVAL > 0.", "In_Parameter", io::message::warning);

    if (param.pcouple.scale == math::pcouple_full_anisotropic && do_ls) {
      io::messages.add("NONBONDED block: Full ansiotropic pressure scaling "
              "could not be tested yet.", "In_Parameter", io::message::warning);
    }

    if (param.nonbonded.influence_function_rms_force_error <= 0.0)
      io::messages.add("NONBONDED block: Illegal value for FACCUR (>0.0)",
            "In_Parameter", io::message::error);

    if (param.nonbonded.influence_function_read ||
            param.nonbonded.influence_function_write)
      io::messages.add("NONBONDED block: Influence function IO not implemented."
            " Set NRDGRD and NWRGRD to 0.",
            "In_Parameter", io::message::error);

    if (param.nonbonded.lj_correction)
      io::messages.add("NONBONDED block: LJ long range correction not implemented."
            " Set NLRLJ to 0.",
            "In_Parameter", io::message::error);

    if (param.nonbonded.lj_solvent_density <= 0.0)
      io::messages.add("NONBONDED block: Illegal value for SLVDNS (>0.0)",
            "In_Parameter", io::message::error);

  } else {
    io::messages.add("no NONBONDED block", "In_Parameter", io::message::error);
    return;
  }
}

/**
 * @section sasa SASA block
 * @verbatim
 SASA
 # NTSASA
 # 0 : not used (default)
 # 1 : use SASA implicit solvent model
 # NTVOL
 # 0 : not used (default)
 # 1 : use VOLUME correction to SASA implicit solvent model (requires NTSASA = 1)
 # P_12 > 0, < 1 pair parameter for SASA reduction for first neighbours
 # P_13 > 0, < 1 pair parameter for SASA reduction for second neighbours
 # P_1X > 0, < 1 pair parameter for SASA reduction for third and higher neighbours
 # SIGMAV scaling parameter for volume energy term (kJ.mol^-1.nm^-3)
 # RSOLV > 0 radius of solvent molecule for SASA calculation (nm)
 # AS1 > 0 an atom with SASA below this contributes to the VOLUME correction (nm^2)
 # AS2 > 0 an atom with SASA above this is not considered for the VOLUME correction (nm^2)
 # atoms with AS1 < SASA < AS2 have a partial contribution determined by a switching function
 #   NTSASA      NTVOL       P_12      P_13     P_1X   SIGMAV  RSOlV    AS1    AS2
          1          1     0.8875    0.3516   0.3516     -100   0.14   0.01   0.02
 END
 @endverbatim
 */
void io::In_Parameter::read_SASA(simulation::Parameter & param, std::ostream & os) {
  DEBUG(8, "read SASA");

  std::vector<std::string> buffer;

  buffer = m_block["SASA"];

  // if there is no SASA term
  if (!buffer.size()) {
    return;
  }

  block_read.insert("SASA");

  _lineStream.clear();
  std::string s;
  _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

  int switch_sasa, switch_volume;
  _lineStream >> switch_sasa
              >> switch_volume
              >> param.sasa.p_12
              >> param.sasa.p_13
              >> param.sasa.p_1x
              >> param.sasa.sigma_v
              >> param.sasa.r_solv
              >> param.sasa.min_cut
              >> param.sasa.max_cut;

  if (_lineStream.fail()) {
    io::messages.add("bad line in SASA block", "In_Parameter", io::message::error);
  }

  switch(switch_sasa) {
    case 0:
      param.sasa.switch_sasa = false;
      break;
    case 1:
      param.sasa.switch_sasa = true;
      break;
    default:
      io::messages.add("SASA block: NTSASA has to be 0,1","In_Parameter", io::message::error);
      param.sasa.switch_sasa = false;
  }

  switch(switch_volume) {
    case 0:
      param.sasa.switch_volume = false;
      break;
    case 1:
      param.sasa.switch_volume = true;
      break;
    default:
      io::messages.add("SASA block: NTSASAVOL has to be 0,1","In_Parameter", io::message::error);
      param.sasa.switch_volume = false;
  }

  // check that vol not used without sasa
  if (!param.sasa.switch_sasa && param.sasa.switch_volume) {
    io::messages.add("SASA block: Cannot have NTSASAVOL without NTSASA",
            "In_Parameter", io::message::error);
  }

  if (param.sasa.p_12 < 0.0 || param.sasa.p_12 > 1.0 ||
          param.sasa.p_13 < 0.0 || param.sasa.p_13 > 1.0 ||
          param.sasa.p_1x < 0.0 || param.sasa.p_1x > 1.0) {
    io::messages.add("SASA block: Probabilities have to be 0.0..1.0",
            "In_Parameter", io::message::error);
  }

  if (param.sasa.r_solv <= 0.0) {
    io::messages.add("SASA block: RSOLV has to be >0.0",
            "In_Parameter", io::message::error);
  }

  // compute and store upper - lower cutoffs
  param.sasa.cut_diff = param.sasa.max_cut - param.sasa.min_cut;
  if (param.sasa.cut_diff < 0.0 ||
          param.sasa.max_cut < 0.0 ||
          param.sasa.min_cut < 0.0) {
    io::messages.add("SASA block: Cutoffs are negative or min. is larger than max.",
            "In_Parameter", io::message::error);
  }

}

/**
 * @section localelev LOCALELEV block
 * @verbatim
LOCALELEV
# NTLES 0,1 controls the use of local elevation.
#    0 : not used [default]
#    1 : local elevation is applied
# NLEPOT >= 0 number of umbrella potentials applied
# NTLESA 0..2 controls the reading of the potential definition
#    0 : read from input file (not supported)
#    1 : read from startup file
#    2 : read from special file (@lud)
# NTWLE >= 0 write umbrellas to trajectory evert NTWLEth step
# NLEPID[1..NLEPOT] IDs of the umbrella potentials
# NLEPFR[1..NLEPOT] 0,1 freeze the umbrella potential
#    0 : build up
#    1 : freeze
# NTLES  NLEPOT  NTLESA  NTWLE
      1       2       1      0
# NLEPID NLEPFR
       1      0
       2      1
END
@endverbatim
 */
void io::In_Parameter::read_LOCALELEV(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read LOCALELEV");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["LOCALELEV"];

  if (buffer.size()) {
    block_read.insert("LOCALELEV");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    int onoff, num, read;
    _lineStream >> onoff >> num >> read >> param.localelev.write;

    if (_lineStream.fail()) {
      io::messages.add("bad line in LOCALELEV block",
              "In_Parameter", io::message::error);
    }

    switch (onoff) {
      case 0:
        param.localelev.localelev = simulation::localelev_off;
        break;
      case 1:
        param.localelev.localelev = simulation::localelev_on;
        break;
      default:
        param.localelev.localelev = simulation::localelev_off;
        io::messages.add("LOCALELEV block: Bad value for NTLES (0,1)",
                "In_Parameter", io::message::error);
    }

    if (num < 0) {
      io::messages.add("LOCALELEV block: Bad value for NLEPOT (>=0)",
              "In_Parameter", io::message::error);
      return;
    }

    switch (read) {
      case 1:
        param.localelev.read = false;
        break;
      case 2:
        param.localelev.read = true;
        break;
      default:
        param.localelev.read = false;
        io::messages.add("LOCALELEV block: Bad value for NTLESA (1,2)",
                "In_Parameter", io::message::error);
    }

    if (param.localelev.write < 0) {
      io::messages.add("LOCALELEV block: Bad value for NTWLE (>=0)",
              "In_Parameter", io::message::error);
      param.localelev.write = 0;
    }

    // read the umbrellas
    for (int i = 0; i < num; ++i) {
      int id, f;
      _lineStream >> id;
      if (_lineStream.fail()) {
        std::ostringstream msg;
        msg << "LOCALELEV block: Bad value for NLEPID[" << i + 1 << "]";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        return;
      }
      _lineStream >> f;
      if (_lineStream.fail() || (f != 0 && f != 1)) {
        std::ostringstream msg;
        msg << "LOCALELEV block: Bad value for NLEPFR[" << i + 1 << "] (0,1)";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        return;
      }

      bool freeze = (f == 1) ? true : false;
      if (param.localelev.umbrellas.find(id) == param.localelev.umbrellas.end()) {
        param.localelev.umbrellas[id] = !freeze;
      } else {
        std::ostringstream msg;
        msg << "LOCALELEV block: duplicated umbrella potential ID (" << id << ")";
        io::messages.add(msg.str(), "In_Parameter", io::message::error);
        return;
      }
    } // for umbrellas
  } // if block
}
/**
 * @section bsleusparam BSLEUS block
 * @verbatim
BSLEUS
#
# The general settings for the B&S-LEUS algorithm
# BSLEUS:   Dow we use B&S-LEUS?
#   0:          Don'use it
#   1:          Use it
# BUILD:    Are we building?
#   0:          No
#   1:          Yes
# WRITE:    >= 0 Do we write the energies and forces of the Umbrella?
#   == 0:          No
#   > 0:           Every nth step
# 
# BSLEUS    BUILD   WRITE
  1         1       0
END
@endverbatim
 */
void io::In_Parameter::read_BSLEUS(simulation::Parameter& param, std::ostream& os)
{
  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["BSLEUS"];

  if (buffer.size()) {
    block_read.insert("BSLEUS");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    
    int use_bsleus, build, write;
    _lineStream >> use_bsleus >> build >> write;
    if (_lineStream.fail()){
      io::messages.add("Bad BSLEUS block!", "In_Parameter", io::message::error);
      return;
    }
    
    if (use_bsleus == 0 || use_bsleus == 1){
      param.bsleus.bsleus = use_bsleus ? simulation::bsleus_on : 
                                         simulation::bsleus_off;
    }
    else {
      io::messages.add("BSLEUS: Bad value for BSLEUS", "In_Parameter", 
              io::message::error);
      return;
    }
    if (build == 0 || build == 1){
      param.bsleus.building = build ? true : false;
    }
    else {
      io::messages.add("BSLEUS: Bad value for BUILD", "In_Parameter", 
              io::message::error);
      return;
    }
    /*if (forceConstIncr > 0.0){
      param.bsleus.forceConstantIncrement = forceConstIncr;
    }
    else {
      io::messages.add("BSLEUS: Bad value for MEMKLE", "In_Parameter", 
              io::message::error);
      return;
    }*/
    if (write >= 0){
      param.bsleus.write = write;
    }
    else {
      io::messages.add("BSLEUS: Bad value for WRITE", "In_Parameter", 
              io::message::error);
      return;
    }
  }
}
/**
 * @section electric ELECTRIC block
 * @verbatim
 ELECTRIC
# FIELD 0,1 controls the use of applied electric field.
#    0 : not used [default]
#    1 : electric field is applied
# DIPOLE 0,1 controls the calculation of box dipole.
#    0 : not used [default]
#    1 : box dipole is calculated and written to special trajectory
# CURRENT 0,1 controls the calculation of electric currents.
#    0 : not used [default]
#    1 : electric current is calculated and written to special trajectory
# ELECTRIC FIELD COMPONENTS (EF_x, EF_y, EF_z)
# 0.0 0.0 0.0
# DIPGRP 0..2 controls the groups considered for box dipole calculation
#    0 : solute only
#    1 : solvent only
#    2 : all
# NTWDIP >= 0 write dipole box every NTWDIPth step
# NCURGRP >=0 number of current groups
# CURGRP [1..NCURGRP] last atom of the group
#  FIELD  DIPOLE CURRENT
       1       1       1
#   EF_x    EF_y    EF_z
     0.0     0.0     0.0
# DIPGRP  NTWDIP
       0       1
# NTWCUR  NCURGRP   CURGRP[1]   CURGRP[2]
       1       2        100        1000
END
@endverbatim
 */
void io::In_Parameter::read_ELECTRIC(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read ELECTRIC");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["ELECTRIC"];

  if (buffer.size()) {
    block_read.insert("ELECTRIC");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int field, dipole, current, dipgrp, ncurgrp;

    _lineStream >> field >> dipole >> current >> param.electric.Ef_x
            >> param.electric.Ef_y >> param.electric.Ef_z
            >> dipgrp >> param.electric.dip_write >> param.electric.cur_write
            >> ncurgrp;

    if (_lineStream.fail()) {
      io::messages.add("bad line in ELECTRIC block",
              "In_Parameter", io::message::error);
      return;
    }

    switch (field) {
      case 0:
        param.electric.electric = simulation::electric_off;
        param.electric.Ef_x = 0.0;
        param.electric.Ef_y = 0.0;
        param.electric.Ef_z = 0.0;
        break;
      case 1:
      {
        param.electric.electric = simulation::electric_on;
        if (param.electric.Ef_x == param.electric.Ef_y &&
            param.electric.Ef_y == param.electric.Ef_z &&
            param.electric.Ef_z == 0.0)
          io::messages.add("Electric field enabled, but all components are zero",
                "In_Parameter", io::message::error);
        if (param.nonbonded.method != simulation::el_reaction_field &&
                param.nonbonded.rf_epsilon != 1.0)
          io::messages.add("To use electric field together with Ewald or P3M, eps_rf must be 1.0",
                "In_Parameter", io::message::error);
        break;
      }
      default:
        param.electric.electric = simulation::electric_off;
        param.electric.Ef_x = 0.0;
        param.electric.Ef_y = 0.0;
        param.electric.Ef_z = 0.0;
        io::messages.add("ELECTRIC block: Bad value for FIELD (0,1)",
                "In_Parameter", io::message::error);
    }

    switch (dipole) {
      case 0:
      {
        param.electric.dipole = false;
        param.electric.dip_write = 0;
        break;
      }
      case 1:
      {
        param.electric.dipole = true;
        if (dipgrp < 0 || dipgrp > 2)
          io::messages.add("ELECTRIC block: DIPGRP should be within 0..2",
                "In_Parameter", io::message::error);
        break;
      }
      default:
      {
        param.electric.dipole = false;
        param.electric.dip_write = 0;
        io::messages.add("ELECTRIC block: Bad value for DIPOLE (0,1)",
                "In_Parameter", io::message::error);
      }
    }

    switch (dipgrp) {
      case 0:
      {
        param.electric.dip_groups = 0;
        break;
      }
      case 1:
      {
        param.electric.dip_groups = 1;
        break;
      }
      case 2:
      {
        param.electric.dip_groups = 2;
        break;
      }
      default:
      {
        io::messages.add("ELECTRIC block: DIPGRP should be within 0..2",
                "In_Parameter", io::message::error);
        break;
      }
    }

    switch (current) {
      case 0:
      {
        param.electric.current = false;
        param.electric.cur_write = 0;
        break;
      }
      case 1:
      {
        param.electric.current = true;
        if (ncurgrp < 0)
          io::messages.add("ELECTRIC block: CURGRP should be >0",
                "In_Parameter", io::message::error);
        param.electric.cur_groups = ncurgrp;
        break;
      }
      default:
      {
        param.electric.current = false;
        param.electric.cur_write = 0;
        io::messages.add("ELECTRIC block: Bad value for CURRENT (0,1)",
                "In_Parameter", io::message::error);
      }
    }

    if (param.electric.current != false) {
      // TO READ THE ELECTRIC GROUPS
      if (_lineStream.fail() || param.electric.cur_groups == (unsigned int) 0)
        io::messages.add("CURRENT enabled, but number of CURRENT (NCURGRP) groups is zero",
              "In_Parameter", io::message::error);
      else {
        unsigned int temp;
        for (unsigned int i = 0; i < param.electric.cur_groups; i++) {
          _lineStream >> temp;
          if (_lineStream.fail()) {
            io::messages.add("CURRENT enabled, but CURGRP[i] < or > NCURGRP",
                    "In_Parameter", io::message::error);
            break;
          }
          param.electric.current_group.push_back(temp);
        }
      }
    }
  } // if block
}

/**
 * @section nemd NEMD block
 * @verbatim
 NEMD
# NEMD 0,1 controls the use of non-equilibrium molecular dynamics.
#    0 : not used [default]
#    1 : nemd is used
# PROPERTY 0- select property to calculate
#    0 : viscosity      
# METHOD 0- select method of NEMD.
#    0 : periodic perturbation method (PPM)
#    1 : internal reservoir method (IRM)
# SLABNUM >=1 number of slabs used in the discretization along z-direction.
#             the effective number is 2xSLABNUM due to periodicity
# PERTFRQ >=1 perturbation frequency: apply perturbation every PERTFRQth timestep
#             [this flag is ignored by the PPM method, but a value must be provided]
# AMPLI   >=0 amplitude of applied field
#             [this flag is ignored by the IRM method, but a value must be provided]
# STDYAFT >=0 first STDYAFTth steps do not contribute for accumulated averages
# WRITE >=1 write flux and average velocities to special trajectory every WRITEth timestep
# NEMD     PROPERTY  METHOD
    1         0        0
# SLABNUM  PERTFRQ    AMPLI   STDYAFT   WRITE
     10       20       10      1000     200
END
@endverbatim
 */
void io::In_Parameter::read_NEMD(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read NEMD");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["NEMD"];

  if (buffer.size()) {
    block_read.insert("NEMD");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int nemd, property, method, slabnum, pertfrq, stdyaft, write;
    double ampli;

    _lineStream >> nemd >> property >> method >> slabnum >> pertfrq >> ampli >> stdyaft >> write;

    if (_lineStream.fail()) {
      io::messages.add("bad line in NEMD block",
              "In_Parameter", io::message::error);
      return;
    }



    switch (nemd) {
      case 0:
        param.nemd.nemd = simulation::nemd_off;
        break;
      case 1:
      {
        param.nemd.nemd = simulation::nemd_on;
        break;
      }
      default:
        param.nemd.nemd = simulation::nemd_off;
        io::messages.add("NEMD block: Bad value for NEMD (0,1)",
                "In_Parameter", io::message::error);
    }

    switch (property){
      case 0:
        param.nemd.property = 0;
        break;
      default:
        io::messages.add("NEMD block: Bad value for PROPERTY, for now only viscosity (0) is available",
                "In_Parameter", io::message::error);
    }


    switch (method) {
      case 0:
      {
        param.nemd.method = 0;
        if (slabnum <=0 || ampli <=0)
          io::messages.add("NEMD block: PPM method used, but found invalid values for SLABNUM and AMPLI",
                "In_Parameter", io::message::error);
        break;
      }
      case 1:
      {
        param.nemd.method = 1;
        if (slabnum <=0 || pertfrq <=0)
          io::messages.add("NEMD block: IRM method used, but found invalid values for SLABNUM and PERTFRQ",
                "In_Parameter", io::message::error);
        break;
      }
      default:
      {
        io::messages.add("NEMD block: Not a valid method(0)",
                "In_Parameter", io::message::error);

      }
    }
 if (write <0)
      io::messages.add("NEMD block: invalid value for WRITE (>=0)",
                "In_Parameter", io::message::error);
    
    param.nemd.slabnum = slabnum;
    param.nemd.pertfrq = pertfrq;
    param.nemd.ampbath = ampli;
    param.nemd.stdyaft = stdyaft;
    param.nemd.write = write;
    

   

    

 
  } // if block
}


/**
 * @section multigradient MULTIGRADIENT block
 * @verbatim
MULTIGRADIENT
# NTMGRE 0,1 enable multiple gradients
#    0: disable gradients
#    1: enable gradients
# NTMGRP 0..3 print of curves
#    0: don't print
#    1: plot the curves
#    2: print that values of the curves
#    3: plot and print the curves
# NTMGRN >= 0 number of gradients
# MGRVAR: vairable name to affect
# MGRFRM: functional form of the curve
#    0: linear interpolation between control points
#    1: cubic spline interpolation between control points
#    2: Bezier curve
#    3: Oscillation: A sin[2Pi/T (t - dt)] + b
#       Note: MGRNCP is 2. A = MGRCPT[1] T = MGRCPV[1] dt = MGRCPT[2] b = MGRCPV[2]
# MGRNCP >= 2: number of control points
# MGRCPT >= 0: time of the control point
# MGRCPV: value of the control point
#
# NTMGRE NTMGRP
       1      1
# NTMGRN
       2
# MGRVAR MGRFRM MGRNCP
  TEMP[0]     0      2
# MGRCPT MGRCPV
  0.0    60.0
  80.0   300.0
# MGRVAR MGRFRM MGRNCP
  CPOR        2      4
# MGRCPT MGRCPV
  0.0    2.5E5
  0.0    2.5E1
 20.0    0.0
 80.0    0.0
END
@endverbatim
 */
void io::In_Parameter::read_MULTIGRADIENT(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read MULTIGRADIENT");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["MULTIGRADIENT"];

  if (buffer.size()) {
    block_read.insert("MULTIGRADIENT");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int enable, plot, num;
    _lineStream >> enable >> plot >> num;
    if (_lineStream.fail()) {
      io::messages.add("Bad line in MULTIGRADIENT block.",
                    "In_Parameter", io::message::error);
      return;
    }

    switch(enable) {
      case 0:
        param.multigradient.multigradient = false;
        break;
      case 1:
        param.multigradient.multigradient = true;
        break;
      default:
        param.multigradient.multigradient = false;
        io::messages.add("MULTIGRADIENT block: NTMGRE must be 0 or 1.",
                    "In_Parameter", io::message::error);
        return;
    }

    switch(plot) {
      case 0:
        param.multigradient.print_graph = false;
        param.multigradient.print_curve = false;
        break;
      case 1:
        param.multigradient.print_graph = true;
        param.multigradient.print_curve = false;
        break;
      case 2:
        param.multigradient.print_graph = false;
        param.multigradient.print_curve = true;
        break;
      case 3:
        param.multigradient.print_graph = true;
        param.multigradient.print_curve = true;
        break;
      default:
        param.multigradient.print_graph = false;
        param.multigradient.print_curve = false;
        io::messages.add("MULTIGRADIENT block: NTMGRP must be 0..3.",
                    "In_Parameter", io::message::error);
        return;
    }

    if (num < 0) {
      io::messages.add("MULTIGRADIENT block: NTMGRN must be >= 0",
                    "In_Parameter", io::message::error);
      return;
    }

    if (num == 0 && param.multigradient.multigradient) {
      io::messages.add("MULTIGRADIENT block: NTMGRN must be > 0 for NTMGRE=1",
                    "In_Parameter", io::message::error);
      return;
    }

    // read the gradient
    for(int i = 0; i < num; ++i) {
      int funct_form, num_p;
      std::string var;
      _lineStream >> var >> funct_form >> num_p;
      if (_lineStream.fail()) {
        io::messages.add("Bad line in MULTIGRADIENT block.",
                "In_Parameter", io::message::error);
        return;
      }


      if (num_p < 2) {
        io::messages.add("MULTIGRADIENT block: MGRNCP must be >= 2.",
                    "In_Parameter", io::message::error);
        return;
      }
      std::vector<std::pair<double, double> > points;
      for(int p = 0; p < num_p; ++p) {
        double t, v;
        _lineStream >> t >> v;
        if (_lineStream.fail()) {
          io::messages.add("Bad line in MULTIGRADIENT block.",
                  "In_Parameter", io::message::error);
          return;
        }
        points.push_back(std::pair<double, double>(t,v));
      }

      param.multigradient.variable.push_back(var);
      param.multigradient.functional_form.push_back(funct_form);
      param.multigradient.control_points.push_back(points);
    } // for gradients
  } // if block
}

/**
 * @section addecouple ADDECOUPLE block
 * @verbatim
ADDECOUPLE
# ADGR    >= 0 number of addiabatic decoupling groups
# ADSTART first atom of the addiabatic decoupling group
# ADEND   last atom of the addiabatic decoupling group
# SM      scaling factor mass
# SV      scaling factor potential energy function
# ST      scaling factor temperature
# TIR     which temperature bath to scale
#  1      translational
#  2      internal-rotatinal
#  3      both
# TMF     tau for calculating mean field
# STAD    printing average to special trajectory # ADGR
      2
# ADSTART ADEND SM SV  ST TIR
  1       1500 10   1  0  1
  1501    3000  1  10  1  3 
# TMF STAD
  0.1 1000
END
@endverbatim
 */
void io::In_Parameter::read_ADDECOUPLE(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read ADDECOUPLE");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["ADDECOUPLE"];

  if (buffer.size()) {
    block_read.insert("ADDECOUPLE");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));
    unsigned int adgr, adstart, eg, tg=0, adend, write, tir;
    double sm, sv, st, tmf;
    _lineStream >> adgr;

    param.addecouple.adgr = adgr;

    if (_lineStream.fail()) {
      io::messages.add("bad line in ADDECOUPLE block",
              "In_Parameter", io::message::error);
    }
    for (unsigned int i = 0; i < adgr; ++i) {
      _lineStream >> adstart >> adend >> sm >> sv >> st >> tir;
      if (adstart < 0.0 || adend < adstart) {
        io::messages.add("ADDECOUPLE block: illegal value for adstart or adend",
                "In_Parameter", io::message::error);
      }
      if (tir != 1 && tir != 2 && tir != 3) {
        io::messages.add("ADDECOUPLE block: illegal value for tir (not 1,2 or 3)",
                "In_Parameter", io::message::error);
      }
      if (st != 1 && param.multibath.couple == false)
        io::messages.add("ADDECOUPLE block: sT bigger 1, but no temperature scaling",
              "In_Parameter", io::message::error);
      //check whether the group is also a temperature group
      if (param.multibath.couple) {
        int addc_bath_index;
        if (param.multibath.multibath.bath_index().size() < adgr && st != 1)
          io::messages.add("ADDECOUPLE block: sT bigger 1, but temperature group and adiabatic decouling group not equivalent",
                "In_Parameter", io::message::error);
        else {
          int check_group = -1;
          for (unsigned int bath_i = 0; bath_i < param.multibath.multibath.bath_index().size(); ++bath_i) {
            if (bath_i > 0
                    && adend - 1 == param.multibath.multibath.bath_index()[bath_i].last_atom
                    && adstart - 2 == param.multibath.multibath.bath_index()[bath_i - 1].last_atom) {
              check_group = 1;
              addc_bath_index = bath_i;
              tg=bath_i;
            } else if (bath_i == 0
                    && adend - 1 == param.multibath.multibath.bath_index()[0].last_atom
                    && adstart == 1) {
              check_group = 1;
              addc_bath_index = bath_i;
              tg=bath_i;
            }
          }
          if (st == 1)
            check_group = 1;
          if (check_group == -1)
            io::messages.add("ADDECOUPLE block: sT bigger 1, but temperature group and adiabatic decouling group not equivalent",
                  "In_Parameter", io::message::error);
            //check whether com and ir are handled "correctly"
          else {
            if (param.multibath.multibath.bath_index()[addc_bath_index].com_bath
                    == param.multibath.multibath.bath_index()[addc_bath_index].ir_bath
                    && tir != 3)
              io::messages.add("ADDECOUPLE block: seperate scaling for this temperature group not possible",
                    "In_Parameter", io::message::error);
            else if (st != 1) {
              int com_bath = param.multibath.multibath.bath_index()[addc_bath_index].com_bath;
              int ir_bath = param.multibath.multibath.bath_index()[addc_bath_index].ir_bath; //scale temperature
              if (com_bath == ir_bath && tir == 3)
                param.multibath.multibath[ir_bath].temperature *= st;
              else if (tir == 1)
                param.multibath.multibath[com_bath].temperature *= st;
              else if (tir == 2)
                param.multibath.multibath[ir_bath].temperature *= st;
              else if (tir == 3) {
                param.multibath.multibath[com_bath].temperature *= st;
                param.multibath.multibath[ir_bath].temperature *= st;
              }
              else
                io::messages.add("ADDECOUPLE block: scaling for this temperature group not possible",
                      "In_Parameter", io::message::error);

            }
          }
        }

      }
      //check whether group is also an energy group
      int check_group = -1;
      if (adstart == 1 && (param.force.energy_group[0] + 1) == adend){
        check_group = 1;
        eg=0;
      }
      else
        for (unsigned int energy_i = 0; energy_i < param.force.energy_group.size() - 1; ++energy_i) {
          if (param.force.energy_group[energy_i] + 2 == adstart 
                  && param.force.energy_group[energy_i + 1] + 1 == adend){
            check_group = 1;
            eg=i;
          }
        }
      if (check_group == -1)
        io::messages.add("ADDECOUPLE block: energy and adiabatic groups are not identical", "In_Parameter", io::message::error);
      param.addecouple.add_adc(adstart - 1, adend - 1, sm, sv, st, tir, eg, tg);
    }
    if (adgr > 0) {
      _lineStream >> tmf >> write;
      if (tmf < 0.0 || write < 0.0) {
        io::messages.add("ADDECOUPLE block: illegal value for TMF or STAD",
                "In_Parameter", io::message::error);
      }

      param.addecouple.tmf = tmf;
      param.addecouple.write = write;
    }

    if (_lineStream.fail()) {
      io::messages.add("bad line in ADDECOUPLE block",
              "In_Parameter", io::message::error);
    }
  }
}



// two helper data types to simply unsupported block handling

enum unsupported_block_type {
  ub_unknown, // I know it and know that I don't use it but I have no idea why
  ub_renamed, // it was renamed. e.g. from previous versions
  ub_promd, // it is a PROMD block. Tell alternative if there is any
  ub_g96 // it is a G96 block. Tell alternative if there is any
};

// give a valid block name as an alternative and it will tell the user to
// use it.

struct unsupported_block {

  unsupported_block() :
  alternative(""), type(ub_unknown) {
  }

  unsupported_block(std::string a, unsupported_block_type t) :
  alternative(a), type(t) {
  }

  std::string alternative;
  unsupported_block_type type;
};

void io::In_Parameter::read_known_unsupported_blocks() {
  std::map<std::string, unsupported_block> ub;
  // add all those unknown blocks
  ub["ANATRAJ"] = unsupported_block("READTRAJ", ub_renamed);
  ub["MINIMISE"] = unsupported_block("ENERGYMIN", ub_renamed);
  ub["STOCHASTIC"] = unsupported_block("STOCHDYN", ub_renamed);
  ub["BOUNDARY"] = unsupported_block("BOUNDCOND", ub_renamed);
  ub["THERMOSTAT"] = unsupported_block("MULTIBATH", ub_promd);
  ub["TCOUPLE"] = unsupported_block("MULTIBATH", ub_g96);
  ub["BAROSTAT"] = unsupported_block("PRESSURESCALE", ub_promd);
  ub["VIRIAL"] = unsupported_block("PRESSURESCALE", ub_promd);
  ub["PCOUPLE"] = unsupported_block("PRESSURESCALE", ub_g96);
  ub["PCOUPLE03"] = unsupported_block("PRESSURESCALE", ub_renamed);
  ub["GEOMCONSTRAINT"] = unsupported_block("CONSTRAINT", ub_promd);
  ub["SHAKE"] = unsupported_block("CONSTRAINT", ub_g96);
  ub["PATHINT"] = unsupported_block("", ub_promd);
  ub["NEIGHBOURLIST"] = unsupported_block("PAIRLIST", ub_promd);
  ub["PLIST"] = unsupported_block("PAIRLIST", ub_g96);
  ub["PLIST03"] = unsupported_block("PAIRLIST", ub_renamed);
  ub["LONGRANGE"] = unsupported_block("NONBONDED", ub_g96);
  ub["START"] = unsupported_block("INITIALISE", ub_g96);
  ub["OVERALLTRANSROT"] = unsupported_block("COMTRANSROT", ub_promd);
  ub["CENTREOFMASS"] = unsupported_block("COMTRANSROT", ub_g96);
  ub["POSREST"] = unsupported_block("POSITIONRES", ub_g96);
  ub["DISTREST"] = unsupported_block("DISTANCERES", ub_g96);
  ub["DIHEREST"] = unsupported_block("DIHEDRALRES", ub_g96);
  ub["J-VAL"] = unsupported_block("JVALUERES", ub_g96);
  ub["J-VAL03"] = unsupported_block("JVALUERES", ub_renamed);
  ub["PERTURB"] = unsupported_block("PERTURBATION", ub_g96);
  ub["PERTURB03"] = unsupported_block("PERTURBATION", ub_renamed);
  ub["UMBRELLA"] = unsupported_block("", ub_promd);
  ub["PRINT"] = unsupported_block("PRINTOUT", ub_g96);
  ub["WRITE"] = unsupported_block("WRITETRAJ", ub_g96);
#ifdef NDEBUG
  ub["DEBUG"] = unsupported_block("--enable-debug at compile time and "
          "the @verb argument", ub_promd);
#else
  ub["DEBUG"] = unsupported_block("the @verb argument", ub_promd);
#endif
  ub["FOURDIM"] = unsupported_block("", ub_g96);
  ub["LOCALELEVATION"] = unsupported_block("LOCALELEV", ub_g96);
  ub["SUBMOLECULES"] = unsupported_block("SOLUTEMOLECULES and moved to "
          "the topology", ub_renamed);
  ub["FORCEFIELD"] = unsupported_block("COVALENTFORM", ub_renamed);
  ub["PSCALE"] = unsupported_block("PERSCALE", ub_renamed);
  ub["REPLICA03"] = unsupported_block("REPLICA", ub_renamed);

  std::map<std::string, unsupported_block>::const_iterator
  it = ub.begin(),
          to = ub.end();

  // loop over unsupported blocks;
  for (; it != to; ++it) {
    // if it is present
    if (m_block[it->first].size()) {
      block_read.insert(it->first);

      std::ostringstream msg;
      msg << it->first << " block";

      switch (it->second.type) {
        case ub_renamed:
          msg << " was renamed to " << it->second.alternative;
          break;
        case ub_promd:
          msg << " is PROMD specific.";
          if (it->second.alternative != "")
            msg << " Use " << it->second.alternative << " instead.";
          break;
        case ub_g96:
          msg << " is GROMOS96 specific.";
          if (it->second.alternative != "")
            msg << " Use " << it->second.alternative << " instead.";
          break;
        default: // don't know what to do.
          msg << " is known to be not supported.";
      }

      io::messages.add(msg.str(), "In_Parameter", io::message::error);
    }
  }
}

/**
 * @section qmmmb QMMM block
 * @verbatim
QMMM
# NTQMMM 0,1 apply QM/MM
#    0: do not apply QM/MM
#    1: apply QM/MM
# NTQMSW 0 QM software package to use
#    0: MNDO
#    1: Turbomole
# RCUTQ >= 0.0 cutoff for inclusion of MM charge groups
#     0.0: include all atoms
#    >0.0: include atoms of charge groups closer than RCUTQ
#          to QM zone.
# NTWQMMM >= 0 write QM/MM related data to special trajectory
#    0: do not write
#   >0: write every NTWQMMMth step
#
# NTQMMM  NTQMSW  RCUTQ  NTWQMMM
       1       0    0.0        0
END

@endverbatim
 */
void io::In_Parameter::read_QMMM(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read QMMM");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["QMMM"];

  if (buffer.size()) {
    block_read.insert("QMMM");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int enable, software, write;
    double cutoff;
    _lineStream >> enable >> software >> cutoff >> write;
    if (_lineStream.fail()) {
      io::messages.add("Bad line in QMMM block.",
                    "In_Parameter", io::message::error);
      return;
    }

    switch(enable) {
      case 0:
        param.qmmm.qmmm = simulation::qmmm_off;
        break;
      case 1:
        param.qmmm.qmmm = simulation::qmmm_on;
        break;
      default:
        param.qmmm.qmmm = simulation::qmmm_off;
        io::messages.add("QMMM block: NTQMMM must be 0 or 1.",
                    "In_Parameter", io::message::error);
        return;
    }

    switch (software) {
      case 0:
        param.qmmm.software = simulation::qmmm_software_mndo;
        break;
      case 1:
        param.qmmm.software = simulation::qmmm_software_turbomole;
        break;
      default:
        param.qmmm.software = simulation::qmmm_software_mndo;
        io::messages.add("QMMM block: NTQMSW invalid choice of software",
                    "In_Parameter", io::message::error);
        return;
    }
    
    if (cutoff < 0.0) {
      io::messages.add("QMMM block: RCUTQ must be >= 0.0",
                    "In_Parameter", io::message::error);
      return;
    }
    param.qmmm.cutoff = cutoff;

    if (write < 0) {
      io::messages.add("QMMM block: NTWQMMM must be >= 0",
                    "In_Parameter", io::message::error);
      return;
    }
    param.qmmm.write = write;
  } // if block
}

/**
 * @section symres SYMRES block
 * @verbatim
SYMRES
# NTSYM 0..2 apply symmetry restraints
#    0: do not apply symmetry restraints
#    1: apply symmetry restraints
#    2: apply symmetry constraints
# CSYM >= 0.0 force constants
#
# NTSYM     CSYM  
       1     0.0 
END

@endverbatim
 */
void io::In_Parameter::read_SYMRES(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read SYMRES");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["SYMRES"];

  if (buffer.size()) {
    block_read.insert("SYMRES");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int enable;
    double fc;
    _lineStream >> enable >> fc;
    if (_lineStream.fail()) {
      io::messages.add("Bad line in SYMRES block.",
                    "In_Parameter", io::message::error);
      return;
    }

    switch(enable) {
      case 0:
        param.symrest.symrest = simulation::xray_symrest_off;
        break;
      case 1:
        param.symrest.symrest = simulation::xray_symrest_ind;
        break;
      case 2:
        param.symrest.symrest = simulation::xray_symrest_constr;
        break;
      default:
        param.symrest.symrest = simulation::xray_symrest_off;
        io::messages.add("SYMRES block: NTSYM must be 0..2.",
                    "In_Parameter", io::message::error);
        return;
    }

    
    if (fc < 0.0) {
      io::messages.add("SYMRES block: CSYM must be >= 0.0",
                    "In_Parameter", io::message::error);
      return;
    }
    param.symrest.force_constant = fc;

  } // if block
}

