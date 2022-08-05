/**
 * @file in_parameter.cc
 * implements methods of In_Parameter
 */
#include <sstream>
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
  read_STOCHDYN(param); // has to be read in before REPLICA_EDS
  read_REPLICA_EDS(param); // has to be read in after MULTIBATH
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
  read_ANGLERES(param);
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
  read_EWARN(param);
  read_MULTISTEP(param);
  read_CHEMICALMONTECARLO(param);
  read_POLARISE(param);
  read_RANDOMNUMBERS(param);
  read_EDS(param);
  read_AEDS(param); // needs to be called after EDS
  read_LAMBDAS(param); // needs to be called after FORCE
  read_PRECALCLAM(param); // ANITA
  read_LOCALELEV(param);
  read_BSLEUS(param);
  read_ELECTRIC(param);
  read_SASA(param);
  read_ADDECOUPLE(param); // needs to be called after MULTIBATH and FORCE
  read_NEMD(param);
  read_MULTIGRADIENT(param);
  read_QMMM(param);
  read_SYMRES(param);
  read_AMBER(param);

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
 * @snippet snippets/snippets.cc SYSTEM
 */
void io::In_Parameter::read_SYSTEM(simulation::Parameter &param,
                                   std::ostream & os) {
    DEBUG(8, "reading SYSTEM");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "SYSTEM\n";
    exampleblock << "# NPM : protein molecules (0, 1)\n";
    exampleblock << "# NSM : solvent molecules (>= 0)\n";
    exampleblock << "#      NPM      NSM\n";
    exampleblock << "         1        0\n";
    exampleblock << "END\n";

    std::string blockname = "SYSTEM";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert("SYSTEM");

        block.get_next_parameter("NPM", param.system.npm, "", "0,1");
        block.get_next_parameter("NSM", param.system.nsm, ">=0", "");

        if (param.system.npm == 0 && param.system.nsm == 0)
            io::messages.add("SYSTEM block: no molecules: both NSM and NPM are 0",
                             "io::In_Parameter", io::message::error);

        block.get_final_messages();
    } else {
        param.system.nsm = 0;
        param.system.npm = 0;
        return;
    }
}

/**
 * @section energymin ENERGYMIN block
 * @snippet snippets/snippets.cc ENERGYMIN
 */
void io::In_Parameter::read_ENERGYMIN(simulation::Parameter &param,
                                      std::ostream & os) {
    DEBUG(8, "reading ENERGYMIN");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "ENERGYMIN\n";
    exampleblock << "# NTEM: 0..3 controls energy minimisation mode\n";
    exampleblock << "#       0: do not do energy minimisation (default)\n";
    exampleblock << "#       1: steepest-descent minimisation\n";
    exampleblock << "#       2: Fletcher-Reeves conjugate-gradient minimisation\n";
    exampleblock << "#       3: Polak-Ribiere conjugate-gradient minimisation\n";
    exampleblock << "# NCYC: >0 number of steps before resetting the conjugate-gradient search direction\n";
    exampleblock << "#       =0 reset only if the energy grows in the search direction\n";
    exampleblock << "# DELE: >0.0 energy threshold for convergence\n";
    exampleblock << "#       >0.0 (conjugate-gradient) RMS force threshold for convergence\n";
    exampleblock << "# DX0: >0.0 initial step size\n";
    exampleblock << "# DXM: >0.0 maximum step size\n";
    exampleblock << "# NMIN >0 minimum number of minimisation steps\n";
    exampleblock << "# FLIM >=0.0 limit force to maximum value (FLIM > 0.0 is not recommended)\n";
    exampleblock << "# CGIM >0 (conjugate-gradient) maximum number of cubic interpolations per step\n";
    exampleblock << "# CGIC >0.0 (conjugate-gradient) displacement threshold after interpolation\n";
    exampleblock << "#     NTEM    NCYC    DELE    DX0     DXM    NMIN    FLIM\n";
    exampleblock << "         1       0     0.1   0.01    0.05     100     0.0\n";
    exampleblock << "# ---- OR: example for NTEM > 1:\n";
    exampleblock << "#     NTEM    NCYC    DELE    DX0     DXM    NMIN    FLIM    CGIM    CGIC\n";
    exampleblock << "         3       0    1e-3   5e-6    5e-4     100     0.0       3    1e-4\n";
    exampleblock << "END\n";

    std::string blockname = "ENERGYMIN";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert("ENERGYMIN");

        block.get_next_parameter("NTEM", param.minimise.ntem, "", "0,1,2,3");
        block.get_next_parameter("NCYC", param.minimise.ncyc, ">=0", "");
        block.get_next_parameter("DELE", param.minimise.dele, ">0", "");
        block.get_next_parameter("DX0", param.minimise.dx0, ">0", "");
        std::string str_dx0=io::to_string(param.minimise.dx0);
        block.get_next_parameter("DXM", param.minimise.dxm, ">="+str_dx0, "");
        block.get_next_parameter("NMIN", param.minimise.nmin, ">0", "");
        block.get_next_parameter("FLIM", param.minimise.flim, ">=0", "");
        if (param.minimise.ntem >= 2) {
            block.get_next_parameter("CGIM", param.minimise.cgim, ">0", "");
            block.get_next_parameter("CGIC", param.minimise.cgic, ">0", "");
        }

        if (param.minimise.ntem == 1 && param.minimise.ncyc > 0)
           io::messages.add("ENERGYMIN block: NCYC > 0 has no effect for steepest descent",
               "io::In_Parameter",
               io::message::warning);

        if (param.minimise.flim > 0)
            io::messages.add("ENERGYMIN: FLIM > 0 may result in "
                             "failure of the minimisation procedure."
                             " Only to be used in special cases.",
                             "io::In_Parameter",
                             io::message::warning);

        block.get_final_messages();
    }
} // ENERGYMIN

/**
 * @section step STEP block
 * @snippet snippets/snippets.cc STEP
 */
void io::In_Parameter::read_STEP(simulation::Parameter &param,
                                 std::ostream & os) {
    DEBUG(8, "reading STEP");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "STEP\n";
    exampleblock << "#   NSTLIM  >0 number of steps\n";
    exampleblock << "#   T       >=0 initial time\n";
    exampleblock << "#           -1  read time from configuration file\n";
    exampleblock << "#   DT      >0 time step\n";
    exampleblock << "#\n";
    exampleblock << "#   NSTLIM         T        DT\n";
    exampleblock << "       100       0.0     0.005\n";
    exampleblock << "END\n";


    std::string blockname = "STEP";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("NSTLIM", param.step.number_of_steps, ">0", "");
        block.get_next_parameter("T", param.step.t0, ">=0", "-1");
        block.get_next_parameter("DT", param.step.dt, ">0", "");

        block.get_final_messages();
    }
}

/**
 * @section constraint CONSTRAINT block
 * @snippet snippets/snippets.cc CONSTRAINT
 */
void io::In_Parameter::read_CONSTRAINT(simulation::Parameter &param,
                                       std::ostream & os) {
    DEBUG(8, "reading CONSTRAINT");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "CONSTRAINT\n";
    exampleblock << "#	NTC\n";
    exampleblock << "#		1	solvent    solvent only\n";
    exampleblock << "#		2	hydrogen   solvent and solute bonds containing hydrogens and \n";
    exampleblock << "#		               constraints in the topology's CONSTRAINTS block\n";
    exampleblock << "#		3	all        solvent and solute all bonds\n";
    exampleblock << "#		4	specified  solvent and constraints in the topology's CONSTRAINTS block\n";
    exampleblock << "    3\n";
    exampleblock << "#       NTCP: solute algorithm\n";
    exampleblock << "#               1        shake\n";
    exampleblock << "#               2        lincs\n";
    exampleblock << "#               3        flexshake\n";
    exampleblock << "#       NTCP\n";
    exampleblock << "        1\n";
    exampleblock << "#       NTCP0(1)..(3): algorithm options\n";
    exampleblock << "#         - shake: tolerance\n";
    exampleblock << "#         - lincs: order\n";
    exampleblock << "#         - flexshake: tolerance, readin, order\n";
    exampleblock << "#       NTCP0(1)   NTCP0(2)   NTCP0(3)\n";
    exampleblock << "        0.0001\n";
    exampleblock << "#	NTCS: solvent algorithm\n";
    exampleblock << "#               1        shake\n";
    exampleblock << "#               2        lincs\n";
    exampleblock << "#               3        flexshake\n";
    exampleblock << "#               4        settle\n";
    exampleblock << "#               5        m_shake (only implemented for water and methanol!)\n";
    exampleblock << "#               6        gpu_shake\n";
    exampleblock << "#       NTCS\n";
    exampleblock << "        1\n";
    exampleblock << "#       NTCS0(1):  algorithm options\n";
    exampleblock << "#         - shake: tolerance\n";
    exampleblock << "#         - lincs: order\n";
    exampleblock << "#         - flexshake: tolerance, readin, order\n";
    exampleblock << "#         - settle: no arguments\n";
    exampleblock << "#         - m_shake: tolerance\n";
    exampleblock << "#       NTCS0(1)\n";
    exampleblock << "        0.0001\n";
    exampleblock << "#       NTCG: Number of GPUs\n";
    exampleblock << "#       NTCD: Device number of the GPU\n";
    exampleblock << "#       NTCG          NTCD\n";
    exampleblock << "        1             0\n";
    exampleblock << "END\n";


    std::string blockname = "CONSTRAINT";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        std::string sntc;
        block.get_next_parameter("NTC", sntc, "", "off, 0, solvent, 1, hydrogen, 2, all, 3, specified, 4");

        if (sntc == "off" || sntc == "0") param.constraint.ntc = 0;
        else if (sntc == "solvent" || sntc == "1") param.constraint.ntc = 1;
        else if (sntc == "hydrogen" || sntc == "2") param.constraint.ntc = 2;
        else if (sntc == "all" || sntc == "3") param.constraint.ntc = 3;
        else if (sntc == "specified" || sntc == "4") param.constraint.ntc = 4;


        if (block.error()) {
            param.constraint.solute.algorithm = simulation::constr_off;
            param.constraint.solvent.algorithm = simulation::constr_off;
            block.get_final_messages();
            return;
        }

        // SOLUTE
        std::string salg;
        block.get_next_parameter("NTCP", salg, "", "shake, 1, lincs, 2, flexshake, 3");

        DEBUG(8, "Constraints (solute): " << salg);

        if (salg == "shake" || salg == "1") {
            DEBUG(9, "constraints solute shake");

            if (param.constraint.ntc > 1)
                param.constraint.solute.algorithm = simulation::constr_shake;
            else param.constraint.solute.algorithm = simulation::constr_off;

            block.get_next_parameter("NTCP[0]", param.constraint.solute.shake_tolerance, ">0", "");

        } else if (salg == "flexshake" || salg == "3") {
            DEBUG(9, "constraints solute flexshake");

            if (param.constraint.ntc > 1)
                param.constraint.solute.algorithm = simulation::constr_flexshake;
            else param.constraint.solute.algorithm = simulation::constr_off;

            block.get_next_parameter("NTCP[0]", param.constraint.solute.shake_tolerance, ">0", "");
            block.get_next_parameter("NTCP[1]", param.constraint.solute.flexshake_readin, "", "");
            block.get_next_parameter("NTCP[2]", param.constraint.solute.flexshake_mode, ">=0 && <=3", "");

        } else if (salg == "lincs" || salg == "2") {
            DEBUG(9, "constraints solute lincs");

            if (param.constraint.ntc > 1)
                param.constraint.solute.algorithm = simulation::constr_lincs;
            else
                param.constraint.solute.algorithm = simulation::constr_off;

            block.get_next_parameter("NTCP[0]", param.constraint.solute.lincs_order, ">=1", "");

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
        block.get_next_parameter("NTCS", salg, "", "shake, 1, lincs, 2, flexshake, 3, settle, 4, m_shake, 5, gpu_shake, 6");

        DEBUG(8, "constraints solvent: " << salg);

        if (salg == "shake" || salg == "1") {
            DEBUG(9, "constraints solvent shake");

            param.constraint.solvent.algorithm = simulation::constr_shake;
            block.get_next_parameter("NTCS[0]", param.constraint.solvent.shake_tolerance, ">0", "");

        } else if (salg == "flexshake" || salg == "3") {
            DEBUG(9, "constraints solvent flexshake");

            param.constraint.solvent.algorithm = simulation::constr_flexshake;
            block.get_next_parameter("NTCS[0]", param.constraint.solvent.shake_tolerance, ">0", "");

        } else if (salg == "lincs" || salg == "2") {
            DEBUG(9, "constraints solvent lincs");

            param.constraint.solvent.algorithm = simulation::constr_lincs;
            block.get_next_parameter("NTCS[0]", param.constraint.solvent.lincs_order, ">=1", "");

        } else if (salg == "settle" || salg == "4") {
            DEBUG(9, "constraints solvent settle");

            param.constraint.solvent.algorithm = simulation::constr_settle;
        } else if (salg == "m_shake" || salg == "5") {
            DEBUG(9, "constraints solvent m_shake");

            param.constraint.solvent.algorithm = simulation::constr_m_shake;
            block.get_next_parameter("NTCS[0]", param.constraint.solvent.shake_tolerance, ">0", "");

        } else if (salg == "gpu_shake" || salg == "6") {
            DEBUG(9, "constraints solvent gpu_shake");

            param.constraint.solvent.algorithm = simulation::constr_gpu_shake;
            block.get_next_parameter("NTCS[0]", param.constraint.solvent.shake_tolerance, ">0", "");

            // Get number of GPUs and their IDs
            block.get_next_parameter("NTCG", param.constraint.solvent.number_gpus, ">0", "");
            int temp = 0;
            bool fail = false;
            for (unsigned int i = 0; i < param.constraint.solvent.number_gpus; i++) {
                std::string idx=io::to_string(i);
                if (block.get_next_parameter("NTCD["+idx+"]", temp, ">=-1", "", true)) {
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
                block.reset_error();
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

        block.get_final_messages();
    } else {
        param.constraint.ntc = 1;
        param.constraint.solute.algorithm = simulation::constr_off;
        param.constraint.solvent.algorithm = simulation::constr_shake;

        io::messages.add("no CONSTRAINT block", "In_Parameter",
                         io::message::error);
        return;
    }
} //CONSTRAINT

/**
 * @section printout PRINTOUT block
 * @snippet snippets/snippets.cc PRINTOUT
 */
void io::In_Parameter::read_PRINTOUT(simulation::Parameter &param,
                                     std::ostream & os) {
    DEBUG(8, "reading PRINTOUT");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "PRINTOUT\n";
    exampleblock << "#  NTPR: print out energies, etc. every NTPR steps\n";
    exampleblock << "#  NTPP: =1 perform dihedral angle transition monitoring\n";
    exampleblock << "#     NTPR      NTPP\n";
    exampleblock << "         0         0\n";
    exampleblock << "END\n";


    std::string blockname = "PRINTOUT";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("NTPR", param.print.stepblock, ">=0", "");
        block.get_next_parameter("NTPP", param.print.monitor_dihedrals, "", "0,1");

        block.get_final_messages();
    }
}

/**
 * @section writetraj WRITETRAJ block
 * @snippet snippets/snippets.cc WRITETRAJ
 */
void io::In_Parameter::read_WRITETRAJ(simulation::Parameter &param,
                                      std::ostream & os) {
    DEBUG(8, "reading WRITETRAJ");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "WRITETRAJ\n";
    exampleblock << "# NTWX       controls writing of coordinate trajectory\n";
    exampleblock << "#       0: no coordinate trajectory is written (default)\n";
    exampleblock << "#      >0: write solute and solvent coordinates every NTWX steps\n";
    exampleblock << "#      <0: write solute coordinates every |NTWX| steps\n";
    exampleblock << "# NTWSE >= 0 selection criteria for coordinate trajectory writing\n";
    exampleblock << "#       0: write normal coordinate trajectory\n";
    exampleblock << "#      >0: write minimum-energy coordinate and energy trajectory (based on the\n";
    exampleblock << "#          energy entry selected by NTWSE and as blocks of length NTWX)\n";
    exampleblock << "#          (see configuration/energy.cc or ene_ana library for indices)\n";
    exampleblock << "# NTWV       controls writing of velocity trajectory\n";
    exampleblock << "#       0: no velocity trajectory is written (default)\n";
    exampleblock << "#      >0: write solute and solvent velocities every NTWV steps\n";
    exampleblock << "#      <0: write solute velocities every |NTWV| steps\n";
    exampleblock << "# NTWF       controls writing of force trajectory\n";
    exampleblock << "#       0: no force trajectory is written (default)\n";
    exampleblock << "#      >0: write solute and solvent forces every NTWF steps\n";
    exampleblock << "#      <0: write solute forces every |NTWF| steps\n";
    exampleblock << "# NTWE >= 0 controls writing of energy trajectory\n";
    exampleblock << "#       0: no energy trajectory is written (default)\n";
    exampleblock << "#      >0: write energy trajectory every NTWE steps\n";
    exampleblock << "# NTWG >= 0 controls writing of free energy trajectory\n";
    exampleblock << "#       0: no free energy trajectory is written (default)\n";
    exampleblock << "#      >0: write free energy trajectory every NTWG steps\n";
    exampleblock << "# NTWB >= 0 controls writing of block-averaged energy trajectory\n";
    exampleblock << "#       0: no block averaged energy trajectory is written (default)\n";
    exampleblock << "#      >0: write block-averaged energy variables every |NTWB| steps\n";
    exampleblock << "#          (and free energies if NTWG > 0) trajectory\n";
    exampleblock << "#\n";
    exampleblock << "#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB\n";
    exampleblock << "      100          0         0         0       100         0       100\n";
    exampleblock << "END\n";


    std::string blockname = "WRITETRAJ";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("NTWX", param.write.position, "", "");
        std::string str_max_ene = io::to_string(configuration::Energy::MAX_ENERGY_INDEX);
        block.get_next_parameter("NTWSE", param.write.energy_index, ">=0 && <="+str_max_ene, "");
        block.get_next_parameter("NTWV", param.write.velocity, "", "");
        block.get_next_parameter("NTWF", param.write.force, "", "");
        block.get_next_parameter("NTWE", param.write.energy, ">=0", "");
        block.get_next_parameter("NTWG", param.write.free_energy, ">=0", "");
        block.get_next_parameter("NTWB", param.write.block_average, ">=0", "");

        if (param.write.energy_index > 0) {
            if (param.write.position == 0) {
                io::messages.add("WRITETRAJ block: NTWX must be a block size > 0 for "
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

        block.get_final_messages();
    }
} // WRITETRAJ

/**
 * @section pressurescale PRESSURESCALE block
 * @snippet snippets/snippets.cc PRESSURESCALE
 */
void io::In_Parameter::read_PRESSURESCALE(simulation::Parameter &param,
                                          std::ostream & os) {
    DEBUG(8, "reading PRESSURESCALE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "PRESSURESCALE\n";
    exampleblock << "#	COUPLE:	off(0), calc(1), scale(2)\n";
    exampleblock << "#	SCALE:  off(0), iso(1), aniso(2), full(3), semianiso(4)\n";
    exampleblock << "#	VIRIAL: none(0), atomic(1), group(2)\n";
    exampleblock << "#\n";
    exampleblock << "#   COUPLE  SCALE   COMP        TAUP    VIRIAL\n";
    exampleblock << "    calc    iso     4.575E-4    0.5     atomic\n";
    exampleblock << "#   SEMI (semianisotropic couplings: X, Y, Z)\n";
    exampleblock << "#       e.g. 1 1 2: x and y jointly coupled and z separately coupled\n";
    exampleblock << "#       e.g. 0 0 1: constant area (xy-plane) and z coupled to a bath\n";
    exampleblock << "    1 1 2\n";
    exampleblock << "#   reference pressure\n";
    exampleblock << "    0.06102     0.00000     0.00000\n";
    exampleblock << "    0.00000     0.06102     0.00000\n";
    exampleblock << "    0.00000     0.00000     0.06102\n";
    exampleblock << "END\n";


    std::string blockname = "PRESSURESCALE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        std::string couple, scale, vir;
        int xs = 0, ys = 0, zs = 0;
        block.get_next_parameter("COUPLE", couple, "", "off, 0, calc, 1, scale, 2");
        block.get_next_parameter("SCALE", scale, "", "off, 0, iso, 1, aniso, 2, full, 3, semianiso, 4");
        block.get_next_parameter("COMP", param.pcouple.compressibility, ">0", "");
        block.get_next_parameter("TAUP", param.pcouple.tau, ">0", "");
        block.get_next_parameter("VIRIAL", vir, "", "none, 0, atomic, 1, group, molecular, 2");
        block.get_next_parameter("SEMIX", xs, ">=0 && <=2", "");
        block.get_next_parameter("SEMIY", ys, ">=0 && <=2", "");
        block.get_next_parameter("SEMIZ", zs, ">=0 && <=2", "");

        for (int i = 0; i < 3; ++i) {
            std::string idx = io::to_string(i);
            for (int j = 0; j < 3; ++j) {
                std::string jdx = io::to_string(j);
                block.get_next_parameter("PRES["+idx+"]["+jdx+"]", param.pcouple.pres0(i, j), ">=0", "");
            }
        }

        if (couple == "off" || couple == "0") {
            param.pcouple.calculate = false;
            param.pcouple.scale = math::pcouple_off;
        } else if (couple == "calc" || couple == "1") {
            param.pcouple.calculate = true;
            param.pcouple.scale = math::pcouple_off;
        } else if (couple == "scale" || couple == "2") {
            param.pcouple.calculate = true;

            if (scale == "off" || scale == "0") {
                io::messages.add("PRESSURESCALE block: requesting scaling but SCALE set to OFF",
                                 "In_Parameter", io::message::error);
                param.pcouple.scale = math::pcouple_off;
            } else if (scale == "iso" || scale == "1")
                param.pcouple.scale = math::pcouple_isotropic;
            else if (scale == "aniso" || scale == "2")
                param.pcouple.scale = math::pcouple_anisotropic;
            else if (scale == "full" || scale == "3")
                param.pcouple.scale = math::pcouple_full_anisotropic;
            else if (scale == "semianiso" || scale == "4") {
                param.pcouple.scale = math::pcouple_semi_anisotropic;
                if (xs == ys && ys == zs)
                    io::messages.add("PRESSURESCALE block: x_semi = y_semi = z_semi in SEMI, maybe you want isotropic pressure scaling?",
                                     "In_Parameter", io::message::error);
                param.pcouple.x_semi = xs;
                param.pcouple.y_semi = ys;
                param.pcouple.z_semi = zs;
            } else {
                param.pcouple.scale = math::pcouple_off;
            }

        } else {
            param.pcouple.calculate = false;
        }

        if (param.pcouple.calculate) {
            if (vir == "none" || vir == "0") {
                param.pcouple.virial = math::no_virial;
            } else if (vir == "atomic" || vir == "1")
                param.pcouple.virial = math::atomic_virial;
            else if (vir == "molecular" || vir == "group" || vir == "2")
                param.pcouple.virial = math::molecular_virial;
            else {
                param.pcouple.virial = math::no_virial;
            }
        } else
            param.pcouple.virial = math::no_virial;

        if (param.pcouple.calculate == false && param.pcouple.scale != math::pcouple_off)
            io::messages.add("PRESSURESCALE block: pressure coupling activated but "
                             "not calculating pressure",
                             "In_Parameter",
                             io::message::error);
        if (param.pcouple.calculate && param.pcouple.virial == math::no_virial)
            io::messages.add("requesting pressure calculation but "
                             "no virial specified",
                             "In_Parameter", io::message::error);

        block.get_final_messages();
    }
} // PRESSURESCALE block

/**
 * @section boundcond BOUNDCOND block
 * @snippet snippets/snippets.cc BOUNDCOND
 */
void io::In_Parameter::read_BOUNDCOND(simulation::Parameter &param,
                                      std::ostream & os) {
    DEBUG(8, "reading BOUNDCOND");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "BOUNDCOND\n";
    exampleblock << "#  NTB: boundary conditions\n";
    exampleblock << "#       -1 : truncated octahedral\n";
    exampleblock << "#        0 : vacuum\n";
    exampleblock << "#        1 : rectangular\n";
    exampleblock << "#        2 : triclinic\n";
    exampleblock << "#  NDFMIN: number of degrees of freedom subtracted for temperature\n";
    exampleblock << "#\n";
    exampleblock << "#         NTB    NDFMIN\n";
    exampleblock << "            1         0\n";
    exampleblock << "END\n";


    std::string blockname = "BOUNDCOND";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert(blockname);

        int ntb = 0;
        block.get_next_parameter("NTB", ntb, "", "-1, 0, 1, 2");
        block.get_next_parameter("NDFMIN", param.boundary.dof_to_subtract, ">=0", "");

        if (ntb == 0) {
            param.boundary.boundary = math::vacuum;
        }
        else if (ntb == 1) param.boundary.boundary = math::rectangular;
        else if (ntb == 2) param.boundary.boundary = math::triclinic;
        else if (ntb == -1) {
            param.boundary.boundary = math::truncoct;
            io::messages.add("Truncated octahedral: the box is converted to triclinic boundary conditions internally",
                             "In_Parameter", io::message::notice);
        } else {
            std::ostringstream msg;
            param.boundary.boundary = math::vacuum;
        }

        block.get_final_messages();

    } else {
        param.boundary.boundary = math::vacuum;
        return;
    }

}

/**
 * @section perturbation PERTURBATION block
 * @snippet snippets/snippets.cc PERTURBATION
 */
void io::In_Parameter::read_PERTURBATION(simulation::Parameter &param,
                                         std::ostream & os) {
    DEBUG(8, "reading PERTURBATION");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "PERTURBATION\n";
    exampleblock << "#    NTG: 0..1 controls use of free-energy calculation.\n";
    exampleblock << "#         0: no free-energy calculation (default)\n";
    exampleblock << "#         1: calculate dH/dRLAM\n";
    exampleblock << "#  NRDGL: 0,1 controls reading of initial value for RLAM.\n";
    exampleblock << "#         0: use initial RLAM parameter from PERTURBATION block\n";
    exampleblock << "#         1: read from configuration\n";
    exampleblock << "#   RLAM: 0.0..1.0 initial value for lambda\n";
    exampleblock << "#  DLAMT: >= 0.0 rate of lambda increase in time.\n";
    exampleblock << "# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter\n";
    exampleblock << "#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter\n";
    exampleblock << "#   NLAM: > 0 power dependence of lambda coupling\n";
    exampleblock << "# NSCALE: 0..2 controls use of interaction scaling\n";
    exampleblock << "#         0: no interaction scaling\n";
    exampleblock << "#         1: interaction scaling\n";
    exampleblock << "#         2: perturbation for all atom pairs with scaled\n";
    exampleblock << "#            interactions. No perturbation for others.\n";
    exampleblock << "#\n";
    exampleblock << "#     NTG   NRDGL    RLAM   DLAMT\n";
    exampleblock << "        0       0     0.0     0.0\n";
    exampleblock << "#  ALPHLJ   ALPHC    NLAM  NSCALE\n";
    exampleblock << "      0.0     0.0       1       0\n";
    exampleblock << "END\n";


    std::string blockname = "PERTURBATION";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        std::string b, s1, s2;
        int ntg = 0, scale = 0, nrdgl = 0;
        block.get_next_parameter("NTG", ntg, "", "0,1");
        block.get_next_parameter("NRDGL", nrdgl, "", "0,1");
        block.get_next_parameter("RLAM", param.perturbation.lambda, ">=0 && <=1", "");
        block.get_next_parameter("DLAMT", param.perturbation.dlamt, ">=0", "");
        block.get_next_parameter("ALPHLJ", param.perturbation.soft_vdw, ">=0", "");
        block.get_next_parameter("ALPHC", param.perturbation.soft_crf, ">=0", "");
        block.get_next_parameter("NLAM", param.perturbation.lambda_exponent, ">0", "");
        block.get_next_parameter("NSCALE", scale, "", "0,1,2");

        switch (ntg) {
            case 0:
                param.perturbation.perturbation = false;
                break;
            case 1:
                param.perturbation.perturbation = true;
                break;
            default:
                break;
        }

        switch (nrdgl) {
            case 0:     // use from input file
                param.perturbation.read_initial = false;
                break;
            case 1:     // use from configuration
                param.perturbation.read_initial = true;
                break;
            default:
                break;
        }

        switch (scale) {
            case 0:     // no scaling
                param.perturbation.scaling = false;
                param.perturbation.scaled_only = false;
                break;
            case 1:     // scaling on
                param.perturbation.scaling = true;
                param.perturbation.scaled_only = false;
                break;
            case 2:     // scaled only
                param.perturbation.scaling = true;
                param.perturbation.scaled_only = true;
                break;
            default:
                break;
        }
        block.get_final_messages();
    }

} //PERTURBATION

/**
 * @section force FORCE block
 * @snippet snippets/snippets.cc FORCE
 */
void io::In_Parameter::read_FORCE(simulation::Parameter &param,
                                  std::ostream & os) {
    DEBUG(8, "reading FORCE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "FORCE\n";
    exampleblock << "# NTF(1..6): 0,1 determines terms used in force calculation\n";
    exampleblock << "#             0: do not include terms\n";
    exampleblock << "#             1: include terms\n";
    exampleblock << "# NEGR: ABS(NEGR): number of energy groups\n";
    exampleblock << "#             > 0: use energy groups\n";
    exampleblock << "#             < 0: use energy and force groups\n";
    exampleblock << "# NRE(1..NEGR): >= 1 last atom in each energy group\n";
    exampleblock << "# NTF(1)    NTF(2)    NTF(3)    NTF(4)    NTF(5)        NTF(6)\n";
    exampleblock << "# bonds     angles    improper  dihedral  electrostatic vdW\n";
    exampleblock << "  0         1         1         1         1             1\n";
    exampleblock << "# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)\n";
    exampleblock << "     1        60\n";
    exampleblock << "END\n";


    std::string blockname = "FORCE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("NTF(1)", param.force.bond, "", "0,1");
        block.get_next_parameter("NTF(2)", param.force.angle, "", "0,1");
        block.get_next_parameter("NTF(3)", param.force.improper, "", "0,1");
        block.get_next_parameter("NTF(4)", param.force.dihedral, "", "0,1");
        block.get_next_parameter("NTF(5)", param.force.nonbonded_crf, "", "0,1");
        block.get_next_parameter("NTF(6)", param.force.nonbonded_vdw, "", "0,1");

        int snum = 0;
        unsigned int num = 0, e = 0, old_e = 0;
        block.get_next_parameter("NEGR", snum, "", "");

        if (snum < 0) {
#ifdef XXFORCEGROUPS
            param.force.force_groups = true;
#else
            io::messages.add("Force groups requested but not compiled with support for them."
                             "Use --enable-forcegroups for configure.", "In_Parameter",
                             io::message::error);
#endif
        }
        num = abs(snum);

        for (unsigned int i = 0; i < num; ++i) {
            std::string idx=io::to_string(i);
            block.get_next_parameter("NRE["+idx+"]", e, ">0", "");
            DEBUG(10, "\tadding energy group " << e - 1);
            param.force.energy_group.push_back(e - 1);
            if (e <= old_e) {
                DEBUG(10, "energy groups not in order...");
                io::messages.add("FORCE block: energy groups are not in order",
                                 "In_Parameter", io::message::error);
            }
            old_e = e;
        }

        DEBUG(10, "number of energy groups: " << param.force.energy_group.size());

        if ((!param.force.nonbonded_crf) && param.force.nonbonded_vdw)
            io::messages.add("FORCE block: setting charges to zero",
                             "In_Parameter", io::message::notice);

        if (param.force.nonbonded_crf && (!param.force.nonbonded_vdw))
            io::messages.add("FORCE block: setting atom types to dummy",
                             "In_Parameter", io::message::notice);

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

        block.get_final_messages();
    }
} //FORCE

/**
 * @section covalentform COVALENTFORM block
 * @snippet snippets/snippets.cc COVALENTFORM
 */
void io::In_Parameter::read_COVALENTFORM(simulation::Parameter &param,
                                         std::ostream & os) {
    DEBUG(8, "reading COVALENTFORM");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "COVALENTFORM\n";
    exampleblock << "# NTBBH: 0,1 controls bond-stretching potential\n";
    exampleblock << "#        0: quartic form (default)\n";
    exampleblock << "#        1: harmonic form\n";
    exampleblock << "# NTBAH: 0,1 controls bond-angle bending potential\n";
    exampleblock << "#        0: cosine-harmonic (default)\n";
    exampleblock << "#        1: harmonic\n";
    exampleblock << "# NTBDN: 0,1 controls torsional dihedral potential\n";
    exampleblock << "#        0: arbitrary phase shifts (default)\n";
    exampleblock << "#        1: phase shifts limited to 0 and 180 degrees.\n";
    exampleblock << "#   NTBBH    NTBAH    NTBDN\n";
    exampleblock << "        0        0        0\n";
    exampleblock << "END\n";


    std::string blockname = "COVALENTFORM";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int bond = 0, angle = 0, dihedral = 0;
        block.get_next_parameter("NTBBH", bond, "", "0,1");
        block.get_next_parameter("NTBAH", angle, "", "0,1");
        block.get_next_parameter("NTBDN", dihedral, "", "0,1");

        if (param.force.bond != 0) {
            switch (bond) {
                case 1: param.force.bond = 2;
                    break;
                case 0:
                default: param.force.bond = 1;
            }
        }

        if (param.force.angle != 0) {
            switch (angle) {
                case 1: param.force.angle = 2;
                    break;
                case 0:
                default: param.force.angle = 1;
            }
        }

        if (param.force.dihedral != 0) {
            switch (dihedral) {
                case 1: param.force.dihedral = 2;
                    break;
                case 0:
                default: param.force.dihedral = 1;
            }
        }

        block.get_final_messages();
    }
}

/**
 * @section initialise INITIALISE block
 * @snippet snippets/snippets.cc INITIALISE
 */
void io::In_Parameter::read_INITIALISE(simulation::Parameter &param,
                                       std::ostream & os) {
    DEBUG(8, "reading INITIALISE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "INITIALISE\n";
    exampleblock << "# NTIVEL: 0,1 controls generation of initial velocities.\n";
    exampleblock << "#         0: read from configuration (default)\n";
    exampleblock << "#         1: generate from Maxell distribution at temperature TEMPI\n";
    exampleblock << "# NTISHK: 0..3 controls shaking of initial configuration\n";
    exampleblock << "#         0: no intial SHAKE (default)\n";
    exampleblock << "#         1: initial SHAKE on coordinates only\n";
    exampleblock << "#         2: initial SHAKE on velocities only\n";
    exampleblock << "#         3: initial SHAKE on coordinates and velocities\n";
    exampleblock << "# NTINHT: 0,1 controls generation of initial Nose-Hoover chain variables\n";
    exampleblock << "#         0: read from configuration (default)\n";
    exampleblock << "#         1: reset variables to zero.\n";
    exampleblock << "# NTINHB: 0,1 controls generation of initial Nose-Hoover (chain) barostat\n";
    exampleblock << "#             variables\n";
    exampleblock << "#         0: read from strartup file (if applicable) (default)\n";
    exampleblock << "#         1: reset variables to zero\n";
    exampleblock << "# NTISHI: 0,1 controls initial setting for lattice shift vectors\n";
    exampleblock << "#         0: read from configuration (default)\n";
    exampleblock << "#         1: reset shifts to zero.\n";
    exampleblock << "# NTIRTC: 0,1 controls initial setting of positions and orientations for\n";
    exampleblock << "#             roto-translational constraints\n";
    exampleblock << "#         0: read from configuration (default)\n";
    exampleblock << "#         1: reset based on initial configuraion of startup file\n";
    exampleblock << "# NTICOM: 0,1 controls initial removal of COM motion\n";
    exampleblock << "#         0: no initial system COM motion removal (default)\n";
    exampleblock << "#         1: initial COM translation is removed\n";
    exampleblock << "#         2: initial COM rotation is removed\n";
    exampleblock << "# NTISTI: 0,1 controls generation of stochastic integrals\n";
    exampleblock << "#         0: read stochastic integrals and IG from configuration (default)\n";
    exampleblock << "#         1: set stochastic integrals to zero and use IG from here.\n";
    exampleblock << "# IG:     random number generator seed\n";
    exampleblock << "# TEMPI:  initial temperature\n";
    exampleblock << "#\n";
    exampleblock << "#  NTIVEL  NTISHK  NTINHT  NTINHB\n";
    exampleblock << "        0       0       0       0\n";
    exampleblock << "#  NTISHI  NTIRTC  NTICOM\n";
    exampleblock << "        0       0       0\n";
    exampleblock << "#  NTISTI\n";
    exampleblock << "        0\n";
    exampleblock << "#      IG   TEMPI\n";
    exampleblock << "        0     0.0\n";
    exampleblock << "END\n";


    std::string blockname = "INITIALISE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert(blockname);

        int ntivel = 0, ntishk = 0, ntinht = 0, ntinhb = 0, ntishi = 0, ntirtc = 0, nticom = 0, ntisti = 0;
        block.get_next_parameter("NTIVEL", ntivel, "", "0,1");
        block.get_next_parameter("NTISHK", ntishk, "", "0,1,2,3");
        block.get_next_parameter("NTINHT", ntinht, "", "0,1");
        block.get_next_parameter("NTINHB", ntinhb, "", "0,1");
        block.get_next_parameter("NTISHI", ntishi, "", "0,1");
        block.get_next_parameter("NTIRTC", ntirtc, "", "0,1");
        block.get_next_parameter("NTICOM", nticom, "", "0,1,2");
        block.get_next_parameter("NTISTI", ntisti, "", "0,1");
        block.get_next_parameter("IG", param.start.ig, "", "");
        block.get_next_parameter("TEMPI", param.start.tempi, ">=0", "");

        // generation of initial velocities
        switch (ntivel) {
            case 0: param.start.generate_velocities = false;
                break;
            case 1: param.start.generate_velocities = true;
                break;
            default:
                break;
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
            default:
              break;
        }

        // controls reading of Nose-Hoover chain variables.
        switch (ntinht) {
            case 0: param.start.read_nosehoover_chains = true;
              break;
            case 1: param.start.read_nosehoover_chains = false;
              break;     // reset them
            default:
              break;
        }

        // controls reading of Nose-Hoover chain barostat variables: not implemented.
        switch (ntinhb) {
            case 0: param.start.read_nosehoover_barostat = true;
              break;
            case 1: param.start.read_nosehoover_barostat = false;
              break;     // reset them
            default:
              break;
        }

        switch (ntishi) {
            case 0: param.start.read_lattice_shifts = true;
              break;
            case 1: param.start.read_lattice_shifts = false;
              break;
            default:
              break;
        }

        // controls reading of restart data for roto-translational constraints
        switch (ntirtc) {
            case 0: param.start.read_rottrans = true;
              break;
            case 1: param.start.read_rottrans = false;
              break;
            default:
              break;
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
            default:
              break;
        }

        // controls reading of stochastic integrals
        switch (ntisti) {
            case 0:
              param.stochastic.generate_integral = false;
              break;
          case 1:
              param.stochastic.generate_integral = true;
              break;
          default:
              break;
        }

    block.get_final_messages();
  }
}

/**
 * @section comtransrot COMTRANSROT block
 * @snippet snippets/snippets.cc COMTRANSROT
 */
void io::In_Parameter::read_COMTRANSROT(simulation::Parameter &param,
                                        std::ostream & os) {
    DEBUG(8, "reading COMTRANSROT");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "COMTRANSROT\n";
    exampleblock << "#    NSCM : controls system centre-of-mass (com) motion removal\n";
    exampleblock << "#           0: no com motion removal (default)\n";
    exampleblock << "#         < 0: com translation and rotation are removed every abs(NSCM)\n";
    exampleblock << "#              steps.\n";
    exampleblock << "#         > 0: com translation is removed every NSCM steps.\n";
    exampleblock << "#     NSCM\n";
    exampleblock << "         0\n";
    exampleblock << "END\n";


    std::string blockname = "COMTRANSROT";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int nscm = 0;
        block.get_next_parameter("NSCM", nscm, "", "");

        if (nscm > 0) {
            param.centreofmass.skip_step = nscm;
            param.centreofmass.remove_rot = false;
            param.centreofmass.remove_trans = true;
        } else if (nscm < 0) {
            param.centreofmass.skip_step = -nscm;
            param.centreofmass.remove_rot = true;
            param.centreofmass.remove_trans = true;
        } else {         // nscm == 0;
            param.centreofmass.skip_step = 0;
            param.centreofmass.remove_rot = false;
            param.centreofmass.remove_trans = false;
        }

        block.get_final_messages();
    }
}

/**
 * @section hoomd HOOMD block (optional)
 * @snippet snippets/snippets.cc HOOMD
 */
void io::In_Parameter::read_HOOMD(simulation::Parameter &param,
                                  std::ostream & os) {
    DEBUG(8, "reading HOOMD");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "HOOMD\n";
    exampleblock << "#       PROCESSOR: cpu gpus\n";
    exampleblock << "#\n";
    exampleblock << "#       PROCESSOR\n";
    exampleblock << "    gpus\n";
    exampleblock << "END\n";

    std::string blockname = "HOOMD";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);
        param.setDevelop("HOOMD is under development.");

        DEBUG(10, "hoomd block");

        // try a HOOMD

        std::string s1;
        block.get_next_parameter("PROCESSOR", s1, "", "cpu, gpus");

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
        block.get_final_messages();
#endif
    }
}

/**
 * @section pairlist PAIRLIST block
 * @snippet snippets/snippets.cc PAIRLIST
 */
void io::In_Parameter::read_PAIRLIST(simulation::Parameter &param,
                                     std::ostream & os) {
    DEBUG(8, "reading PAIRLIST");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "PAIRLIST\n";
    exampleblock << "#       ALGORITHM  standard(0) (gromos96 like pairlist)\n";
    exampleblock << "#                  grid(1) (md++ grid pairlist)\n";
    exampleblock << "#                  grid_cell(2) (creates a mask)\n";
    exampleblock << "#       NSNB  >0    frequency (number of steps) a pairlist is constructed\n";
    exampleblock << "#       RCUTP >0.0  short-range cut-off in twin-range\n";
    exampleblock << "#       RCUTL >0.0  intermediate-range cut-off in twin-range\n";
    exampleblock << "#       SIZE  >0    grid cell size (or auto = 0.5 * RCUTP)\n";
    exampleblock << "#       TYPE  chargegoup(0) (chargegroup based cutoff)\n";
    exampleblock << "#             atomic(1)     (atom based cutoff)\n";
    exampleblock << "#\n";
    exampleblock << "#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE\n";
    exampleblock << "        grid            5       0.8     1.4     auto    chargegroup\n";
    exampleblock << "#\n";
    exampleblock << "END\n";


    std::string blockname = "PAIRLIST";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert(blockname);

        std::string s1, s2, s3;
        block.get_next_parameter("ALGORITHM", s1, "", "standard, 0, grid, 1, grid_cell, 2");
        block.get_next_parameter("NSNB", param.pairlist.skip_step, ">0", "");
        block.get_next_parameter("RCUTP", param.pairlist.cutoff_short, ">=0", "");
        std::string str_rcutp=io::to_string(param.pairlist.cutoff_short);
        block.get_next_parameter("RCUTL", param.pairlist.cutoff_long, ">="+str_rcutp, "");
        block.get_next_parameter("SIZE", s2, "", "");
        block.get_next_parameter("TYPE", s3, "", "chargegroup, 0, atomic, 1");

        if (s1 == "standard" || s1 == "0") param.pairlist.grid = 0;
        else if (s1 == "grid" || s1 == "1") param.pairlist.grid = 1;
        else if (s1 == "grid_cell" || s1 == "2") param.pairlist.grid = 2;
        else param.pairlist.grid = 0;

        if (param.pairlist.grid) {
            if (s2 == "auto")
                param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
            else {
                std::istringstream css;
                css.str(s2);
                css >> param.pairlist.grid_size;
                if (!param.pairlist.grid_size)
                    param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
                if (css.fail()) {
                    io::messages.add("PAIRLIST block: wrong pairlist grid size chosen",
                                     "In_Parameter", io::message::error);
                    param.pairlist.grid_size = 0.5 * param.pairlist.cutoff_short;
                }
            }
            if (param.pairlist.grid_size <= 0)
                io::messages.add("PAIRLIST block: Illegal value for grid size (>0)",
                                 "In_Parameter", io::message::error);
        } else param.pairlist.grid_size = 0;

        if (s3 == "atomic" || s3 == "1") param.pairlist.atomic_cutoff = true;
        else if (s3 == "chargegroup" || s3 == "0") param.pairlist.atomic_cutoff = false;
        else param.pairlist.atomic_cutoff = false;
        block.get_final_messages();
    }
} // PAIRLIST

/**
 * @section cgrain CGRAIN block
 * @snippet snippets/snippets.cc CGRAIN
 */
void io::In_Parameter::read_CGRAIN(simulation::Parameter &param,
                                   std::ostream & os) {
    DEBUG(8, "reading CGRAIN");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "CGRAIN\n";
    exampleblock << "# NTCGRAN 0..3 coarse grain selection\n";
    exampleblock << "#         0: atomistic (off)\n";
    exampleblock << "#         1: coarse-grained using MARTINI model (on)\n";
    exampleblock << "#         2: coarse-grained using GROMOS model (on)\n";
    exampleblock << "#         3: mixed-grained using GROMOS model (on)\n";
    exampleblock << "#     EPS >= 0.0 dielectric constant for coarse grained coulombic interaction\n";
    exampleblock << "#    EPSM >= 0.0 dielectric constant for mixed CG-FG coulombic interaction\n";
    exampleblock << "# NTCGRAN     EPS     EPSM\n";
    exampleblock << "        1      20        1\n";
    exampleblock << "END\n";


    std::string blockname = "CGRAIN";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int ntcgran = 0;
        block.get_next_parameter("NTCGRAN", ntcgran, "", "0,1,2,3");
        block.get_next_parameter("EPS", param.cgrain.EPS, ">=0", "");
        block.get_next_parameter("EPSM", param.cgrain.EPSM, ">=0", "");

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

        block.get_final_messages();
    }
}

/**
 * @section multibath MULTIBATH block
 * @snippet snippets/snippets.cc MULTIBATH
 */
void io::In_Parameter::read_MULTIBATH(simulation::Parameter &param,
                                      std::ostream & os) {

    DEBUG(8, "reading MULTIBATH");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "MULTIBATH\n";
    exampleblock << "# NTBTYP: temperature coupling algorithm\n";
    exampleblock << "#		weak-coupling(0)\n";
    exampleblock << "#		nose-hoover(1)\n";
    exampleblock << "#		nose-hoover-chains(2)	num\n";
    exampleblock << "#		(where num is the number of chains to use)\n";
    exampleblock << "#   NTBTYP           NUM\n";
    exampleblock << "    nose-hoover-chains	3\n";
    exampleblock << "#   NBATHS: number of baths\n";
    exampleblock << "    2\n";
    exampleblock << "#   TEMP0  TAU\n";
    exampleblock << "    300    0.10\n";
    exampleblock << "    300    0.10\n";
    exampleblock << "#   DOFSET: number of different couplings\n";
    exampleblock << "    1\n";
    exampleblock << "#   LAST   COM-BATH  IR-BATH\n";
    exampleblock << "    60     1         2\n";
    exampleblock << "#   (this couples the first 60 atoms com motion to bath 1 and\n";
    exampleblock << "#    the internal / rotational motion to bath 2)\n";
    exampleblock << "END\n";

    std::string blockname = "MULTIBATH";
    Block block(blockname, exampleblock.str());

    param.multibath.couple = false;

    if (!block.read_buffer(m_block[blockname], false)) {
        block_read.insert(blockname);

        param.multibath.found_multibath = true;
        param.multibath.found_tcouple = false;

        // the algorithm
        std::string alg;
        block.get_next_parameter("NTBTYP", alg, "", "weak-coupling, 0, nose-hoover, 1, nose-hoover-chains, 2");

        if (alg == "weak-coupling" || alg == "0")
            param.multibath.algorithm = 0;
        else if (alg == "nose-hoover" || alg == "1")
            param.multibath.algorithm = 1;
        else if (alg == "nose-hoover-chains" || alg == "2") {
            param.multibath.algorithm = 2;
            // read in the number of chains => overwrites algorithm number!!
            block.get_next_parameter("NUM", param.multibath.algorithm, ">=2", "");
        }

        int num_baths = 0, num_dof = 0;
        unsigned int last = 0, com_bath = 0, ir_bath = 0;
        double temp = 0.0, tau = 0.0;

        // the baths
        block.get_next_parameter("NBATHS", num_baths, ">=0", "");

        if (block.error()) {
          block.get_final_messages();
          return;
        }

        for (int i = 0; i < num_baths; ++i) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("TEMP["+idx+"]", temp, ">=0.0", "");
            block.get_next_parameter("TAU["+idx+"]", tau, ">=0.0", "-1");

            param.multibath.multibath.add_bath(temp, tau);
            if (tau != -1) param.multibath.couple = true;
        }

        // now the DOF sets
        block.get_next_parameter("DOFSET", num_dof, ">=0", "");

        if (block.error()) {
          block.get_final_messages();
          return;
        }

        if (param.multibath.multibath.size() == 0) {
            if (num_dof != 0) {
                io::messages.add("MULTIBATH block: no baths but coupling groups specified",
                                 "In_Parameter", io::message::error);
                num_dof = 0;
            } else {
                io::messages.add("MULTIBATH block: no baths and no coupling groups => no temperature coupling",
                                 "In_Parameter", io::message::warning);
                param.multibath.couple = false;
                return;
            }

        }
        std::string str_bathsize = io::to_string(param.multibath.multibath.size());

        for (int i = 0; i < num_dof; ++i) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("LAST["+idx+"]", last, ">=1", "");
            block.get_next_parameter("COM-BATH["+idx+"]", com_bath, ">=1 && <="+ str_bathsize, "");
            block.get_next_parameter("IR-BATH["+idx+"]", ir_bath, ">=1 && <="+ str_bathsize, "");

            if (last < 1) last = 1;
            if (com_bath < 1) com_bath = 1;
            if (ir_bath < 1) ir_bath = 1;

            if (com_bath > param.multibath.multibath.size()) com_bath = param.multibath.multibath.size();
            if (ir_bath > param.multibath.multibath.size()) ir_bath = param.multibath.multibath.size();

            param.multibath.multibath.add_bath_index(last - 1, 0, com_bath - 1, ir_bath - 1);
        }

        block.get_final_messages();
    }
} // MULTIBATH

/**
 * @section positionres POSITIONRES block
 * @snippet snippets/snippets.cc POSITIONRES
 */
void io::In_Parameter::read_POSITIONRES(simulation::Parameter &param,
                                        std::ostream & os) {
    DEBUG(8, "reading POSITIONRES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "POSITIONRES\n";
    exampleblock << "#    NTPOR 0..3 controls atom positions re(con)straining.\n";
    exampleblock << "#          0: no position re(con)straints (default)\n";
    exampleblock << "#          1: restraining with force constant CPOR\n";
    exampleblock << "#          2: restraining with force constant CPOR weighted by\n";
    exampleblock << "#             atomic B-factors\n";
    exampleblock << "#          3: constraining\n";
    exampleblock << "#   NTPORB 0,1 controls reading of reference positions and\n";
    exampleblock << "#              B-factors\n";
    exampleblock << "#          0: read reference positions from startup file.\n";
    exampleblock << "#          1: read reference positions and B-factors from\n";
    exampleblock << "#             special file\n";
    exampleblock << "#   NTPORS 0,1 controls scaling of reference positions upon\n";
    exampleblock << "#              pressure scaling\n";
    exampleblock << "#          0: do not scale reference positions\n";
    exampleblock << "#          1: scale reference positions\n";
    exampleblock << "#     CPOR >= 0 position restraining force constant\n";
    exampleblock << "#\n";
    exampleblock << "#   NTPOR  NTPORB  NTPORS    CPOR\n";
    exampleblock << "        0       0       0   2.5E4\n";
    exampleblock << "END\n";


    std::string blockname = "POSITIONRES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int ntpor = 0, ntpors = 0, read = 0;
        block.get_next_parameter("NTPOR", ntpor, "", "0,1,2,3");
        block.get_next_parameter("NTPORB", read, "", "0,1");
        block.get_next_parameter("NTPORS", ntpors, "", "0,1");
        block.get_next_parameter("CPOR", param.posrest.force_constant, ">=0", "");

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
                break;
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
        }

        block.get_final_messages();
    }
} // POSITIONRES

/**
 * @section xrayres XRAYRES block
 * @snippet snippets/snippets.cc XRAYRES
 */
void io::In_Parameter::read_XRAYRES(simulation::Parameter &param,
                                    std::ostream & os) {
    DEBUG(8, "reading XRAYRES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "XRAYRES\n";
    exampleblock << "#    NTXR   -2: time-averaged electron density restraints\n";
    exampleblock << "#           -1: instantaneous electron density restraints\n";
    exampleblock << "#            0: no xray restraints.\n";
    exampleblock << "#            1: instantaneous structure factor restraints\n";
    exampleblock << "#            2: time-averaged structure factor restraints\n";
    exampleblock << "#            3: biquadratic/timeaveraged structure factor restraints\n";
    exampleblock << "#    NTXLE   0: do not perform local elevation\n";
    exampleblock << "#            1: do perform local elevation\n";
    exampleblock << "#    CXR     >= 0 xray restraining force constant\n";
    exampleblock << "#    NTWXR   >= 0 write xray data to output file\n";
    exampleblock << "#            0: don't write xray data\n";
    exampleblock << "#            > 0 write every NTPXRth step\n";
    exampleblock << "#    NTWDE   0..3 write density-maps\n";
    exampleblock << "#            0: write nothing\n";
    exampleblock << "#            1: write electron densitiy map\n";
    exampleblock << "#            2: write asymmetric-unit-only electron densitiy map\n";
    exampleblock << "#            3: write both\n";
    exampleblock << "#    NTWXM   >= 0 write every NTWXMth step electron density map(s) to external file\n";
    exampleblock << "#    CXTAU   >=0 xray time-average restraining memory-constant\n";
    exampleblock << "#    RDAVG   0/1 read sf-timeaverages (from job to job)\n";
    exampleblock << "#\n";
    exampleblock << "#   NTXR   NTXLE     CXR   NTWXR   NTWDE   NTWXM   CXTAU   RDAVG\n";
    exampleblock << "       0       0     0.0       0       0      0      0.0       0\n";
    exampleblock << "END\n";


    std::string blockname = "XRAYRES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);
        param.setDevelop("XRAY restraining is under development.");

        int ntxr = 0, ntxle = 0;
        block.get_next_parameter("NTXR", ntxr, "", "-2,-1,0,1,2,3");
        block.get_next_parameter("NTXLE", ntxle, "", "0,1");
        block.get_next_parameter("CXR", param.xrayrest.force_constant, ">=0", "");
        block.get_next_parameter("NTWXR", param.xrayrest.write, ">=0", "");
        block.get_next_parameter("NTWDE", param.xrayrest.writedensity, "", "0,1,2,3");
        block.get_next_parameter("NTWXM", param.xrayrest.writexmap, ">=0", "");
        block.get_next_parameter("CXTAU", param.xrayrest.tau, ">0", "");
        block.get_next_parameter("RDAVG", param.xrayrest.readavg, "", "0,1");

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
                break;
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
                break;
        }

        if (param.xrayrest.local_elevation) {
            if (param.xrayrest.mode != simulation::xrayrest_mode_electron_density ||
                param.xrayrest.xrayrest == simulation::xrayrest_biq) {
                io::messages.add("XRAYRES block: NTXLE 1 requires NTXR -2 or -1.",
                                 "In_Parameter", io::message::error);
            }
        }

        block.get_final_messages();
    }
} // XRAYRES

/**
 * @section distanceres DISTANCERES block
 * @snippet snippets/snippets.cc DISTANCERES
 */
void io::In_Parameter::read_DISTANCERES(simulation::Parameter &param,
                                        std::ostream & os) {
    DEBUG(8, "reading DISTANCERES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "DISTANCERES\n";
    exampleblock << "#   NTDIR -2..2 controls distance restraining\n";
    exampleblock << "#         0: no distance restraining (default)\n";
    exampleblock << "#         1: instantaneous, using force constant CDIR\n";
    exampleblock << "#         2: instantaneous, using force constant CDIR x W0\n";
    exampleblock << "#        -1: time-averaged, using force constant CDIR\n";
    exampleblock << "#        -2: time-averaged, using force constant CDIR x W0\n";
    exampleblock << "#  NTDIRA 0,1 controls values for initial distance averages\n";
    exampleblock << "#         0: generate initial averages\n";
    exampleblock << "#         1: read from configuration\n";
    exampleblock << "#    CDIR >= 0.0 force constant for distance restraining\n";
    exampleblock << "#    DIR0 > 0.0 distance offset in restraining function\n";
    exampleblock << "#  TAUDIR > 0.0 coupling time for time averaging\n";
    exampleblock << "# FORCESCALE 0..2 controls approximation of force scaling\n";
    exampleblock << "#         0: approximate d<r>/dr = 1\n";
    exampleblock << "#         1: approximate d<r>/dr = (1.0 - exp(-Dt/tau))\n";
    exampleblock << "#         2: use d<r>/dr = (1.0 - exp(-Dt/tau))*(<r>/r)^4\n";
    exampleblock << "#    VDIR 0,1 controls contribution to virial\n";
    exampleblock << "#         0: no contribution\n";
    exampleblock << "#         1: distance restraints contribute to virial\n";
    exampleblock << "#  NTWDIR >= 0 write every NTWDIRth step dist. restr. information to external file\n";
    exampleblock << "#   NTDIR  NTDIRA    CDIR    DIR0  TAUDIR  FORCESCALE  VDIR  NTWDIR\n";
    exampleblock << "        0       0     0.0     1.0     0.0           0     0       0\n";
    exampleblock << "END\n";


    std::string blockname = "DISTANCERES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int ntdira = 0;
        block.get_next_parameter("NTDIR", param.distanceres.distanceres, "", "0, 1, -1, 2, -2");
        block.get_next_parameter("NTDIRA", ntdira, "", "0,1");
        block.get_next_parameter("CDIR", param.distanceres.K, ">=0", "");
        block.get_next_parameter("DIR0", param.distanceres.r_linear, ">=0", "");
        block.get_next_parameter("TAUDIR", param.distanceres.tau, ">0", "");
        block.get_next_parameter("FORCESCALE", param.distanceres.forcescale, "", "0, 1, 2");
        block.get_next_parameter("VDIR", param.distanceres.virial, "", "0,1");
        block.get_next_parameter("NTWDIR", param.distanceres.write, ">=0", "");

        switch (ntdira) {
            case 0: param.distanceres.read = false;
                break;
            case 1: param.distanceres.read = true;
                break;
            default: param.distanceres.read = false;
                break;
        }

        if (param.distanceres.read && param.distanceres.distanceres > 0) {
            io::messages.add("DISTANCERES block: NTDIRA=1 but NTDIR > 0 - DISRESEXPAVE ignored",
                             "In_Parameter", io::message::warning);
        }

        block.get_final_messages();
    }
} // DISTANCERES

/**
 * @section distancefield DISTANCEFIELD block
 * @snippet snippets/snippets.cc DISTANCEFIELD
 */
void io::In_Parameter::read_DISTANCEFIELD(simulation::Parameter &param,
                                          std::ostream & os) {
    DEBUG(8, "reading DISTANCEFIELD");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "DISTANCEFIELD\n";
    exampleblock << "#   NTDFR 0,1         controls distance field restraining\n";
    exampleblock << "#         0: no distance field restraining (default)\n";
    exampleblock << "#         1: apply distance field restraining\n";
    exampleblock << "#   GRID  > 0.0       grid size for distance field\n";
    exampleblock << "#   PROTEINOFFSET >= 0 penalty for distances through the host\n";
    exampleblock << "#   PROTEINCUTOFF >= 0 distance to protein atoms to be considered inside\n";
    exampleblock << "#   PROTECT >= 0      protect grid points within this radius around the zero-distance\n";
    exampleblock << "#                     point from being flagged as protein\n";
    exampleblock << "#   UPDATE > 0        update frequency for grid\n";
    exampleblock << "#   RL >= 0           linearize forces for distances larger than RL\n";
    exampleblock << "#   SMOOTH >= 0       smoothen the protein boundary after grid construction\n";
    exampleblock << "#                     by SMOOTH layers\n";
    exampleblock << "#   NTWDF >= 0        write every NTWDF step disfield information to external file\n";
    exampleblock << "#   PRINTGRID = 0,1   write grid to final configuration file\n";
    exampleblock << "#\n";
    exampleblock << "#   NTDFR\n";
    exampleblock << "        1\n";
    exampleblock << "#    GRID   PROTEINOFFSET  PROTEINCUTOFF  PROTECT\n";
    exampleblock << "      0.2   15             0.2            0\n";
    exampleblock << "#  UPDATE   SMOOTH   RL    NTWDF   PRINTGRID\n";
    exampleblock << "      100   1        1.0      50           0\n";
    exampleblock << "END\n";


    std::string blockname = "DISTANCEFIELD";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int printgrid = -1;

        block.get_next_parameter("NTDFR", param.distancefield.distancefield, "", "0,1");
        block.get_next_parameter("GRID", param.distancefield.grid, ">0.0", "");
        block.get_next_parameter("PROTEINOFFSET", param.distancefield.proteinoffset, ">=0.0", "");
        block.get_next_parameter("PROTEINCUTOFF", param.distancefield.proteincutoff, ">=0.0", "");
        block.get_next_parameter("PROTECT", param.distancefield.protect, ">=0", "");
        block.get_next_parameter("UPDATE", param.distancefield.update, ">0", "");
        block.get_next_parameter("SMOOTH", param.distancefield.smooth, ">=0", "");
        block.get_next_parameter("RL", param.distancefield.r_l, ">=0", "");
        block.get_next_parameter("NTWDF", param.distancefield.write, ">=0", "");
        block.get_next_parameter("PRINTGRID", printgrid, "", "0,1");

        if(printgrid==1) param.distancefield.printgrid = true;
        else param.distancefield.printgrid = false;

        block.get_final_messages();
    }
} // DISTANCEFIELD


/**
 * @section angleres ANGLERES block
 * @snippet snippets/snippets.cc ANGLERES
 */
void io::In_Parameter::read_ANGLERES(simulation::Parameter &param,
                                        std::ostream & os) {
    DEBUG(8, "reading ANGLERES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "ANGLERES\n";
    exampleblock << "# NTALR   0...3 controls angle restraining and constraining\n";
    exampleblock << "#         0:    off [default]\n";
    exampleblock << "#         1:    angle restraining using CALR\n";
    exampleblock << "#         2:    angle restraining using CALR * WALR\n";
    exampleblock << "#         3:    angle constraining\n";
    exampleblock << "#\n";
    exampleblock << "# CALR    >=0.0 force constant for angle restraining [kJ/mol/degree^2]\n";
    exampleblock << "# VARES 0,1 controls contribution to virial\n";
    exampleblock << "#         0: no contribution\n";
    exampleblock << "#         1: angle restraints contribute to virial\n";
    exampleblock << "# NTWALR  >=0   write every NTWALR step angle restraint information to external file\n";
    exampleblock << "# TOLBAC  >0    tolerance for constraint deviation (in degrees)\n";
    exampleblock << "#\n";
    exampleblock << "# NTALR  CALR  VARES    NTWALR TOLBAC\n";
    exampleblock << "  1      1.0      0       100    0.01\n";
    exampleblock << "END\n";


    std::string blockname = "ANGLERES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        double K = 0.0, tolerance = 0.0;
        int angrest = 0;
        block.get_next_parameter("NTALR", angrest, "", "0,1,2,3");
        block.get_next_parameter("CALR", K, ">=0", "");
        block.get_next_parameter("VARES", param.angrest.virial, "", "0,1");
        block.get_next_parameter("NTWALR", param.angrest.write, ">=0", "");
        block.get_next_parameter("TOLBAC", tolerance, ">=0", "");

        switch (angrest) {
            case 0:
                param.angrest.angrest = simulation::angle_restr_off;
                break;
            case 1:
                param.angrest.angrest = simulation::angle_restr_inst;
                break;
            case 2:
                param.angrest.angrest = simulation::angle_restr_inst_weighted;
                break;
            case 3:
                param.angrest.angrest = simulation::angle_constr;
                break;
            default:
                break;
        }

        param.angrest.K = K*180*180 / (math::Pi * math::Pi);
        param.angrest.tolerance = tolerance * math::Pi / 180;

        if (param.angrest.angrest == simulation::angle_constr) {
            if (param.constraint.ntc == 1 && param.constraint.solute.algorithm == simulation::constr_off)
                param.constraint.solute.algorithm = simulation::constr_shake;

            if (param.constraint.solute.algorithm != simulation::constr_shake) {
                io::messages.add("ANGLERES block: needs SHAKE as (solute) constraints algorithm",
                                 "In_Parameter",
                                 io::message::error);
            }
        }
        block.get_final_messages();
    }
} // ANGLERES

/**
 * @section dihedralres DIHEDRALRES block
 * @snippet snippets/snippets.cc DIHEDRALRES
 */
void io::In_Parameter::read_DIHEDRALRES(simulation::Parameter &param,
                                        std::ostream & os) {
    DEBUG(8, "reading DIHEDRALRES");
    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "DIHEDRALRES\n";
    exampleblock << "# NTDLR   0...3 controls dihedral-angle restraining and constraining\n";
    exampleblock << "#         0:    off [default]\n";
    exampleblock << "#         1:    dihedral restraining using CDLR\n";
    exampleblock << "#         2:    dihedral restraining using CDLR * WDLR\n";
    exampleblock << "#         3:    dihedral constraining (using sine and cosine)\n";
    exampleblock << "#\n";
    exampleblock << "# CDLR    >=0.0 force constant for dihedral restraining [kJ/mol/degree^2]\n";
    exampleblock << "# PHILIN  >0.0  deviation after which the potential energy function is linearized\n";
    exampleblock << "# VDIH 0,1 controls contribution to virial\n";
    exampleblock << "#         0: no contribution\n";
    exampleblock << "#         1: dihedral restraints contribute to virial\n";
    exampleblock << "# NTWDLR  >=0   write every NTWDLR step dihedral information to external file\n";
    exampleblock << "# TOLDAC  >0    tolerance for constraint deviation (in degrees)\n";
    exampleblock << "#\n";
    exampleblock << "# NTDLR  CDLR      PHILIN  VDIH  NTWDLR  TOLDAC\n";
    exampleblock << "  1      100.0     180.0     0   100       0.01\n";
    exampleblock << "END\n";


    std::string blockname = "DIHEDRALRES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        double phi_lin = 0.0, K = 0.0, tolerance = 0.0;
        int dihrest = 0;
        block.get_next_parameter("NTDLR", dihrest, "", "0,1,2,3");
        block.get_next_parameter("CDLR", K, ">=0", "");
        block.get_next_parameter("PHILIN", phi_lin, ">0", "");
        block.get_next_parameter("VDIH", param.dihrest.virial, "", "0,1");
        block.get_next_parameter("NTWDLR", param.dihrest.write, ">=0", "");
        block.get_next_parameter("TOLDAC", tolerance, ">=0", "");

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
                break;
            default:
                break;
        }

        param.dihrest.K = K*180*180 / (math::Pi * math::Pi);
        param.dihrest.phi_lin = phi_lin * 2 * math::Pi / 360;
        param.dihrest.tolerance = tolerance * math::Pi / 180;

        if (param.dihrest.dihrest == simulation::dihedral_constr) {
            if (param.constraint.ntc == 1 && param.constraint.solute.algorithm == simulation::constr_off)
                param.constraint.solute.algorithm = simulation::constr_shake;

            if (param.constraint.solute.algorithm != simulation::constr_shake) {
                io::messages.add("DIHEDRALRES block: needs SHAKE as (solute) constraints algorithm",
                                 "In_Parameter",
                                 io::message::error);
            }
        }
        block.get_final_messages();
    }
} // DIHEDRALRES

/**
 * @section jval JVALUERES block
 * @snippet snippets/snippets.cc JVALUERES
 */
void io::In_Parameter::read_JVALUERES(simulation::Parameter &param,
                                      std::ostream & os) {
    DEBUG(8, "reading JVALUERES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "JVALUERES\n";
    exampleblock << "# NTJVR    -3..2\n";
    exampleblock << "#          -3                biquadratic using CJVR * WJVR\n";
    exampleblock << "#          -2                time-averaged using CJVR * WJVR\n";
    exampleblock << "#          -1                time-avaraged using CJVR\n";
    exampleblock << "#           0                no J-value restraints [default]\n";
    exampleblock << "#           1                instantaneous using CJVR\n";
    exampleblock << "#           2                instantaneous using CJVR * WJVR\n";
    exampleblock << "# NTJVRA    0                controls reading of averages from startup file\n";
    exampleblock << "#           0                start from initial values of J0 [default]\n";
    exampleblock << "#           1                read time averages from startup file (for continuation time-averaged run)\n";
    exampleblock << "# CJVR   >= 0                J-value restraining force constant\n";
    exampleblock << "#                            (weighted by individual WJVR)\n";
    exampleblock << "# TAUJVR >= 0                coupling time for time-averaging\n";
    exampleblock << "# NJVRTARS  0,1              omits or includes force scaling by memory decay factor in case of time-averaging\n";
    exampleblock << "#           0                omit factor (set (1 - exp(-Dt/tau)) = 1)\n";
    exampleblock << "#           1                scale force by (1 - exp(-Dt/tau))\n";
    exampleblock << "# NJVRBIQW  0..2             controls weights (X in Eq. 19 of MD98.17) of the two terms in biquadratic restraining\n";
    exampleblock << "#           0                X = 1\n";
    exampleblock << "#           1                X = (1 - exp(-Dt/tau))\n";
    exampleblock << "#           2                X = 0\n";
    exampleblock << "# LE        0,1              local elevation restraining\n";
    exampleblock << "#           0                local elevation off [default]\n";
    exampleblock << "#           1                local elevation on\n";
    exampleblock << "# NGRID   > 1                number of grid points in local elevation restraining\n";
    exampleblock << "# DELTA  >= 0.0              no elevation of potential if J is within DELTA of J0\n";
    exampleblock << "# NTWJV  >= 0                write J-value averages and LE grid to special trajectory\n";
    exampleblock << "#           0                don't write [default]\n";
    exampleblock << "#         > 0                write every NTWJVth step\n";
    exampleblock << "#\n";
    exampleblock << "#       NTJVR  NTJVRA  CJVR   TAUJVR  NJVRTARS   NJVRBIQW   LE    NGRID   DELTA  NTWJV\n";
    exampleblock << "           -3  0       10.0      5.0     0          0       1       16     0.5      0\n";
    exampleblock << "END\n";


    std::string blockname = "JVALUERES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int ntjvr = 0;
        block.get_next_parameter("NTJVR", ntjvr, "", "-3,-2,-1,0,1,2");
        block.get_next_parameter("NTJVRA", param.jvalue.read_av, "", "0,1");
        block.get_next_parameter("CJVR", param.jvalue.K, ">=0", "");
        block.get_next_parameter("TAUJVR", param.jvalue.tau, ">0", "");
        block.get_next_parameter("NJVRTARS", param.jvalue.tarfscale, "", "0,1");
        block.get_next_parameter("NJVRBIQW", param.jvalue.biqweight, "", "0,1,2");
        block.get_next_parameter("LE", param.jvalue.le, "", "0,1");
        block.get_next_parameter("NGRID", param.jvalue.ngrid, ">1", "");
        block.get_next_parameter("DELTA", param.jvalue.delta, ">=0.0", "");
        block.get_next_parameter("NTWJV", param.jvalue.write, ">=0", "");

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
                break;
        }

        if (param.jvalue.read_av && (ntjvr >= 0 && !param.jvalue.le)) {
            io::messages.add("JVALUERES block: NTJVRA: Continuation only needed "
                             "with time-averaging or LE.",
                             "In_Parameter", io::message::error);
        }
        block.get_final_messages();
    }
} // JVALUERES

/**
 * @section orderparamres ORDERPARAMRES block
 * @snippet snippets/snippets.cc ORDERPARAMRES
 */
void io::In_Parameter::read_ORDERPARAMRES(simulation::Parameter &param,
                                          std::ostream & os) {
    DEBUG(8, "reading ORDERPARAMRES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "ORDERPARAMRES\n";
    exampleblock << "# NTOPR    -2..2\n";
    exampleblock << "#          -2                time-averaged using COPR * WOPR\n";
    exampleblock << "#          -1                time-averaged using COPR\n";
    exampleblock << "#           0                no order-parameter restraints [default]\n";
    exampleblock << "#           1                window-averaged using COPR\n";
    exampleblock << "#           2                window-averaged using COPR * WOPR\n";
    exampleblock << "# NTOPRA    0                controls reading of averages from startup file\n";
    exampleblock << "#           0                start from initial values of S0 [default]\n";
    exampleblock << "#           1                read time averages from startup file (for continuation time-averaged run)\n";
    exampleblock << "# COPR   >= 0.0              order-parameter restraining force constant\n";
    exampleblock << "#                            (weighted by individual WOPR)\n";
    exampleblock << "# TAUOPR >= 0.0              coupling time for time-averaging\n";
    exampleblock << "# UPDOPR  > 0                update average every UPDOPRth step\n";
    exampleblock << "# NTWOP  >= 0                write order-parameter to special trajectory\n";
    exampleblock << "#           0                don't write [default]\n";
    exampleblock << "#         > 0                write every NTWOP step\n";
    exampleblock << "#\n";
    exampleblock << "#       NTOPR  NTOPRA  COPR   TAUOPR   UPDOPR  NTWOP\n";
    exampleblock << "           -2  0       10.0      5.0        1      0\n";
    exampleblock << "END\n";

    std::string blockname = "ORDERPARAMRES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int ntopr = 0;

        block.get_next_parameter("NTOPR", ntopr, "", "-2,-1,0,1,2");
        block.get_next_parameter("NTOPRA", param.orderparamrest.read, "", "0,1");
        block.get_next_parameter("COPR", param.orderparamrest.K, ">=", "");
        block.get_next_parameter("TAUOPR", param.orderparamrest.tau, ">=", "");
        block.get_next_parameter("UPDOPR", param.orderparamrest.update_step, ">0", "");
        block.get_next_parameter("NTWOP", param.orderparamrest.write, ">=0", "");


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
                break;
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

        block.get_final_messages();
    }
} // ORDERPARAMRES


/**
 * @section rdcres RDCRES block
 * @snippet snippets/snippets.cc RDCRES
 */

void io::In_Parameter::read_RDCRES(simulation::Parameter &param,
                                   std::ostream & os)
{
    DEBUG(8, "reading RDCRES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "RDCRES\n";
    exampleblock << "# NTRDCR -4..2                 RDC restraining\n";
    exampleblock << "#        -4:                   biquadratic using CRDCR * WRDCR\n";
    exampleblock << "#        -3:                   biquadratic using CRDCR\n";
    exampleblock << "#        -2:                   time averaged using CRDCR * WRDCR\n";
    exampleblock << "#        -1:                   time averaged using CRDCR\n";
    exampleblock << "#         0:                   no RDC restraints [default]\n";
    exampleblock << "#         1:                   instantaneous using CRDCR\n";
    exampleblock << "#         2:                   instantaneous using CRDCR * WRDCR\n";
    exampleblock << "# NTRDCRA 0,1                  controls reading of average RDCs\n";
    exampleblock << "#         0:                   take initial values RDC0 from the RDC restraint file\n";
    exampleblock << "#         1:                   read time averages from initial coordinate file\n";
    exampleblock << "#                              (for continuation run)\n";
    exampleblock << "#\n";
    exampleblock << "# NTRDCT  0..2                 Type of alignment representation\n";
    exampleblock << "#         0:                   cartesian magnetic field vectors\n";
    exampleblock << "#         1:                   alignment tensor\n";
    exampleblock << "#         2:                   spherical harmonics\n";
    exampleblock << "# NTALR   0,1                  controls reading of values in the chosen representation\n";
    exampleblock << "#         0:                   start from values given in RDC restraint file\n";
    exampleblock << "#         1:                   read values from initial coordinate file (for continuation run)\n";
    exampleblock << "#\n";
    exampleblock << "# METHOD  0..2                 Method of updating the magnetic field vectors\n";
    exampleblock << "#         0:                   Energy minimisation\n";
    exampleblock << "#         1:                   Stochastic dynamics\n";
    exampleblock << "#         2:                   Molecular dynamics\n";
    exampleblock << "# EMGRAD  > 0.0                (METHOD = 0, EM) stop minimisation if gradient is below EMGRAD\n";
    exampleblock << "# EMDX0   > 0.0                (METHOD = 0, EM) initial step size\n";
    exampleblock << "# EMNMAX  > 0                  (METHOD = 0, EM) maximum number of minimisation steps\n";
    exampleblock << "# SDCFRIC >= 0.0               (METHOD = 1, SD) global friction coefficient gamma\n";
    exampleblock << "# TEMP  >= 0.0                 temperature of stochastic bath (SD) and temperature used for initial velocities (MD, SD)\n";
    exampleblock << "# DELTA   >= 0                 the flatbottom potential is 2 DELTA wide [ps^-1]\n";
    exampleblock << "# CRDCR   >= 0                 RDC restraining force constant [kJ*ps^2]\n";
    exampleblock << "#                              (weighted by individual WRDCR)\n";
    exampleblock << "# TAU     >= 0                 coupling time for time averaging [ps]\n";
    exampleblock << "# NRDCRTARS 0,1                omits or includes force scaling by memory decay factor in case of time-averaging\n";
    exampleblock << "#           0                  omit factor (set (1-exp(-dt/tau))=1 )\n";
    exampleblock << "#           1                  scale force by (1-exp(-dt/tau))\n";
    exampleblock << "# NRDCRBIQW 0..2               controls weights of the two terms in biquadratic restraining\n";
    exampleblock << "#           0                  X = 1\n";
    exampleblock << "#           1                  X = (1 - exp(-dt/tau))\n";
    exampleblock << "#           2                  X = 0\n";
    exampleblock << "# NTWRDC   >= 0                write output to special trajectory\n";
    exampleblock << "#           0:                 don't write\n";
    exampleblock << "#          >0:                 write every NTWRDCth step.\n";
    exampleblock << "#\n";
    exampleblock << "#      NTRDCR  NTRDCRA  NTRDCT  NTALR  METHOD\n";
    exampleblock << "            2        0       0      0       0\n";
    exampleblock << "#      EMGRAD  EMDX0  EMNMAX  SDCFRIC    TEMP    DELTA  CRDCR  TAU   NRDCRTARS NRDCRBIQW   NTWRDC   NTWRDC\n";
    exampleblock << "        0.001   0.01    1000       20     300      0      1     1           0    0          0        10000\n";
    exampleblock << "END\n";


    std::string blockname = "RDCRES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        //set at develop flag
        param.setDevelop("RDC restraining is under development!");

        int ntrdcr = 0, ntrdct = 0, method = 0;
        block.get_next_parameter("NTRDCR", ntrdcr, "", "-4,-3,-2,-1,0,1,2");
        block.get_next_parameter("NTRDCRA", param.rdc.read_av, "", "0,1");
        block.get_next_parameter("NTRDCT", ntrdct, "", "0,1,2");
        block.get_next_parameter("NTALR", param.rdc.read_align, "", "0,1");
        block.get_next_parameter("METHOD", method, "", "0,1,2");
        block.get_next_parameter("EMGRAD", param.rdc.emgradient, ">0", "");
        block.get_next_parameter("EMDX0", param.rdc.emstepsize, ">0", "");
        block.get_next_parameter("EMNMAX", param.rdc.emmaxiter, ">0", "");
        block.get_next_parameter("SDCFRIC", param.rdc.sdfric, ">=0", "");
        block.get_next_parameter("TEMP", param.rdc.temp, ">=0", "");
        block.get_next_parameter("DELTA", param.rdc.delta, ">=0", "");
        block.get_next_parameter("CRDCR", param.rdc.K, ">=0", "");
        block.get_next_parameter("TAU", param.rdc.tau, ">=0", "");
        block.get_next_parameter("NRDCRTARS", param.rdc.tAVfactor, "", "0,1");
        block.get_next_parameter("NRDCRBIQW", param.rdc.biqfactor, "", "0,1,2");
        block.get_next_parameter("NTWRDC", param.rdc.write, ">=0", "");

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
                break;
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
                break;
        }

        if(method != simulation::rdc_em) {
            io::messages.add("RDCRES block: Your choice of METHOD is legal. However, be aware, that rdc_em is faster and always a better choice except when you are experimenting.",
                             "In_Parameter", io::message::warning);
        }
        if(ntrdct != simulation::rdc_t) {
            io::messages.add("RDCRES block: Your choice of NTRDCT is legal. However, be aware, that rdc_t is faster and always a better choice except when you are experimenting.",
                             "In_Parameter", io::message::warning);
        }

        switch (ntrdcr) {
            case -4:
                param.rdc.mode = simulation::rdc_restr_biq_weighted;
                param.rdc.K *= 10e20;         // multiply with 10e20 ps^2  to roughly compensate for second quadratic term
                io::messages.add("Multyplying RDC force constant K with 10e20 ps^2 to compensate for second quadratic term", "In_Parameter", io::message::notice);
                break;
            case -3:
                param.rdc.mode = simulation::rdc_restr_biq;
                param.rdc.K *= 10e20;         // multiply with 10e20 ps^2  to roughly compensate for second quadratic term
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
                break;
        }

        if (param.rdc.mode != simulation::rdc_restr_off) {
            if (param.rdc.tau <= 0.0 && ( param.rdc.mode != simulation::rdc_restr_inst && param.rdc.mode != simulation::rdc_restr_inst_weighted)) {
                io::messages.add("RDCRES block: bad value for TAU, should be > 0.0",
                                 "In_Parameter", io::message::error);
            }
            if (param.rdc.K < 0.0) {
                io::messages.add("RDCRES block: bad value for CRDCR(K) in RDCRES block, should be > 0.0",
                                 "In_Parameter", io::message::error);
            }
            if (param.rdc.read_av && (param.rdc.mode != simulation::rdc_restr_av
                                      && param.rdc.mode != simulation::rdc_restr_av_weighted
                                      && param.rdc.mode != simulation::rdc_restr_biq
                                      && param.rdc.mode != simulation::rdc_restr_biq_weighted )) {
                io::messages.add("RDCRES block: Continuation only needed with averaging.",
                                 "In_Parameter", io::message::error);
            }
        }         // rdcs not switched off

        if(!(param.rdc.tAVfactor == 0 || param.rdc.tAVfactor == 1 ))
            io::messages.add("RDCRES block: The only possible values for NRDCRTARS are 0 and 1.", "In_Parameter", io::message::error);

        if(!(param.rdc.biqfactor == 0 || param.rdc.biqfactor == 1 || param.rdc.biqfactor == 2))
            io::messages.add("RDCRES block: The only possible values for NRDCRBIQW are 0, 1 and 2.", "In_Parameter", io::message::error);

        block.get_final_messages();
    }
} // RDCRES



/**
 * @section perscale PERSCALE block
 * @snippet snippets/snippets.cc PERSCALE
 */
void io::In_Parameter::read_PERSCALE(simulation::Parameter &param,
                                     std::ostream & os) {
    DEBUG(8, "reading PERSCALE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "PERSCALE\n";
    exampleblock << "# RESTYPE		special energy term to which periodic scaling should\n";
    exampleblock << "#            be applied\n";
    exampleblock << "#	   0		don't apply periodic scaling\n";
    exampleblock << "#	   1		apply periodic scaling to J-value restraints\n";
    exampleblock << "#\n";
    exampleblock << "# parameters for periodic scaling of J-value restraints\n";
    exampleblock << "# KDIH	>= 0		maximum scaling factor for dihedral angle potential\n";
    exampleblock << "# KJ	>= 0		maximum scaling factor for J-Value restraint potential\n";
    exampleblock << "# T	>= 0		period of cosine scaling function\n";
    exampleblock << "# DIFF	>= 0		minimum deviation from target value at which to start\n";
    exampleblock << "#                scaling period\n";
    exampleblock << "# RATIO	>= 0		minimum fraction of T that needs to be passed before\n";
    exampleblock << "#                starting a new scaling period\n";
    exampleblock << "# READ	   0,1		controls reading of scaling parameters\n";
    exampleblock << "#          0		reset scaling parameters\n";
    exampleblock << "#          1		read from configuration\n";
    exampleblock << "#\n";
    exampleblock << "# RESTYPE\n";
    exampleblock << "      0\n";
    exampleblock << "#    KDIH      KJ       T   DIFF    RATIO    READ\n";
    exampleblock << "      0.1     0.1     0.2    0.7      1.0       0\n";
    exampleblock << "END\n";


    std::string blockname = "PERSCALE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        std::string s1;
        int read = 0;
        block.get_next_parameter("RESTYPE", s1, "", "off, 0, jrest, 1");
        block.get_next_parameter("KDIH", param.pscale.KDIH, ">=0", "");
        block.get_next_parameter("KJ", param.pscale.KJ, ">=0", "");
        block.get_next_parameter("T", param.pscale.T, ">=0", "");
        block.get_next_parameter("DIFF", param.pscale.difference, ">=0", "");
        block.get_next_parameter("RATIO", param.pscale.ratio, "", "");
        block.get_next_parameter("READ", read, "", "");

        if (s1 == "jrest" || s1 == "1") {
            param.pscale.jrest = true;
        } else if (s1 == "off" || s1 == "0") {
            param.pscale.jrest = false;
        }

        switch (read) {
            case 0:
                param.pscale.read_data = false;
                break;
            case 1:
                param.pscale.read_data = true;
                break;
            default:
                break;
        }
        block.get_final_messages();
    }
} // PERSCALE

/**
 * @section rottrans ROTTRANS block
 * @snippet snippets/snippets.cc ROTTRANS
 */
void io::In_Parameter::read_ROTTRANS(simulation::Parameter &param,
                                     std::ostream & os) {
    DEBUG(8, "reading ROTTRANS");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "ROTTRANS\n";
    exampleblock << "# roto-translational constraints\n";
    exampleblock << "# use either centre of mass removal or roto-translational constraints\n";
    exampleblock << "# not both!\n";
    exampleblock << "#\n";
    exampleblock << "#     RTC: 0,1 controls roto-translational constraints\n";
    exampleblock << "#          0 don't use roto-translational constraints (default)\n";
    exampleblock << "#          1 use roto-translational constraints\n";
    exampleblock << "# RTCLAST: last atom to be roto-translationally constrained\n";
    exampleblock << "#     RTC  RTCLAST\n";
    exampleblock << "        1     1155\n";
    exampleblock << "END\n";


    std::string blockname = "ROTTRANS";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int rtc = 0;
        block.get_next_parameter("RTC", rtc, "", "0,1");
        block.get_next_parameter("RTCLAST", param.rottrans.last, ">0", "");

        switch (rtc) {
            case 0:
                param.rottrans.rottrans = false;
                break;
            case 1:
                param.rottrans.rottrans = true;
                break;
            default:
                break;
        }

        block.get_final_messages();
    }
}

/**
 * @section spc_loops INNERLOOP block
 * @snippet snippets/snippets.cc INNERLOOP
 */
void io::In_Parameter::read_INNERLOOP(simulation::Parameter &param,
                                      std::ostream & os) {
    DEBUG(8, "reading INNERLOOP");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "INNERLOOP\n";
    exampleblock << "# NTILM: 0..4, acceleration method used\n";
    exampleblock << "#        0: use standard solvent loops [default]\n";
    exampleblock << "#        1: use fast generic solvent loops\n";
    exampleblock << "#        2: use solvent loops with hardcoded parameters\n";
    exampleblock << "#        3: use solvent loops with tabulated forces and energies\n";
    exampleblock << "#        4: use solvent loops with CUDA library\n";
    exampleblock << "# NTILS: 0..1, solvent used\n";
    exampleblock << "#        0: use topology [default]\n";
    exampleblock << "#        1: use SPC\n";
    exampleblock << "# NGPUS: number of GPUs to use\n";
    exampleblock << "# NDEVG: Which GPU device number to use. If not given driver will determine.\n";
    exampleblock << "# NTILM NTILS NGPUS NDEVG\n";
    exampleblock << "      4     0     2   0 1\n";
    exampleblock << "END\n";


    std::string blockname = "INNERLOOP";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int method = 0, solvent = 0;
        block.get_next_parameter("NTILM", method, "", "0,1,2,3,4");
        block.get_next_parameter("NTILS", solvent, "", "0,1");

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
#ifdef HAVE_LIBCUDART
                // cuda library
                param.innerloop.method = simulation::sla_cuda;
#else
                param.innerloop.method = simulation::sla_off;
                io::messages.add("INNERLOOP block: CUDA solvent loops are not available "
                                 "in your compilation. Use --with-cuda for compiling.",
                                 "In_Parameter", io::message::error);
#endif
                break;
            }
            default:
            {
                param.innerloop.method = simulation::sla_off;
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
            }
        }

        // Get the number of gpus and their device number
        if (param.innerloop.method == simulation::sla_cuda) {
            block.get_next_parameter("NGPUS", param.innerloop.number_gpus, ">0", "");

            int temp = 0;
            bool fail = false;
            for (unsigned int i = 0; i < param.innerloop.number_gpus; i++) {
                std::string idx=io::to_string(i);
                if (block.get_next_parameter("NDEVG["+idx+"]", temp, ">=-1", "", true)) {
                    fail = true;
                    break;
                }
                param.innerloop.gpu_device_number.push_back(temp);
            }
            if (fail) {
                // if not enough device numbers are given, set all numibers to -1
                // and do not report this as an error
                param.innerloop.gpu_device_number.clear();
                param.innerloop.gpu_device_number.resize(param.innerloop.number_gpus, -1);
                io::messages.add("CUDA driver will determine devices for nonbonded interaction evaluation.",
                                    "In_Parameter", io::message::notice);
                block.reset_error();
            }
        }
        block.get_final_messages();
    }
} //INNERLOOP

/**
 * @section replica REPLICA block
 * @snippet snippets/snippets.cc REPLICA
 */
void io::In_Parameter::read_REPLICA(simulation::Parameter &param,
                                    std::ostream & os) {
    DEBUG(8, "reading REPLICA");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "REPLICA\n";
    exampleblock << "#    RETL >= 0   : turn off REplica exchange - Temperature and/or Lambda Coupled";
    exampleblock << "#             1   : turn on  ";
    exampleblock << "#     NRET >= 1 number of replica exchange temperatures\n";
    exampleblock << "#    RET() >= 0.0 temperature for each replica\n";
    exampleblock << "# LRESCALE 0,1 controls temperature scaling\n";
    exampleblock << "#          0 don't scale temperatures after exchange trial\n";
    exampleblock << "#          1 scale temperatures after exchange trial\n";
    exampleblock << "#   NRELAM >= 1 number of replica exchange lambda values\n";
    exampleblock << "#  RELAM() >= 0.0 lambda value for each lambda-replica\n";
    exampleblock << "#   RETS() >= 0.0 timestep of each lambda-replica\n";
    exampleblock << "# NRETRIAL >= 0 number of overall exchange trials\n";
    exampleblock << "#  NREQUIL >= 0 number of exchange periods to equilibrate\n";
    exampleblock << "#               (disallow switches)\n";
    exampleblock << "#     CONT >= 0 continuation run\n";
    exampleblock << "#             0 start from one configuration file\n";
    exampleblock << "#             1 start from multiple configuration files\n";
    exampleblock << "#\n";
    exampleblock << "# RETL\n";
    exampleblock << "  1\n";
    exampleblock << "# NRET\n";
    exampleblock << "  10\n";
    exampleblock << "# RET(1..NRET)\n";
    exampleblock << "  300.0  320.0  340.0 360.0 380.0\n";
    exampleblock << "  400.0  420.0  440.0 460.0 480.0\n";
    exampleblock << "# LRESCALE\n";
    exampleblock << "  1\n";
    exampleblock << "# NRELAM\n";
    exampleblock << "  10\n";
    exampleblock << "# RELAM(1..NRELAM)\n";
    exampleblock << "  0.0    0.1    0.2   0.3   0.4\n";
    exampleblock << "  0.5    0.6    0.7   0.8   0.9\n";
    exampleblock << "# RETS(1..NRELAM)\n";
    exampleblock << "  0.002  0.001  0.001 0.001 0.002\n";
    exampleblock << "  0.003  0.002  0.001 0.001 0.002\n";
    exampleblock << "# NERTRIAL\n";
    exampleblock << "  100\n";
    exampleblock << "# NREQUIL\n";
    exampleblock << "  10\n";
    exampleblock << "# CONT\n";
    exampleblock << "  0\n";
    exampleblock << "END\n";

    DEBUG(1, "REPLICA BLOCK\t START");

    std::string blockname = "REPLICA";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("RETL", param.replica.retl, "", "0,1");

        block.get_next_parameter("NRET", param.replica.num_T, ">=1", "");

        if (block.error()) {
          block.get_final_messages();
          return;
        }
        param.replica.temperature.resize(param.replica.num_T, 0.0);

        for (int i = 0; i < param.replica.num_T; ++i) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("RET["+idx+"]", param.replica.temperature[i], ">=0", "");
        }

        int scale = 0;
        block.get_next_parameter("LRESCALE", scale, "", "0,1");
        switch (scale) {
            case 0:
                param.replica.scale = false;
                break;
            case 1:
                param.replica.scale = true;
                break;
            default:
                break;
        }
        if (param.replica.scale && param.replica.num_T<2)
            io::messages.add("REPLICA block: NRESCALE=1 requires more than one replica exchange temperature.",
                             "In_Parameter", io::message::error);

        block.get_next_parameter("NRELAM", param.replica.num_l, ">=1", "");

        if (block.error()) {
          block.get_final_messages();
          return;
        }

        param.replica.lambda.resize(param.replica.num_l, 0.0);
        param.replica.dt.resize(param.replica.num_l, 0.0);

        for (int i = 0; i < param.replica.num_l; ++i) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("RELAM["+idx+"]", param.replica.lambda[i], ">=0 && <=1", "");
        }
        for (int i = 0; i < param.replica.num_l; ++i) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("RETS["+idx+"]", param.replica.dt[i], ">=0", "");
        }

        if (block.error()) {
          block.get_final_messages();
          return;
        }

        block.get_next_parameter("NRETRIAL", param.replica.trials, ">=0", "");
        block.get_next_parameter("NREQUIL", param.replica.equilibrate, ">=0", "");

        // do continuation run
        block.get_next_parameter("", param.replica.cont, "", "0,1");


        if (block.error()) {
            param.replica.num_T = 0;
            param.replica.num_l = 0;

            param.replica.temperature.clear();
            param.replica.lambda.clear();
            param.replica.dt.clear();
        }
        block.get_final_messages();
    }
    DEBUG(1, "REPLICA BLOCK\t DONE\n");

}

/**
 * @section replica REPLICA_EDS block
 * @verbatim
REPLICA_EDS
#    REEDS >= 0   : turn off Reeds
#             1   : turn on
#             2   : turn on 2D Reeds (s & Eoff)
#    NRES >= number of replica exchange eds smoothing values
#    NUMSTATES >= 2 Number of states
#    RES > 0 for each replica smoothing value
#    EIR (NUMSTATES X NRES): energy offsets for states and replicas
#    NRETRIAL >= 0 number of overall exchange trials
#    NREQUIL >= 0 number of exchange periods to equilibrate
#               (disallow switches)
#    CONT >= 0 continuation run
#             0 start from one configuration file
#             1 start from multiple configuration files
#    EDS_STAT_OUT >= 0     creates output files for each replica, which contains for each exchange trial
#                          the potential energies with the given coordinates for all s value. This data
#                           can be used to optimize the s distribution.
#                 0 eds stat turned off
#                 1 eds stat turned on
#   PERIODIC >= 0 2D periodic boundary (Eoff only)
#               0 periodic boundary off
#               1 periodic boundary on
#
#   REEDS
    1
#  NRES NUMSTATES NUMEOFF
    12  5 12
# RES(1 ... NRES)
  1.0 0.7 0.5 0.3 0.1 0.07 0.05 0.03 0.01 0.007 0.005 0.003
# EIR (NUMSTATES x NRES)
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
# NRETRIAL   NREQUIL    CONT    EDS_STAT_OUT PERIODIC
       10         0         1           1      0
END
@endverbatim
 */

void io::In_Parameter::read_REPLICA_EDS(simulation::Parameter &param, std::ostream & os){
    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "REPLICA_EDS                                                    \n";
    exampleblock << "#    REEDS >= 0   : turn off Reeds                             \n";
    exampleblock << "#             1   : turn on 1D Reeds (s)                       \n";
    exampleblock << "#             2   : turn on 2D Reeds (s & Eoff)                \n";
    exampleblock << "#             3   : turn on 1D Simulated Annealing-eds (s)                \n";
    exampleblock << "#    NRES >= number of replica exchange eds smoothing values   \n";
    exampleblock << "#    NEOFF >= number of replica exchange eds energy Offset vectors (only used for 2D REEDS - still needs to be present ;))   \n";
    exampleblock << "#    NUMSTATES >= 2 Number of states                              \n";
    exampleblock << "#    RES > 0 for each replica smoothing value                 \n";
    exampleblock << "#    EIR (NUMSTATES X NRES): energy offsets for states and replicas   \n";
    exampleblock << "#    NRETRIAL >= 0 number of overall exchange trials                  \n";
    exampleblock << "#    NREQUIL >= 0 number of exchange periods to equilibrate          \n";
    exampleblock << "#               (disallow switches)                                \n";
    exampleblock << "#     CONT >= 0 continuation run                                   \n";
    exampleblock << "#             0 start from one configuration file                  \n";
    exampleblock << "#             1 start from multiple configuration files            \n";
    exampleblock << "# EDS_STAT_OUT >= 0     creates output files for each replica, which contains for each exchange trial  \n";
    exampleblock << "#                       the potential energies with the given coordinates for all s value. This data   \n";
    exampleblock << "#                       can be used to optimize the s distribution.                                    \n";
    exampleblock << "#             0 eds stat turned off                                                                    \n";
    exampleblock << "#             1 eds stat turned on                                                                     \n";
    exampleblock << "# PERIODIC >= 0 2D periodic boundary (Eoff only)               \n";  //REMOVE THIS PART @bschroed
    exampleblock << "#             0 periodic boundary off                          \n";  //REMOVE THIS PART @bschroed
    exampleblock << "#             1 periodic boundary on                           \n";  //REMOVE THIS PART @bschroed
    exampleblock << "#          \n";
    exampleblock << "#  REEDS    \n";
    exampleblock << "   1       \n";
    exampleblock << "#  NRES  NUMSTATES  NEOFF  \n";
    exampleblock << "   12    5   12  \n";
    exampleblock << "# RES(1 ... NRES)  \n";
    exampleblock << "  1.0 0.7 0.5 0.3 0.1 0.07 0.05 0.03 0.01 0.007 0.005 0.003    \n";
    exampleblock << "# EIR (NUMSTATES x NEOFF)   \n";
    exampleblock << "  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  \n";
    exampleblock << "  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  \n";
    exampleblock << "  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  \n";
    exampleblock << "  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  \n";
    exampleblock << "  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  \n";
    exampleblock << "# NRETRIAL   NREQUIL    CONT    EDS_STAT_OUT    PERIODIC       \n";
    exampleblock << "       10         0         1           1           0          \n";
    exampleblock << "END\n";

    DEBUG(1, "REPLICA_EDS BLOCK\t START");
    //check that EDS Block was not read in before,because either REPLICA_EDS or EDS possible
    if (param.eds.eds) {
        std::ostringstream msg;
        msg << "REPLICA_EDS block cannot be used at the same time with EDS block";
          io::messages.add(msg.str(), "In_Parameter", io::message::error);
          return;
        }

    std::string blockname = "REPLICA_EDS";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        DEBUG(2, "REPLICA_EDS BLOCK: reading Block an translating to vars");
        //init_vars
        unsigned int reeds_control = 0, num_s = 0, num_states = 0, num_eoff=0;
        unsigned int ntrials = 0, nEquilibrate = 0, cont_run = 0, eds_stat_out=0;
        bool periodic=1;

        // GET BLOCKVARS
        //SYS Settings
        block.get_next_parameter("REEDS", reeds_control, "", "0,1,2,3");
        DEBUG(3, "REPLICA_EDS BLOCK: reeds_control= " << reeds_control);
        block.get_next_parameter("NRES", num_s, ">0", "");
        block.get_next_parameter("NUMSTATES", num_states, ">0", "");
        block.get_next_parameter("NEOFF", num_eoff, ">0", "");
        DEBUG(3, "REPLICA_EDS BLOCK: NEOFF= " << num_eoff);

        //get RES-Vector
        std::vector<double> s_vals(num_s, 0.0);
        for (unsigned int i = 0; i < num_s; i++) {
          std::string idx = io::to_string(i);
          block.get_next_parameter("RES[" + idx + "]", s_vals[i], "", "");
        }

        //get EIR-Matrix
        if(reeds_control == 1 || reeds_control == 3){ // in 1D REEDS case, NEOFFS must be == num_s
            num_eoff = num_s;
        }
        std::vector<std::vector<float>> eir(num_eoff);
        for(unsigned int replicaJ=0; replicaJ<num_eoff; replicaJ++){//init eir vectors
          std::vector<float> eir_vector_J(num_states, 0.0);
          eir[replicaJ] = eir_vector_J;
        }
        for (unsigned int stateI = 0; stateI < num_states; stateI++) {
          std::string stateI_idx = io::to_string(stateI);
          for (unsigned int replicaJ=0; replicaJ< num_eoff; replicaJ++){
            std::string replicaJ_idx = io::to_string(replicaJ);
            block.get_next_parameter("EIR[" + replicaJ_idx+ "]["+stateI_idx+"]", eir[replicaJ][stateI], "", "");    //Comment "this function only reads line by line! doesn't matter the indices in the string
          }
        }


        // general Settings
        block.get_next_parameter("NRETRIAL", ntrials, ">=0", "");
        block.get_next_parameter("NREQUIL", nEquilibrate, ">=0", "");
        block.get_next_parameter("CONT", cont_run, "", "0,1");
        block.get_next_parameter("EDS_STAT_OUT", eds_stat_out, "", "0,1");
        block.get_next_parameter("PERIODIC", periodic, "", "0,1");  //REMOVE THIS PART @bschroed
        DEBUG(3, "REPLICA_EDS BLOCK: PERIODIC= " << periodic);  //REMOVE THIS PART @bschroed


        // SET SETTINGS
        DEBUG(2, "REPLICA_EDS BLOCK: Set settings for sim.");
        // READ:REEDS control
        param.reeds.reeds = reeds_control;
        DEBUG(3, "REPLICA_EDS BLOCK: reeds_control= " << param.reeds.reeds);

        param.reeds.num_s = num_s;
        param.reeds.num_states = num_states;
        param.reeds.num_eoff = num_eoff;
        DEBUG(3, "REPLICA_EDS BLOCK: NEOFF= " << param.reeds.num_eoff);
        param.reeds.trials = ntrials;

        //param.reeds.svals.resize();
        param.reeds.svals=s_vals;

        //param.reeds.svals=s_vals;
        param.reeds.equilibrate = nEquilibrate;
        param.reeds.cont =cont_run;
        param.reeds.eds_stat_out = eds_stat_out;
        param.reeds.periodic = periodic;
        DEBUG(3, "REPLICA_EDS BLOCK: PERIODIC= " << periodic);

        // Replica temperatures - has to be the same for each replica // Not sure if this is optimal? bschroed
        // handle case that we have a stochastic dynamics simulation
        if(!param.stochastic.sd)
          param.reeds.temperature = param.multibath.multibath.bath(0).temperature;
        else
          param.reeds.temperature = param.stochastic.temp;

        DEBUG(2, "REPLICA_EDS BLOCK: assigned all reeds params");
        //set size of vectors in param.reeds
        switch(reeds_control) {
              case 0: case 1: case 3:
                  param.reeds.eds_para.resize(param.reeds.num_s);
                  param.reeds.dt.resize(param.reeds.num_s);
                  param.reeds.svals.resize(param.reeds.num_s);
                  break;
              case 2:
                  param.reeds.eds_para.resize(param.reeds.num_s * param.reeds.num_eoff);
                  param.reeds.dt.resize(param.reeds.num_s * param.reeds.num_eoff);
                  param.reeds.svals.resize(param.reeds.num_s * param.reeds.num_eoff);
                  break;
        }


        //Loop over all replicas in order to initialize complete eds_struct for each replica
        //initvars
        std::vector<double> dtV;    //is necessary to give replicas the paramesters
        std::vector<double> temperatureV;
        switch(reeds_control){
          case 0: case 1: case 3:
              for (int i = 0; i < param.reeds.num_s; ++i) {
                  dtV.push_back(param.step.dt);
                  temperatureV.push_back(param.reeds.temperature);

                  //READ:NUMSTATES
                  param.reeds.eds_para[i].eds = true;
                  param.reeds.eds_para[i].numstates=param.reeds.num_states;
                  //num_eoff not used in eds_struct only in reeds_struct
                  //param.reeds.eds_para[i].num_eoff=param.reeds.num_eoff;

                  //indicate only one parameter s used for reference state hamiltonian
                  param.reeds.eds_para[i].form = simulation::single_s;

                  //RES - give s_values
                  param.reeds.eds_para[i].s.resize(1);//only one parameter s per replica
                  param.reeds.eds_para[i].s[0]=s_vals[i];

                  //initialize size of EIR
                  param.reeds.eds_para[i].eir.resize(param.reeds.eds_para[i].numstates);

                  //init:
                  param.reeds.eds_para[i].visitedstates.resize(param.eds.numstates, false);
                  param.reeds.eds_para[i].visitcounts.resize(param.eds.numstates, 0);
                  param.reeds.eds_para[i].avgenergy.resize(param.eds.numstates, 0.0);

                  if(param.reeds.eds_para[i].s[0] < 0.0){
                      std::ostringstream msg;
                      msg << "REPLICA_EDS block: RES(" << i + 1 << ") must be >= 0.0";
                      io::messages.add(msg.str(), "In_Parameter", io::message::error);
                  }

                  DEBUG(2, "REPLICA_EDS BLOCK: assign all eds params - EIR");
		  //TODO: assertion such that NEOFF = NRES must hold!!
                  for(unsigned int j = 0; j < param.reeds.eds_para[0].numstates; ++j){
                      param.reeds.eds_para[i].eir[j] = eir[i][j];
                      DEBUG(3, "REPLICA_EDS BLOCK: eir[i][j]= " << param.reeds.eds_para[i].eir[j]);

                  }
              }
              break;
          case 2:
              for (int i = 0; i < param.reeds.num_s * param.reeds.num_eoff; ++i) {
                  dtV.push_back(param.step.dt);
                  temperatureV.push_back(param.reeds.temperature);

                  //READ:NUMSTATES
                  param.reeds.eds_para[i].eds = true;
                  param.reeds.eds_para[i].numstates=param.reeds.num_states;
                  //num_eoff not used in eds_struct only in reeds_struct
                  //param.reeds.eds_para[i].num_eoff=param.reeds.num_eoff;

                  //indicate only one parameter s used for reference state hamiltonian
                  param.reeds.eds_para[i].form = simulation::single_s;

                  //RES - give s_values
                  param.reeds.eds_para[i].s.resize(1);//only one parameter s per replica
                  param.reeds.eds_para[i].s[0]=s_vals[i%num_s];
                  DEBUG(3, "REPLICA_EDS BLOCK: s[i]= " << i << "\t" << param.reeds.eds_para[i].s[0] << "\n");

                  //initialize size of EIR
                  param.reeds.eds_para[i].eir.resize(param.reeds.eds_para[i].numstates);

                  //init:
                  param.reeds.eds_para[i].visitedstates.resize(param.eds.numstates, false);
                  param.reeds.eds_para[i].visitcounts.resize(param.eds.numstates, 0);
                  param.reeds.eds_para[i].avgenergy.resize(param.eds.numstates, 0.0);

                  if(param.reeds.eds_para[i].s[0] < 0.0){
                      std::ostringstream msg;
                      msg << "REPLICA_EDS block: RES(" << i + 1 << ") must be >= 0.0";
                      io::messages.add(msg.str(), "In_Parameter", io::message::error);
                  }
              }
              DEBUG(2, "REPLICA_EDS BLOCK: assign all eds params - EIR");
              int count = 0;
              for(int i = 0; i < param.reeds.num_eoff; ++i){
                for(int k = 0; k < param.reeds.num_s; ++k){
                  for(unsigned int j = 0; j < param.reeds.eds_para[0].numstates; ++j){
                      param.reeds.eds_para[count].eir[j] = eir[i][j];
                      DEBUG(3, "REPLICA_EDS BLOCK: eir[i][j]= " << count << "\t" << param.reeds.eds_para[count].eir[j]);

                  }
                  ++count;
                }
              }
              break;
        }

        DEBUG(2, "REPLICA_EDS BLOCK: assigned all eds params");

        // turn on eds for pertubation reading - Overwrite:
        param.eds.eds = true;
        param.eds.numstates = param.reeds.num_states;

        // turn not on Pertubation block!: hope thats ok -> killed warning
        param.perturbation.perturbation = false;

        //REPLICA Set replica settings:
        DEBUG(2, "REPLICA_EDS BLOCK: assign all replicas param:");

        // check whether all baths have the same temperature (unambiguous kT)
        param.replica.num_T =1;
        param.replica.temperature = temperatureV;
        param.replica.lambda = param.reeds.svals;
        param.replica.dt = dtV;
        param.replica.num_l = param.reeds.num_s ;
        param.replica.trials = param.reeds.trials;
        param.replica.equilibrate = param.reeds.equilibrate;
        param.replica.cont = param.reeds.cont;

        DEBUG(2, "REPLICA_EDS BLOCK: assigned all replicas param");

        //CHECK SETTINGS
        DEBUG(2, "REPLICA_EDS BLOCK: Check Settings:");

        if(!param.stochastic.sd){ // only check the multibath, if we don't have stochastic dynamics
          for (unsigned int i = 1; i < param.multibath.multibath.size(); i++) {
            if (param.multibath.multibath.bath(i).temperature !=
                    param.multibath.multibath.bath(0).temperature) {
              io::messages.add("Error in RE_EDS block: all baths must have the same temperature.",
                      "In_Parameter", io::message::error);
            }
          }
        }
        DEBUG(2, "REPLICA_EDS BLOCK: Checked Settings");
        block.get_final_messages();
    }
    DEBUG(1, "REPLICA_EDS BLOCK\t DONE\n");
}

/**
 * @section multicell MULTICELL block
 * @snippet snippets/snippets.cc MULTICELL
 */
void io::In_Parameter::read_MULTICELL(simulation::Parameter & param,
                                      std::ostream & os) {
    DEBUG(8, "reading MULTICELL");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "MULTICELL\n";
    exampleblock << "#  NTM: 0,1 switch for multiple-unit-cell simulation.\n";
    exampleblock << "#       0 : single-unit-cell simulation [default]\n";
    exampleblock << "#       1 : multiple-unit-cell simulation\n";
    exampleblock << "#         NTM\n";
    exampleblock << "            0\n";
    exampleblock << "#  number of subdivisions along axis\n";
    exampleblock << "#   NCELLA    NCELLB    NCELLC\n";
    exampleblock << "         1         1         1\n";
    exampleblock << "#  periodicity checks (relative tolerance)\n";
    exampleblock << "#  not available in md++ -> 0.0\n";
    exampleblock << "#    TOLPX     TOLPV     TOLPF    TOLPFW\n";
    exampleblock << "       0.0       0.0       0.0       0.0\n";
    exampleblock << "END\n";


    std::string blockname = "MULTICELL";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        param.setDevelop("MULTICELL is under development.");

        int ntm = 0;
        double tolpx = 0.0, tolpv = 0.0, tolpf = 0.0, tolpfw = 0.0;
        block.get_next_parameter("NTM", ntm, "", "0,1");
        block.get_next_parameter("NCELLA", param.multicell.x, ">=1", "");
        block.get_next_parameter("NCELLB", param.multicell.y, ">=1", "");
        block.get_next_parameter("NCELLC", param.multicell.z, ">=1", "");

        block.get_next_parameter("TOLPX", tolpx, "", "");
        block.get_next_parameter("TOLPV", tolpv, "", "");
        block.get_next_parameter("TOLPF", tolpf, "", "");
        block.get_next_parameter("TOLPW", tolpfw, "", "");

        switch (ntm) {
            case 1:
                param.multicell.multicell = true;
                break;
            case 0:
                param.multicell.multicell = false;
                break;
            default:
                param.multicell.multicell = false;
        }

        if (param.multicell.multicell) {
            /*
               // disable because broken
               param.multicell.multicell = false;
               io::messages.add("MULTICELL simulations are broken in MD++",
                               "In_Parameter", io::message::error);      */

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

        block.get_final_messages();
    }
}

/*
 * @section readtraj READTRAJ block
 * @snippet snippets/snippets.cc READTRAJ
 */
void io::In_Parameter::read_READTRAJ(simulation::Parameter & param,
                                     std::ostream & os) {
    DEBUG(8, "reading READTRAJ");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "READTRAJ\n";
    exampleblock << "# NTRD  0,1 controls trajectory-reevaluation mode\n";
    exampleblock << "#       0: do not use trajectory-reevaluation mode (default)\n";
    exampleblock << "#       1: use trajectory-reevaluation mode\n";
    exampleblock << "# NTSTR stride: should be the NTWX used to produce the analyzed trajectory\n";
    exampleblock << "# NTRB  read box (must be 1)\n";
    exampleblock << "# NTSHK 0,1 controls application of constraints\n";
    exampleblock << "#       0 apply constraints with respect to previous coordinates\n";
    exampleblock << "#       1 apply constraints with respect to current coordinates\n";
    exampleblock << "#       2 do not apply constraints (neither solute nor solvent)\n";
    exampleblock << "#\n";
    exampleblock << "#   NTRD   NTSTR    NTRB   NTSHK\n";
    exampleblock << "       0       0       1       0\n";
    exampleblock << "END\n";


    std::string blockname = "READTRAJ";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int ntrd = 0, ntrb = 0, ntshk = 0;
        block.get_next_parameter("NTRD", ntrd, "", "0,1");
        block.get_next_parameter("NTSTR", param.analyze.stride, "", "");
        block.get_next_parameter("NTRB", ntrb, "", "1");
        block.get_next_parameter("NTSHK", ntshk, "", "0,1,2");

        if (block.error()) {
          block.get_final_messages();
          return;
        }

        switch (ntrd) {
            case 1:
                param.analyze.analyze = true;
                io::messages.add("READTRAJ block: make sure NTSTR is set to the value of NTWX used for writing the trajectory to be analyzed!", "In_Parameter",
                             io::message::notice);
                break;
            case 0:
                param.analyze.analyze = false;
                break;
            default:
                break;
        }

        if (ntrb != 1)
            io::messages.add("READTRAJ block: NTRB must be 1.", "In_Parameter",
                             io::message::error);

        switch (ntshk) {
            case 1:
                param.analyze.copy_pos = true;
                param.analyze.no_constraints = false;
                break;
            case 0:
                param.analyze.copy_pos = false;
                param.analyze.no_constraints = false;
                break;
            case 2:
                param.analyze.copy_pos = false;
                param.analyze.no_constraints = true;
                break;
            default:
                break;
        }

        block.get_final_messages();
    }
} //READTRAJ

/**
 * @section integrate INTEGRATE block
 * @snippet snippets/snippets.cc INTEGRATE
 */
void io::In_Parameter::read_INTEGRATE(simulation::Parameter & param,
                                      std::ostream & os) {
    DEBUG(8, "reading INTEGRATE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "INTEGRATE\n";
    exampleblock << "#  NINT 0..1 selects integration method\n";
    exampleblock << "#	0: no integration performed\n";
    exampleblock << "#	1: leap-frog integration scheme performed (default)\n";
    exampleblock << "#\n";
    exampleblock << "#    NINT\n";
    exampleblock << "        1\n";
    exampleblock << "END\n";

    std::string blockname = "INTEGRATE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int nint = 0;
        block.get_next_parameter("NINT", nint, "", "0,1");

        switch (nint) {
            case 0:
                param.integrate.method = simulation::integrate_off;
                break;
            case 1:
                param.integrate.method = simulation::integrate_leap_frog;
                break;
            default:
                break;
        }
        block.get_final_messages();
    }
}

/**
 * @section stochdyn STOCHDYN block
 * @snippet snippets/snippets.cc STOCHDYN
 */
void io::In_Parameter::read_STOCHDYN(simulation::Parameter & param,
                                     std::ostream & os) {
    DEBUG(8, "reading STOCHDYN");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "STOCHDYN\n";
    exampleblock << "# NTSD    0,1 controls stochastic dynamics mode\n";
    exampleblock << "#         0: do not do stochastic dynamics (default)\n";
    exampleblock << "#         1: do stochastic dynamics\n";
    exampleblock << "# NTFR    0..3 defines atomic friction coefficients gamma\n";
    exampleblock << "#         0: set gamma to 0.0 (default)\n";
    exampleblock << "#         1: set gamma to CFRIC\n";
    exampleblock << "#         2: set gamma to CFRIC*GAM0\n";
    exampleblock << "#         3: set gamma to CFRIC*w where w approximates the solvent-accessible \n";
    exampleblock << "#            surface area as described in the Stochastic Dynamics Chapter in Vol.2 of the manual \n";
    exampleblock << "# NSFR    > 0 recalculate gamma every NSFR steps\n";
    exampleblock << "# NBREF   > 0 threshold number of neighbour atoms for a buried atom\n";
    exampleblock << "# RCUTF   >= 0.0 interatomic distance considered when calculating gamma\n";
    exampleblock << "# CFRIC   >= 0.0 global weighting for gamma\n";
    exampleblock << "# TEMPSD  >= 0.0 temperature of stochastic bath\n";
    exampleblock << "#\n";
    exampleblock << "#     NTSD     NTFR     NSFR   NBREF  RCUTF    CFRIC    TEMPSD\n";
    exampleblock << "         0        1        0       6    0.3     91.0     300.0\n";
    exampleblock << "END\n";


    std::string blockname = "STOCHDYN";
    Block block(blockname, exampleblock.str());

    param.stochastic.sd=0;

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("NTSD", param.stochastic.sd, "", "0,1");
        block.get_next_parameter("NTFR", param.stochastic.ntfr, "", "0,1,2,3");
        block.get_next_parameter("NSFR", param.stochastic.nsfr, ">0", "");
        block.get_next_parameter("NBREF", param.stochastic.nbref, ">0", "");
        block.get_next_parameter("RCUTF", param.stochastic.rcutf, ">=0.0", "");
        block.get_next_parameter("CFRIC", param.stochastic.cfric, ">=0.0", "");
        block.get_next_parameter("TEMPSD", param.stochastic.temp, ">=0.0", "");

        block.get_final_messages();
    }
}

/**
 * @section ewarn EWARN block
 * @snippet snippets/snippets.cc EWARN
 */
void io::In_Parameter::read_EWARN(simulation::Parameter & param,
                                  std::ostream & os) {
    DEBUG(8, "reading EWARN");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "EWARN\n";
    exampleblock << "# MAXENER issue a warning if total energy is larger than this value\n";
    exampleblock << "#\n";
    exampleblock << "# MAXENER\n";
    exampleblock << "   100000\n";
    exampleblock << "END\n";

    std::string blockname = "EWARN";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("MAXENER", param.ewarn.limit, "", "");

        block.get_final_messages();
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
 * @snippet snippets/snippets.cc MULTISTEP
 */
void io::In_Parameter::read_MULTISTEP(simulation::Parameter & param,
                                      std::ostream & os) {
    DEBUG(8, "reading MULTISTEP");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "MULTISTEP\n";
    exampleblock << "#   STEPS calculate non-bonded every STEPSth step.\n";
    exampleblock << "#   BOOST 0,1\n";
    exampleblock << "#         0: stored forces of STEPSth step are added every step\n";
    exampleblock << "#         1: stored forces of STEPSth setp are multiplied by STEPS\n";
    exampleblock << "#            and added every STEPSth step (default)\n";
    exampleblock << "#\n";
    exampleblock << "#   STEPS   BOOST\n";
    exampleblock << "        0       0\n";
    exampleblock << "END\n";

    std::string blockname = "MULTISTEP";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int boost = 0;
        block.get_next_parameter("STEPS", param.multistep.steps, ">=0", "");
        block.get_next_parameter("BOOST", boost, "", "0,1");

        switch (boost) {
            case 0:
                param.multistep.boost = false;
                break;
            case 1:
                param.multistep.boost = true;
                break;
            default:
                break;
        }
        block.get_final_messages();
    }
}

/**
 * @section montecarlo CHEMICALMONTECARLO block
 * @snippet snippets/snippets.cc CHEMICALMONTECARLO
 */
void io::In_Parameter::read_CHEMICALMONTECARLO(simulation::Parameter & param,
                                       std::ostream & os) {
    DEBUG(8, "reading CHEMICALMONTECARLO");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "CHEMICALMONTECARLO\n";
    exampleblock << "#\n";
    exampleblock << "#     MC  MCSTEPS   MCDLAM\n";
    exampleblock << "       0        1      0.5\n";
    exampleblock << "END\n";

    std::string blockname = "CHEMICALMONTECARLO";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int mc = 0;
        block.get_next_parameter("MC", mc, "", "0,1");
        block.get_next_parameter("MCSTEPS", param.montecarlo.steps, ">=0", "");
        block.get_next_parameter("MCDLAM", param.montecarlo.dlambda, ">=0", "");

        switch (mc) {
            case 0: param.montecarlo.mc = mc;
                break;
            case 1: param.montecarlo.mc = mc;
                break;
            default:
                break;
        }

        block.get_final_messages();
    }
}

/**
 * @section polarise POLARISE block
 * @snippet snippets/snippets.cc POLARISE
 */
void io::In_Parameter::read_POLARISE(simulation::Parameter & param,
                                     std::ostream & os) {
    DEBUG(8, "reading POLARISE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "POLARISE\n";
    exampleblock << "# COS      0,1,2 use polarisation\n";
    exampleblock << "#          0: don't use polarisation (default)\n";
    exampleblock << "#          1: use charge-on-spring model for dipolar polarisation\n";
    exampleblock << "#          2: use charge-on-spring model for dipolar polarisation with off atom site\n";
    exampleblock << "# EFIELD   0,1 controls evaluation site for electric field\n";
    exampleblock << "#          0: evaluate at atom position\n";
    exampleblock << "#          1: evaluate at cos position\n";
    exampleblock << "# MINFIELD >0.0 convergence criterion in terms of the el. field (dU/(|qO|*dOH))\n";
    exampleblock << "#          where dU  .. maximum change in energy due to the change in field\n";
    exampleblock << "#                       [typically 2.5 kJ/mol]\n";
    exampleblock << "#                qO  .. charge of an oxygen atom\n";
    exampleblock << "#                dOH .. length of OH bond\n";
    exampleblock << "# DAMP     0,1 controls polarisability damping\n";
    exampleblock << "#          0: don't damp polarisability\n";
    exampleblock << "#          1: damp polarisability (with paramaters from topology)\n";
    exampleblock << "# WRITE    > 0 write COS positions to special trajectory\n";
    exampleblock << "#          0: don't write\n";
    exampleblock << "#         >0: write COS positions every WRITEth step\n";
    exampleblock << "#\n";
    exampleblock << "#     COS  EFIELD MINFIELD    DAMP  WRITE\n";
    exampleblock << "        0       0      2.5       0      0\n";
    exampleblock << "END\n";


    std::string blockname = "POLARISE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int cos = 0, damp = 0, efield = 0;
        block.get_next_parameter("COS", cos, "", "0,1,2");
        block.get_next_parameter("EFIELD", efield, "", "0,1");
        block.get_next_parameter("MINFIELD", param.polarise.minfield, ">0.0", "");
        block.get_next_parameter("DAMP", damp, "", "0,1");
        block.get_next_parameter("WRITE", param.polarise.write, ">=0", "");

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
                break;
        }

        switch (efield) {
            case 0: param.polarise.efield_site = simulation::ef_atom;
                break;
            case 1: param.polarise.efield_site = simulation::ef_cos;
                break;
            default:
                break;
        }

        switch (damp) {
            case 0: param.polarise.damp = false;
                break;
            case 1: param.polarise.damp = true;
                break;
            default:
                break;
        }

        if (!param.polarise.cos)
            param.polarise.write = 0;

        if (param.polarise.damp && !param.polarise.cos) {
            io::messages.add("POLARISE block: DAMP is ignored if no polarisation is used",
                             "In_Parameter", io::message::warning);
        }

        block.get_final_messages();
    }
} //POLARISE

/**
 * @section randomnumbers RANDOMNUMBERS block
 * @snippet snippets/snippets.cc RANDOMNUMBERS
 */
void io::In_Parameter::read_RANDOMNUMBERS(simulation::Parameter & param,
                                          std::ostream & os) {
    DEBUG(8, "reading RANDOMNUMBERS");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "RANDOMNUMBERS\n";
    exampleblock << "# NTRNG 0,1 random number generator\n";
    exampleblock << "#         0 use G96 algorithm (default)\n";
    exampleblock << "#         1 use GSL library\n";
    exampleblock << "# NTGSL -1.. GSL random number generation algorithm\n";
    exampleblock << "#         -1: use default algorithm (mt19937)\n";
    exampleblock << "#       >=0 : run contrib/rng_gsl for a list of possible arguments\n";
    exampleblock << "#\n";
    exampleblock << "#   NTRNG   NTGSL\n";
    exampleblock << "        1      -1\n";
    exampleblock << "END\n";


    std::string blockname = "RANDOMNUMBERS";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int rng = 0;
        block.get_next_parameter("NTRNG", rng, "", "0,1");
        block.get_next_parameter("NTGSL", param.rng.gsl_rng, ">=0", "-1");

        switch (rng) {
            case 0:
                param.rng.rng = simulation::random_g96;
                break;
            case 1:
                param.rng.rng = simulation::random_gsl;
                break;
            default:
                break;
        }

        math::RandomGenerator::check(param);
        block.get_final_messages();
    }
} // RANDOMNUMBERS

/**
 * @section EDS EDS block
 * @snippet snippets/snippets.cc EDS

 */
void io::In_Parameter::read_EDS(simulation::Parameter & param,
                                std::ostream & os) {
    DEBUG(8, "reading EDS");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line should be the blockname and is used as snippet tag
    exampleblock << "EDS\n";
    exampleblock << "# EDS        0,1\n";
    exampleblock << "#              0: no enveloping distribution sampling (EDS) [default]\n";
    exampleblock << "#              1: enveloping distribution sampling\n";
    exampleblock << "# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter\n";
    exampleblock << "#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter\n";
    exampleblock << "# FORM       1-3\n";
    exampleblock << "#              1: Single s Hamiltonian\n";
    exampleblock << "#              2: Hamiltonian with NUMSTATES*(NUMSTATES-1)/2 (pairwise) S parameters\n";
    exampleblock << "#              3: Hamiltonian with (NUMSTATES-1) S parameters\n";
    exampleblock << "# NUMSTATES >1  : number of states\n";
    exampleblock << "# if NUMSTATES != 3:\n";
    exampleblock << "# S         >0.0: smoothness parameter(s)\n";
    exampleblock << "# if NUMSTATES == 3:\n";
    exampleblock << "# i   j   S     : state pair i j and associated s parameter\n";
    exampleblock << "# EIR           : energy offsets for states\n";
    exampleblock << "#\n";
    exampleblock << "# EDS\n";
    exampleblock << "  1\n";
    exampleblock << "# ALPHLJ  ALPHC  FORM  NUMSTATES\n";
    exampleblock << "  0.0     0.0       2          3\n";
    exampleblock << "# S\n";
    exampleblock << "  0.2  0.01 0.1\n";
    exampleblock << "# EIR\n";
    exampleblock << "  0   20   40\n";
    exampleblock << "#\n";
    exampleblock << "# ---- OR: example for FORM = 3:\n";
    exampleblock << "#\n";
    exampleblock << "# EDS\n";
    exampleblock << "  1\n";
    exampleblock << "# ALPHLJ  ALPHC  FORM  NUMSTATES\n";
    exampleblock << "  0.0     0.0       3          3\n";
    exampleblock << "# i  j  S\n";
    exampleblock << "  1  2  0.1\n";
    exampleblock << "  2  3  0.5\n";
    exampleblock << "# EIR\n";
    exampleblock << "  0   20   40\n";
    exampleblock << "END\n";


    std::string blockname = "EDS";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int eds = 0, form = 0;
        double soft_lj = 0.0, soft_crf = 0.0;
        block.get_next_parameter("EDS", eds, "", "0,1");
        block.get_next_parameter("ALPHLJ", soft_lj, ">=0", "");
        block.get_next_parameter("ALPHC", soft_crf, ">=0", "");
        block.get_next_parameter("FORM", form, "", "1,2,3");
        block.get_next_parameter("NUMSTATES", param.eds.numstates, ">=2", "");

        switch (eds) {
            case 0:
                param.eds.eds = 0;
                param.eds.numstates = 0;
                break;
            case 1:
                param.eds.eds = 1;
                break;
            default:
                break;
        }

        if (!param.eds.eds) {
          block.get_final_messages();
          return;
        }

        param.eds.soft_vdw = soft_lj;
        param.eds.soft_crf = soft_crf;
        if (soft_lj > 0.0 || soft_crf > 0.0)
            param.setDevelop("Soft-core EDS is under development.");

        switch (form) {
            case 1: {
                param.eds.form = simulation::single_s;
                // read in 1 S value
                param.eds.s.resize(1, 1.0);
                block.get_next_parameter("S[0]", param.eds.s[0], ">0", "");
                break;
            }
            case 2: {
                param.eds.form = simulation::multi_s;
                const unsigned int n = param.eds.numstates;
                param.eds.s.resize((n * (n - 1)) / 2, 1.0);
                for (unsigned int pair = 0; pair < param.eds.s.size(); pair++) {
                    std::string idx = io::to_string(pair);
                    block.get_next_parameter("S["+idx+"]", param.eds.s[pair], ">0", "");
                }
                break;
            }
            case 3: {
                param.eds.form = simulation::pair_s;
                const unsigned int n = param.eds.numstates;
                param.eds.s.resize(n - 1, 1.0);
                param.eds.pairs.resize(n - 1);
                for (unsigned int pair = 0; pair < param.eds.s.size(); pair++) {
                    std::string idx = io::to_string(pair);
                    block.get_next_parameter("i["+idx+"]", param.eds.pairs[pair].i, "", "");
                    block.get_next_parameter("j["+idx+"]", param.eds.pairs[pair].j, "", "");
                    block.get_next_parameter("S["+idx+"]", param.eds.s[pair], ">0", "");
                }
                break;
            }
            default:
                break;
        }

        param.eds.eir.resize(param.eds.numstates, 0.0);
        for (unsigned int i = 0; i < param.eds.numstates; i++) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("EIR["+idx+"]", param.eds.eir[i], "", "");
        }

        block.get_final_messages();
    }
}

/**
* @section AEDS AEDS block
* @snippet snippets/snippets.cc AEDS

*/
void io::In_Parameter::read_AEDS(simulation::Parameter & param,
  std::ostream & os) {
  DEBUG(8, "reading AEDS");

  std::stringstream exampleblock;
  // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
  // will be used to generate snippets that can be included in the doxygen doc;
  // the first line should be the blockname and is used as snippet tag
  exampleblock << "AEDS\n";
  exampleblock << "# AEDS       0,1\n";
  exampleblock << "#              0: no accelerated enveloping distribution sampling (A-EDS) [default]\n";
  exampleblock << "#              1: accelerated enveloping distribution sampling\n";
  exampleblock << "# ALPHLJ: >= 0.0 Lennard-Jones soft-core parameter\n";
  exampleblock << "#  ALPHC: >= 0.0 Coulomb-RF soft-core parameter\n";
  exampleblock << "# FORM       1-4\n";
  exampleblock << "#              1: A-EDS with fixed parameters\n";
  exampleblock << "#              2: fixed Emax and Emin parameters, search for offset parameters\n";
  exampleblock << "#              3: search for Emax and Emin parameters, fixed offset parameters\n";
  exampleblock << "#              4: search for Emax, Emin and offset parameters\n";
  exampleblock << "# NUMSTATES >1  : number of states\n";
  exampleblock << "# EMAX          : A-EDS parameter Emax\n";
  exampleblock << "# EMIN          : A-EDS parameter Emin\n";
  exampleblock << "# EIR           : energy offsets for states\n";
  exampleblock << "# NTIAEDSS   0,1\n";
  exampleblock << "#              0: read A-EDS parameter search configuration from input configuration\n";
  exampleblock << "#              1: initialize A-EDS parameter search\n";
  exampleblock << "# RESTREMIN  0,1\n";
  exampleblock << "#              0: do not restrict Emin >= minimum average end-state energy\n";
  exampleblock << "#              1: restrict Emin >= minimum average end-state energy before all states have been visited at least once\n";
  exampleblock << "# BMAXTYPE   1,2\n";
  exampleblock << "#              1: absolute maximum energy barrier between the states in energy units\n";
  exampleblock << "#              2: multiples of the standard deviation of the energy of the end-state with the lowest average energy\n";
  exampleblock << "# BMAX          : maximum energy barrier parameter\n";
  exampleblock << "# ASTEPS        : have-life in simulation steps of the exponential averaged energy difference between the end-states at the begining of the run\n";
  exampleblock << "# BSTEPS        : have-life in simulation steps of the exponential averaged energy difference between the end-states at the end of the run\n";
  exampleblock << "#\n";
  exampleblock << "# AEDS\n";
  exampleblock << "  1\n";
  exampleblock << "# ALPHLJ  ALPHC  FORM  NUMSTATES\n";
  exampleblock << "  0.0     0.0       4          5\n";
  exampleblock << "# EMAX  EMIN\n";
  exampleblock << "  10    -50\n";
  exampleblock << "# EIR\n";
  exampleblock << "  0   -5   -140   -560   -74\n";
  exampleblock << "# NTIAEDSS  RESTREMIN  BMAXTYPE  BMAX  ASTEPS  BSTEPS\n";
  exampleblock << "  1         1          2         3     500     50000\n";
  exampleblock << "END\n";



  std::string blockname = "AEDS";
  Block block(blockname, exampleblock.str());

  if (block.read_buffer(m_block[blockname], false) == 0) {
    block_read.insert(blockname);

    int aeds = 0, form = 0;
    double soft_lj = 0.0, soft_crf = 0.0;
    block.get_next_parameter("AEDS", aeds, "", "0,1");
    block.get_next_parameter("ALPHLJ", soft_lj, ">=0", "");
    block.get_next_parameter("ALPHC", soft_crf, ">=0", "");
    block.get_next_parameter("FORM", form, "", "1,2,3,4");
    block.get_next_parameter("NUMSTATES", param.eds.numstates, ">=2", "");

    if (param.eds.eds != 1) {
      switch (aeds) {
      case 0:
        param.eds.eds = 0;
        param.eds.numstates = 0;
        break;
      case 1:
        param.eds.eds = 2;
        break;
      default:
        break;
      }
    }
    else {
      io::messages.add("AEDS block: AEDS cannot be used in combination with EDS.",
                                 "In_Parameter", io::message::error);
    }

    if (!param.eds.eds) {
      block.get_final_messages();
      return;
    }

    param.eds.soft_vdw = soft_lj;
    param.eds.soft_crf = soft_crf;
    if (soft_lj > 0.0 || soft_crf > 0.0)
      param.setDevelop("Soft-core EDS is under development.");

    switch (form) {
    case 1: {
      param.eds.form = simulation::aeds;
      break;
    }
    case 2: {
      param.eds.form = simulation::aeds_search_eir;
      break;
    }
    case 3: {
      param.eds.form = simulation::aeds_search_emax_emin;
      break;
    }
    case 4: {
      param.eds.form = simulation::aeds_search_all;
      break;
    }
    default:
      break;
    }

    block.get_next_parameter("EMAX", param.eds.emax, "", "");
    block.get_next_parameter("EMIN", param.eds.emin, "", "");

    if (param.eds.emin > param.eds.emax) {
      io::messages.add("AEDS paramater EMIN is larger than EMAX",
        "In_Parameter", io::message::warning);
      return;
    }

    param.eds.eir.resize(param.eds.numstates, 0.0);
    for (unsigned int i = 0; i < param.eds.numstates; i++) {
      std::string idx = io::to_string(i);
      block.get_next_parameter("EIR[" + idx + "]", param.eds.eir[i], "", "");
    }

    int ntia = 0, restremin = 0;
    block.get_next_parameter("NTIAEDSS",ntia, "", "0,1");
    switch (ntia) {
    case 0:
      param.eds.initaedssearch = false;
      break;
    case 1:
      param.eds.initaedssearch = true;
      break;
    default:
      break;
    }
    block.get_next_parameter("RESTREMIN", restremin, "", "0,1");
    switch (restremin) {
    case 0:
      param.eds.fullemin = true;
      break;
    case 1:
      param.eds.fullemin = false;
      break;
    default:
      break;
    }
    block.get_next_parameter("BMAXTYPE", param.eds.bmaxtype, "", "1,2");
    block.get_next_parameter("BMAX", param.eds.setbmax, ">0", "");
    block.get_next_parameter("ASTEPS", param.eds.asteps, ">0", "");
    block.get_next_parameter("BSTEPS", param.eds.bsteps, ">0", "");

    param.eds.searchemax = 0.0;
    param.eds.emaxcounts = 0;
    param.eds.oldstate = 0;

    param.eds.lnexpde.resize(param.eds.numstates, 0.0);
    param.eds.statefren.resize(param.eds.numstates, 0.0);
    param.eds.visitedstates.resize(param.eds.numstates, false);
    param.eds.visitcounts.resize(param.eds.numstates, 0);
    param.eds.avgenergy.resize(param.eds.numstates, 0.0);
    param.eds.eiravgenergy.resize(param.eds.numstates, 0.0);
    param.eds.bigs.resize(param.eds.numstates, 0.0);
    param.eds.stdevenergy.resize(param.eds.numstates, 0.0);

    block.get_final_messages();
  }
}

/**
 * @section LAMBDAS LAMBDAS block
 * @snippet snippets/snippets.cc LAMBDAS
 */
void io::In_Parameter::read_LAMBDAS(simulation::Parameter & param,
                                    std::ostream & os) {
    DEBUG(8, "reading LAMBDAS");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "LAMBDAS\n";
    exampleblock << "# NTIL    off(0), on(1)\n";
    exampleblock << "#         0: no special treatment of interactions with individual lambda-values\n";
    exampleblock << "#         1: interactions are treated with special individual lambda-values\n";
    exampleblock << "# NTLI(1..)  interaction type to treat with individual lambda:\n";
    exampleblock << "#            bond(1), angle(2), dihedral(3), improper(4), vdw(5), vdw_soft(6),\n";
    exampleblock << "#            crf(7), crf_soft(8), distanceres(9), distancefield(10),\n";
    exampleblock << "#            dihedralres(11), mass(12), angleres(13)\n";
    exampleblock << "# NILG1, NILG2 energy groups of interactions that are treated with individual\n";
    exampleblock << "#              lambda values\n";
    exampleblock << "# ALI, BLI, CLI, DLI, ELI polynomial coefficients linking the individual lambda-\n";
    exampleblock << "#                         values to the overall lambda-value\n";
    exampleblock << "# NTIL\n";
    exampleblock << "   1\n";
    exampleblock << "# NTLI NILG1  NILG2  ALI   BLI   CLI   DLI   ELI\n";
    exampleblock << "    7      1      3    0     0     1     0     0\n";
    exampleblock << "END\n";


    std::string blockname = "LAMBDAS";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int num = block.numlines()-3;

        std::string nm;
        simulation::interaction_lambda_enum j;
        int n1 = 0, n2 = 0;
        double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
        block.get_next_parameter("NTIL", nm, "", "on,off,1,0");
        DEBUG(10, "read NTIL " << nm);

        if (nm == "on" || nm == "1")
            param.lambdas.individual_lambdas = true;
        else if (nm == "off" || nm == "0") {
            param.lambdas.individual_lambdas = false;
            return;
        }

        if (param.perturbation.perturbation == false) {
            io::messages.add("LAMBDAS block without perturbation is ignored",
                             "In_Parameter", io::message::warning);
            return;
        }

        // the full matrices for the energy groups were created and
        // filled with a, b, c, e = 0 and d=1 when the FORCE block was read in
        // that way, they are also defined if you do not give the LAMBDAS block
        // we cannot initialize them earlier, because they depend on the
        // energy groups

        int maxnilg = param.force.energy_group.size();
        for (int i = 0; i < num; ++i) {
            std::string idx = io::to_string(i);
            block.get_next_parameter("NTLI["+idx+"]", nm, "", "1, bond, 2, angle, 3, dihedral, 4, improper, 5, vdw, 6, vdw_soft, 7, crf, 8, crf_soft, 9, distanceres, 10, distancefield, 11, dihedralres, 12, mass");
            block.get_next_parameter("NILG1["+idx+"]", n1, ">0", "");
            block.get_next_parameter("NILG2["+idx+"]", n2, ">0", "");
            block.get_next_parameter("ALI["+idx+"]", a, "", "");
            block.get_next_parameter("BLI["+idx+"]", b, "", "");
            block.get_next_parameter("CLI["+idx+"]", c, "", "");
            block.get_next_parameter("DLI["+idx+"]", d, "", "");
            block.get_next_parameter("ELI["+idx+"]", e, "", "");
            DEBUG(10, "read : " << nm << n1 << n2 << a << b << c << d << e);

            if (n1 > maxnilg || n2 > maxnilg) {
                io::messages.add("LAMBDAS block: NILG1 and NILG2 need to be smaller than the number of energy groups",
                                 "In_Parameter", io::message::error);
                return;
            }

            if (n2 < n1) {
                io::messages.add("only give NILG2 >= NILG1 in LAMBDAS BLOCK",
                                 "In_Parameter", io::message::error);
            }

            if (block.error()){
              block.get_final_messages();
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
        block.get_final_messages();
    }
}

/**
 * @section precalclam PRECALCLAM block
 * @snippet snippets/snippets.cc PRECALCLAM
 */
void io::In_Parameter::read_PRECALCLAM(simulation::Parameter & param,
        std::ostream & os) {
    DEBUG(8, "read PRECALCLAM");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "PRECALCLAM\n";
    exampleblock << "# NRLAM   0  : off\n";
    exampleblock << "#         >1 : precalculating energies for NRLAM extra lambda values\n";
    exampleblock << "# MINLAM  between 0 and 1: minimum lambda value to precalculate energies\n";
    exampleblock << "# MAXLAM  between MINLAM and 1: maximum lambda value to precalculate energies\n";
    exampleblock << "# NRLAM	  MINLAM   MAXLAM\n";
    exampleblock << "   100      0.0        1.0\n";
    exampleblock << "END\n";

    std::string blockname = "PRECALCLAM";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        block.get_next_parameter("NRLAM", param.precalclam.nr_lambdas, ">=0", "");
        block.get_next_parameter("MINLAM", param.precalclam.min_lam, ">=0 && <=1", "");
        block.get_next_parameter("MAXLAM", param.precalclam.max_lam, ">=0 && <=1", "");

        if (param.precalclam.min_lam >= param.precalclam.max_lam)
          io::messages.add("PRECALCLAM block: MINLAM should be smaller than MAXLAM",
            "In_Parameter", io::message::error);

    block.get_final_messages();

  }

} // PRECALCLAM

/**
 * @section nonbonded NONBONDED block
 * @snippet snippets/snippets.cc NONBONDED
 */
void io::In_Parameter::read_NONBONDED(simulation::Parameter & param,
                                      std::ostream & os) {
    DEBUG(8, "reading NONBONDED");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "NONBONDED\n";
    exampleblock << "# NLRELE    1-3 method to handle electrostatic interactions\n";
    exampleblock << "#    -1 : reaction-field (LSERF compatibility mode)\n";
    exampleblock << "#     0 : no electrostatic interactions\n";
    exampleblock << "#     1 : reaction-field\n";
    exampleblock << "#     2 : Ewald method\n";
    exampleblock << "#     3 : P3M method\n";
    exampleblock << "# APPAK     >= 0.0 reaction-field inverse Debye screening length\n";
    exampleblock << "# RCRF      >= 0.0 reaction-field radius\n";
    exampleblock << "#   0.0 : set to infinity\n";
    exampleblock << "# EPSRF     = 0.0 || > 1.0 reaction-field permittivity\n";
    exampleblock << "#   0.0 : set to infinity\n";
    exampleblock << "# NSLFEXCL  0,1 contribution of excluded atoms to reaction field\n";
    exampleblock << "#     0 : contribution turned off\n";
    exampleblock << "#     1 : contribution considered (default)\n";
    exampleblock << "# NSHAPE    -1..10 lattice sum charge-shaping function\n";
    exampleblock << "#    -1 : gaussian\n";
    exampleblock << "# 0..10 : polynomial\n";
    exampleblock << "# ASHAPE    > 0.0 width of the lattice sum charge-shaping function\n";
    exampleblock << "# NA2CALC   0..4 controls evaluation of lattice sum A2 term\n";
    exampleblock << "#     0 : A2 = A2~ = 0\n";
    exampleblock << "#     1 : A2~ exact, A2 = A2~\n";
    exampleblock << "#     2 : A2 numerical, A2~ = A2\n";
    exampleblock << "#     3 : A2~ exact from Ewald or from mesh and atom coords, A2 numerical\n";
    exampleblock << "#     4 : A2~ averaged from mesh only, A2 numerical\n";
    exampleblock << "# TOLA2     > 0.0 tolerance for numerical A2 evaluation\n";
    exampleblock << "# EPSLS      = 0.0 || > 1.0 lattice sum permittivity (0.0 = tinfoil)\n";
    exampleblock << "# NKX, NKY, NKZ > 0 maximum absolute Ewald k-vector components\n";
    exampleblock << "# KCUT       > 0.0 Ewald k-space cutoff\n";
    exampleblock << "# NGX, NGY, NGZ > 0 P3M number of grid points\n";
    exampleblock << "# NASORD    1..5 order of mesh charge assignment function\n";
    exampleblock << "# NFDORD    0..5 order of the mesh finite difference operator\n";
    exampleblock << "#     0 : ik - differentiation\n";
    exampleblock << "#  1..5 : finite differentiation\n";
    exampleblock << "# NALIAS    > 0 number of mesh alias vectors considered\n";
    exampleblock << "# NSPORD        order of SPME B-spline functions (not available)\n";
    exampleblock << "# NQEVAL    >= 0 controls accuracy reevaluation\n";
    exampleblock << "#     0 : do not reevaluate\n";
    exampleblock << "#   > 0 : evaluate every NQEVAL steps\n";
    exampleblock << "# FACCUR    > 0.0 rms force error threshold to recompute influence function\n";
    exampleblock << "# NRDGRD    0,1 read influence function\n";
    exampleblock << "#     0 : calculate influence function at simulation start up\n";
    exampleblock << "#     1 : read influence function from file (not yet implemented)\n";
    exampleblock << "# NWRGRD    0,1 write influence function\n";
    exampleblock << "#     0 : do not write\n";
    exampleblock << "#     1 : write at the end of the simulation (not yet implemented)\n";
    exampleblock << "# NLRLJ     0,1 controls long-range Lennard-Jones corrections\n";
    exampleblock << "#     0 : no corrections\n";
    exampleblock << "#     1 : do corrections (not yet implemented)\n";
    exampleblock << "# SLVDNS    > 0.0 average solvent density for long-range LJ correction (ignored)\n";
    exampleblock << "#\n";
    exampleblock << "#   NLRELE\n";
    exampleblock << "         1\n";
    exampleblock << "#    APPAK      RCRF     EPSRF   NSLFEXCL\n";
    exampleblock << "       0.0       1.4      61.0          1\n";
    exampleblock << "#   NSHAPE    ASHAPE    NA2CLC      TOLA2    EPSLS\n";
    exampleblock << "        -1       1.4         2     0.1E-9      0.0\n";
    exampleblock << "#      NKX       NKY       NKZ       KCUT\n";
    exampleblock << "        10        10        10      100.0\n";
    exampleblock << "#      NGX       NGY       NGZ    NASORD    NFDORD    NALIAS    NSPORD\n";
    exampleblock << "        32        32        32         3         2         3         4\n";
    exampleblock << "#   NQEVAL    FACCUR    NRDGRD    NWRGDR\n";
    exampleblock << "    100000       1.6         0         0\n";
    exampleblock << "#    NLRLJ    SLVDNS\n";
    exampleblock << "         0      33.3\n";
    exampleblock << "END\n";


    std::string blockname = "NONBONDED";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], true) == 0) {
        block_read.insert(blockname);

        int method = 0, ls_calculate_a2 = 0;

        block.get_next_parameter("NLRELE", method, "", "-1,0,1,2,3");
        block.get_next_parameter("APPAK", param.nonbonded.rf_kappa, ">=0", "");
        block.get_next_parameter("RCRF", param.nonbonded.rf_cutoff, ">=0", "");
        block.get_next_parameter("EPSRF", param.nonbonded.rf_epsilon, ">=1", "0.0");
        block.get_next_parameter("NSLFEXCL", param.nonbonded.selfterm_excluded_atoms, "", "0,1");
        block.get_next_parameter("NSHAPE", param.nonbonded.ls_charge_shape, ">=-1 && <=10", "");
        block.get_next_parameter("ASHAPE", param.nonbonded.ls_charge_shape_width, ">0", "");
        block.get_next_parameter("NA2CLC", ls_calculate_a2, "", "0,1,2,3,4");
        block.get_next_parameter("TOLA2", param.nonbonded.ls_a2_tolerance, ">0.0", "");
        block.get_next_parameter("EPSLS", param.nonbonded.ls_epsilon, ">=1.0", "0.0");
        block.get_next_parameter("NKX", param.nonbonded.ewald_max_k_x, ">0", "");
        block.get_next_parameter("NKY", param.nonbonded.ewald_max_k_y, ">0", "");
        block.get_next_parameter("NKZ", param.nonbonded.ewald_max_k_z, ">0", "");
        block.get_next_parameter("KCUT", param.nonbonded.ewald_kspace_cutoff, ">0.0", "");
        block.get_next_parameter("NGX", param.nonbonded.p3m_grid_points_x, ">0", "");
        block.get_next_parameter("NGY", param.nonbonded.p3m_grid_points_y, ">0", "");
        block.get_next_parameter("NGZ", param.nonbonded.p3m_grid_points_z, ">0", "");
        block.get_next_parameter("NASORD", param.nonbonded.p3m_charge_assignment, "", "1,2,3,4,5");
        block.get_next_parameter("NFDORD", param.nonbonded.p3m_finite_differences_operator, "", "0,1,2,3,4,5");
        block.get_next_parameter("NALIAS", param.nonbonded.p3m_mesh_alias, ">0", "");
        block.get_next_parameter("NSPORD", param.nonbonded.spme_bspline, "", "");
        block.get_next_parameter("NQEVAL", param.nonbonded.accuracy_evaluation, ">=0", "");
        block.get_next_parameter("FACCUR", param.nonbonded.influence_function_rms_force_error, ">0.0", "");
        block.get_next_parameter("NRDGRD", param.nonbonded.influence_function_read, "", "0,1");
        block.get_next_parameter("NWRGRD", param.nonbonded.influence_function_write, "", "0,1");
        block.get_next_parameter("NLRLJ", param.nonbonded.lj_correction, "", "0,1");
        block.get_next_parameter("SLVDNS", param.nonbonded.lj_solvent_density, ">0.0", "");

        if (block.error()) {
          block.get_final_messages();
          return;
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
                break;
        }

        // switch to turn of the contribution of the excluded atoms to the reaction field
        if (param.nonbonded.method == simulation::el_reaction_field
            && param.nonbonded.selfterm_excluded_atoms == 0) {
            param.nonbonded.rf_excluded = false;
            io::messages.add("NONBONDED block: contribution of excluded atoms to the reaction field turned off!",
                             "In_Parameter", io::message::warning);
        }

        if (param.nonbonded.method != simulation::el_reaction_field)
            param.force.interaction_function = simulation::lj_ls_func;

        if (do_ls && param.nonbonded.ls_charge_shape_width > param.pairlist.cutoff_short &&
            fabs(param.nonbonded.ls_charge_shape_width - param.pairlist.cutoff_short) > math::epsilon)
            io::messages.add("NONBONDED block: charge width greater than cutoff! (ASHAPE > RCUTP)",
                             "In_Parameter", io::message::warning);


        switch (ls_calculate_a2) {
            case 0:
                param.nonbonded.ls_calculate_a2 = simulation::ls_a2_zero;
                if (param.nonbonded.method == simulation::el_p3m || param.nonbonded.method == simulation::el_ewald)
                    io::messages.add("NONBONDED block: you are using p3m or ewald but NA2CLC=0.",
                                     "In_Parameter", io::message::warning);
                break;
            case 1:
                param.nonbonded.ls_calculate_a2 = simulation::ls_a2t_exact;
                if (param.nonbonded.method != simulation::el_ewald) {
                    io::messages.add("NONBONDED block: exact A2~ calculation needs Ewald.",
                                     "In_Parameter", io::message::error);
                }
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

        if (param.nonbonded.p3m_grid_points_x % 2 != 0 ||
            param.nonbonded.p3m_grid_points_y % 2 != 0 ||
            param.nonbonded.p3m_grid_points_z % 2 != 0)
            io::messages.add("NONBONDED block: Illegal value for NGA, NGB or NGC (even)",
                             "In_Parameter", io::message::error);

        if (param.nonbonded.influence_function_read ||
            param.nonbonded.influence_function_write)
            io::messages.add("NONBONDED block: Influence function IO not implemented."
                             " Set NRDGRD and NWRGRD to 0.",
                             "In_Parameter", io::message::error);

        if (param.nonbonded.lj_correction)
            io::messages.add("NONBONDED block: LJ long range correction Switched on."
                             "In_Parameter", io::message::warning);

        block.get_final_messages();
    }
} // NONBONDED

/**
 * @section sasa SASA block
 * @snippet snippets/snippets.cc SASA
 */
void io::In_Parameter::read_SASA(simulation::Parameter & param, std::ostream & os) {
    DEBUG(8, "reading SASA");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "SASA\n";
    exampleblock << "# NTSASA\n";
    exampleblock << "# 0 : not used (default)\n";
    exampleblock << "# 1 : use SASA implicit solvent model\n";
    exampleblock << "# NTVOL\n";
    exampleblock << "# 0 : not used (default)\n";
    exampleblock << "# 1 : use VOLUME correction to SASA implicit solvent model (requires NTSASA = 1)\n";
    exampleblock << "# P_12 >= 0, <= 1 pair parameter for SASA reduction for first neighbours\n";
    exampleblock << "# P_13 >= 0, <= 1 pair parameter for SASA reduction for second neighbours\n";
    exampleblock << "# P_1X >= 0, <= 1 pair parameter for SASA reduction for third and higher neighbours\n";
    exampleblock << "# SIGMAV >0 scaling parameter for volume energy term (kJ.mol^-1.nm^-3)\n";
    exampleblock << "# RSOLV > 0 radius of solvent molecule for SASA calculation (nm)\n";
    exampleblock << "# AS1 > 0 an atom with SASA below this contributes to the VOLUME correction (nm^2)\n";
    exampleblock << "# AS2 > 0 an atom with SASA above this is not considered for the VOLUME correction (nm^2)\n";
    exampleblock << "# atoms with AS1 < SASA < AS2 have a partial contribution determined by a switching function\n";
    exampleblock << "#   NTSASA      NTVOL       P_12      P_13     P_1X   SIGMAV  RSOlV    AS1    AS2\n";
    exampleblock << "         1          1     0.8875    0.3516   0.3516     -100   0.14   0.01   0.02\n";
    exampleblock << "END\n";


    std::string blockname = "SASA";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int switch_sasa = 0, switch_volume = 0;
        block.get_next_parameter("NTSASA", switch_sasa, "", "0,1");
        block.get_next_parameter("NTVOL", switch_volume, "", "0,1");
        block.get_next_parameter("P_12", param.sasa.p_12, ">=0.0 && <=1.0", "");
        block.get_next_parameter("P_13", param.sasa.p_13, ">=0.0 && <=1.0", "");
        block.get_next_parameter("P_1X", param.sasa.p_1x, ">=0.0 && <=1.0", "");
        block.get_next_parameter("SIGMAV", param.sasa.sigma_v, ">0", "");
        block.get_next_parameter("RSOLV", param.sasa.r_solv, ">0", "");
        block.get_next_parameter("AS1", param.sasa.min_cut, ">0", "");
        block.get_next_parameter("AS2", param.sasa.max_cut, ">0", "");

        switch(switch_sasa) {
            case 0:
                param.sasa.switch_sasa = false;
                break;
            case 1:
                param.sasa.switch_sasa = true;
                break;
            default:
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
                param.sasa.switch_volume = false;
        }

        // check that vol not used without sasa
        if (!param.sasa.switch_sasa && param.sasa.switch_volume) {
            io::messages.add("SASA block: Cannot have NTSASAVOL without NTSASA",
                             "In_Parameter", io::message::error);
        }

        // compute and store upper - lower cutoffs
        param.sasa.cut_diff = param.sasa.max_cut - param.sasa.min_cut;
        if (param.sasa.cut_diff < 0.0 ) {
            io::messages.add("SASA block: Cutoffs: min. is larger than max.",
                             "In_Parameter", io::message::error);
        }
        block.get_final_messages();
    }
}

/**
 * @section localelev LOCALELEV block
 * @snippet snippets/snippets.cc LOCALELEV
 */
void io::In_Parameter::read_LOCALELEV(simulation::Parameter & param,
                                      std::ostream & os) {
    DEBUG(8, "reading LOCALELEV");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "LOCALELEV\n";
    exampleblock << "# NTLES 0,1 controls the use of local elevation.\n";
    exampleblock << "#    0 : not used [default]\n";
    exampleblock << "#    1 : local elevation is applied\n";
    exampleblock << "# NLEPOT >= 0 number of umbrella potentials applied\n";
    exampleblock << "# NTLESA 1..2 controls the reading of the potential definition\n";
    exampleblock << "#    1 : read from startup file\n";
    exampleblock << "#    2 : read from special file (@lud)\n";
    exampleblock << "# NTWLE >= 0 write umbrellas to trajectory every NTWLEth step\n";
    exampleblock << "# NLEPID[1..NLEPOT] IDs of the umbrella potentials\n";
    exampleblock << "# NTLEPFR[1..NLEPOT] 0,1 freeze the umbrella potential\n";
    exampleblock << "#    0 : build up\n";
    exampleblock << "#    1 : freeze\n";
    exampleblock << "# NTLES  NLEPOT  NTLESA  NTWLE\n";
    exampleblock << "      1       2       1      0\n";
    exampleblock << "# NLEPID NLEPFR\n";
    exampleblock << "       1      0\n";
    exampleblock << "       2      1\n";
    exampleblock << "END\n";


    std::string blockname = "LOCALELEV";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int onoff = 0, num = 0, read = 0;
        block.get_next_parameter("NTLES", onoff, "", "0,1");
        block.get_next_parameter("NLEPOT", num, ">=0", "");
        block.get_next_parameter("NTLESA", read, "", "1,2");
        block.get_next_parameter("NTWLE", param.localelev.write, ">=0", "");

        switch (onoff) {
            case 0:
                param.localelev.localelev = simulation::localelev_off;
                break;
            case 1:
                param.localelev.localelev = simulation::localelev_on;
                break;
            default:
                param.localelev.localelev = simulation::localelev_off;
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
        }

        // read the umbrellas
        for (int i = 0; i < num; ++i) {
            int id = 0, f = 0;
            block.get_next_parameter("NLEPID", id, ">=1", "");
            block.get_next_parameter("NTLEPFR", f, "", "0,1");
            if (block.error()) {
                block.get_final_messages();
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
        }         // for umbrellas
        block.get_final_messages();
    }     // if block
}  //LOCALELEV

/**
 * @section bsleusparam BSLEUS block
 * @snippet snippets/snippets.cc BSLEUS
 */
void io::In_Parameter::read_BSLEUS(simulation::Parameter& param, std::ostream& os)
{
    DEBUG(8, "reading BSLEUS");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "BSLEUS\n";
    exampleblock << "#\n";
    exampleblock << "# The general settings for the B&S-LEUS algorithm\n";
    exampleblock << "# BSLEUS:   Dow we use B&S-LEUS?\n";
    exampleblock << "#   0:          Don'use it (default)\n";
    exampleblock << "#   1:          Use it\n";
    exampleblock << "# BUILD:    Are we building?\n";
    exampleblock << "#   0:          No\n";
    exampleblock << "#   1:          Yes\n";
    exampleblock << "# WRITE:    >= 0 Do we write the energies and forces of the Umbrella?\n";
    exampleblock << "#   == 0:          No\n";
    exampleblock << "#   > 0:           Every nth step\n";
    exampleblock << "#\n";
    exampleblock << "# BSLEUS    BUILD   WRITE\n";
    exampleblock << "  1         1       0\n";
    exampleblock << "END\n";

    std::string blockname = "BSLEUS";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int use_bsleus = 0, build = 0;
        block.get_next_parameter("BSLEUS", use_bsleus, "", "0,1");
        block.get_next_parameter("BUILD", build, "", "0,1");
        block.get_next_parameter("WRITE", param.bsleus.write, ">=0", "");

        if (block.error()) {
          block.get_final_messages();
          return;
        }


        param.bsleus.bsleus = use_bsleus ? simulation::bsleus_on :
                              simulation::bsleus_off;
        param.bsleus.building = build ? true : false;

        block.get_final_messages();
    }
}
/**
 * @section electric ELECTRIC block
 * @snippet snippets/snippets.cc ELECTRIC
 */
void io::In_Parameter::read_ELECTRIC(simulation::Parameter & param,
                                     std::ostream & os) {
    DEBUG(8, "reading ELECTRIC");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "ELECTRIC\n";
    exampleblock << "# FIELD 0,1 controls the use of applied electric field.\n";
    exampleblock << "#    0 : not used [default]\n";
    exampleblock << "#    1 : electric field is applied\n";
    exampleblock << "# DIPOLE 0,1 controls the calculation of box dipole.\n";
    exampleblock << "#    0 : not used [default]\n";
    exampleblock << "#    1 : box dipole is calculated and written to special trajectory\n";
    exampleblock << "# CURRENT 0,1 controls the calculation of electric currents.\n";
    exampleblock << "#    0 : not used [default]\n";
    exampleblock << "#    1 : electric current is calculated and written to special trajectory\n";
    exampleblock << "# ELECTRIC FIELD COMPONENTS (EF_x, EF_y, EF_z)\n";
    exampleblock << "# 0.0 0.0 0.0\n";
    exampleblock << "# DIPGRP 0..2 controls the groups considered for box dipole calculation\n";
    exampleblock << "#    0 : solute only\n";
    exampleblock << "#    1 : solvent only\n";
    exampleblock << "#    2 : all\n";
    exampleblock << "# NTWDIP >= 0 write dipole box every NTWDIPth step\n";
    exampleblock << "# NTWCUR >= 0 write current every NTWDIPth step\n";
    exampleblock << "# NCURGRP >=0 number of current groups\n";
    exampleblock << "# CURGRP [1..NCURGRP] last atom of the group\n";
    exampleblock << "#  FIELD  DIPOLE CURRENT\n";
    exampleblock << "       1       1       1\n";
    exampleblock << "#   EF_x    EF_y    EF_z\n";
    exampleblock << "     0.0     0.0     0.0\n";
    exampleblock << "# DIPGRP  NTWDIP\n";
    exampleblock << "       0       1\n";
    exampleblock << "# NTWCUR  NCURGRP   CURGRP[1]   CURGRP[2]\n";
    exampleblock << "       1       2        100        1000\n";
    exampleblock << "END\n";


    std::string blockname = "ELECTRIC";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int field = 0, dipole = 0, current = 0, ncurgrp = 0;
        block.get_next_parameter("FIELD", field, "", "0,1");
        block.get_next_parameter("DIPOLE", dipole, "", "0,1");
        block.get_next_parameter("CURRENT", current, "", "0,1");
        block.get_next_parameter("EF_x", param.electric.Ef_x, "", "");
        block.get_next_parameter("EF_y", param.electric.Ef_y, "", "");
        block.get_next_parameter("EF_z", param.electric.Ef_z, "", "");
        block.get_next_parameter("DIPGRP", param.electric.dip_groups, "", "0,1,2");
        block.get_next_parameter("NTWDIP", param.electric.dip_write, ">=0", "");
        block.get_next_parameter("NTWCUR", param.electric.cur_write, ">=0", "");
        block.get_next_parameter("NCURGRP", ncurgrp, ">=0", "");

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
                break;
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
                break;
            }
            default:
                break;
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
                param.electric.cur_groups = ncurgrp;
                break;
            }
            default:
                break;
        }

        if (param.electric.current != false) {
            // TO READ THE ELECTRIC GROUPS
            unsigned int temp = 0;
            for (unsigned int i = 0; i < param.electric.cur_groups; i++) {
                std::string idx = io::to_string(i);
                block.get_next_parameter("CURGRP["+idx+"]", temp, ">0", "");
                if (block.error()) break;
                param.electric.current_group.push_back(temp);
            }
        }
        block.get_final_messages();
    }     // if block
} //ELECTRIC

/**
 * @section nemd NEMD block
 * @snippet snippets/snippets.cc NEMD
 */
void io::In_Parameter::read_NEMD(simulation::Parameter & param,
                                 std::ostream & os) {
    DEBUG(8, "reading NEMD");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "NEMD\n";
    exampleblock << "# NEMD 0,1 controls the use of non-equilibrium molecular dynamics.\n";
    exampleblock << "#    0 : not used [default]\n";
    exampleblock << "#    1 : nemd is used\n";
    exampleblock << "# PROPERTY 0- select property to calculate\n";
    exampleblock << "#    0 : viscosity\n";
    exampleblock << "# METHOD 0- select method of NEMD.\n";
    exampleblock << "#    0 : periodic perturbation method (PPM)\n";
    exampleblock << "#    1 : internal reservoir method (IRM)\n";
    exampleblock << "# SLABNUM >=1 number of slabs used in the discretization along z-direction.\n";
    exampleblock << "#             the effective number is 2xSLABNUM due to periodicity\n";
    exampleblock << "# PERTFRQ >=1 perturbation frequency: apply perturbation every PERTFRQth timestep\n";
    exampleblock << "#             [this flag is ignored by the PPM method, but a value must be provided]\n";
    exampleblock << "# AMPLI   >=0 amplitude of applied field\n";
    exampleblock << "#             [this flag is ignored by the IRM method, but a value must be provided]\n";
    exampleblock << "# STDYAFT >=0 first STDYAFTth steps do not contribute for accumulated averages\n";
    exampleblock << "# WRITE >=1 write flux and average velocities to special trajectory every WRITEth timestep\n";
    exampleblock << "# NEMD     PROPERTY  METHOD\n";
    exampleblock << "    1         0        0\n";
    exampleblock << "# SLABNUM  PERTFRQ    AMPLI   STDYAFT   WRITE\n";
    exampleblock << "     10       20       10      1000     200\n";
    exampleblock << "END\n";

    std::string blockname = "NEMD";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);
        param.setDevelop("NEMD is under development.");

        int nemd = 0, method = 0;
        block.get_next_parameter("NEMD", nemd, "", "0,1");
        block.get_next_parameter("PROPERTY", param.nemd.property, "", "0");
        block.get_next_parameter("METHOD", method, "", "0,1");
        block.get_next_parameter("SLABNUM", param.nemd.slabnum, ">=0", "");
        block.get_next_parameter("PERTFRQ", param.nemd.pertfrq, ">=0", "");
        block.get_next_parameter("AMPLI", param.nemd.ampbath, ">=0", "");
        block.get_next_parameter("STDYAFT",  param.nemd.stdyaft, ">=0", "");
        block.get_next_parameter("WRITE", param.nemd.write, ">=0", "");


        if (block.error()) {
          block.get_final_messages();
          return;
        }

        param.nemd.nemd = nemd ? simulation::nemd_on : simulation::nemd_off;

        switch (method) {
            case 0:
            {
                param.nemd.method = 0;
                if (param.nemd.slabnum <=0 || param.nemd.ampbath <=0)
                    io::messages.add("NEMD block: PPM method used, but found invalid values for SLABNUM and AMPLI",
                                     "In_Parameter", io::message::error);
                break;
            }
            case 1:
            {
                param.nemd.method = 1;
                if (param.nemd.slabnum <=0 || param.nemd.pertfrq <=0)
                    io::messages.add("NEMD block: IRM method used, but found invalid values for SLABNUM and PERTFRQ",
                                     "In_Parameter", io::message::error);
                break;
            }
            default:
                break;
        }
        block.get_final_messages();
    }     // if block
} // NEMD

/**
 * @section multigradient MULTIGRADIENT block
 * @snippet snippets/snippets.cc MULTIGRADIENT
 */
void io::In_Parameter::read_MULTIGRADIENT(simulation::Parameter & param,
                                          std::ostream & os) {
    DEBUG(8, "reading MULTIGRADIENT");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "MULTIGRADIENT\n";
    exampleblock << "# NTMGRE 0,1 enable multiple gradients\n";
    exampleblock << "#    0: disable gradients (default)\n";
    exampleblock << "#    1: enable gradients\n";
    exampleblock << "# NTMGRP 0..3 print of curves\n";
    exampleblock << "#    0: don't print\n";
    exampleblock << "#    1: plot the curves\n";
    exampleblock << "#    2: print that values of the curves\n";
    exampleblock << "#    3: plot and print the curves\n";
    exampleblock << "# NTMGRN >= 0 number of gradients\n";
    exampleblock << "# MGRVAR: variable name to affect, available are:\n";
    exampleblock << "    TEMP0, CPOR, CDIR, RESO, CXR, COPR\n";
    exampleblock << "# MGRFRM: functional form of the curve\n";
    exampleblock << "#    0: linear interpolation between control points\n";
    exampleblock << "#    1: cubic spline interpolation between control points\n";
    exampleblock << "#    2: Bezier curve\n";
    exampleblock << "#    3: Oscillation: A sin[2Pi/T (t - dt)] + b\n";
    exampleblock << "#       Note: MGRNCP is 2. A = MGRCPT[1] T = MGRCPV[1] dt = MGRCPT[2] b = MGRCPV[2]\n";
    exampleblock << "# MGRNCP >= 2: number of control points\n";
    exampleblock << "# MGRCPT >= 0: time of the control point\n";
    exampleblock << "# MGRCPV: value of the control point\n";
    exampleblock << "#\n";
    exampleblock << "# NTMGRE NTMGRP\n";
    exampleblock << "       1      1\n";
    exampleblock << "# NTMGRN\n";
    exampleblock << "       2\n";
    exampleblock << "# MGRVAR MGRFRM MGRNCP\n";
    exampleblock << "  TEMP0[0]     0      2\n";
    exampleblock << "# MGRCPT MGRCPV\n";
    exampleblock << "  0.0    60.0\n";
    exampleblock << "  80.0   300.0\n";
    exampleblock << "# MGRVAR MGRFRM MGRNCP\n";
    exampleblock << "  CPOR        2      4\n";
    exampleblock << "# MGRCPT MGRCPV\n";
    exampleblock << "  0.0    2.5E5\n";
    exampleblock << "  0.0    2.5E1\n";
    exampleblock << " 20.0    0.0\n";
    exampleblock << " 80.0    0.0\n";
    exampleblock << "END\n";


    std::string blockname = "MULTIGRADIENT";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int enable = 0, plot = 0, num = 0;
        block.get_next_parameter("NTMGRE", enable, "", "0,1");
        block.get_next_parameter("NTMGRP", plot, "", "0,1,2,3");
        block.get_next_parameter("NTMGRN", num, ">=0", "");

        if (block.error()) {
          block.get_final_messages();
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
                break;
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
                break;
        }

        // read the gradient
        for(int i = 0; i < num; ++i) {
            int funct_form = 0, num_p = 0;
            std::string var;
            block.get_next_parameter("MGRVAR", var, "", "");
            block.get_next_parameter("MGRFRM", funct_form, "", "0,1,2,3");
            block.get_next_parameter("MGRNCP", num_p, ">=2", "");

            if (block.error()) {
              block.get_final_messages();
              return;
            }

            std::vector<std::pair<double, double> > points;
            for(int p = 0; p < num_p; ++p) {
                double t = 0.0, v = 0.0;
                block.get_next_parameter("MGRCPT", t, ">=0", "");
                block.get_next_parameter("MGRCPV", v, "", "");
                if (block.error()) return;
                points.push_back(std::pair<double, double>(t,v));
            }

            param.multigradient.variable.push_back(var);
            param.multigradient.functional_form.push_back(funct_form);
            param.multigradient.control_points.push_back(points);
        } // for gradients
        block.get_final_messages();
    }  // if block
} // MULTIGRADIENT

/**
 * @section addecouple ADDECOUPLE block
 * @snippet snippets/snippets.cc ADDECOUPLE
 */
void io::In_Parameter::read_ADDECOUPLE(simulation::Parameter & param,
                                       std::ostream & os) {
    DEBUG(8, "reading ADDECOUPLE");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "ADDECOUPLE\n";
    exampleblock << "# ADGR    >= 0 number of adiabatic decoupling groups\n";
    exampleblock << "# ADSTART first atom of the adiabatic decoupling group\n";
    exampleblock << "# ADEND   last atom of the adiabatic decoupling group\n";
    exampleblock << "# SM      scaling factor mass\n";
    exampleblock << "# SV      scaling factor potential energy function\n";
    exampleblock << "# ST      scaling factor temperature\n";
    exampleblock << "# TIR     which temperature bath to scale\n";
    exampleblock << "#  1      translational\n";
    exampleblock << "#  2      internal-rotatinal\n";
    exampleblock << "#  3      both\n";
    exampleblock << "# TMF     tau for calculating mean field\n";
    exampleblock << "# STAD    printing average to special trajectory\n";
    exampleblock << "# ADGR\n";
    exampleblock << "      2\n";
    exampleblock << "# ADSTART ADEND SM SV  ST TIR\n";
    exampleblock << "  1       1500 10   1  0  1\n";
    exampleblock << "  1501    3000  1  10  1  3\n";
    exampleblock << "# TMF STAD\n";
    exampleblock << "  0.1 1000\n";
    exampleblock << "END\n";

    std::string blockname = "ADDECOUPLE";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        unsigned int adstart = 0, eg = 0, tg = 0, adend = 0, tir = 0;
        double sm = 0.0, sv = 0.0, st = 0.0;

        block.get_next_parameter("ADGR", param.addecouple.adgr, "=>0", "");

        if (block.error()) {
          block.get_final_messages();
          return;
        }

        for (unsigned int i = 0; i < param.addecouple.adgr; ++i) {
            block.get_next_parameter("ADSTART", adstart, ">=0", "");
            std::string str_adstart  = io::to_string(adstart);
            block.get_next_parameter("ADEND", adend, ">="+str_adstart, "");
            block.get_next_parameter("SM", sm, ">=0", "");
            block.get_next_parameter("SV", sv, "", "");
            block.get_next_parameter("ST", st, ">=0", "");
            block.get_next_parameter("TIR", tir, "", "1,2,3");

            if (st != 1 && param.multibath.couple == false)
                io::messages.add("ADDECOUPLE block: ST>1, but no temperature scaling",
                                 "In_Parameter", io::message::error);

            //check whether the group is also a temperature group
            if (param.multibath.couple) {
                int addc_bath_index = 0;
                if (param.multibath.multibath.bath_index().size() < param.addecouple.adgr && st != 1)
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
                        io::messages.add("ADDECOUPLE block: ST bigger 1, but temperature group and adiabatic decouling group not equivalent",
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
                            int ir_bath = param.multibath.multibath.bath_index()[addc_bath_index].ir_bath;                             //scale temperature
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
            if (adstart == 1 && (param.force.energy_group[0] + 1) == adend) {
                check_group = 1;
                eg=0;
            }
            else
                for (unsigned int energy_i = 0; energy_i < param.force.energy_group.size() - 1; ++energy_i) {
                    if (param.force.energy_group[energy_i] + 2 == adstart
                        && param.force.energy_group[energy_i + 1] + 1 == adend) {
                        check_group = 1;
                        eg=i;
                    }
                }
            if (check_group == -1)
                io::messages.add("ADDECOUPLE block: energy and adiabatic groups are not identical", "In_Parameter", io::message::error);
            param.addecouple.add_adc(adstart - 1, adend - 1, sm, sv, st, tir, eg, tg);
        }
        if (param.addecouple.adgr > 0) {
            block.get_next_parameter("TMF", param.addecouple.tmf, ">=0", "");
            block.get_next_parameter("STAD", param.addecouple.write, ">=0", "");
        }

        block.get_final_messages();
    }
} // ADDECOUPLE

/**
 * @section qmmm QMMM block
 * @snippet snippets/snippets.cc QMMM

 */
void io::In_Parameter::read_QMMM(simulation::Parameter & param,
                                 std::ostream & os) {
    DEBUG(8, "reading QMMM");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "QMMM\n";
    exampleblock << "# NTQMMM -1..3 apply QM/MM\n";
    exampleblock << "#    0: do not apply QM/MM\n";
    exampleblock << "#   -1: apply mechanical embedding scheme with constant QM charges\n";
    exampleblock << "#    1: apply mechanical embedding scheme with dynamic QM charges\n";
    exampleblock << "#    2: apply electrostatic embedding scheme\n";
    exampleblock << "#    3: apply polarisable embedding scheme\n";
    exampleblock << "# NTQMSW 0..5 QM software package to use\n";
    exampleblock << "#    0: MNDO\n";
    exampleblock << "#    1: Turbomole\n";
    exampleblock << "#    2: DFTB\n";
    exampleblock << "#    3: MOPAC\n";
    exampleblock << "#    4: Gaussian\n";
    exampleblock << "#    5: Schnetpack NN\n";
    exampleblock << "#    6: Orca\n";
    exampleblock << "#    7: XTB\n";
    exampleblock << "# RCUTQM: ABS(RCUTQM): cutoff for inclusion of MM atoms in QM calculation\n";
    exampleblock << "#         (ignored for NTQMMM = 1)\n";
    exampleblock << "#     0.0: include all atoms\n";
    exampleblock << "#    >0.0: chargegroup-based cutoff\n";
    exampleblock << "#    <0.0: atom-based cutoff\n";
    exampleblock << "# NTWQMMM >= 0 write QM/MM related data to special trajectory\n";
    exampleblock << "#    0: do not write\n";
    exampleblock << "#   >0: write every NTWQMMMth step\n";
    exampleblock << "# QMLJ 0, 1 apply LJ between QM atoms \n";
    exampleblock << "#    0: do not apply\n";
    exampleblock << "#    1: apply\n";
    exampleblock << "# QMCON 0, 1 keep distance constraints in QM zone \n";
    exampleblock << "#    0: remove\n";
    exampleblock << "#    1: keep\n";
    exampleblock << "# MMSCALE scale mm-charges with (2/pi)*atan(x*(r_{qm}-r_{mm})) (optional) \n";
    exampleblock << "#     > 0.0: scaling-factor x\n";
    exampleblock << "#     < 0.0: don't scale (default)\n";
    exampleblock << "#\n";
    exampleblock << "# NTQMMM  NTQMSW  RCUTQM  NTWQMMM QMLJ QMCON MMSCALE\n";
    exampleblock << "       1       0     0.0        0    0     0    -1.0\n";
    exampleblock << "END\n";


    std::string blockname = "QMMM";
    Block block(blockname, exampleblock.str());
    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);


    int enable = 0,software = 0,write = 0,qmlj = 0,qmcon = 0;
    double mm_scale = -1.;
    double cutoff = 0.0;
    block.get_next_parameter("NTQMMM", enable, "", "-1,0,1,2,3");
    block.get_next_parameter("NTQMSW", software, "", "0,1,2,3,4,5,6,7");
    block.get_next_parameter("RCUTQM", cutoff, "", "");
    block.get_next_parameter("NTWQMMM", write, ">=0", "");
    block.get_next_parameter("QMLJ", qmlj, "", "0,1");
    block.get_next_parameter("QMCON", qmcon, "", "0,1");
    block.get_next_parameter("MMSCALE", mm_scale, ">=0.0 || <0.0", "", true);

    if (block.error()) {
        block.get_final_messages();
        return;
    }
    switch(enable) {
        case 0:
            param.qmmm.qmmm = simulation::qmmm_off;
            break;
        case -1:
            param.qmmm.qmmm = simulation::qmmm_mechanical;
            param.qmmm.qm_ch = simulation::qm_ch_constant;
            break;
        case 1:
            param.qmmm.qmmm = simulation::qmmm_mechanical;
            param.qmmm.qm_ch = simulation::qm_ch_dynamic;
            break;
        case 2:
            param.qmmm.qmmm = simulation::qmmm_electrostatic;
            break;
        case 3:
            param.qmmm.qmmm = simulation::qmmm_polarisable;
            break;
        default:
            break;
    }

    switch (software) {
        case 0:
            param.qmmm.software = simulation::qm_mndo;
            break;
        case 1:
            param.qmmm.software = simulation::qm_turbomole;
            break;
        case 2:
            param.qmmm.software = simulation::qm_dftb;
            break;
        case 3:
            param.qmmm.software = simulation::qm_mopac;
            break;
        case 4:
            param.qmmm.software = simulation::qm_gaussian;
            break;
        case 5:
#ifdef HAVE_PYBIND11
            param.qmmm.software = simulation::qm_nn;
#else
            io::messages.add("QMMM block: Schnetpack NN interface is not available "
                                "in your compilation. Use --enable-schnetpack for compiling.",
                                "In_Parameter", io::message::error);
#endif
            break;
        case 6:
            param.qmmm.software = simulation::qm_orca;
            break;
        case 7:
#ifdef WITH_XTB
            param.qmmm.software = simulation::qm_xtb;
#else       
            io::messages.add("QMMM block: XTB interface is not available "
                                "in your compilation. Use --with-xtb=PATH/TO/XTB for compiling.",
                                "In_Parameter", io::message::error);   
#endif
            break;
        default:
            break;
    }

    switch (qmlj) {
        case 0:
            param.qmmm.qm_lj = simulation::qm_lj_off;
            break;
        case 1:
            param.qmmm.qm_lj = simulation::qm_lj_on;
            break;
        default:
            break;
    }

    switch (qmcon) {
        case 0:
            param.qmmm.qm_constraint = simulation::qm_constr_off;
            break;
        case 1:
            param.qmmm.qm_constraint = simulation::qm_constr_on;
            break;
        default:
            break;
    }

    param.qmmm.mm_scale = mm_scale;
    if (cutoff < 0.0)
        param.qmmm.atomic_cutoff = true;
    param.qmmm.cutoff = fabs(cutoff);
    param.qmmm.write = write;
    if (param.qmmm.qmmm != simulation::qmmm_mechanical && param.qmmm.software == simulation::qm_nn)
        io::messages.add("QMMM block: Schnetpack NN works only with mechanical embedding scheme",
            "io::In_Parameter",
            io::message::error);
    if (param.qmmm.qmmm == simulation::qmmm_mechanical && param.qmmm.cutoff != 0.0)
        io::messages.add("QMMM block: RCUTQM > 0.0 has no effect for mechanical embedding scheme",
            "In_Parameter",
            io::message::warning);

    block.get_final_messages();
    }     // if block
} // QMMM

// two helper data types to simply unsupported block handling

enum unsupported_block_type {
    ub_unknown,     // I know it and know that I don't use it but I have no idea why
    ub_renamed,     // it was renamed. e.g. from previous versions
    ub_promd,     // it is a PROMD block. Tell alternative if there is any
    ub_g96     // it is a G96 block. Tell alternative if there is any
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
                default:         // don't know what to do.
                    msg << " is known to be not supported.";
            }

            io::messages.add(msg.str(), "In_Parameter", io::message::error);
        }
    }
}

/**
 * @section symres SYMRES block
 * @snippet snippets/snippets.cc SYMRES
 */
void io::In_Parameter::read_SYMRES(simulation::Parameter & param,
                                   std::ostream & os) {
    DEBUG(8, "reading SYMRES");

    std::stringstream exampleblock;
    // lines starting with 'exampleblock<<"' and ending with '\n";' (spaces don't matter)
    // will be used to generate snippets that can be included in the doxygen doc;
    // the first line is the tag
    exampleblock << "SYMRES\n";
    exampleblock << "# NTSYM 0..2 apply symmetry restraints\n";
    exampleblock << "#    0: do not apply symmetry restraints (default)\n";
    exampleblock << "#    1: apply symmetry restraints\n";
    exampleblock << "#    2: apply symmetry constraints\n";
    exampleblock << "# CSYM >= 0.0 force constants\n";
    exampleblock << "#\n";
    exampleblock << "# NTSYM     CSYM\n";
    exampleblock << "       1     0.0\n";
    exampleblock << "END\n";

    std::string blockname = "SYMRES";
    Block block(blockname, exampleblock.str());

    if (block.read_buffer(m_block[blockname], false) == 0) {
        block_read.insert(blockname);

        int enable = 0;
        block.get_next_parameter("NTSYM", enable, "", "0,1,2");
        block.get_next_parameter("CSYM", param.symrest.force_constant, ">=0", "");

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
                break;
        }
        block.get_final_messages();
    }     // if block
} // SYMRES

/**
 * @section amber AMBER block
 * @verbatim
AMBER
# AMBER 0..1 use AMBER topology
#    0: GROMOS topology
#    1: AMBER topology
# AMBSCAL >= 0.0 scaling factor for 1,4-electrostatic-interactions
#                (will be directly inverted when reading in)
#                Default: 1.2
#
# AMBER     AMBSCAL
      1     1.2
END

@endverbatim
 */
void io::In_Parameter::read_AMBER(simulation::Parameter & param,
        std::ostream & os) {
  DEBUG(8, "read AMBER");

  std::vector<std::string> buffer;
  std::string s;

  buffer = m_block["AMBER"];

  if (buffer.size()) {
    block_read.insert("AMBER");
    _lineStream.clear();
    _lineStream.str(concatenate(buffer.begin() + 1, buffer.end() - 1, s));

    int amber = 0;
    double factor = 0.0;
    _lineStream >> amber >> factor;

    if (_lineStream.fail()) {
      io::messages.add("Bad line in AMBER block.",
                    "In_Parameter", io::message::error);
      return;
    }

    if (amber != 1 && amber != 0) {
      io::messages.add("AMBER block: AMBER must be 0 or 1",
                    "In_Parameter", io::message::error);
    }
    param.amber.amber = (bool)amber;

    if (factor < 0.0) {
      io::messages.add("AMBER block: AMBSCAL must be >= 0.0",
                    "In_Parameter", io::message::error);
      return;
    } else if (factor > 0.0) {
      param.amber.coulomb_scaling = 1.0 / factor;
    } else {
      param.amber.coulomb_scaling = 1.0;
    }
  } // if block

} // AMBER
