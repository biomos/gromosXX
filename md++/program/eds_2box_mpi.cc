/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file eds_2box_mpi.cc
 * the main md program for 2-box eds simulations (MPI version)
 */

/**
 * @page programs Program Documentation
 *
 * @anchor eds_2box_mpi
 * @section eds_2box_mpi parallel molecular dynamics with MPI
 * @date 16.09.2011
 *
 * Program eds_2box_mpi is used to run parallel eds simulations with 2 boxes.
 * In order to use this program one has to compile GROMOS with MPI support:
 *
 * See @ref eds_2box for the documentation of the command line arguments
 *
*/

#ifdef XXMPI
#include <mpi.h>
#endif
#include <math/fft.h>
#include <stdheader.h>
#include <numeric>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <math/periodicity.h>
#include <algorithm/constraints/shake.h>
#include <algorithm/constraints/m_shake.h>
#include <algorithm/integration/monte_carlo.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>
//#include <signal.h>
#include <algorithm/integration/eds.h>

//bool exit_md;

//void signal_handler(int sig) {
//  std::cerr << "\nThe program will exit as soon as the MD step has finished. Press CTRL-C to force quit." << std::endl;
//  exit_md = true;
//  signal(SIGINT, SIG_DFL);
//}


int main(int argc, char *argv[]){
#ifdef XXMPI

    const double start = util::now();
    //  exit_md = false;
    //  signal(SIGINT, signal_handler);

    util::Known knowns;
    knowns << "topo1" << "topo2" << "conf1" << "conf2" << "input1" << "input2"
         << "verb" << "pttopo1" << "pttopo2" << "trc1" << "trc2" << "fin1"
         << "fin2" << "trv1" << "trv2" << "trf1" << "trf2" << "trs1" << "trs2"
         << "tre1" << "tre2" << "trg1" << "trg2" << "bae" << "bag"
         << "posresspec" << "refpos" << "distrest1" << "distrest2" << "dihrest"
         << "jval" << "xray" << "order" << "rdc" << "tfrdc" << "zanglerest"
         << "lud" << "led" << "anatrj" << "print" << "friction" << "qmmm"
         << "version";

    std::string usage;
    util::get_usage(knowns, usage, argv[0]);
    usage += "#\n\n";

    // master or slave : that's the question
    MPI_Init(&argc, &argv);
    FFTW3(mpi_init());

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create an output file (for the slaves)
    std::ostringstream oss;
    oss << "slave_" << rank << ".out";
    std::ofstream ofs(oss.str().c_str());

    bool quiet = false;
    std::ostream * os;
    if (rank == 0){
        os = &std::cout;
    }
    else{
        os = &ofs;
        quiet = true;
    }

    io::Argument args;

    if (args.parse(argc, argv, knowns)){
        if (rank == 0){
            std::cerr << usage << std::endl;
        }
        FFTW3(mpi_cleanup());
        MPI_Finalize();
        return 1;
    }

    util::print_title(true, *os);
    if (args.count("version") >= 0){
        MPI_Finalize();
        FFTW3(mpi_cleanup());
        return 0;
    }

    // parse the verbosity flag and set debug levels
    if (util::parse_verbosity(args)){
        if (rank == 0){
            std::cerr << "could not parse verbosity argument" << std::endl;
        }
        FFTW3(mpi_cleanup());
        MPI_Finalize();
        return 1;
    }

    // create the simulation classes
    topology::Topology topo1, topo2;
    configuration::Configuration conf1, conf2;
    simulation::Simulation sim1, sim2;
    algorithm::Algorithm_Sequence md1, md2;
    algorithm::Algorithm_Sequence md1_first(false), md2_first(false),
          md1_second(false), md2_second(false);

    // enable mpi for nonbonded terms
    sim1.mpi = true;
    sim2.mpi = true;

    //Build Attributes:
    ////GENERATE SIM SPECIFIC SIMULATION COMM
    MPI_Comm simulationCOMM;
    MPI_Comm_split(MPI_COMM_WORLD, sim1.mpiControl().mpiColor, rank, &simulationCOMM);
    MPI_Barrier(simulationCOMM);
    ////RANK and Size of thread in the simulationCOMM
    int simulation_rank, simulation_size;
    MPI_Comm_rank(simulationCOMM, &simulation_rank);
    MPI_Comm_size(simulationCOMM, &simulation_size);
    ////build up vector containing the replica_ids owned by thread
    std::vector<unsigned int> simulationOwnedThreads(size);
    std::iota(std::begin(simulationOwnedThreads), std::end(simulationOwnedThreads), 0);

    //Pass the information to MPI_control Struct
    sim1.mpiControl().numberOfThreads = simulation_size;
    sim1.mpiControl().threadID = simulation_rank;
    sim1.mpiControl().simulationOwnedThreads = simulationOwnedThreads;
    sim1.mpiControl().comm = simulationCOMM;

    //Control Output
    if (sim1.mpiControl().threadID == sim1.mpiControl().masterID) {
        std::cout << "MPI INFORMATION \n" <<
              "\tSimulation ID " << sim1.mpiControl().simulationID << "\n" <<
              "\tSimulation Master Thread ID " << sim1.mpiControl().masterID << "\n" <<
              "\tSimulation Total Num Threads " << sim1.mpiControl().numberOfThreads << "\n" <<
              "\tSimulation All Thread IDS:  " << "\n\t\t";
        for (auto x : sim1.mpiControl().simulationOwnedThreads) {
            std::cout << x << "\t ";
        }
        std::cout << "\n";
    }

    std::cout.flush();
    MPI_Barrier(sim1.mpiControl().comm);
    std::cout << "Simulation This Thread ID " << sim1.mpiControl().threadID << "\n";
    std::cout.flush();
    MPI_Barrier(sim1.mpiControl().comm);
    if (sim1.mpiControl().threadID == sim1.mpiControl().masterID) {
        std::cout << "END\n";
    }
    sim2.mpiControl() = sim1.mpiControl();


    io::Argument args1, args2;
    for(io::Argument::const_iterator it = args.begin(), to = args.end();
          it != to; ++it) {
        const std::string & argname = it->first;
        if (argname[argname.size() - 1] == '1') {
            std::string newargname(argname.substr(0, argname.size() - 1));
            args1.insert(std::pair<std::string, std::string>(newargname, it->second));
        }
    }

    for(io::Argument::const_iterator it = args.begin(), to = args.end();
          it != to; ++it) {
        const std::string & argname = it->first;
        if (argname[argname.size() - 1] == '2') {
            std::string newargname(argname.substr(0, argname.size() - 1));
            args2.insert(std::pair<std::string, std::string>(newargname, it->second));
        }
    }

    if (io::read_input(args1, topo1, conf1, sim1,  md1, *os, quiet)){
        io::messages.display(std::cout);
        std::cout << "\nErrors during initialization of box 1!\n" << std::endl;
        FFTW3(mpi_cleanup());
        MPI_Finalize();
        return 1;
    }

    if (io::read_input(args2, topo2, conf2, sim2,  md2, *os, quiet)){
        io::messages.display(std::cout);
        std::cout << "\nErrors during initialization of box 2!\n" << std::endl;
        FFTW3(mpi_cleanup());
        MPI_Finalize();
        return 1;
    }

    algorithm::EDS * eds1 = (algorithm::EDS *) md1.algorithm("EDS");
    if (eds1 != NULL) {
        eds1->set_conf2(conf2);
    }
    algorithm::EDS * eds2 = (algorithm::EDS *) md2.algorithm("EDS");
    if (eds2 != NULL) {
        eds2->set_conf2(conf1);
    }

    // initialises all algorithms (and therefore also the forcefield)
    md1.init(topo1, conf1, sim1, *os, quiet);
    md2.init(topo2, conf2, sim2, *os, quiet);

    bool second = false;
    for(algorithm::Algorithm_Sequence::const_iterator it = md1.begin(),
          to = md1.end(); it != to; ++it) {
        if (*it == eds1) {
          second = true;
          continue;
        }

        if (second)
          md1_second.push_back(*it);
        else
          md1_first.push_back(*it);
        }
        second = false;
        for(algorithm::Algorithm_Sequence::const_iterator it = md2.begin(),
              to = md2.end(); it != to; ++it) {
            if (*it == eds2) {
                second = true;
                continue;
            }

            if (second){
                md2_second.push_back(*it);
            }
            else{
                md2_first.push_back(*it);
            }
            }

    //For Debugging the Algorithm CHAIN
    //std::cout << "\nAlgorithm Chain Splitting \n";f
    //std::cout << "MD1 - first: \n";
    //md1_first.printSequence();
    //std::cout << "MD1 - second: \n";
    //md1_second.printSequence();
    //std::cout << "\n";

    //std::cout << "MD2 - first: \n";
    //md2_first.printSequence();
    //std::cout << "MD2 - second: \n";
    //md2_second.printSequence();
    //std::cout << "\n\n";


    int error1, error2;

    if(rank == 0){


        io::Out_Configuration traj1(GROMOSXX "\n"), traj2(GROMOSXX "\n");

        traj1.title(GROMOSXX "\n" + sim1.param().title);
        traj2.title(GROMOSXX "\n" + sim2.param().title);

        // create output files...
        traj1.init(args1, sim1.param());
        traj2.init(args2, sim2.param());

        std::cout << "\nMESSAGES FROM INITIALISATION\n";
        if (io::messages.display(std::cout) >= io::message::error){
            // exit
            std::cout << "\nErrors during initialisation!\n" << std::endl;
            MPI_Finalize();
            return 1;
        }

        io::messages.clear();

        std::cout.precision(5);
        std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

        std::cout << "\nenter the next level of molecular "
              << "dynamics simulations\n" << std::endl;

        MPI_Barrier(sim1.mpiControl().comm);

        int percent = 0;
        std::cout << "==================================================\n"
              << " MAIN MD LOOP\n"
              << "==================================================\n"
              << std::endl;

        std::cout << "MPI master node(of " << size << "nodes)" << std::endl;

        const double init_time = util::now() - start;

        //int error1, error2;

        int next_step = 1;

        while(int(sim1.steps()) < sim1.param().step.number_of_steps &&
              int(sim2.steps()) < sim2.param().step.number_of_steps ){

            traj1.write(conf1, topo1, sim1, io::reduced);
            traj2.write(conf2, topo2, sim2, io::reduced);

            // run a step
            //std::cout << "MASTER: STEP - "<< int(sim1.steps()) << " \n";
            error1 = md1_first.run(topo1, conf1, sim1);
            //std::cout << "\tMASTER: first md1Eval \n";
            error2 = md2_first.run(topo2, conf2, sim2);
            //std::cout << "\tMASTER: first md2Eval \n";
            error1 += eds1->apply(topo1, conf1, sim1);
            //std::cout << "\tMASTER: first eds1Eval \n";
            error2 += eds2->apply(topo2, conf2, sim2);
            //std::cout << "\tMASTER: first eds2Eval \n";
            error1 += md1_second.run(topo1, conf1, sim1);
            //std::cout << "\tMASTER: first md1REST \n";
            error2 += md2_second.run(topo2, conf2, sim2);
            //std::cout << "\tMASTER: first md2REST \n";

            if ((error1 || error2)){
                if (error1 == E_MINIMUM_REACHED){
                    conf1.old().energies.calculate_totals();
                    traj1.print_timestep(sim1, traj1.output());
                    io::print_ENERGY(traj1.output(), conf1.old().energies,
                         topo1.energy_groups(),
                         "MINIMUM ENERGY", "EMIN_");

                    error1 = 0; // clear error condition
                    break;
                }
                else {
                    // try to print energies anyway
                    // if (error == E_NAN){
                    io::print_ENERGY(traj1.output(), conf1.current().energies,
                         topo1.energy_groups(),
                         "OLDERROR", "OLDERR_");

                      io::print_ENERGY(traj2.output(), conf2.current().energies,
                         topo2.energy_groups(),
                         "OLDERROR", "OLDERR_");

                    io::print_ENERGY(traj1.output(), conf1.old().energies,
                         topo1.energy_groups(),
                         "ERROR", "ERR_");

                    io::print_ENERGY(traj2.output(), conf2.old().energies,
                         topo2.energy_groups(),
                         "ERROR", "ERR_");

                }

                std::cout << "\nError during MD run!\n" << std::endl;
                // send error status to slaves
                next_step = 0;
                std::cout << "Telling slaves to quit." << std::endl;
                MPI_Bcast(&next_step, 1, MPI_INT, sim1.mpiControl().masterID, sim1.mpiControl().comm);
                MPI_Bcast(&next_step, 1, MPI_INT, sim2.mpiControl().masterID, sim2.mpiControl().comm);
                //std::cout << "\tat step " << sim1.steps() << " (time " << sim1.time() << ")\n" << std::endl;
                // try to save the final structures...
                break;
            }

            //if (exit_md)
            //std::cout << "\nMD run terminated by SIGINT (CTRL-C)." << std::endl;

            // tell the slaves to continue
            MPI_Bcast(&next_step, 1, MPI_INT, sim1.mpiControl().masterID, sim1.mpiControl().comm);
            MPI_Bcast(&next_step, 1, MPI_INT, sim2.mpiControl().masterID, sim2.mpiControl().comm);
            traj1.print(topo1, conf1, sim1);
            traj2.print(topo2, conf2, sim2);

            ++sim1.steps();
            ++sim2.steps();
            sim1.time() = sim1.param().step.t0 + sim1.steps()*sim1.time_step_size();
            sim2.time() = sim2.param().step.t0 + sim2.steps()*sim2.time_step_size();

            if ((sim1.param().step.number_of_steps / 10 > 0) &&
            (sim1.steps() % (sim1.param().step.number_of_steps / 10) == 0)){
            ++percent;
            const double spent = util::now() - start;
            const int hh = int(spent / 3600);
            const int mm = int((spent - hh * 3600) / 60);
            const int ss = int(spent - hh * 3600 - mm * 60);

            std::cerr << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
            std::cout << "MD:       " << std::setw(4) << percent * 10 << "% done..." << std::endl;
            std::cerr << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;
            std::cout << "MD: spent " << hh << ":" << mm << ":" << ss << std::endl;

            const double eta_spent = spent / sim1.steps() * sim1.param().step.number_of_steps - spent;
            const int eta_hh = int(eta_spent / 3600);
            const int eta_mm = int((eta_spent - eta_hh * 3600) / 60);
            const int eta_ss = int(eta_spent - eta_hh * 3600 - eta_mm * 60);

            std::cerr << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
            std::cout << "MD: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
            }
        } // main md loop

        std::cout << "writing final configuration" << std::endl;

        traj1.write(conf1, topo1, sim1, io::final);
        traj1.print_final(topo1, conf1, sim1);

        traj2.write(conf2, topo2, sim2, io::final);
        traj2.print_final(topo2, conf2, sim2);

        const double sim_time = (util::now() - start) - init_time;

        std::cout << "\nMESSAGES FROM SIMULATION\n";
        io::message::severity_enum err_msg = io::messages.display(std::cout);

        std::cout << "\n\n";

        md1.print_timing(std::cout);

        std::setprecision(5);
        std::cout << std::endl;    
        std::cout << "Wall time initialisation (s):  " << std::setw(10) << init_time << std::endl;
        std::cout << "Wall time simulation (s):      " << std::setw(10) << sim_time << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Wall time total (s):           " << std::setw(10) << (init_time + sim_time) << std::endl;
        std::cout << std::endl;
        
        const double ns_calculated = (sim1.param().step.number_of_steps * sim1.param().step.dt)/1000.0;
        //std::cout << "Simulated period (ns):         " << std::setw(10) << ns_calculated << std::endl; 
        
        const double performance_per_day = ns_calculated / (sim_time/(60.0*60.0*24.0));    
        std::cout << "Performance (ns/day):          " << std::setw(10) << performance_per_day << std::endl; 
        std::cout << std::endl;
        std::cout << std::endl;

        const time_t time_now = time_t(util::now());
        std::cout << ctime(&time_now) << "\n\n";

        if (error1 || error2)
          std::cout << "\nErrors encountered during run - check above!\n" << std::endl;
        else if(err_msg > io::message::notice){
          std::cout << "\n" GROMOSXX " finished."
                    << "Check the messages for possible problems during the run."
                    << std::endl;
        }
        else{
          std::cout << "\n" GROMOSXX " finished successfully\n" << std::endl;
        }
    }

  ////////////////////////////////////////////////////////////////////////////////
  // MPI Slave
  ////////////////////////////////////////////////////////////////////////////////

    else{

        (*os) << "MPI slave " << rank << " of " << size-1 << std::endl;

        // check whether we have nonbonded interactions
        bool do_nonbonded1 = (sim1.param().force.nonbonded_crf ||
                              sim1.param().force.nonbonded_vdw);

        bool do_nonbonded2 = (sim2.param().force.nonbonded_crf ||
                              sim2.param().force.nonbonded_vdw);

        // let's get the forcefield
        interaction::Forcefield * ff1 =
          dynamic_cast<interaction::Forcefield *>(md1.algorithm("Forcefield"));

        interaction::Forcefield * ff2 =
          dynamic_cast<interaction::Forcefield *>(md2.algorithm("Forcefield"));

        if (do_nonbonded1 && ff1 == NULL){
          std::cerr << "MPI slave: could not access forcefield 1\n\t(internal error)" << std::endl;
          MPI_Finalize();
          return 1;
        }

        if (do_nonbonded2 && ff2 == NULL){
          std::cerr << "MPI slave: could not access forcefield 2\n\t(internal error)" << std::endl;
          MPI_Finalize();
          return 1;
        }

        interaction::Interaction * nb1 = ff1->interaction("NonBonded");
        if (do_nonbonded1 && nb1 == NULL){
          std::cerr << "MPI slave: could not get NonBonded interactions from forcefield 1"
            << "\n\t(internal error)"
            << std::endl;
          MPI_Finalize();
          return 1;
        }

        interaction::Interaction * nb2 = ff2->interaction("NonBonded");
        if (do_nonbonded2 && nb2 == NULL){
          std::cerr << "MPI slave: could not get NonBonded interactions from forcefield 2"
            << "\n\t(internal error)"
            << std::endl;
          MPI_Finalize();
          return 1;
        }

        // get shake and check whether we do it for solvent
        bool do_shake1 = sim1.param().constraint.solvent.algorithm == simulation::constr_shake; //sim1.param().system.nsm && <-- doesn't make sense here, constraints are on whole system @bschroed

        bool do_shake2 = sim2.param().constraint.solvent.algorithm == simulation::constr_shake; //sim2.param().system.nsm && <-- doesn't make sense here, constraints are on whole system @bschroed


        algorithm::Shake * shake1 =
          dynamic_cast<algorithm::Shake *>(md1.algorithm("Shake"));
        if (do_shake1 && shake1 == NULL) {
            std::cerr << "MPI slave: could not get Shake algorithm from MD1 sequence."
                    << "\n\t(internal error)"
                    << std::endl;
          MPI_Finalize();
          return 1;
        }

        algorithm::Shake * shake2 =
                dynamic_cast<algorithm::Shake *>(md2.algorithm("Shake"));
        if (do_shake2 && shake2 == NULL) {
            std::cerr << "MPI slave: could not get Shake algorithm from MD2 sequence."
                    << "\n\t(internal error)"
                    << std::endl;
          MPI_Finalize();
          return 1;
        }

        // get m_shake and check whether we do it for solvent
        bool do_m_shake1 = sim1.param().system.nsm &&
                sim1.param().constraint.solvent.algorithm == simulation::constr_m_shake;

        bool do_m_shake2 = sim2.param().system.nsm &&
                sim2.param().constraint.solvent.algorithm == simulation::constr_m_shake;

        //std::cout <<" SLave SHAKER: "<< do_m_shake1 << "  "<< do_m_shake2 << "\n";

        algorithm::M_Shake * m_shake1 =
          dynamic_cast<algorithm::M_Shake *>(md1.algorithm("M_Shake"));
        if (do_m_shake1 && m_shake1 == NULL) {
            std::cerr << "MPI slave: could not get M_Shake algorithm from MD1 sequence."
                    << "\n\t(internal error)"
                    << std::endl;
          MPI_Finalize();
          return 1;
        }

        algorithm::M_Shake * m_shake2 =
          dynamic_cast<algorithm::M_Shake *>(md2.algorithm("M_Shake"));
        if (do_m_shake2 && m_shake2 == NULL) {
            std::cerr << "MPI slave: could not get M_Shake algorithm from MD2 sequence."
                    << "\n\t(internal error)"
                    << std::endl;
          MPI_Finalize();
          return 1;
        }

        // get chemical monte carlo
        bool do_cmc1 = sim1.param().montecarlo.mc;
        bool do_cmc2 = sim2.param().montecarlo.mc;

        algorithm::Monte_Carlo * monte_carlo1 =
                dynamic_cast<algorithm::Monte_Carlo *>(md1.algorithm("MonteCarlo"));
        if(do_cmc1 && monte_carlo1 == NULL){
          std::cerr << "MPI slave: could not get Monte Carlo algorithm from MD1 sequence."
                  << "\n\t(internal error)"
                  << std:: endl;
          MPI_Finalize();
          return 1;
        }

        algorithm::Monte_Carlo * monte_carlo2 =
                dynamic_cast<algorithm::Monte_Carlo *>(md2.algorithm("MonteCarlo"));
        if(do_cmc2 && monte_carlo2 == NULL){
          std::cerr << "MPI slave: could not get Monte Carlo algorithm from MD2 sequence."
                  << "\n\t(internal error)"
                  << std:: endl;
          MPI_Finalize();
          return 1;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////
        // run the simulation
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        MPI_Barrier(sim1.mpiControl().comm);

        const double init_time = util::now() - start;
        int next_step = 0 ;

        while(int(sim1.steps()) < sim1.param().step.number_of_steps &&
              int(sim2.steps()) < sim2.param().step.number_of_steps){
          // run a step
          // (*os) << "waiting for master (nonbonded interaction)" << std::endl;
          // DEBUG(10, "slave " << rank << " waiting for master");
          if (do_nonbonded1 && (error1 = nb1->calculate_interactions(topo1, conf1, sim1)) != 0){
        std::cout << "MPI slave " << rank << ": error in nonbonded calculation box 1!\n" << std::endl;
          }
          //std::cout << "\t SLAVE: first md1Eval \n";

          if (do_nonbonded2 && (error2 = nb2->calculate_interactions(topo2, conf2, sim2)) != 0){
        std::cout << "MPI slave " << rank << ": error in nonbonded calculation box 2!\n" << std::endl;
          }
          //std::cout << "\t SLAVE: first md2Eval \n";

          // DEBUG(10, "slave " << rank << " step done");
          // (*os) << "step done (it really worked?)" << std::endl;
          if (do_cmc1 && (error1 = monte_carlo1->apply(topo1, conf1, sim1)) != 0){
            std::cout << "MPI slave " << rank << ": error in Monte Carlo algorithm box 1!\n"
                    << std::endl;
          }

          if (do_cmc2 && (error2 = monte_carlo2->apply(topo2, conf2, sim2)) != 0){
            std::cout << "MPI slave " << rank << ": error in Monte Carlo algorithm box 2!\n"
                    << std::endl;
          }

          if (do_shake1 && (error1 = shake1->apply(topo1, conf1, sim1)) != 0) {
            std::cout << "MPI slave " << rank << ": error in Shake algorithm box 1!\n" << std::endl;
          }
          //std::cout << "\t SLAVE: first shake1Eval \n";

          if (do_shake2 && (error2 = shake2->apply(topo2, conf2, sim2)) != 0) {
            std::cout << "MPI slave " << rank << ": error in Shake algorithm box 2!\n" << std::endl;
          }
          //std::cout << "\t SLAVE: first shake2Eval \n";

          //std::cout << "\t SLAVE: BC Exchange - START\n";

          MPI_Bcast(&next_step, 1, MPI_INT, sim1.mpiControl().masterID, sim1.mpiControl().comm);
          MPI_Bcast(&next_step, 1, MPI_INT, sim2.mpiControl().masterID, sim2.mpiControl().comm);

          //std::cout << "\t SLAVE: BC Exchange - DONE\n";

          if (!next_step) {
            (*os) << "There was an error in the master. Check output file for details." << std::endl
                  << "Exiting from MD main loop." << std::endl;
            error1 = 1;
            break;
          }

          ++sim1.steps();
          ++sim2.steps();
          sim1.time() = sim1.param().step.t0 + sim1.steps()*sim1.time_step_size();
          sim2.time() = sim2.param().step.t0 + sim2.steps()*sim2.time_step_size();
          //std::cout << "\t SLAVE: STEPDone\n";
        }

        const double sim_time = (util::now() - start) - init_time;

        (*os) << "\nMESSAGES FROM SIMULATION\n";
        io::message::severity_enum err_msg = io::messages.display(*os);

        (*os) << "\n\n";

        md1.print_timing(*os);
        md2.print_timing(*os);

        std::setprecision(5);
        (*os) << std::endl;    
        (*os) << "Wall time initialisation (s):  " << std::setw(10) << init_time << std::endl;
        (*os) << "Wall time simulation (s):      " << std::setw(10) << sim_time << std::endl;
        (*os) << "-----------------------------------------" << std::endl;
        (*os) << "Wall time total (s):           " << std::setw(10) << (init_time + sim_time) << std::endl;
        (*os) << std::endl;
        
        const double ns_calculated = (sim1.param().step.number_of_steps * sim1.param().step.dt)/1000.0;
        //(*os) << "Simulated period (ns):         " << std::setw(10) << ns_calculated << std::endl; 
        
        const double performance_per_day = ns_calculated / (sim_time/(60.0*60.0*24.0));    
        (*os) << "Performance (ns/day):          " << std::setw(10) << performance_per_day << std::endl; 
        (*os) << std::endl;
        (*os) << std::endl;

        const time_t time_now = time_t(util::now());
        (*os) << ctime(&time_now) << "\n\n";

        if (error1 || error2)
          (*os) << "\nErrors encountered during run - check above!\n" << std::endl;
        else if(err_msg > io::message::notice){
          (*os) << "\n" GROMOSXX " MPI slave " << rank << " finished. "
            << "Check the messages for possible problems during the run."
            << std::endl;
        }
        else{
          (*os) << "\n" GROMOSXX " MPI slave " << rank << " finished successfully\n" << std::endl;
        }

    } // end of slave

  os = NULL;

  ofs.flush();
  ofs.close();

  // and exit...
  FFTW3(mpi_cleanup());
  MPI_Finalize();
  return error1;

#else
  std::cout << argv[0] << " needs MPI to run\n\tuse --enable-mpi at configure and appropriate compilers\n" << std::endl;
  return 1;
#endif

}
