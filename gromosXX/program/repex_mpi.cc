/**
 * @file repex_mpi.cc
 * the main md program for replica exchange simulations using the message passing interface
 */
/**
 * @page programs Program Documentation
 *
 * @anchor repex_mpi
 * @section repex_mpi replica exchange
 * @date 20.08.2019
 *
 * Program repex_mpi is used to run replica exchange simulations.
 *
 * See @ref md for the documentation of the command line arguments.
 * Addition command line arguments are:
 * <table border=0 cellpadding=0>
 * <tr><td> \@repdat</td><td>&lt;name of the replica exchange data file&gt; </td><td style="color:#088A08">in</td></tr>
 * <tr><td> \@repout</td><td>&lt;name of the replica exchange output files&gt; </td><td style="color:#088A08">in</td></tr>
 * </table>
 */
#include <stdheader.h>
#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/parameter/check_parameter.h>
#include <simulation/parameter.h>

#include <io/print_block.h>

#include <thread>
#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>
#include <math/gmath.h>

#ifdef XXMPI
    #include <mpi.h>
#endif

#include <util/replicaExchange/replica_mpi_tools.h>
#include <util/replicaExchange/replica/replica.h>

#include <util/replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_master.h>
#include <util/replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_slave.h>
#include <util/replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_master_eds.h>
#include <util/replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_slave_eds.h>
#include <util/replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_master_2d_s_eoff_eds.h>
#include <util/replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_slave_2d_s_eoff_eds.h>

#include <util/replicaExchange/repex_mpi.h>
#include <util/replicaExchange/replica_mpi_tools.h>

#include "util/replicaExchange/replica_graph_control.h"

//Debug Instructions
#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

int main(int argc, char *argv[]) {
#ifdef XXMPI
    //initializing MPI
    MPI_Init(&argc, &argv);
    const double start = util::now();

    int ttotalNumberOfThreads;
    int tglobalThreadID;

    MPI_Comm_size(MPI_COMM_WORLD, &ttotalNumberOfThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &tglobalThreadID);

    unsigned int totalNumberOfThreads = ttotalNumberOfThreads;
    unsigned int globalThreadID = tglobalThreadID;

    if (globalThreadID == 0) {
        std::string msg("\n==================================================\n\tGROMOS Replica Exchange:\n==================================================\n");
        std::cout << msg;
        std::cerr << msg;
    }


    // reading arguments
    util::Known knowns;
    knowns << "topo" << "conf" << "input" << "verb" << "pttopo"
            << "trc" << "fin" << "trv" << "trf" << "trs" << "tre" << "trg"
            << "bae" << "bag" << "posresspec" << "refpos" << "distrest" << "dihrest"
            << "jval" << "rdc" << "xray" << "lud" << "led" << "print" << "friction" // << "anatrj"
            << "version" << "repdat" << "repout";

    std::string usage;
    util::get_usage(knowns, usage, argv[0]);
    usage += "#\n\n";


    MPI_DEBUG(1, "RANK: "<<globalThreadID<<" Parse ARGs \n");

    // Parse command line arguments
    io::Argument args;
    try {
        if (args.parse(argc, argv, knowns)) {
            std::cerr << usage << std::endl;
            MPI_Finalize();
            return 1;
        }

        if (args.count("version") >= 0) {
            MPI_Finalize();
            return 0;
        }

        // parse the verbosity flag and set debug levels
        if (util::parse_verbosity(args)) {
            std::cerr << "could not parse verbosity argument" << std::endl;
            MPI_Finalize();
            return 1;
        }

    } catch (const std::exception &e) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\t Exception was thrown in Commandline Parsing!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an Exception in initial Args parse!\n\n";
        std::cout << msg << e.what() << std::endl;
        std::cerr << msg << e.what() << std::endl;
        MPI_Finalize();
        MPI_Abort(MPI_COMM_WORLD, E_USAGE);
        return -1;
    } catch (...) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\t Exception was thrown in Commandline Parsing!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial Args parse!!\n Hit the developers!\n\n";
        std::cout << msg << std::endl;
        std::cerr << msg << std::endl;
        MPI_Finalize();
        return -1;
    }


    MPI_DEBUG(1, "RANK: "<<globalThreadID<<" totalRankNum: "<< totalNumberOfThreads<<": hello world!\n");


    /**
     * GLOBAL SETTING VARIABLES
     */


    int reedsSim;
    unsigned int equil_runs;
    unsigned int sim_runs;

    unsigned int numReplicas;
    unsigned int numEDSstates;

    unsigned int numAtoms;
    unsigned int numEoff;
    unsigned int numSVals;
    unsigned int cont;
    bool periodic;

    //simulation dependend MPI Vars
    unsigned  simulationID; //MPI thread belongs to replica Simulation:
    unsigned int threadsPerReplicaSimulation; //remove this param - also from exchangers
    unsigned int simulationTrheadID;

    std::map<unsigned int, unsigned int> thread_id_replica_map; // where is which replica
    std::vector<std::vector<unsigned int> > replica_owned_threads; // set IDs for each replica


    try {
            simulation::Simulation sim;
            topology::Topology topo;
            configuration::Configuration conf;
            algorithm::Algorithm_Sequence md;

            // read in parameters
            bool quiet = true;
            if (globalThreadID == 0) {
                quiet = true;
            }

            // read in (test and get initial inf)
            if (io::read_input_repex(args, topo, conf, sim, md,  globalThreadID, totalNumberOfThreads, std::cout, quiet)) {
                if (globalThreadID == 0) {
                    std::cerr << "\n\t########################################################\n"
                            << "\n\t\tErrors during initial Parameter reading!\n"
                            << "\n\t########################################################\n";
                    io::messages.display(std::cout);
                    io::messages.display(std::cerr);
                    MPI_Abort(MPI_COMM_WORLD, E_USAGE);
                }
                MPI_Finalize();
                return 1;
            }

            MPI_DEBUG(2, "RANK: "<<globalThreadID<<" done with repex_in\n");

            //set global parameters
            cont = sim.param().replica.cont;
            equil_runs = sim.param().replica.equilibrate;
            sim_runs = sim.param().replica.trials;
            numAtoms = topo.num_atoms();
            reedsSim = sim.param().reeds.reeds;

            /**
             * HERE STARTS taking mapping of available threads to the replicas
             */

            numEoff = sim.param().reeds.num_eoff;
            numSVals = sim.param().reeds.num_l;
            periodic = sim.param().reeds.periodic;

            switch(reedsSim) {
                  case 0:
                      numReplicas = sim.param().replica.num_T * sim.param().replica.num_l;
                      numEDSstates = 0;
                      DEBUG(3, "numReps & numEDSstates: " << numReplicas << ", " << numEDSstates << "\n");
                      break;
                  case 1:
                      numReplicas = numSVals;
                      numEDSstates = sim.param().reeds.eds_para[0].numstates;
                      DEBUG(3, "numReps & numEDSstates: " << numReplicas << ", " << numEDSstates << "\n");
                      break;
                  case 2:
                      numReplicas = numSVals * numEoff;
                      numEDSstates = sim.param().reeds.eds_para[0].numstates;
                      DEBUG(3, "numReps & numEDSstates: " << numReplicas << ", " << numEDSstates << "\n");
                      break;
              }


            //MPI THREAD SIMULATION SPLITTING
            replica_owned_threads.resize(numReplicas);
            // counts through every replica and assigns it to respective node
            // starts at beginning if numReplicas > size
            // could be optimized by putting neighboring replicas on same node; less communication...
            //MPI THREAD SPLITING ONTO Simulation - REPLICAS
            threadsPerReplicaSimulation = totalNumberOfThreads / numReplicas;
            unsigned int leftOverThreads = totalNumberOfThreads % numReplicas;

            DEBUG(3, "threadsperSim & leftOverThreads: " << threadsPerReplicaSimulation << ", " << leftOverThreads << "\n");

            unsigned int threadID =0;
            int replica_offset = 0;
            for (unsigned int replicaSimulationID = 0; replicaSimulationID < numReplicas; replicaSimulationID++) {
                for (unsigned int replicaSimulationSubThread = 0; replicaSimulationSubThread < threadsPerReplicaSimulation; replicaSimulationSubThread++) {
                    threadID = replicaSimulationSubThread + replicaSimulationID*threadsPerReplicaSimulation+replica_offset;
                    thread_id_replica_map.insert(std::pair<unsigned int, unsigned int>(threadID, replicaSimulationID));
                    replica_owned_threads[replicaSimulationID].push_back(threadID);

                }
                if(leftOverThreads>0 and replicaSimulationID < leftOverThreads){    //left over threads are evenly distirbuted.
                    threadID = threadsPerReplicaSimulation+replicaSimulationID*threadsPerReplicaSimulation+replica_offset;
                    thread_id_replica_map.insert(std::pair<unsigned int, unsigned int>(threadID, replicaSimulationID));
                    replica_owned_threads[replicaSimulationID].push_back(threadID);
                    replica_offset++;
                }
            }

            simulationID = thread_id_replica_map[globalThreadID];

            int counter = 0;
            for(unsigned int x : replica_owned_threads[simulationID]){
                if(x==globalThreadID){
                    simulationTrheadID = counter;
                    break;
                }
                counter++;
            }



    } catch (const std::exception &e) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\t Exception was thrown in initial test parsing!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an Exception in initial test Parsing of the Gromos Parameters!\n\nMessage:\t";
        std::cout << msg << e.what() << std::endl;
        std::cerr << msg << e.what() << std::endl;
        MPI_Finalize();
        return -1;
    } catch (...) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\t Exception was thrown!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial test Parsing of the  Gromos Parameters!\n Hit the developers!\n\nMessage:\t";
        std::cout << msg << std::endl;
        std::cerr << msg << std::endl;
        MPI_Finalize();
        return -1;
    }

    MPI_DEBUG(1, "RANK: "<<globalThreadID<<" Done with init parse\n");    //TODO: lower verb level
    io::messages.clear();
    MPI_DEBUG(1, "REPLICA_ID \t " << globalThreadID << "\t Simulation_ID\t"<< simulationID << "\t SimulationThread\t"<<simulationTrheadID<<"\n")


    //////////////////////////////////////
    //MPI SETTINGS!
    //////////////////////////////////////

    //simulation parallelization
    if(totalNumberOfThreads%numReplicas != 0){
        std::cerr << "\n\t########################################################\n"
                << "\n\t\t ERROR please make sure giving the same number of threads as replicas or a multiple of these!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR please make sure giving the same number of threads as replicas or a multiple of these !\n";
        std::cout << msg  << "( "<<numReplicas<<" replica * x)\n\n" << std::endl;
        std::cerr << msg  << "( "<<numReplicas<<" replica * x)\n\n" << std::endl;
        MPI_Finalize();
        MPI_Abort(MPI_COMM_WORLD, E_USAGE);
        return -1;    }
    ////GENERATE SIM SPECIFIC SIMULATION COMM
    MPI_Comm simulationCOMM = MPI_Comm();    //this is used for different replica simulation parallelisations
    MPI_Comm_split(MPI_COMM_WORLD, simulationID, globalThreadID, &simulationCOMM);

    MPI_DEBUG(1, "REPLICA_ID \t " << globalThreadID << "\t Simulation_ID\t"<< simulationID << "\t SIMCOMM ESTABLISHED\n");
    //get comm specifics
    int simulation_rank,simulation_size;
    MPI_Comm_rank(simulationCOMM, &simulation_rank);
    MPI_Comm_size(simulationCOMM, &simulation_size);
    simulation::MpiControl replica_mpi_control = simulation::MpiControl(true, simulationID,
            replica_owned_threads[simulationID].size(), // thread per Simulation
            0, //Master thread
            simulation_rank,     //this thread ID in the simulation  thread_id_replica_map
            simulationID,           // MPI-Color
            replica_owned_threads[simulationID], //Threads belonging to the simulatoin
            simulationCOMM );   //simComm

    MPI_Barrier(replica_mpi_control.comm);    //wait for all threads to register!
    // Simulation Control DONE

    //////////////////////////////////////
    //for RE-Graph
    //////////////////////////////////////
    // communication structure build up for Replica Exchange
    MPI_Comm  replicaGraphCOMM = MPI_Comm();    //this is used for different replica GRAPH parallelisations
    int graphMPIColor = 9999;
    int graphThreadID = -1;
    int graphSize = -1;

    if(replica_mpi_control.masterID == replica_mpi_control.threadID){
        MPI_Comm_split(MPI_COMM_WORLD, graphMPIColor, globalThreadID, &replicaGraphCOMM);
        MPI_Barrier(replicaGraphCOMM);    //wait for all threads to register!
        MPI_Comm_rank(replicaGraphCOMM, &graphThreadID);
        MPI_Comm_size(replicaGraphCOMM, &graphSize);
        std::cout << graphThreadID << "\n";
    }
    else{
        MPI_Comm_split(MPI_COMM_WORLD,graphMPIColor+3, globalThreadID, &replicaGraphCOMM);
    }



    util::replica_graph_control reGMPI(0,   //GraphID
                                       0,   // MasterID
                                graphThreadID,  // This Thread
                                totalNumberOfThreads,   //total number of threads
                                replica_owned_threads.size(),  //replicas in the graph
                                graphMPIColor,   //MPI-color
                                replicaGraphCOMM);  //Com of this graph

    //
    std::vector<unsigned int> replica_master_IDS = std::vector<unsigned int>();
    for (auto replicaThread : replica_owned_threads){
        replica_master_IDS.push_back(replicaThread[0]);
    }

    reGMPI.replicaMasterIDs = replica_master_IDS;
    reGMPI.replicaThreads = replica_owned_threads;
    reGMPI.threadReplicaMap = thread_id_replica_map;

    //replica Graph out:
    if(globalThreadID == 0){
        std::cerr << "num of threads: "  << totalNumberOfThreads << "\n";
        std::cerr << "threads per rep: " << threadsPerReplicaSimulation << "\n";
        std::cout << "Masters of "<<reGMPI.numberOfReplicas<< "replicas \t";
    }
    std::cout << "\n";
    MPI_DEBUG(1, "REPLICA_ID \t " << globalThreadID << "\t Simulation_ID\t"<< simulationID << "\t RE_GRAPH COMM ESTABLISHED\n");
    MPI_Barrier(MPI_COMM_WORLD);    //wait for all threads to register!

    MPI_Datatype MPI_VEC;

    try { // SUBCOMS //Exchange structures
        // Vector
        MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_VEC);
        MPI_Type_commit(&MPI_VEC);

        // Box
        MPI_Type_contiguous(3, MPI_VEC, &MPI_BOX);
        MPI_Type_commit(&MPI_BOX);

        // VArray with size of system
        MPI_Type_contiguous(numAtoms, MPI_VEC, &MPI_VARRAY);
        MPI_Type_commit(&MPI_VARRAY);

        // defining struct with non static replica information
        int blocklen[] = {3, 3};
        MPI_Datatype typ[] = {MPI_INT, MPI_DOUBLE};
        MPI_Aint intext;
        MPI_Aint lb; //not used
        MPI_Type_get_extent(MPI_INT, &lb, &intext);

        MPI_Aint disps[] = {(MPI_Aint) 0, 4 * intext};
        MPI_Type_create_struct(2, blocklen, disps, typ, &MPI_REPINFO);
        MPI_Type_commit(&MPI_REPINFO);

        if (reedsSim > 0) {
            MPI_Type_contiguous(numEDSstates, MPI_DOUBLE, &MPI_EDSINFO);
            MPI_Type_commit(&MPI_EDSINFO);
        }
    } catch (const std::exception &e) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\tErrors during MPI-Initialisation!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an Exception in MPI-Initialisation!\n\nMessage:\t";
        std::cout << msg << e.what() << std::endl;
        std::cerr << msg << e.what() << std::endl;
        MPI_Comm_free(&replica_mpi_control.comm); //Clean up
        MPI_Finalize();

        return 1;
    } catch (...) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\tErrors during MPI-Initialisation!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial test Parsing of the MPI Parameters!\n Hit the developers!\n\nMessage:\t";
        std::cout << msg << std::endl;
        std::cerr << msg << std::endl;
        MPI_Comm_free(&replica_mpi_control.comm); //Clean up
        MPI_Finalize();
        return -1;
    }


    if (globalThreadID == 0) { //Nice message
        std::ostringstream msg;
        msg.clear();
        msg << "\t Short RE-Setting Overview:\n";
        msg << "\t continuation:\t" << cont << "\n";
        msg << "\t equilibration runs:\t" << equil_runs << "\n";
        msg << "\t Exchange Trials runs:\t" << sim_runs << "\n";

        msg << "\t numReplicas:\t" << numReplicas << "\n";

        if (reedsSim > 0) {
            msg << "\n\t RE-EDS:\n";
            msg << "\t reeds_control:\t" << reedsSim << "\n";
            msg << "\t numStates:\t" << numEDSstates << "\n";
            msg << "\t numSValues:\t" << numSVals << "\n";
            msg << "\t numEoffs:\t" << numEoff << "\n";
            if (reedsSim == 2) {
              msg << "\t periodic:\t" << periodic << "\n";
            }
        }

        msg << "\n\t Mpi Settings:\t\n";


        // make sure all nodes have initialized everything
        msg << "\t repIDS:\t\n";
        for(unsigned int x=0; x < replica_owned_threads.size(); x++){
            msg << "\t " << x << ". Replica: ";
            for(int threadID :  replica_owned_threads[x]){
                msg << "\t" << threadID;
            }
            msg<< "\n";
        }
        msg<< "\n";


        /*msg<< "\t repMap:\n";
        for(std::pair<unsigned int, unsigned int> p : thread_id_replica_map){
            msg << "\t " << p.first << "\t" << p.second <<"\n";
        }
        msg << "\n";
         */

        msg << "\n==================================================\n\tFinished Setting Up Simulation\n\n==================================================\n";
        std::cout << msg.str();
        std::cerr << msg.str();
    }

    MPI_DEBUG(1, "FINISHED MPI MAPPING!");


    //////////////////////////////
    /// Starting master-slave Pattern
    //////////////////////////////

    if (globalThreadID == 0) { //MASTER

        //nice messages
        std::cout << "\n==================================================\n\tStart REPLICA EXCHANGE SIMULATION:\n\n==================================================\n";
        std::cout << "numreplicas:\t " << numReplicas << "\n";
        std::cout << "num Slaves:\t " << numReplicas - 1 << "\n";
        std::cerr << "\n==================================================\n\tStart REPLICA EXCHANGE SIMULATION:\n\n==================================================\n";
        std::cerr << "numreplicas:\t " << numReplicas << "\n";
        std::cerr << "num Slaves:\t " << numReplicas - 1 << "\n";
        MPI_DEBUG(1, "Master \t " << globalThreadID<< "simulation: " << simulationID);

        //CONSTRUCT
        // Select repex Implementation - Polymorphism
        MPI_DEBUG(1, "MASTER " << globalThreadID << "::Constructor: START ");
        util::replica_exchange_master_interface * Master;

        switch(reedsSim) {
              case 0:
                  DEBUG(1, "Master \t Constructor\n");
                  std::cerr << "Constructor:\t Master\n";
                  Master = new util::replica_exchange_master(args, cont, globalThreadID, reGMPI, replica_mpi_control);
                  break;
              case 1:
                  DEBUG(1, "Master_eds \t Constructor\n")
                  std::cerr << "Constructor:\t Master_eds\n";
                  Master = new util::replica_exchange_master_eds(args, cont, globalThreadID, reGMPI, replica_mpi_control);
                  break;
              case 2:
                  DEBUG(1, "Master_2d_s_eoff_eds \t Constructor\n")
                  std::cerr << "Constructor:\t Master_2d_s_eoff_eds\n";
                  Master = new util::replica_exchange_master_2d_s_eoff_eds(args, cont, globalThreadID, reGMPI, replica_mpi_control);
                  break;
          }
        MPI_DEBUG(1, "MASTER " << globalThreadID << "::Constructor: DONE ");


        MPI_DEBUG(1, "Master \t INIT START");

        Master->init();

        Master->init_repOut_stat_file();

        MPI_DEBUG(1, "Master \t INIT DONE")


        //do md:
        unsigned int trial = 0;
        MPI_DEBUG(1, "Master \t Equil: " << equil_runs<< " steps")
        for (; trial < equil_runs; ++trial) { // for equilibrations
            Master->run_MD();
        }
        //MPI_Finalize();
        //return 0;

        //Vars for timing
        int hh, mm, ss, eta_hh, eta_mm, eta_ss = 0;
        double eta_spent, percent, spent = 0.0;
        trial = 0; //reset trials
        MPI_DEBUG(1, "Master \t MD: " << sim_runs<< " steps")
        for (; trial < sim_runs; ++trial) { //for repex execution
            MPI_DEBUG(2, "Master " << globalThreadID << " \t MD trial: " << trial << "\n");
            MPI_DEBUG(2, "Master " << globalThreadID << " \t run_MD START " << trial << "\n");
            Master->run_MD();
            MPI_DEBUG(2, "Master " << globalThreadID << " \t swap START " << trial << "\n")
            Master->swap();
            MPI_DEBUG(2, "Master " << globalThreadID << " \t receive START " << trial << "\n");
            Master->receive_from_all_slaves();
            MPI_DEBUG(2, "Master " << globalThreadID << " \t write START " << trial << "\n");
            Master->write();


            if ((sim_runs / 10 > 0) && (trial % (sim_runs / 10) == 0)) { //Timer
                percent = double(trial) / double(sim_runs);

                spent = util::now() - start;
                hh = int(spent / 3600);
                mm = int((spent - hh * 3600) / 60);
                ss = int(spent - hh * 3600 - mm * 60);

                if(trial > 10){
                    eta_spent = spent / trial * sim_runs - spent;
                    eta_hh = int(eta_spent / 3600);
                    eta_mm = int((eta_spent - eta_hh * 3600) / 60);
                    eta_ss = int(eta_spent - eta_hh * 3600 - eta_mm * 60);
                }

                std::cout << "\nREPEX:       " << std::setw(2) << percent * 100 << "% done..." << std::endl;
                std::cout << "REPEX: spent " << hh << ":" << mm << ":" << ss << std::endl;
                std::cout << "REPEX: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
                std::cerr << "\nREPEX:       " << std::setw(2) << percent * 100 << "% done..." << std::endl;
                std::cerr << "REPEX: spent " << hh << ":" << mm << ":" << ss << std::endl;
                std::cerr << "REPEX: ETA   " << eta_hh << ":" << eta_mm << ":" << eta_ss << std::endl;
            }
        }
        MPI_DEBUG(1, "Master \t \t finalize ")

        Master->write_final_conf();
        std::cout << "\n=================== Master Node " << globalThreadID << "  finished successfully!\n";

        } else { //SLAVES
        MPI_DEBUG(1, "Slave " << globalThreadID << "simulation: " << simulationID);


        // Select repex Implementation - Polymorphism
        MPI_DEBUG(1, "Slave " << globalThreadID << "::Constructor: START ");
        util::replica_exchange_slave_interface* Slave;

        switch(reedsSim) {
              case 0:
                  Slave = new util::replica_exchange_slave(args, cont, globalThreadID, reGMPI, replica_mpi_control);
                  break;
              case 1:
                  Slave = new util::replica_exchange_slave_eds(args, cont, globalThreadID, reGMPI, replica_mpi_control);
                  break;
              case 2:
                  Slave = new util::replica_exchange_slave_2d_s_eoff_eds(args, cont, globalThreadID, reGMPI, replica_mpi_control);
                  break;
          }
        MPI_DEBUG(1, "Slave " << globalThreadID << "::Constructor: DONE ")


        //INIT
        MPI_DEBUG(1, "Slave " << globalThreadID << " \t INIT START");
        Slave->init();

        MPI_DEBUG(1, "Slave " << globalThreadID << " \t INIT Done")




        //do md:
        unsigned int trial = 0;
        MPI_DEBUG(1, "Slave \t Equil: " << equil_runs<< " steps")
        for (; trial < equil_runs; ++trial) { // for equilibrations
            Slave->run_MD();
        }

        //MPI_Finalize();
        //return 0;
        trial = 0; //reset trials
        MPI_DEBUG(1, "Slave " << globalThreadID << " \t MD " << sim_runs << " steps")
        for (; trial < sim_runs; ++trial) { //for repex execution
            MPI_DEBUG(2, "Slave " << globalThreadID << " \t MD trial: " << trial << "\n")
            MPI_DEBUG(2, "Slave " << globalThreadID << " \t run_MD START " << trial << "\n")
            Slave->run_MD();
            MPI_DEBUG(2, "Slave " << globalThreadID << " \t swap START " << trial << "\n");
            Slave->swap();
            MPI_DEBUG(2, "Slave " << globalThreadID << " \t send START " << trial << "\n");
            Slave->send_to_master();
            MPI_DEBUG(2, "Slave " << globalThreadID << " \t send Done " << trial << "\n");
        }
        MPI_DEBUG(1, "Slave " << globalThreadID << " \t Finalize");

        Slave->write_final_conf();
        std::cout << "\n=================== Slave Node " << globalThreadID << "  finished successfully!\n";
    }
    MPI_DEBUG(1, "Lets end this!\n")

    MPI_Barrier(MPI_COMM_WORLD); //Make sure all processes finished.

    //last MASTER OUTPUT
    if (globalThreadID == 0) {
        //FINAL OUTPUT - Time used:
        double end = util::now();
        double duration = end - start;
        int durationMin = int(duration / 60);
        int durationHour = int(std::floor(durationMin / 60));
        int durationMinlHour = int(std::floor(durationMin - durationHour * 60));
        int durationSlMin = int(std::floor(duration - (durationMinlHour + durationHour * 60)*60));

        //Todo: if (cond){finished succ}else{not} bschroed
        std::string msg("\n====================================================================================================\n\tREPLICA EXCHANGE SIMULATION finished!\n====================================================================================================\n");
        std::cout << msg;
        std::cerr << msg;

        std::cout << "TOTAL TIME USED: \n\th:min:s\t\tseconds\n"
                << "\t" << durationHour << ":" << durationMinlHour << ":" << durationSlMin << "\t\t" << duration << "\n";
        std::cerr << "TOTAL TIME USED: \n\th:min:s\t\tseconds\n"
                << "\t" << durationHour << ":" << durationMinlHour << ":" << durationSlMin << "\t\t" << duration << "\n";
    }

    MPI_DEBUG(1, "Before Finalize!\n")
    MPI_Comm_free(&reGMPI.comm); //Clean up
    MPI_Comm_free(&replica_mpi_control.comm); //Clean up

    MPI_Finalize();
    return 0;
#else
    std::cout << "\n==================================================\n\tGROMOS Replica Exchange:\n==================================================\n";
    std::cout << argv[0] << " needs MPI to run\n\tuse --enable-mpi at configure and appropriate compilers\n" << std::endl;
    return 1;
#endif
}
