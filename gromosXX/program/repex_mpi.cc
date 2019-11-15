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
#include <io/print_block.h>

#include <thread>
#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>
#include <math/gmath.h>


#ifdef XXMPI
#include <mpi.h>
#endif

#include <util/replicaExchange/mpi_tools.h>
#include <util/replicaExchange/replica/replica.h>
#include <util/replicaExchange/replica/replica_reeds.h>

#include <util/replicaExchange/replica_exchangers/replica_exchange_master.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_slave.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_master_eds.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_slave_eds.h>

#include <util/replicaExchange/repex_mpi.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange
#define XXMPI

int main(int argc, char *argv[]) {
#ifdef XXMPI
    //initializing MPI
    MPI_Init(&argc, &argv);
    const double start = MPI_Wtime();
    int totalNumberOfThreads;
    int uniqueThreadID;

    MPI_Comm_size(MPI_COMM_WORLD, &totalNumberOfThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &uniqueThreadID);

    if (uniqueThreadID == 0) {
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
    MPI_DEBUG(1, "RANK: "<<uniqueThreadID<<" totalRankNum: "<< totalNumberOfThreads<<": hello world!\n");

    //GLOBAL SETTING VARIABLES
    bool reedsSim;

    unsigned int equil_runs;
    unsigned int sim_runs;
    unsigned int total_runs;
    unsigned int numAtoms;
    int numReplicas;
    unsigned int numEDSstates;
    int cont;


    int subThreadOfSimulation; //MPI thread belongs to replica Simulation:
    int threadInSimulation;
    std::map<unsigned int, unsigned int> repMap; // where is which replica
    std::vector< std::vector<int> > repIDs; // set IDs for each replica
    int threadsPerReplicaSimulation; //remove this param - also from exchangers
    
    try {
        {
            topology::Topology topo;
            configuration::Configuration conf;
            algorithm::Algorithm_Sequence md;
            simulation::Simulation sim;

            // read in parameters
            bool quiet = true;
            if (uniqueThreadID == 0) {
                quiet = true;
            }

            if (io::read_parameter(args, sim, std::cout, true)) {
                if (uniqueThreadID == 0) {
                    std::cerr << "\n\t########################################################\n"
                            << "\n\t\tErrors during read Parameters reading!\n"
                            << "\n\t########################################################\n";
                    io::messages.display(std::cout);
                    io::messages.display(std::cerr);

                }
                MPI_Finalize();
                return 1;
            }
            
            MPI_DEBUG(2, "RANK: "<<uniqueThreadID<<" done with params\n");    //TODO: lower verb level
            
            //set global parameters
            cont = sim.param().replica.cont;
            equil_runs = sim.param().replica.equilibrate;
            sim_runs = sim.param().replica.trials;
            total_runs = sim.param().replica.trials + equil_runs;
            numAtoms = topo.num_atoms();
            reedsSim = sim.param().reeds.reeds;
            
            if (reedsSim) {
                numReplicas = sim.param().reeds.num_l;
                numEDSstates = sim.param().reeds.eds_para[0].numstates;
            } else {
                numReplicas = sim.param().replica.num_T * sim.param().replica.num_l;
                numEDSstates = 0;
            }
            
            //MPI THREAD SIMULATION SPLITTING
            //needed to be calculated here
            //repIDs every node gets one element of that vector
            repIDs.resize(numReplicas);

            // counts through every replica and assigns it to respective node
            // starts at beginning if numReplicas > size
            // could be optimized by putting neighboring replicas on same node; less communication...
            //TODO: not sure if this code really does what it should.! (idea would)

            //MPI THREAD SPLITING ONTO Simulation - REPLICAS
            threadsPerReplicaSimulation = totalNumberOfThreads / numReplicas;
            int leftOverThreads = totalNumberOfThreads % numReplicas;
            
            unsigned int threadID =0;
            int replica_offset = 0;
            std::cout << "left_over "<< leftOverThreads << "\n";
            for (int replicaSimulationID = 0; replicaSimulationID < numReplicas; replicaSimulationID++) {

                for (int replicaSimulationSubThread = 0; replicaSimulationSubThread < threadsPerReplicaSimulation; replicaSimulationSubThread++) {

                    threadID = replicaSimulationSubThread + replicaSimulationID*threadsPerReplicaSimulation+replica_offset;
                    repMap.insert(std::pair<unsigned int, unsigned int>(threadID, replicaSimulationID));
                    repIDs[replicaSimulationID].push_back(threadID);

                }
                if(leftOverThreads>0 and replicaSimulationID < leftOverThreads){    //left over threads are evenly distirbuted.

                    threadID = threadsPerReplicaSimulation + replicaSimulationID*totalNumberOfThreads+replica_offset;
                    repMap.insert(std::pair<unsigned int, unsigned int>(threadID, replicaSimulationID));
                    repIDs[replicaSimulationID].push_back(threadID);
                    replica_offset++;

                }    
            }            
            
            if(uniqueThreadID==0){  //if DEBUG TODO: remove bschroed
                std::cout<< "repIDS:\t";
                for(unsigned int x=0; x < repIDs.size(); x++){
                    std::cout << "\t" << x << "G: ";
                    for(unsigned int y=0; y< repIDs[x].size(); y++){
                        std::cout << "\t" << repIDs[x][y];
                    }
                    std::cout<< "\n";
                }
                std::cout<< "\n";

                std::cout<< "repMap:\n";
                for(std::pair<unsigned int, unsigned int> p : repMap){
                    std::cout << "\t " << p.first << "\t" << p.second <<"\n";
                }
                std::cout<< "\n";
            }       //endif DEBUG TODO: remove bschroed
            subThreadOfSimulation = repMap[uniqueThreadID];
            
            int counter = 0;
            for(int x : repIDs[subThreadOfSimulation]){
                if(x==uniqueThreadID){
                    threadInSimulation = counter;
                    break;
                }
                counter++;
            } 
            
            // read in the rest
            if (io::read_input_repex(args, topo, conf, sim, md, subThreadOfSimulation, uniqueThreadID, std::cout, quiet)) {
                if (uniqueThreadID == 0) {
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
            MPI_DEBUG(2, "RANK: "<<uniqueThreadID<<" done with repex_in\n");    //TODO: lower verb level
            
            //if any replica Ex block - present   
            if (sim.param().reeds.reeds == false && sim.param().replica.retl == false) {
                if (uniqueThreadID == 0) {
                    std::cerr << "\n\t########################################################\n"
                            << "\n\t\tErrors during initial Parameter reading! "
                            << "\n\t\t    No repex block was satisfied!\n"
                            << "\n\t########################################################\n";
                    std::cerr << "\n Please add one RE-block (e.g.:REPLICA or REEDS) to the imd file.\n";
                    std::cout << "\n Please add one RE-block (e.g.:REPLICA or REEDS)  to the imd file.\n";
                }
                MPI_Finalize();
                return 1;
            }
            MPI_DEBUG(2, "RANK: "<<uniqueThreadID<<" done with replica Ex check\n");    //TODO: lower verb level
            
            if (io::check_parameter(sim)) {
                if (uniqueThreadID == 0) {
                    std::cerr << "\n\t########################################################\n"
                            << "\n\t\tErrors during initial Parameter reading!\n"
                            << "\n\t########################################################\n";
                    io::messages.display(std::cout);
                    io::messages.display(std::cerr);
                }
                MPI_Finalize();
                return 1;
            }
            MPI_DEBUG(2, "RANK: "<<uniqueThreadID<<" done with param check\n");    //TODO: lower verb level
            
  
            //SOME additional Checks
            //if enough threads avail
            if (totalNumberOfThreads < numReplicas) {
                if (uniqueThreadID == 0) {
                    std::cerr << "\n\t########################################################\n"
                            << "\n\t\tErrors during initial Parameter reading!\n"
                            << "\n\t########################################################\n";
                    std::cerr << "\n There were not enough MPI threads assigned to this run!\n"
                            << "FOUND THREADS: " << totalNumberOfThreads << "\tNEED: " << numReplicas << "\n";
                    std::cout << "\n There were not enough MPI thread assigned to this run!\n"
                            << "FOUND THREADS: " << totalNumberOfThreads << "\tNEED: " << numReplicas << "\n";
                    MPI_Finalize();

                }
                return -1;
            }

            if (uniqueThreadID == 0) { //Nice message
                std::ostringstream msg;
                msg << "\t Short RE-Setting Overview:\n";
                msg << "\t continuation:\t" << cont << "\n";
                msg << "\t equilibration runs:\t" << equil_runs << "\n";
                msg << "\t Exchange Trials runs:\t" << sim_runs << "\n";
                msg << "\t Simulation Steps Between Trials:\t" << sim.param().step.number_of_steps << "\n";
                msg << "\t Total Simulation Time:\t" << sim.param().step.number_of_steps * sim_runs * sim.param().step.dt << "ps\n";
                msg << "\t numReplicas:\t" << numReplicas << "\n";

                if (reedsSim) {
                    msg << "\n\t RE-EDS:\n";
                    msg << "\t numSValues:\t" << sim.param().reeds.num_l << "\n";
                    msg << "\t numStates:\t" << numEDSstates << "\n";
                }

                msg << "\n==================================================\n\tFinished Initial Parsing\n\n==================================================\n";
                std::cout << msg.str();
                std::cerr << msg.str();

            }
        }
    } catch (const std::exception &e) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\t Exception was thrown in initial parsing!\n"
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
    MPI_DEBUG(1, "RANK: "<<uniqueThreadID<<" Done with init parse\n");    //TODO: lower verb level
   
    io::messages.clear();
    MPI_DEBUG(1, "REPLICA_ID \t " << uniqueThreadID << "\t Simulation_ID\t"<< subThreadOfSimulation << "\t SimulationThread\t"<<threadInSimulation<<"\n")

    //////////////////////////
    // defining MPI Datatypes
    //////////////////////////
    //GLOBAL MPI VARIABLES:
    MPI_Datatype MPI_VEC;
    
    try { 

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
        MPI_Type_extent(MPI_INT, &intext);

        MPI_Aint disps[] = {(MPI_Aint) 0, 4 * intext};
        MPI_Type_create_struct(2, blocklen, disps, typ, &MPI_REPINFO);
        MPI_Type_commit(&MPI_REPINFO);

        if (reedsSim) {
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
        return 1;
    } catch (...) {
        std::cerr << "\n\t########################################################\n"
                << "\n\t\tErrors during MPI-Initialisation!\n"
                << "\n\t########################################################\n";
        std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial test Parsing of the MPI Parameters!\n Hit the developers!\n\nMessage:\t";
        std::cout << msg << std::endl;
        std::cerr << msg << std::endl;
        MPI_Finalize();
        return -1;
    }

    // make sure all nodes have initialized everything
    MPI_Barrier(MPI_COMM_WORLD);


    //////////////////////////////
    /// Starting master-slave Pattern
    //////////////////////////////

    if (uniqueThreadID == 0) { //MASTER
        //std::cout << "RANK: "<< uniqueThreadID <<"\tSLEPPING Master\n";
        //MPI_Barrier(MPI_COMM_WORLD);
                
        //nice messages
        std::cout << "\n==================================================\n\tStart REPLICA EXCHANGE SIMULATION:\n\n==================================================\n";
        std::cout << "numreplicas:\t " << numReplicas << "\n";
        std::cout << "num Slaves:\t " << numReplicas - 1 << "\n";
        std::cerr << "\n==================================================\n\tStart REPLICA EXCHANGE SIMULATION:\n\n==================================================\n";
        std::cerr << "numreplicas:\t " << numReplicas << "\n";
        std::cerr << "num Slaves:\t " << numReplicas - 1 << "\n";
        MPI_DEBUG(1, "Master \t " << uniqueThreadID<< "simulation: " << subThreadOfSimulation)

        // Select repex Implementation - Polymorphism
        util::replica_exchange_master* Master;
        try {
            if (reedsSim) {
                DEBUG(1, "Master_eds \t Constructor")    
                Master = new util::replica_exchange_master_eds(args, cont, uniqueThreadID, threadInSimulation, subThreadOfSimulation, threadsPerReplicaSimulation, totalNumberOfThreads, numReplicas , repIDs[uniqueThreadID], repMap);
            } else {
                DEBUG(1, "Master \t Constructor")
                Master = new util::replica_exchange_master(args, cont, uniqueThreadID, threadInSimulation, subThreadOfSimulation, threadsPerReplicaSimulation,  totalNumberOfThreads, numReplicas, repIDs[uniqueThreadID], repMap);
            }
        } catch (...) {
            std::cerr << "\n\t########################################################\n"
                    << "\n\t\tErrors during Master Class Init!\n"
                    << "\n\t########################################################\n";
            std::cerr << "\nREPEX:       Failed in initialization of Replica rank: " << uniqueThreadID << "" << std::endl;
            std::cout << "\nREPEX:       Failed in initialization of Replica rank: " << uniqueThreadID << "" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, E_USAGE);
            return 1;
        }

        MPI_DEBUG(1, "Master \t INIT")
        Master->init();
        Master->init_repOut_stat_file();

        //do md:
        unsigned int trial = 0;
        MPI_DEBUG(1, "Master \t \t \t Equil: " << equil_runs)
        for (; trial < equil_runs; ++trial) { // for equilibrations
            Master->run_MD();
        }
        MPI_DEBUG(1, "Master \t \t MD: " << total_runs)

                //Vars for timing
                int hh, mm, ss, eta_hh, eta_mm, eta_ss = 0;
        double eta_spent, percent, spent = 0.0;
        trial = 0; //reset trials
        for (; trial < sim_runs; ++trial) { //for repex execution
            DEBUG(2, "Master " << uniqueThreadID << " \t MD trial: " << trial << "\n")\
            DEBUG(2, "Master " << uniqueThreadID << " \t run_MD START " << trial << "\n")
            Master->run_MD();
            DEBUG(2, "Master " << uniqueThreadID << " \t run_MD DONE " << trial << "\n")

            DEBUG(2, "Master " << uniqueThreadID << " \t swap START " << trial << "\n")
            Master->swap();
            DEBUG(2, "Master " << uniqueThreadID << " \t run_MD DONE " << trial << "\n")

            DEBUG(2, "Master " << uniqueThreadID << " \t receive START " << trial << "\n")
            Master->receive_from_all_slaves();
            DEBUG(2, "Master " << uniqueThreadID << " \t write START " << trial << "\n")
            Master->write();

            if ((total_runs / 10 > 0) && (trial % (total_runs / 10) == 0)) { //Timer 
                percent = double(trial) / double(total_runs);
                spent = util::now() - start;
                hh = int(spent / 3600);
                mm = int((spent - hh * 3600) / 60);
                ss = int(spent - hh * 3600 - mm * 60);
                if (trial > 10) {
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
        
        if (reedsSim) {
            DEBUG(1, "Master_eds \t FINAL_OUT")    
            dynamic_cast<util::replica_exchange_master_eds*>(Master)->write_final_conf();
        } else {
            DEBUG(1, "Master \t FINAL_OUT")
            dynamic_cast<util::replica_exchange_master*>(Master)->write_final_conf();
        }
       
    } else { //SLAVES    
        //if(uniqueThreadID < 3){
        //std::cout << "RANK: "<< uniqueThreadID <<"\tSLEPPING SLAVE\n";
        //MPI_Barrier(MPI_COMM_WORLD);
        //}
        MPI_DEBUG(1, "Slave " << uniqueThreadID << "simulation: " << subThreadOfSimulation)

        // Select repex Implementation - Polymorphism
        util::replica_exchange_slave* Slave;
        if (reedsSim) {
            Slave = new util::replica_exchange_slave_eds(args, cont, uniqueThreadID, threadInSimulation, subThreadOfSimulation, threadsPerReplicaSimulation, repIDs[subThreadOfSimulation], repMap);
        } else {
            Slave = new util::replica_exchange_slave(args, cont, uniqueThreadID, threadInSimulation, subThreadOfSimulation, threadsPerReplicaSimulation,repIDs[subThreadOfSimulation], repMap);
        }

        MPI_DEBUG(1, "Slave " << uniqueThreadID << " \t INIT")
        Slave->init();

        //do md:
        unsigned int trial = 0;
        MPI_DEBUG(1, "Slave " << uniqueThreadID << " \t EQUIL " << equil_runs << " steps")
        for (; trial < equil_runs; ++trial) { // for equilibrations
            Slave->run_MD();
        }

        MPI_DEBUG(1, "Slave " << uniqueThreadID << " \t MD " << total_runs << " steps")
        for (; trial < total_runs; ++trial) { //for repex execution
            DEBUG(2, "Slave " << uniqueThreadID << " \t MD trial: " << trial << "\n")
            DEBUG(2, "Slave " << uniqueThreadID << " \t run_MD START " << trial << "\n")
            Slave->run_MD();
            DEBUG(2, "Slave " << uniqueThreadID << " \t swap START " << trial << "\n")
            Slave->swap();
            DEBUG(2, "Slave " << uniqueThreadID << " \t send START " << trial << "\n")
            Slave->send_to_master();
        }

        MPI_DEBUG(1, "Slave " << uniqueThreadID << " \t Finalize")
        
        Slave->write_final_conf();
        
        std::cout << "\n=================== Slave Node " << uniqueThreadID << "  finished successfully!\n";
    }

    MPI_Barrier(MPI_COMM_WORLD); //Make sure all processes finished.

    //last MASTER OUTPUT
    if (uniqueThreadID == 0) {
        //FINAL OUTPUT - Time used:
        double end = MPI_Wtime();
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

    MPI_Finalize();
    return 0;
#else
    std::cout << "\n==================================================\n\tGROMOS Replica Exchange:\n==================================================\n";
    std::cout << argv[0] << " needs MPI to run\n\tuse --enable-mpi at configure and appropriate compilers\n" << std::endl;
    return 1;
#endif
}

