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

#include <util/replicaExchange/replica_exchange_master.h>
#include <util/replicaExchange/replica_exchange_slave.h>
#include <util/replicaExchange/replica_exchange_master_eds.h>
#include <util/replicaExchange/replica_exchange_slave_eds.h>

#include <util/replicaExchange/repex_mpi.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

int main(int argc, char *argv[]) {
#ifdef XXMPI
  //initializing MPI
  MPI_Init(&argc, &argv);
  const double start = MPI_Wtime();
  int size;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if(rank == 0){
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
  try{
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
  }
  catch (const std::exception &e){
      std::cerr <<"\n\t########################################################\n" 
        << "\n\t\t Exception was thrown in Commandline Parsing!\n"
        <<"\n\t########################################################\n" ;
      std::string msg = "ERROR!\n Uh OH! Caught an Exception in initial Args parse!\n\n";
      std::cout << msg << e.what() << std::endl;
      std::cerr << msg << e.what() << std::endl;
      MPI_Finalize();
      MPI_Abort(MPI_COMM_WORLD, E_USAGE);
      return -1;
  }
  catch (...){
      std::cerr <<"\n\t########################################################\n" 
        << "\n\t\t Exception was thrown in Commandline Parsing!\n"
        <<"\n\t########################################################\n" ;
      std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial Args parse!!\n Hit the developers!\n\n";
      std::cout << msg << std::endl;
      std::cerr << msg << std::endl;
      MPI_Finalize();
      return -1;
  }
  
  //GLOBAL SETTING VARIABLES
  bool reedsSim;

  unsigned int equil_runs;
  unsigned int sim_runs;
  unsigned int total_runs;
  unsigned int numAtoms;
  int numReplicas;
  unsigned int numEDSstates;
  int cont;

  try{
    {
        topology::Topology topo;
        configuration::Configuration conf;
        algorithm::Algorithm_Sequence md;
        simulation::Simulation sim;

        // read in parameters
        bool quiet = true;
        if(rank == 0){
              quiet = true;
        }

        if (io::read_parameter(args,sim,std::cout, true)){
              if (rank == 0) {
                  std::cerr <<"\n\t########################################################\n" 
                    << "\n\t\tErrors during read Parameters reading!\n"
                    <<"\n\t########################################################\n" ;
                  io::messages.display(std::cout);
                  io::messages.display(std::cerr);

              }
              MPI_Finalize();
              return 1;
        }

        // read in the rest
        if(io::read_input_repex(args, topo, conf, sim, md, rank, std::cout, quiet)){
              if (rank == 0) {
                  std::cerr <<"\n\t########################################################\n" 
                          << "\n\t\tErrors during initial Parameter reading!\n"
                          <<"\n\t########################################################\n" ;
                  io::messages.display(std::cout);
                  io::messages.display(std::cerr);
                  MPI_Abort(MPI_COMM_WORLD, E_USAGE);
              }
              MPI_Finalize();
              return 1;
          }

        //if any replica Ex block - present   
        if(sim.param().reeds.reeds == false && sim.param().replica.retl  == false ){
              if(rank == 0){
                  std::cerr <<"\n\t########################################################\n" 
                          << "\n\t\tErrors during initial Parameter reading! "
                          << "\n\t\t    No repex block was satisfied!\n"
                          <<"\n\t########################################################\n" ;
                  std::cerr << "\n Please add one RE-block (e.g.:REPLICA or REEDS) to the imd file.\n";
                  std::cout <<  "\n Please add one RE-block (e.g.:REPLICA or REEDS)  to the imd file.\n";
              }
              MPI_Finalize();
              return 1;
        }
      
        if (io::check_parameter(sim)){
          if (rank == 0) {
                std::cerr <<"\n\t########################################################\n" 
                  << "\n\t\tErrors during initial Parameter reading!\n"
                  <<"\n\t########################################################\n" ;
                io::messages.display(std::cout);
                io::messages.display(std::cerr);
          }
          MPI_Finalize();
          return 1;
        }

        //set global parameters
        cont = sim.param().replica.cont;
        equil_runs = sim.param().replica.equilibrate;
        sim_runs = sim.param().replica.trials;
        total_runs = sim.param().replica.trials + equil_runs;
        numReplicas = sim.param().replica.num_T * sim.param().replica.num_l;
        numAtoms = topo.num_atoms();
        reedsSim = sim.param().reeds.reeds;

        if(reedsSim){
            numEDSstates=sim.param().reeds.eds_para[0].numstates;
        }else{
            numEDSstates=0;
        }

        //SOME additional Checks
        //if enough threads avail
        if( size  < numReplicas){
          if(rank == 0){
              std::cerr <<"\n\t########################################################\n" 
                      << "\n\t\tErrors during initial Parameter reading!\n"
                      <<"\n\t########################################################\n" ;
              std::cerr << "\n There were not enough MPI thread assigned to this run!\n"
                      << "FOUND THREAD: "<<size<< "\tNEED: "<<numReplicas<<"\n";
              std::cout << "\n There were not enough MPI thread assigned to this run!\n"
                        << "FOUND THREAD: "<<size<< "\tNEED: "<<numReplicas<<"\n";
              MPI_Finalize();

          }
          return -1;
        }

        if(rank == 0){  //Nice message
            std::ostringstream msg;
            msg << "\t Short RE-Setting Overview:\n";
            msg << "\t continuation:\t"<<cont<<"\n";
            msg << "\t equilibration runs:\t"<<equil_runs<<"\n";
            msg << "\t Exchange Trials runs:\t"<<sim_runs<<"\n";
            msg << "\t Simulation Steps Between Trials:\t"<<sim.param().step.number_of_steps<<"\n";
            msg << "\t Total Simulation Time:\t"<<sim.param().step.number_of_steps*sim_runs*sim.param().step.dt<<"ps\n";
            msg << "\t numReplicas:\t"<<numReplicas<<"\n";

            if(reedsSim){
                msg << "\n\t RE-EDS:\n";
                msg << "\t numSValues:\t"<<sim.param().reeds.num_l<<"\n";
                msg << "\t numStates:\t"<<numEDSstates<<"\n";
            }

            msg << "\n==================================================\n\tFinished Initial Parsing\n\n==================================================\n";
            std::cout << msg.str();
            std::cerr << msg.str();

        }
    }
  }
  catch (const std::exception &e){
      std::cerr <<"\n\t########################################################\n" 
        << "\n\t\t Exception was thrown in initial parsing!\n"
        <<"\n\t########################################################\n" ;
      std::string msg = "ERROR!\n Uh OH! Caught an Exception in initial test Parsing of the Parameters!\n\n";
      std::cout << msg << e.what() << std::endl;
      std::cerr << msg << e.what() << std::endl;
      MPI_Finalize();
      return -1;
  }
  catch (...){
      std::cerr <<"\n\t########################################################\n" 
        << "\n\t\t Exception was thrown!\n"
        <<"\n\t########################################################\n" ;
      std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial test Parsing of the Parameters!\n Hit the developers!\n\n";
      std::cout << msg << std::endl;
      std::cerr << msg << std::endl;
      MPI_Finalize();
      return -1;
  }


  io::messages.clear();

  //////////////////////////
  // defining MPI Datatypes
  //////////////////////////
  //GLOBAL MPI VARIABLES:
  MPI_Datatype MPI_VEC;
  std::map<unsigned int, unsigned int> repMap;     // where is which replica
  std::vector< std::vector<int> > repIDs;     // set IDs for each replica
  
  try{
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

    if(reedsSim){
      MPI_Type_contiguous(numEDSstates, MPI_DOUBLE, &MPI_EDSINFO);
      MPI_Type_commit(&MPI_EDSINFO);
    }

    assert(numReplicas > 0);// is there a positive nonzero num of Replicas?

    //repIDs every node gets one element of that vector
    repIDs.resize(size);

    // counts through every replica and assigns it to respective node
    // starts at beginning if numReplicas > size
    // could be optimized by putting neighboring replicas on same node; less communication...
    for (int i = 0; i < ceil((float) numReplicas / (float) size); ++i) {
      for (int j = 0; j < size; ++j) {
        unsigned int ID = j + i*size;
        if (ID >= numReplicas)
          break;
        repMap.insert(std::pair<unsigned int, unsigned int>(ID, j));
        repIDs[j].push_back(ID);
      }
    }
  }
  catch (const std::exception &e){
      std::cerr <<"\n\t########################################################\n" 
        << "\n\t\tErrors during MPI-Initialisation!\n"
        <<"\n\t########################################################\n" ;
      std::string msg = "ERROR!\n Uh OH! Caught an Exception in MPI-Initialisation!\n\n";
      std::cout << msg << e.what() << std::endl;
      std::cerr << msg << e.what() << std::endl;
      return 1;
  }
  catch (...){
      std::cerr <<"\n\t########################################################\n" 
        << "\n\t\tErrors during MPI-Initialisation!\n"
        <<"\n\t########################################################\n" ;
      std::string msg = "ERROR!\n Uh OH! Caught an non standard Exception in initial test Parsing of the Parameters!\n Hit the developers!\n\n";
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

  if (rank == 0) {  //MASTER
        //nice messages
        std::cout << "\n==================================================\n\tStart REPLICA EXCHANGE SIMULATION:\n\n==================================================\n";
        std::cout << "numreplicas:\t "<< numReplicas<<"\n";
        std::cout << "num Slaves:\t "<< numReplicas-1<<"\n";
        std::cerr << "\n==================================================\n\tStart REPLICA EXCHANGE SIMULATION:\n\n==================================================\n";
        std::cerr << "numreplicas:\t "<< numReplicas<<"\n";
        std::cerr << "num Slaves:\t "<< numReplicas-1<<"\n";
    DEBUG(1, "Master \t "<< rank)
    
    // Select repex Implementation - Polymorphism
    util::replica_exchange_master* Master;
    try{
        if(reedsSim){
            DEBUG(1, "Master_eds \t Constructor")
            Master = new util::replica_exchange_master_eds(args, cont, rank, size, numReplicas, repIDs[rank], repMap);
          } else{
            DEBUG(1, "Master \t Constructor")
            Master = new util::replica_exchange_master(args, cont, rank, size, numReplicas, repIDs[rank], repMap);
          }
    }
    catch(...){
        std::cerr <<"\n\t########################################################\n" 
            << "\n\t\tErrors during Master Class Init!\n"
            <<"\n\t########################################################\n" ;
        std::cerr << "\nREPEX:       Failed in initialization of Replica rank: "<<rank<<"" << std::endl;
        std::cout << "\nREPEX:       Failed in initialization of Replica rank: "<<rank<<"" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, E_USAGE);
        return 1;
    }
        
    DEBUG(1, "Master \t INIT")
    Master->init();
    Master->init_repOut_stat_file();
    
    //do md:
    unsigned int trial=0;
    DEBUG(1, "Master \t \t \t Equil: "<< equil_runs)
    for( ;trial<equil_runs; ++trial){    // for equilibrations
        Master->run_MD();
    }
    DEBUG(1, "Master \t \t MD: "<< total_runs)
    
    //Vars for timing
    int hh, mm, ss, eta_hh, eta_mm, eta_ss = 0;
    double eta_spent, percent, spent = 0.0;
    trial=0;    //reset trials
    for ( ; trial < sim_runs; ++trial){ //for repex execution
        DEBUG(2, "Master "<< rank <<" \t MD trial: "<< trial << "\n")\
        DEBUG(2, "Master "<< rank <<" \t run_MD START "<<trial<<"\n")  
        Master->run_MD();
        DEBUG(2, "Master "<< rank <<" \t run_MD DONE "<<trial<<"\n")  

        DEBUG(2, "Master " << rank << " \t swap START "<<trial<<"\n")    
        Master->swap();
        DEBUG(2, "Master "<< rank <<" \t run_MD DONE "<<trial<<"\n")  

        DEBUG(2, "Master " << rank << " \t receive START "<<trial<<"\n")    
        Master->receive_from_all_slaves();
        DEBUG(2, "Master " << rank << " \t write START "<<trial<<"\n")    
        Master->write();
        
        if ((total_runs/ 10 > 0) && (trial % (total_runs / 10) == 0)){  //Timer 
          percent = double(trial)/double(total_runs);
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
    DEBUG(1, "Master \t \t finalize ")
            
    Master->write_final_conf();
       
  } else {  //SLAVES    
    DEBUG(1, "Slave " << rank)    
            
    // Select repex Implementation - Polymorphism
    util::replica_exchange_slave* Slave;     
    if(reedsSim){
       Slave = new util::replica_exchange_slave_eds(args, cont, rank, repIDs[rank], repMap);
    } else{
       Slave = new util::replica_exchange_slave(args, cont, rank, repIDs[rank], repMap);
    }
    
    DEBUG(1, "Slave "<< rank <<" \t INIT")    
    Slave->init();
    
    //do md:
    unsigned int trial=0;
    DEBUG(1, "Slave "<< rank <<" \t EQUIL "<< equil_runs << " steps")    
    for( ;trial<equil_runs; ++trial){    // for equilibrations
        Slave->run_MD();
    }
    
    DEBUG(1, "Slave "<< rank <<" \t MD "<< total_runs << " steps")    
    for ( ; trial < total_runs; ++trial){ //for repex execution
      DEBUG(2, "Slave "<< rank <<" \t MD trial: "<< trial << "\n")    
      DEBUG(2, "Slave "<< rank <<" \t run_MD START "<<trial<<"\n")    
      Slave->run_MD();
      DEBUG(2, "Slave "<< rank <<" \t swap START "<<trial<<"\n")    
      Slave->swap();
      DEBUG(2, "Slave "<< rank <<" \t send START "<<trial<<"\n")    
      Slave->send_to_master();
    }
    
    DEBUG(1, "Slave "<< rank <<" \t Finalize") ;   
    Slave->write_final_conf();
    std::cout << "\n=================== Slave Node "<< rank << "  finished successfully!\n";
  }
  
  MPI_Barrier(MPI_COMM_WORLD);  //Make sure all processes finished.

  //last MASTER OUTPUT
  if(rank == 0){
    //FINAL OUTPUT - Time used:
    double end = MPI_Wtime();
    double duration = end - start;
    int durationMin = int(duration/60);
    int durationHour = int(std::floor(durationMin/60));
    int durationMinlHour = int(std::floor(durationMin-durationHour*60));
    int durationSlMin = int(std::floor(duration - (durationMinlHour+durationHour*60)*60));

    //Todo: if (cond){finished succ}else{not} bschroed
    std::string msg("\n====================================================================================================\n\tREPLICA EXCHANGE SIMULATION finished!\n====================================================================================================\n");
    std::cout << msg;
    std::cerr << msg;
 
    std::cout<< "TOTAL TIME USED: \n\th:min:s\t\tseconds\n"
       << "\t" << durationHour << ":"<<durationMinlHour << ":" << durationSlMin << "\t\t" << duration << "\n";
    std::cerr<< "TOTAL TIME USED: \n\th:min:s\t\tseconds\n"
       << "\t" << durationHour << ":"<<durationMinlHour << ":" << durationSlMin << "\t\t" << duration << "\n";

  }

  MPI_Finalize();
  return 0;
#else
  std::cout << "\n==================================================\n\tGROMOS Replica Exchange:\n==================================================\n";
  std::cout << argv[0] << " needs MPI to run\n\tuse --enable-mpi at configure and appropriate compilers\n" << std::endl;
  return 1;
#endif
}

