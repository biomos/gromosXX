/**
 * @file repex_mpi.cc
 * the main md program for replica exchange simulations using the message passing interface
 */
/**
 * @page programs Program Documentation
 *
 * @anchor repex_mpi
 * @section repex_mpi replica exchange
 * @date 13.05.2011
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
#include <string>
#include <sstream>
#include <util/error.h>

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
    std::cout << "START REPEXMPI\n";  
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

  io::Argument args;

  if (args.parse(argc, argv, knowns)) {
    std::cerr << usage << std::endl;
    MPI_Abort(MPI_COMM_WORLD, E_USAGE);
    return 1;
  }

  if (args.count("version") >= 0) {
    MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
    return 0;
  }

  // parse the verbosity flag and set debug levels
  if (util::parse_verbosity(args)) {
    std::cerr << "could not parse verbosity argument" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
    return 1;
  }
  
  bool reedsSim;
  unsigned int equil_runs;
  unsigned int total_runs;
  unsigned int numAtoms;
  unsigned int numReplicas;
  unsigned int numEDSstates;

  int cont;
  
  {
    topology::Topology topo;
    configuration::Configuration conf;
    algorithm::Algorithm_Sequence md;
    simulation::Simulation sim;
    // TODO Bschroed: add nice starting message for rank 0
    // read in parameters
    
    bool quiet = true;
    
    if(rank == 0){
        quiet = false;
    }
    
    if (io::read_parameter(args,sim,std::cout, quiet)){
      if (rank == 0) {
        io::messages.display(std::cout);
        std::cout << "\nErrors in in_parameters!\n" << std::endl;
      }
      return -1;
    }

    // make a copy, don't change the original args
    io::Argument args2(args);
    
    // look for modified coordinate file if continuation run
    if(sim.param().replica.cont == 1){
      std::multimap< std::string, std::string >::iterator it = args2.lower_bound(("conf"));
      size_t pos = (*it).second.find_last_of(".");
      std::stringstream tmp;
      tmp << "_" << 1;
      (*it).second.insert(pos, tmp.str());
    }
    
    // read in the rest
    if(io::read_input_repex(args2, topo, conf, sim, md, rank, std::cout, quiet)){
        std::cerr << "\nErrors during initialization!\n" << std::endl;
        io::messages.display(std::cout);
        MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
        return 1;
    }
    
    if (io::check_parameter(sim) != 0){
        io::messages.display();
        return -1; //reactivated check param at end.  
    }

    cont = sim.param().replica.cont;
    equil_runs = sim.param().replica.equilibrate;
    total_runs = sim.param().replica.trials + equil_runs;
    numReplicas = sim.param().replica.num_T * sim.param().replica.num_l;
    numAtoms = topo.num_atoms();
    reedsSim = sim.param().reeds.reeds;
    if(reedsSim){
        numEDSstates=sim.param().reeds.eds_para[0].numstates;//BSchr: Assume: eds.numstates is constant over all replicas;
    }else{
        numEDSstates=0;
    }
    //Todo bschroed: Nice Messaging
    if(rank == 0){
        std::cout << "done Reading input\n";
        std::cout.flush();
    }
  }
  
  io::messages.clear();

  //////////////////////////
  // defining MPI Datatypes
  //////////////////////////

  // Vector
  MPI_Datatype MPI_VEC;
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
  
  assert(numReplicas > 0);

  // where is which replica
  std::map<unsigned int, unsigned int> repMap;

  // every node gets one element of that vector
  std::vector< std::vector<int> > repIDs;
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

  // make sure all nodes have initialized everything
  MPI_Barrier(MPI_COMM_WORLD);
  
  //TODO integrate reeds with system.param().reeds.reeds
    
  if(rank == 0){
    std::cout << "Went trough MPI fun and parsing!\n"; //todo bschroed: remove!
    std::cout.flush();//todo bschroed: remove!
  }  
  //////////////////////////////
  /// Starting master-slave mode
  //////////////////////////////
  if (rank == 0) {  //MASTER
    //print Initial Master text:
    std::cout  << "\n==================================================\n"
               << "GROMOS REEDS:"
               <<"\n==================================================\n"
               << "Start Master on: " << "Node " << rank << std::endl
               << "numreplicas:\t "<< numReplicas<<std::endl
               << "num Slaves:\t "<< numReplicas-1<<std::endl
               << "reeds:\t "<< reedsSim<<std::endl<<std::endl;
    
    // Select repex Implementation
    util::replica_exchange_master* Master;
    if(reedsSim){
        std::cout <<  "Master REEDS " << rank << std::endl;
        std::cout.flush();
        Master = new util::replica_exchange_master_eds(args, cont, rank, size, numReplicas, repIDs[rank], repMap);

      } else{
        Master = new util::replica_exchange_master(args, cont, rank, size, numReplicas, repIDs[rank], repMap);
      }
    
    Master->init();
    
    //do md:
    std::cout << "\n==================================================\n"
              << " MAIN MD LOOP\n"
              << "==================================================\n\n";
    unsigned int trial;
    for( ;trial<equil_runs; ++trial){    // for equilibrations
        Master->run_MD();
    }
    for ( ; trial < total_runs; ++trial){ //for repex execution
      Master->run_MD();
      Master->swap();
      Master->receive_from_all_slaves();
      Master->write();
    }
    Master->write_final_conf();
    
    //FINAL OUTPUT - Time used:
    double end = MPI_Wtime();
    double duration = end - start;
    double durationMin = duration/60;
    double durationHour = std::floor(durationMin/60);
    double durationMinlHour = std::floor(durationMin-durationHour*60);
    double durationSlMin = std::floor(duration - (durationMinlHour+durationHour*60)*60);
    
    std::cout << "\n==================================================\n"
              << "REPLICA EXCHANGE SIMULATION finished successfully! " << "Node " << rank << " - MASTER\n"
              << "\n==================================================\n"
              << "TOTAL TIME USED: \n\th:min:s\t\tseconds\n"
              << "\t" << durationHour << ":"<<durationMinlHour << ":" << durationSlMin << "\t\t" << duration << "\n";
    MPI_Finalize();

  } else {  //SLAVES
      
        std::cout << "Slave here: "<<rank<<std::endl;      
        
        // Select repex Implementation
        util::replica_exchange_slave* Slave;     
        if(reedsSim){
           std::cout <<  "Slave REEDS " << rank << std::endl;
           Slave = new util::replica_exchange_slave_eds(args, cont, rank, repIDs[rank], repMap);
         } else{
           Slave = new util::replica_exchange_slave(args, cont, rank, repIDs[rank], repMap);
         }

    Slave->init();
    //do md:
    unsigned int trial;
    for( ;trial<equil_runs; ++trial){    // for equilibrations
        Slave->run_MD();
    }
    for ( ; trial < total_runs; ++trial){ //for repex execution
      Slave->run_MD();
      Slave->swap();
      Slave->send_to_master();
    }
    
    Slave->write_final_conf();
    std::cout << "\n=================== REPLICA EXCHANGE SIMULATION finished successfully! " << "Node " << rank << "\n";
  }
  return 0;
  
#else
  std::cout << argv[0] << " needs MPI to run\n\tuse --enable-mpi at configure and appropriate compilers\n" << std::endl;
  return 1;
#endif
}

