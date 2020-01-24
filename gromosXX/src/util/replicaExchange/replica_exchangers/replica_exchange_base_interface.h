/**
 * @file replica_exchange_base.h

 */

#ifndef REPLICA_EXCHANGE_BASE_H
#define	REPLICA_EXCHANGE_BASE_H


#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>


#include <util/replicaExchange/repex_mpi.h>
#include <util/replicaExchange/replica/_replica_Interface.h>
#include <string>
#include <math/random.h>

#ifdef XXMPI
#include <mpi.h>
#endif

namespace util {

  /**
   * @class replica_exchange_base
   * One instance of this class per node managing one or more replicas. Master and
   * slave are derived from this class.
   */
  class replica_exchange_base_interface {
  private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica_exchange_base_interface(const replica_exchange_base_interface &);
  public:
    /**
     * Constructor
     * @param _args io::Argument, passed on to Replica
     * @param rank integer, rank of node
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_base_interface(io::Argument _args, 
                          unsigned int cont, 
                          unsigned int globalThreadID, 
                          std::vector<std::vector<unsigned int> > replica_owned_threads, 
                          std::map<ID_t, rank_t>& thread_id_replica_map);
    /**
     * Destructor
     */
    ~replica_exchange_base_interface();

    //Simulation functions
    /**
     * runs MD simulation for all replicas; one by one
     */
    virtual void run_MD();
    
    /**
     * write coordinates for all replicas to cnf
     */
    virtual void write_final_conf();
    
    /**
     * init MD simulation for all replicas; one by one
     */
    virtual void init();
    /**
     * prints out configuration to a file named \<name\>_\<ID\>.cnf
     * @param name string, name of output file
     */
    void print_coords(std::string name);
    /**
     * Tries a swapping of configuration if possible. Calculates energies, probabilities
     * and sends information via MPI communication if necessary.
     */
    virtual void swap();
    
    /**
     *  GET :
     */
    unsigned int get_stepsPerRun(){
        return stepsPerRun;
    }
  protected:
    /*ATTRIBUTES*/
      //General
    /**
     * input parameters of gromos
     */
    io::Argument args;
    /**
     * the random number generator
     */
    math::RandomGeneratorGSL rng;
    /**
     * continuation? of this class
     */
    unsigned int cont;
    
    //MPI
    //Global MPI
    /**
     * Thead ID of this class in the global context (= MPI-Rank)
     */
    unsigned int globalThreadID;
    /**
     * ID of the Master Thread for this RE-Graph
     */
    unsigned int globalMasterThreadID;
    /**
     *  simulation ID - to which simulation does this thread belong?
     */
    unsigned int numReplicas;

    
    ////Simulation_MPI
    /**
     *  simulation ID - to which simulation does this thread belong?
     */
    unsigned int simulationID;
        /**
     *  simulation ID - to which simulation does this thread belong?
     */
    int simulationMasterID;
    /**
     * How many threads per simulation?
     */
    unsigned int simulationThreads;
    /**
     *  simulation Thread ID - which ID has the thread in the simulation X
     */
    unsigned int simulationThreadID;
    
    ////Helper Constructs
    /**
     * replica IDs
     */
    std::vector<std::vector<unsigned int> > replica_owned_threads;
    /**
     * mapping of IDs and rank to know where to send data to
     */
    const std::map<ID_t, rank_t> thread_id_replica_map;
    
    //REplica Exchange:
    ////Me and my partner
    /**
     * potential energy using current Hamiltonian
     */
    double epot;
    /**
     * ID of partner Replica and master ThreadID
     */
    unsigned int partnerReplicaID;
    unsigned int partnerReplicaMasterThreadID;

    /**
     * potential energy of the partner Hamiltonian (for bookkeeping)
     */
    double epot_partner;
    
    //Exchange?
    /**
     * probability of last switch
     */
    double probability;
    /**
     * switched last time?
     */
    bool switched;
    
    //REPLICA
    /**
     * simulating unit of this Thread
     */
    util::replica_Interface *replica;

    //Simulation Progress
    /**
     * run number of replica
     */
    unsigned int run;

    /**
     * total number of runs per replica
     */
    unsigned int total_runs;

    /**
     * total number of steps simulated
     */
    unsigned int steps;
    /**
     * maximal number of steps to simulate per run
     */
    unsigned int stepsPerRun;
    /**
     * current simulation time
     */
    double time;
        
    /*
     * FUNCTIONS
     */
    //Initialisation Function
    //init Replicas - used in contstructor, initialises the replica objs.
    /**
     * initialise the "replica" in this class with this function
     * @param cont int, continuation of this class?
     * @param rank int, global Thread ID - later remove
     * @return void
     */
    virtual void createReplicas(int cont,  int globalThreadID);
    
    //Replica Exchange Functions
    ////SWAP FUNCTIONS
    /**
     * TODO: BSCHROED - REWRITE FOR PARALLEL REPLICAS
    * Finds partner for current switch
    * @return ID of partner, own ID if no switching in current trial
    */
    virtual int find_partner() const = 0;
    /**
    * finds out if configurations are to be switched; sets switched to true if so
    */
    //virtual void swap_coordinates(const unsigned int partnerReplicaID);
    
    /*
     *TODO: get betteer name!
     */
    virtual void swap_replicas_priv(const unsigned int partnerReplicaID)=0;
    
    /**
     * calculates potential energy for current configuration with current lambda
     */
    virtual double calculate_energy_core() = 0; 
    /**
     * calculates potential energy of current configuration with lambda(Hamiltonian) of partner with the partnerThreadID
     * @param partner ID of partner
     * @return potential energy of configuration with lambda(Hamiltonian) of partner
     */
    virtual double calculate_energy(const unsigned int partnerReplicaID) = 0;
       
    /**
    * switch back the averages
    */
    virtual void exchange_averages();
 
    /**
     * TODO: bring into Exchanger 
     * calculates probability of switch with current partner, may involve MPI communication
     * @param partnerID
     * @param partnerRank
     * @return probability
     */
    virtual double calc_probability(const unsigned int partnerReplicaID) = 0;

    //SWAP COORDINATES FUNCTIONS
    /**TODO: bring into Exchanger 
    * Initiates MPI communication to receive new configuration information
    * @param senderID
    * @param senderRank
    */
   virtual void receive_new_coord(const unsigned int senderReplicaID);
    /**TODO: bring into Exchanger 
     * Initiates MPI communication to send new configuration information
     * @param receiverID
     * @param senderRank
     */
    virtual void send_coord(const unsigned int receiverReplicaID);
    
    /**TODO: bring into Exchanger 
     * Prints information of replica to std::cout for debugging purposes
     */
    virtual void print_info(std::string bla)const = 0;
    
    /**
     * Scales the velocities after an exchange in temperature replica exchange
     */
     virtual void velscale(unsigned int partnerReplica) = 0;

     //setters:
    ////Exchange Param adaptation
     /*
      */
    void updateReplica_params();
  };
}

#endif	/* REPLICA_EXCHANGE_BASE_H */

