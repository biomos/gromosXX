/**
 * @file replica_exchange_master.h
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
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>
#include <util/replicaExchange/repex_mpi.h>
#include <util/replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_base.h>
#include <util/replicaExchange/replica/replica.h>
#include <util/replicaExchange/replica_data.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#ifndef REPLICA_EXCHANGE_MASTER_H
#define	REPLICA_EXCHANGE_MASTER_H

namespace util {

  /**
   * @class replica_exchange_master
   * Additionally to replica_exchange_base: receives and writes data to file.
   */
  class replica_exchange_master : public virtual replica_exchange_base {
  public:
    /**
     * constructor
     * @param args io::Argument, passed on to replica
     * @param rank integer, rank of node
     * @param _size size of mpi comm world
     * @param _numReplicas total number of replicas
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_master(io::Argument & args,
            unsigned int cont,
            unsigned int globalThreadID,
            replica_graph_mpi_control replicaGraphMPIControl,
            simulation::mpi_control_struct replica_mpi_control);
    /**
     * Destructor
     */
    virtual ~replica_exchange_master();
    
    //Simulation functions
    /**
     * receives all information written to output file from the slaves
     */
    virtual void receive_from_all_slaves();
    /**
     * writes data to output file \@repdat
     */
    virtual void write();
    
 
    virtual void init_repOut_stat_file();
    
  protected:


    /**
     * output file stream for repdat output file
     */
    std::ofstream repOut;
    
    /*
     * global Parameters for replica exchange simulation
     * int num_T;
     * int num_l;
     * std::vector<double> temperature;
     * bool scale;
     * std::vector<double> lambda;
     * std::vector<double> dt;
     * int trials;
     * int equilibrate;
     * int slave_runs;
     * int write;
     */
    const simulation::Parameter::replica_struct& repParams;

    /**
     * information of all replicas
     */
    std::vector<util::replica_data> replicaData;
    /**
     * output file Path for repdat output file
     */
    std::string repdatName;
    //virtual void init_repOut_stat_file(std::string repoutPath);

  };
}
#endif	/* REPLICA_EXCHANGE_MASTER_H */

