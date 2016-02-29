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
#include <util/repex_mpi.h>
#include <util/replica_exchange_base.h>
#include <util/replica.h>
#include <util/replica_data.h>

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
  class replica_exchange_master : public replica_exchange_base {
  public:
    /**
     * constructor
     * @param args io::Argument, passed on to replica
     * @param rank integer, rank of node
     * @param _size size of mpi comm world
     * @param _numReplicas total number of replicas
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_master(io::Argument & args,
            int cont,
            int rank,
            int _size,
            int _numReplicas,
            std::vector<int> repIDs,
            std::map<ID_t, rank_t> & repMap);
    /**
     * destructor
     */
    virtual ~replica_exchange_master();
    /**
     * receives all information written to output file from the slaves
     */
    void receive_from_all_slaves();
    /**
     * writes data to output file @repdat
     */
    void write();

  private:
    /**
     * output file stream for output file
     */
    std::ofstream repOut;
    /**
     *  comm world size; number of processors available
     */
    const unsigned int size;
    /**
     * total number of replicas in system
     */
    const unsigned int numReplicas;
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

  };
}
#endif	/* REPLICA_EXCHANGE_MASTER_H */

