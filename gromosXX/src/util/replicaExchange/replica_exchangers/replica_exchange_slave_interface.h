/**
 * @file replica_exchange_slave.h
 * contains replica_exchange_slave class
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
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <util/replicaExchange/replica/replica.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#ifndef REPLICA_EXCHANGE_SLAVE_INTERFACE_H
#define	REPLICA_EXCHANGE_SLAVE_INTERFACE_H

namespace util {

  /**
   * @class replica_exchange_slave
   * 
   */
  class replica_exchange_slave_interface : public virtual replica_exchange_base_interface {
  public:
    /**
     * constructor
     * @param _args io::Argument, passed on to Replica
     * @param rank integer, rank of node
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_slave_interface(io::Argument & _args,
            unsigned int cont,
            unsigned int globalThreadID,
            replica_graph_mpi_control replicaGraphMPIControl,
            simulation::MpiControl replica_mpi_control);
    /**
     * destructor
     */
    virtual ~replica_exchange_slave_interface();
    /**
     * sends information of all replicas to master
     */
    virtual void send_to_master() const;
  protected:
    /**
     * copy constructor
     * not allowed (yet)
     */
    replica_exchange_slave_interface(const replica_exchange_slave_interface &);
    
    /**
     *  Sending data or not?
     */
    bool not_sender;

  };
}
#endif	/* REPLICA_EXCHANGE_SLAVE_H */

