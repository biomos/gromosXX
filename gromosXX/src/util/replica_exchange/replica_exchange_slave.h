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
#include <util/replica_exchange_base.h>
#include <util/repex_mpi.h>
#include <util/replica.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#ifndef REPLICA_EXCHANGE_SLAVE_H
#define	REPLICA_EXCHANGE_SLAVE_H

namespace util {

  /**
   * @class replica_exchange_slave
   * 
   */
  class replica_exchange_slave : public replica_exchange_base {
  public:
    /**
     * constructor
     * @param _args io::Argument, passed on to Replica
     * @param rank integer, rank of node
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_slave(io::Argument & _args,
            int cont,
            int rank,
            std::vector<int> repIDs,
            std::map<ID_t, rank_t> & repMap);
    /**
     * destructor
     */
    virtual ~replica_exchange_slave();
    /**
     * sends information of all replicas to master
     */
    void send_to_master() const;
  private:
    /**
     * copy constructor
     * not allowed (yet)
     */
    replica_exchange_slave(const replica_exchange_slave &);

  };
}
#endif	/* REPLICA_EXCHANGE_SLAVE_H */

