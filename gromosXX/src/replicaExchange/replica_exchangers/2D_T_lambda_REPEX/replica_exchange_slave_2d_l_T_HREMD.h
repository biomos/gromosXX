/**
 * @file replica_exchange_slave_2d_l_T_HREMD.h
 * contains replica_exchange_slave_2d_l_T_HREMD class
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

#include <replicaExchange/replica/replica.h>
#include <replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_base_2d_l_T_HREMD.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#ifndef REPLICA_EXCHANGE_SLAVE_H
#define	REPLICA_EXCHANGE_SLAVE_H

namespace re {

  /**
   * @class replica_exchange_slave_2d_l_T_HREMD
   * 
   */
  class replica_exchange_slave_2d_l_T_HREMD : public virtual replica_exchange_base_2d_l_T_HREMD, public virtual replica_exchange_slave_interface {
  public:
    /**
     * constructor
     * @param _args io::Argument, passed on to Replica
     * @param rank integer, rank of node
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_slave_2d_l_T_HREMD(io::Argument & _args,
                                        unsigned int cont,
                                        unsigned int globalThreadID,
                                        replica_graph_control & replicaGraphMPIControl,
                                        simulation::MpiControl & replica_mpi_control);
    /**
     * destructor
     */
    virtual ~replica_exchange_slave_2d_l_T_HREMD();

  };
}
#endif	/* REPLICA_EXCHANGE_SLAVE_H */

