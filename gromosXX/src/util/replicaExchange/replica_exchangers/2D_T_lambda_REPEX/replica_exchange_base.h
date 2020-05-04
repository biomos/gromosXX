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
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
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
  class replica_exchange_base : public virtual replica_exchange_base_interface {
  private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica_exchange_base(const replica_exchange_base &);
  public:
    /**
     * Constructor
     * @param _args io::Argument, passed on to Replica
     * @param rank integer, rank of node
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_base(io::Argument _args, 
                          unsigned int cont, 
                          unsigned int globalThreadID, 
                          replica_graph_mpi_control replicaGraphMPIControl,
                          simulation::MpiControl replica_mpi_control);
    /**
     * Destructor
     */
    ~replica_exchange_base();
    
    /**
     * Temperature of replica
     */
    double T;
    /**
     * Lambda of replica
     */
    double l;
    
  protected:
       /**
       * Setting RE-Param
       */
      void setParams() override;
  };
}

#endif	/* REPLICA_EXCHANGE_BASE_H */

