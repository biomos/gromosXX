/**
 * @file replica_exchange_base_2d_l_T_HREMD.h

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

#include <replicaExchange/repex_mpi.h>
#include <replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <replicaExchange/replica/_replica_Interface.h>
#include <string>
#include <math/random.h>

#ifdef XXMPI
#include <mpi.h>
#endif

namespace re {

  /**
   * @class replica_exchange_base_2d_l_T_HREMD
   * One instance of this class per node managing one or more replicas. Master and
   * slave are derived from this class.
   */
  class replica_exchange_base_2d_l_T_HREMD : public virtual replica_exchange_base_interface {
  private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica_exchange_base_2d_l_T_HREMD(const replica_exchange_base_2d_l_T_HREMD &);
  public:
    /**
     * Constructor
     * @param _args io::Argument, passed on to Replica
     * @param rank integer, rank of node
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param _repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_base_2d_l_T_HREMD(io::Argument _args,
                                       unsigned int cont,
                                       unsigned int globalThreadID,
                                       replica_graph_control & replicaGraphMPIControl,
                                       simulation::MpiControl & replica_mpi_control);
    /**
     * Destructor
     */
    ~replica_exchange_base_2d_l_T_HREMD();

    
  protected:
       /**
       * Setting RE-Param
       */
      void setParams() override;

        /**
         *
         * @param partnerReplicaID
         * @return
         */
      double calculate_energy(const unsigned int partnerReplicaID) override;

      void calculate_energy_helper(const unsigned int partnerReplicaID) override;


    /**
     * This function should be overriden in subclass
     */
     virtual void determine_switch_probabilities();

     int find_partner() const override;

     void swap_replicas_2D(const unsigned int partnerReplicaMasterThreadID);

     /**
      *    Specifics
      */
     /**
        * Sets lambda parameter to original value of replica
        */
     void set_lambda();
     /**
      * Sets temperature to original value of replica
      */
     void set_temp();
     /**
      * Sets lambda parameter to value of partner (for energy calculations)
      */
     void change_lambda(const unsigned int partnerReplicaID);
     /**
      * Sets temperature to value of partner (for energy calculations)
      */
     void change_temp(const unsigned int partnerReplicaID);


  };
}

#endif	/* REPLICA_EXCHANGE_BASE_H */

