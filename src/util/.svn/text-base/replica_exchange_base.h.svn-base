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

#include <util/repex_mpi.h>
#include <util/replica.h>
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
  class replica_exchange_base {
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
    replica_exchange_base(io::Argument _args, int cont, int rank, std::vector<int> repIDs, 
                          std::map<ID_t, rank_t>& _repMap);
    /**
     * Destructor
     */
    ~replica_exchange_base();

    /**
     * runs MD simulation for all replicas; one by one
     */
    void run_MD();
    
    /**
     * write coordinates for all replicas to cnf
     */
    void write_final_conf();
    
    /**
     * init MD simulation for all replicas; one by one
     */
    void init();
    /**
     * prints out configuration to a file named <name>_<ID>.cnf
     * @param name string, name of output file
     */
    void print_coords(std::string name);
    /**
     * Tries a swapping of configuration if possible. Calculates energies, probabilities
     * and sends information via MPI communication if necessary.
     */
    void swap();

  protected:
    /**
     * Swapping routine if the replicas are on the same node, no MPI communication needed. Always called from swap().
     * @param it repIterator, iterator to replica with lower ID
     * @param partner integer, ID of replica with higher ID
     */
    void swap_on_node(repIterator it, const unsigned int partner);
    /**
     * switches configuration information if replicas are on same node.
     * @param it repIterator, iterator to replica with lower ID
     * @param partner integer, ID of replica with higher ID
     */
    void switch_coords_on_node(repIterator it, const unsigned int partner);
    /**
     * calculates switching probability of two replicas if they are on same node
     * @param rep1 replica*, pointer to replica with lower ID
     * @param rep2 replica*, pointer to replica with higher ID
     * @return double, probability in [0.0,1.0]
     */
    double calc_probability(replica * rep1, replica * rep2);
    /**
     * input parameters
     */
    io::Argument args;
    /**
     * mapping of IDs and rank to know where to send data to
     */
    const std::map<ID_t, rank_t> repMap;
    /**
     * number of replicas in the system
     */
    const unsigned int numReplicas;
    /**
     * all replicas on this node
     */
    std::vector<util::replica *> replicas;
    /**
     * the random number generator
     */
    math::RandomGeneratorGSL rng;


  };
}

#endif	/* REPLICA_EXCHANGE_BASE_H */

