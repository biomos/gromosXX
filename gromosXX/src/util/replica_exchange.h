/**
 * @file replica_exchange.h
 * replica exchange
 */

#ifndef INCLUDED_REPLICA_EXCHANGE_H
#define INCLUDED_REPLICA_EXCHANGE_H

class gsl_rng;

namespace topology
{
  class Topology;
}

namespace simulation
{
  class Simulation;
}

namespace configuration
{
  class Configuration;
}

namespace io
{
  class Argument;
}

namespace util
{
  /**
   * @class Replica_Exchange
   * replica exchange
   */
  class Replica_Exchange
  {
  public:
    /**
     * Constructor
     */
    Replica_Exchange();

    enum state_enum{ waiting=0, ready=1, running=2, terminate=4 };

    /**
     * @struct Replica_Data
     * replica information
     */
    struct Replica_Data
    {
      double     temperature;
      double     lambda;
      double     energy;
      int        switch_replica;
      double     switch_temperature;
      double     switch_lambda;
      double     switch_energy;
      int        TID;
      int        run;
      state_enum state;
    };

    /**
     * @struct Slave_Data
     * slave information
     */
    struct Slave_Data
    {
      /**
       * Default Constructor
       */
      Slave_Data() : state(waiting), replica(-1) 
      {
      }
      /**
       * Constructor
       */
      Slave_Data(state_enum state, int replica) 
	: state(state), replica(replica) 
      {
      }
      
      /**
       * state of the slave
       */
      state_enum state;
      /**
       * assigned replica
       */
      int        replica;
    };

    /**
     * run
     */
    virtual int run(io::Argument & args, int tid, int num_threads) = 0;
    
  protected:

    /**
     * thread ID
     */
    int m_ID;

  };
  
  /**
   * @class Replica_Exchange_Master
   * replica exchange master
   */
  class Replica_Exchange_Master : public Replica_Exchange
  {
  public:
    /**
     * Constructor
     */
    Replica_Exchange_Master();

    /**
     * information of all replicas
     * gets updated from the slaves
     */
    std::vector<Replica_Data> replica_data;
    
    /**
     * information of all slaves
     */
    std::vector<Slave_Data>   slave_data;
    
    /**
     * neighbour list for each replica
     * this should get multi-dimensional later...
     */
    std::vector<int>        neighbour;
    /**
     * position in the neighbour list for each replica
     */
    std::vector<int>        neighbour_pos;

    /**
     * run the thread
     */
    virtual int run(io::Argument & args, int tid, int num_threads);
    
    /**
     * configuration accessor
     */
    configuration::Configuration & conf(int i)
    {
      assert(i < m_conf.size() && i >= 0);
      return m_conf[i];
    }
    
  private:
    /**
     * try a switch between i and j
     */
    int switch_replica(int i, int j);

    /**
     * random number generator
     */
    gsl_rng * m_rng;

    /**
     * (current) configuration of all replicas
     * this is necessary if there are more replicas than slaves
     * (very likely true)
     */
    std::vector<configuration::Configuration> m_conf;

  };

  /**
   * @class Replica_Exchange_Slave
   * replica exchange master
   */
  class Replica_Exchange_Slave : public Replica_Exchange
  {
  public:
    /**
     * Constructor
     */
    Replica_Exchange_Slave();

    /**
     * replica information
     */
    Replica_Data replica_data;
    
    /**
     * slave information
     */
    Slave_Data slave_data;
    
    /**
     * run the slave
     */
    virtual int run(io::Argument & args, int tid, int num_threads);
    
  private:

    /**
     * get thread state from master
     */
    int get_slave_data();
    /**
     * update thread state on master
     */
    int update_slave_data();
    /**
     * get replica data from master
     */
    int get_replica_data();
    /**
     * update replica data on master
     */
    int update_replica_data();
    /**
     * get configuration from master
     */
    int get_configuration(topology::Topology const & topo,
			  configuration::Configuration & conf);
    /**
     * update configuration on master
     */
    int update_configuration(topology::Topology const & topo,
			     configuration::Configuration & conf);
    /**
     * initialise run (coordinates, time, temperature, lambda)
     */
    int init_replica(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim);
    /**
     * run the replica
     */
    int run_md(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim,
	       algorithm::Algorithm_Sequence & md,
	       io::Out_Configuration & traj);

  };

  /**
   * replica exchange accessor
   * for OMP
   */
  extern Replica_Exchange_Master * replica_master;
  
} // util

#endif
