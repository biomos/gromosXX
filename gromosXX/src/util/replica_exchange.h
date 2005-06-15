/**
 * @file replica_exchange.h
 * replica exchange
 */

#ifndef INCLUDED_REPLICA_EXCHANGE_H
#define INCLUDED_REPLICA_EXCHANGE_H

#ifdef XXMPI
#include <mpi.h>
#endif

#include <util/replica_data.h>

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

#ifdef XXMPI

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
    Replica_Exchange() {}
    /**
     * Destructor
     */
    virtual ~Replica_Exchange() {}
    
    /**
     * run
     */
    virtual int run(io::Argument & args) = 0;
    
  protected:

  };
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Master  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

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
     * Destructor
     */
    virtual ~Replica_Exchange_Master()
    {
    }

    /**
     * run the thread
     */
    virtual int run(io::Argument & args);
    
    /**
     * configuration accessor
     */
    configuration::Configuration & conf(int i)
    {
      assert(i>=0 && unsigned(i) < m_conf.size());
      return m_conf[i];
    }
    
  private:
    /**
     * try a switch between i and j
     */
    int switch_replica(int i);

    /**
     * information of all replicas
     * gets updated from the slaves
     */
    std::vector<Replica_Data> replica_data;
    
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Slave  /////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @class Replica_Exchange_Slave
   * replica exchange slave
   */
  class Replica_Exchange_Slave : public Replica_Exchange
  {
  public:
    /**
     * Constructor
     */
    Replica_Exchange_Slave();
    /**
     * Destructor
     */
    virtual ~Replica_Exchange_Slave() {}
    
    /**
     * run the slave
     */
    virtual int run(io::Argument & args);
    
  private:
    /**
     * replica information
     */
    Replica_Data replica_data;

    /**
     * communicator
     */
    MPI_Comm master;

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

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Control  ///////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @class Replica_Exchange_Control
   * replica exchange interactive module
   */
  class Replica_Exchange_Control : public Replica_Exchange
  {
  public:
    /**
     * Constructor
     */
    Replica_Exchange_Control();
    /**
     * Destructor
     */
    virtual ~Replica_Exchange_Control() {}
    
    /**
     * run the slave
     */
    virtual int run(io::Argument & args);
    
  private:
    /**
     * replica information
     */
    std::vector<Replica_Data> replica_data;
    
    /**
     * communicator
     */
    MPI_Comm master;

  };

#endif
  
} // util

#endif
