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
    virtual ~Replica_Exchange_Master() {}

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
     * try to switch replica i
     */
    int switch_replica(int i, simulation::Parameter const & param);

    /**
     * calculate switching probability
     */
    double switch_probability(int i, int j, simulation::Parameter const & param);
    
    /**
     * find the switch partner of replica i
     */
    int find_switch_partner(int i);
    
    /**
     * determine Tj and lj for replica i
     */
    void set_next_switch(int i);

    void print_replica(int r,
		       simulation::Parameter const & param,
		       std::ostream & os);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Data  ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * information of all replicas
     */
    std::vector<Replica_Data> replica_data;
    
    /**
     * random number generator
     */
    gsl_rng * m_rng;

    /**
     * (current) configuration of all replicas
     */
    std::vector<configuration::Configuration> m_conf;

    /**
     * switch temperatures?
     */
    int switch_T;

    /**
     * switch lambdas?
     */
    int switch_l;
    
    std::ofstream rep_out;

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
