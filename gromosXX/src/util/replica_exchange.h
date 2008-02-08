/**
 * @file replica_exchange.h
 * replica exchange
 */

#ifndef INCLUDED_REPLICA_EXCHANGE_H
#define INCLUDED_REPLICA_EXCHANGE_H

#include "replica_data.h"

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

struct addrinfo;

namespace util
{

#ifdef REPEX

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
    Replica_Exchange() : serv_socket(-1), cl_socket(-1), retry(5), timeout(60) {}
    /**
     * Destructor
     */
    virtual ~Replica_Exchange() {}
    
    /**
     * run
     */
    virtual int run(io::Argument & args) = 0;
    
  protected:

    static const bool quiet = true;

    /**
     * get replica data from network
     */
    int get_replica_data(Replica_Data & r, int &i);

    /**
     * get replica data from network
     */
    int get_replica_data(std::vector<Replica_Data> & r, int &i);

    /**
     * put replica data on network
     */
    int put_replica_data(Replica_Data & r);

    /**
     * get configuration from network
     */
    int get_configuration(configuration::Configuration & conf);

    /**
     * put configuration on network
     */
    int put_configuration(configuration::Configuration & conf);

    /**
     * readblock
     */
    ssize_t readblock(char * source, ssize_t size);

    /**
     * writeblock
     */
    ssize_t writeblock(char * dest, ssize_t size);

    /**
     * magic cookie exchange
     * @param master send first, receive after if true
     */
    bool magic_cookie(bool master);

    addrinfo * get_server(io::Argument & args, addrinfo & hints, std::string server_name);

    /**
     * unix network socket
     */
    int serv_socket;
    /**
     * client socket
     */
    int cl_socket;
    /**
     * connecting is re-attempted in case of failure
     */
    int retry;
    /**
     * timeout in [s]
     */
    int timeout;
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
      gsl_rng_free(m_rng);
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

    /**
     * coarse-grained configuration accessor
     */
    configuration::Configuration & cg_conf(int i)
    {
      assert(i>=0 && unsigned(i) < m_cg_conf.size());
      return m_cg_conf[i];
    }
    
  private:
    /**
     * select a job for a client
     */
    int select_job(simulation::Simulation & sim);
    
    /**
     * finish a job of a client
     */
    int finish_job(simulation::Simulation & sim);

    /**
     * interactive session
     */
    int interactive(simulation::Simulation & sim);
    
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

    /**
     * print replica stuff
     */
    void print_replica(int r,
		       simulation::Parameter const & param,
		       std::ostream & os);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Data  ///////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * do multigraining ?
     */
    bool multigraining;

    /**
     * information of all replicas
     */
    std::vector<Replica_Data> replica_data;
    
    /**
     * random number generator
     */
    gsl_rng * m_rng;

    /**
     * (current) (fine-grained) configuration of all replicas
     */
    std::vector<configuration::Configuration> m_conf;

    /**
     * (current) (coarse-grained) configuration of all replicas
     * (only used for multigraining)
     */
    std::vector<configuration::Configuration> m_cg_conf;

    /**
     * switch temperatures?
     */
    int switch_T;

    /**
     * switch lambdas?
     */
    int switch_l;
    
    /**
     * output
     */
    std::ofstream rep_out;
    
    /**
     * trials done
     */
    int trials;
    /**
     * runs completed
     */
    int runs;

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
     * do multigraining
     */
    bool multigraining;

    /**
     * initialise slave
     */
    int init
    (
     io::Argument & args,
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     algorithm::Algorithm_Sequence & md,
     io::Out_Configuration & traj,
     topology::Topology & cg_topo,
     configuration::Configuration & cg_conf,
     simulation::Simulation & cg_sim,
     algorithm::Algorithm_Sequence & cg_md,
     io::Out_Configuration & cg_traj
     );

    /**
     * initialise run (coordinates, time, temperature, lambda)
     */
    int init_replica
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     topology::Topology & cg_topo,
     configuration::Configuration & cg_conf,
     simulation::Simulation & cg_sim
     );
    /**
     * run the replica
     */
    int run_md
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     algorithm::Algorithm_Sequence & md,
     io::Out_Configuration & traj,
     topology::Topology & cg_topo,
     configuration::Configuration & cg_conf,
     simulation::Simulation & cg_sim,
     algorithm::Algorithm_Sequence & cg_md,
     io::Out_Configuration & cg_traj
     );

    /**
     * run the replica
     */
    int recalc_energy
    (
     topology::Topology & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     algorithm::Algorithm_Sequence & md,
     io::Out_Configuration & traj,
     topology::Topology & cg_topo,
     configuration::Configuration & cg_conf,
     simulation::Simulation & cg_sim,
     algorithm::Algorithm_Sequence & cg_md,
     io::Out_Configuration & cg_traj
     );

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
    
  };

#endif
  
} // util

#endif
