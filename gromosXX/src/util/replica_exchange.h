/**
 * @file replica_exchange.h
 * replica exchange
 */

#ifndef INCLUDED_REPLICA_EXCHANGE_H
#define INCLUDED_REPLICA_EXCHANGE_H

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

    /**
     * readblock
     */
    ssize_t readblock(char * source, ssize_t size);

    /**
     * writeblock
     */
    ssize_t writeblock(char * dest, ssize_t size);

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
     * do multigraining
     */
    bool multigraining;

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
		     simulation::Simulation & sim,
		     topology::Topology & cg_topo,
		     simulation::Simulation & cg_sim);
    /**
     * run the replica
     */
    int run_md(topology::Topology & topo,
	       configuration::Configuration & conf,
	       simulation::Simulation & sim,
	       algorithm::Algorithm_Sequence & md,
	       topology::Topology & cg_topo,
	       configuration::Configuration & cg_conf,
	       simulation::Simulation & cg_sim,
	       interaction::Forcefield * cg_ff,
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
    
  };

  inline ssize_t Replica_Exchange::readblock(char * source, ssize_t size)
  {
    ssize_t current;
    ssize_t window = 1024;
    while(size > 0){
      if (size > window)
	current = window;
      else current = size;
      
      ssize_t count;
      if ((count = read(cl_socket, source, current)) == 0){
	std::cerr << "received zero bytes instead of "
		  << current << " !!!" << std::endl;
	throw std::runtime_error("could not read data block");
      }
      
      if (current != count){
	std::cerr << "received only " << count << " bytes..." << std::endl;
	// make window smaller...
	window = count;
      }

      char c = 0;
      if (write(cl_socket, &c, 1) != 1){
	std::cerr << "sending ACK failed" << std::endl;
	throw std::runtime_error("could not send ACK");
      }

      source += count;
      size -= count;
    }
    return 0;
  }
  
  inline ssize_t Replica_Exchange::writeblock(char * dest, ssize_t size)
  {
    ssize_t current;
    ssize_t window = 1024;
    while(size > 0){
      if (size > window)
	current = window;
      else current = size;
      
      ssize_t count;
      if ((count = write(cl_socket, dest, current)) == 0){
	std::cerr << "could not write a single byte!\n"
		  << "tried to send " << current << std::endl;
	throw std::runtime_error("could not write data block");
      }

      if (current != count){
	std::cerr << "sent only " << count << " bytes..." << std::endl;
	// make window smaller...
	window = count;
      }
      
      char c = 0;
      if (read(cl_socket, &c, 1) != 1){
	std::cerr << "getting ACK failed" << std::endl;
	throw std::runtime_error("could not read ACK");
      }
      if (c != 0){
	std::cerr << "wrong ACK received" << std::endl;
	throw std::runtime_error("wrong ACK received");
      }
      
      dest += count;
      size -= count;
    }
    return 0;
  }
  
} // util

#endif
