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
   * replica exchange master
   * and state data
   */
  class Replica_Exchange
  {
  public:
    /**
     * Constructor
     */
    Replica_Exchange();

    enum state_enum{ waiting, ready, running, terminate };

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
    
    std::vector<Replica_Data> replica_data;
    
    std::vector<state_enum> thread_state;
    std::vector<int>        thread_replica;

    std::vector<int>        neighbour;
    std::vector<int>        neighbour_pos;

    int run(io::Argument & args, int tid, int num_threads);
    
  private:
    int switch_replica(int i, int j);
    int synchronise_thread_state(int tid);
    int synchronise_thread_replica(int tid);
    int synchronise_replica(int r);
    
    int update_thread_state(int tid);
    int update_replica(int r);
    
    gsl_rng * m_rng;

    std::vector<configuration::Configuration> m_conf;

  };
  
  /**
   * replica exchange accessor
   * for OMP
   */
  extern Replica_Exchange * replica_master;
  
} // util

#endif
