/**
 * @file src/simulation/simulation/simulation.h
 * the Simulation class
 */

#ifndef INCLUDED_SIMULATION_H
#define INCLUDED_SIMULATION_H

/**
 * @namespace simulation
 * namespace that contains the simulation class
 * and the associated topology, system and parameter clases
 */
namespace simulation
{
  /**
   * @class Simulation
   * holds together the system, the topology and
   * input parameters.
   * provides the simulation properties.
   */
  template<typename t_topo, typename t_system>
  class Simulation
  {
  public:
    typedef t_topo topology_type;
    typedef t_system system_type;
    
    /**
     * Constructor
     */
    explicit Simulation(t_topo &topo, t_system &sys);
    
    /**
     * topology accessor
     */
    t_topo & topology();
    /**
     * system accessor
     */
    t_system & system();
    /**
     * time accessor
     */
    double time();
    /**
     * set (initial) time.
     */
    void time(double t);
    /**
     * steps accessor
     */
    int steps();
    /**
     * pairlist update every n steps.
     */
    void nonbonded_update(int const update_step);
    /**
     * accessor pairlist update.
     */
    int nonbonded_update()const;
    /**
     * set short range cutoff.
     */
    void nonbonded_cutoff_short(double const cutoff_short);
    /**
     * get short range cutoff.
     */
    double nonbonded_cutoff_short()const;
    /**
     * set long range cutoff.
     */
    void nonbonded_cutoff_long(double const cutoff_long);
    /**
     * get long range cutoff.
     */
    double nonbonded_cutoff_long()const;

    /**
     * increase the time by dt.
     */
    void increase_time(double dt);

  private:
    /**
     * the topology.
     */
    topology_type &m_topology;
    /** 
     * the system.
     */
    system_type   &m_system;
    /**
     * the time.
     */
    double m_time;
    /**
     * the number of steps done.
     */
    int m_steps;

    /**
     * nonbonded update.
     */
    int m_nonbonded_update;
    
    /**
     * nonbonded short range cutoff.
     */
    double m_nonbonded_cutoff_short;
    
    /**
     * nonbonded long range cutoff.
     */
    double m_nonbonded_cutoff_long;
    
  }; // class Simulation
  
  
} // namespace simulation

// template definition
#include "simulation.tcc"

#endif
