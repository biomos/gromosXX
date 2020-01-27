/**
 * @file replica.h
 * Holds all the information for a single replica
 */
#include <stdheader.h>
#include "_replica_Interface.h"

#ifndef REPLICA_H
#define REPLICA_H


namespace util {

  /**
   * replica Class
   * All the data members are held public, so Master/Slave classes have full 
   * control over it.
   */
  class replica: public replica_Interface {
  private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica(const replica & r);
  public:
    /**
     * Constructor
     * Reads in all input parameters and assigns rank and ID of the replica
     * @param _args io::Argument, copy of the one initialized in repex_mpi.cc 
     * @param _ID integer, unique identifier of the replica
     * @param _rank integer, for knowing where to send the data
     */
    replica(io::Argument _args, int cont, int globalThreadID, simulation::mpi_control_struct replica_mpi_control);

    /**
     * runs an MD simulation
     */
    void run_MD();
    /**
     * init an MD simulation
     */
    void init();

  };
}

#endif  /* REPLICA_H */
