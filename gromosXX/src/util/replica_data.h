/**
 * @file replica_data.h
 * replica exchange data structure
 */

#ifndef INCLUDED_REPLICA_DATA_H
#define INCLUDED_REPLICA_DATA_H

namespace util
{
  /**
   * @enum state_enum
   * state of a replica
   */
  enum state_enum{ waiting=0, ready=1, running=2, terminate=4, st_error=5 };
  
  /**
   * @struct Replica_Data
   * replica information
   */
  struct Replica_Data
  {
    /**
     * replica ID
     */
    int        ID;
    /**
     * index to temperature array for current run
     */
    int        Ti;
    /**
     * index to lambda array for current run
     */
    int        li;
    /**
     * potential energy using current Hamiltonian
     */
    double     epot_i;
    /**
     * index to temperature array to match switch partner
     * (no influence on switched Hamiltonian
     *  or do you want temperature dependent Hamiltonians???)
     */
    int        Tj;
    /**
     * index to lambda array to match switch partner
     * used to determine the switched Hamiltonian
     */
    int        lj;
    /**
     * potential energy using switched Hamiltonian
     */
    double     epot_j;
    /**
     * run number of replica
     * (only switches between replicas with identical run number)
     */
    int        run;
    /**
     * replica state:
     * - waiting (wait)
     * - ready   (rdy)
     * - running (run)
     * - st_error (err)
     * - terminate (term)
     */
    state_enum state;
    /**
     * probability of last switch
     */
    double     probability;
    /**
     * last switch successfull?
     */
    bool       switched;
  };
}

#endif
