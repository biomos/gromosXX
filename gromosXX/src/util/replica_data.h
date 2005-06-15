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
    int        ID;
    double     temperature;
    double     lambda;
    double     energy;
    double     switch_temperature;
    double     switch_lambda;
    double     switch_energy;
    int        run;
    state_enum state;
    double     probability;
    bool       switched;
  };
}

#endif
