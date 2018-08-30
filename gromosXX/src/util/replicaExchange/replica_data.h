/**
 * @file replica_data.h
 * contains replica_ata struct
 */

#ifndef REPLICA_DATA_H
#define	REPLICA_DATA_H
namespace util
{
  /**
   * @struct replica_data
   * contains all necessary data of a replica that is written to the output file.
   * See util::replica.h for more information.
   */
    struct replica_data
    {
        unsigned int ID;
        double T;
        double l;
        double dt;
        int       switched;
        unsigned int        run;
        unsigned int         partner;
        double     epot;
        double     epot_partner;
        double     probability;
    };
   /**
   * @struct replica_stat_data
   * contains additional data of a replica that is optionally written to the output file during an Replica Exchange EDS run..
   */
    struct Reeds_replica_stat_data
    {
        unsigned int ID;    //Replica ID
        double T;   //Temperature
        double s;   //S-value
        double dt;  // timesteps
        unsigned int        run;    // current run (ammount of trials executed)
        //vector because of eds_stat() fct. if not used, one just stores value at the correspnding replica id position, all others are empty.
        std::vector<double>     epot_vec;
        std::vector<double>     prob_vec;
    };
}
#endif	/* REPLICA_DATA_H */

