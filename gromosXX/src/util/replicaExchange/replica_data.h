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
        std::pair<int, int> pos_info;
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
   * @struct replica_data
   * contains all necessary data of a replica that is written to the output file.
   * See util::replica.h for more information.
   */
    struct reeds_replica_data
    {
        unsigned int ID;
        std::pair<int, int> pos_info;
        double T;
        double l;
        double dt;
        int       switched;
        unsigned int        run;
        unsigned int         partner;
        double     epot;
        double     epot_partner;
        double     probability;

        simulation::Parameter::eds_struct eds_state;
        std::vector<double> Vi; //SID for RE-EDS I want to give out all potential energies for each individual state in the repdat file. //todo_remove
    };
   /**
   * @struct replica_stat_data
   * contains additional data of a replica that is optionally written to the output file.
   */
    struct reeds_replica_stat_data
    {
        unsigned int ID;
        std::pair<int, int> pos_info;
        double T;
        double s;
        double dt;
        unsigned int        run;
        //vector because of eds_stat() fct. if not used, one just stores value at the correspnding replica id position, all others are empty.
        std::vector<double>     epot_vec;
        std::vector<double>     prob_vec;
        simulation::Parameter::eds_struct eds_state;

    };

}
#endif	/* REPLICA_DATA_H */
