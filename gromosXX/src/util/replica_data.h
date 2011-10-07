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
}
#endif	/* REPLICA_DATA_H */

