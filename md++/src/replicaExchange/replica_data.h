/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file replica_data.h
 * contains replica_ata struct
 * Modified June 18, 2021 - bschroed, srieder
 */

#ifndef REPLICA_DATA_H
#define	REPLICA_DATA_H
namespace re
{

  /**
   * @struct replica_data
   * contains all necessary data of a replica that is written to the output file.
   * See re::replica.h for more information.
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
   * @struct replica_data
   * contains all necessary data of a replica that is written to the output file.
   * See re::replica.h for more information.
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
