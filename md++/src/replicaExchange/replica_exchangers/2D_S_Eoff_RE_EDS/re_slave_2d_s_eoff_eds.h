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

/*
 * File:   re_slave_2d_s_eoff_eds.h
 * Author: theosm
 *
 * Created on March 29, 2020, 11:03 AM
 */

#include <replicaExchange/replica_exchangers/replica_exchange_slave_interface.h>
#include <replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/re_base_2d_s_eoff_eds.h>
//for the constructor
#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/usage.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>
#include <unistd.h>

#ifndef REPLICA_EXCHANGE_SLAVE_2D_S_EOFF_EDS_H
#define REPLICA_EXCHANGE_SLAVE_2D_S_EOFF_EDS_H

namespace re{
    class replica_exchange_slave_2d_s_eoff_eds: public replica_exchange_base_2d_s_eoff_eds, public replica_exchange_slave_interface {
    public:
        replica_exchange_slave_2d_s_eoff_eds(io::Argument & _args,
                unsigned int cont,
                unsigned int globalThreadID,
                replica_graph_control &replicaGraphMPIControl,
                simulation::MpiControl &replica_mpi_control);
        /**
        * sends information of all replicas to master
        */
        void send_to_master() const override;

    private:
        //replica_exchange_slave_2d_s_eoff_eds(const replica_exchange_slave_2d_s_eoff_eds& orig);
        virtual ~replica_exchange_slave_2d_s_eoff_eds(){};
    };
}
#endif /* REPLICA_EXCHANGE_SLAVE_2D_S_EOFF_EDS_H */
