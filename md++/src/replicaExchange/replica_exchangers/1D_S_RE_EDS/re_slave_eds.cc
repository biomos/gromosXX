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
 * File:   replica_exchange_slave_eds.cc
 * Author: bschroed
 *
 * Created on August 31, 2018, 10:43 AM
 * Modified June 18, 2021 - bschroed, srieder
 */
#include <replicaExchange/replica_exchangers/1D_S_RE_EDS/re_slave_eds.h>


#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica_exchanger

re::replica_exchange_slave_eds::replica_exchange_slave_eds(io::Argument & _args,
                                                            unsigned int cont,
                                                            unsigned int globalThreadID,
                                                            replica_graph_control & replicaGraphMPIControl,
                                                            simulation::MpiControl & replica_mpi_control):
            replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
            replica_exchange_base_eds(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
            replica_exchange_slave_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control)
{

    DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":Constructor:\t START");

    DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":Constructor:\t DONE"); 
}

void re::replica_exchange_slave_eds::send_to_master() const{
  #ifdef XXMPI
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t START");
  if(replicaInfoSender){
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t Show vPots");

    re::repInfo info;
    std::vector<double> eds_energies;
    info.run = run;
    info.epot = epot;
    info.epot_partner = epot_partner;
    info.partner = partnerReplicaID;
    info.probability = probability;
    info.switched = switched;

    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << "send_to_master:\t epotTot\t "<< info.epot);

    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_REPINFO");
    MPI_Send(&info, 1, MPI_REPINFO, replicaGraphMPIControl().masterID, REPINFO, replicaGraphMPIControl().comm);

    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS");
    
    // If conformations were switched, the values used to determine
    // if we want to exchange are now in conf.old()
    if (switched){
      eds_energies= replica->conf.old().energies.eds_vi;
    } else {
      eds_energies= replica->conf.current().energies.eds_vi;
    }

    MPI_Send(&eds_energies[0], 1, MPI_EDSINFO, replicaGraphMPIControl().masterID, EDSINFO, replicaGraphMPIControl().comm);
    DEBUG(4,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t\t send MPI_EDS \t DONE" );
  }
  DEBUG(2,"replica_exchange_slave_eds " << globalThreadID << ":send_to_master:\t DONE");

  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!";
  #endif
}

