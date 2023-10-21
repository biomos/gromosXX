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
 * File:   replica_exchange_master_2d_l_T_HREMD.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:18 PM
 */

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

#include <io/configuration/out_configuration.h>
#include <replicaExchange/replica/replica.h>
#include <replicaExchange/replica_exchangers/replica_exchange_master_interface.h>
#include <replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_master_2d_l_T_HREMD.h>

#include <string>

#ifdef XXMPI
    #include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica_exchanger

re::replica_exchange_master_2d_l_T_HREMD::replica_exchange_master_2d_l_T_HREMD(io::Argument & args,
                                                                               unsigned int cont,
                                                                               unsigned int globalThreadID,
                                                                               replica_graph_control & replicaGraphMPIControl,
                                                                               simulation::MpiControl & replica_mpi_control) :
        replica_exchange_base_interface(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_base_2d_l_T_HREMD(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_master_interface(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control)
{
#ifdef XXMPI
  DEBUG(2,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t START");

  assert(replicaGraphMPIControl.masterID == replicaGraphMPIControl.threadID);    //TODO: This can be removed in future! bscrhoed
  assert(replicaGraphMPIControl.numberOfReplicas > 0);
  DEBUG(2,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t rep_params THERE?");
  DEBUG(2,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t" << replica->sim.param().replica.num_l);
  DEBUG(2,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t" << replica->sim.param().replica.lambda[0]);

  assert(repParams.num_l > 0);
  
  DEBUG(4,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t Init Replicas \t Next");
  replicaData.resize(replicaGraphMPIControl.numberOfReplicas);
  DEBUG(4,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t Replica_data type \t " << typeid(replicaData).name());

  //initialize data of replicas
  int ID = 0;
  for (int i = 0; i < repParams.num_l; ++i) {
    for (int j = 0; j < repParams.num_T; ++j) {
      replicaData[ID].ID = ID;
      replicaData[ID].T = repParams.temperature[j];
      DEBUG(4,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t Init Replicas \t "<< repParams.temperature[j]);
      replicaData[ID].l = repParams.lambda[i];
      replicaData[ID].dt = repParams.dt[i];
      ++ID;
    }
  }

  // set output file
 DEBUG(2,"replica_exchange_master_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t DONE");
#else
   throw "Cannot initialize replica_exchange_master_2d_l_T_HREMD without MPI!";
#endif
}


re::replica_exchange_master_2d_l_T_HREMD::~replica_exchange_master_2d_l_T_HREMD() {
   repOut.close();
}

