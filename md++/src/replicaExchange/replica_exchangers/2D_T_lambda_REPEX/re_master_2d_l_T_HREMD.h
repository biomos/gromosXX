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
 * @file replica_exchange_master_2d_l_T_HREMD.h
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
#include <replicaExchange/replica_exchangers/replica_exchange_master_interface.h>
#include <replicaExchange/replica_exchangers/2D_T_lambda_REPEX/re_base_2d_l_T_HREMD.h>

#include <replicaExchange/replica/replica.h>
#include <replicaExchange/replica_data.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#ifndef REPLICA_EXCHANGE_MASTER_H
#define	REPLICA_EXCHANGE_MASTER_H

namespace re {

  /**
   * @class replica_exchange_master_2d_l_T_HREMD
   * Additionally to replica_exchange_base_2d_l_T_HREMD: receives and writes data to file.
   */
  class replica_exchange_master_2d_l_T_HREMD : public virtual replica_exchange_base_2d_l_T_HREMD, public virtual replica_exchange_master_interface {
  public:
    /**
     * constructor
     * @param args io::Argument, passed on to replica
     * @param rank integer, rank of node
     * @param _size size of mpi comm world
     * @param _numReplicas total number of replicas
     * @param repIDs std::vector<int>, IDs of replicas the instance has to manage
     * @param repMap std::map<int,int>, maps replica IDs to nodes; needed for communication
     */
    replica_exchange_master_2d_l_T_HREMD(io::Argument & args,
                                         unsigned int cont,
                                         unsigned int globalThreadID,
                                         replica_graph_control & replicaGraphMPIControl,
                                         simulation::MpiControl & replica_mpi_control);
    /**
     * Destructor
     */
    virtual ~replica_exchange_master_2d_l_T_HREMD();

  };
}
#endif	/* REPLICA_EXCHANGE_MASTER_H */

