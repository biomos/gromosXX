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
 * File:   replica_MPI_master.h
 * Author: bschroed
 *
 * Created on November 8, 2019, 12:48 PM
 * Modified June 18, 2021 - bschroed, srieder
 */
#include <stdheader.h>
#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>
#include <io/configuration/out_configuration.h>
#include <io/argument.h>
#include <io/read_input.h>
#include <math/random.h>
#include <replicaExchange/repex_mpi.h>
#include <util/usage.h>
#include "_replica_Interface.h"


#ifndef REPLICA_MPI_MASTER_H
#define REPLICA_MPI_MASTER_H

namespace re{
    class replica_MPI_Master: public virtual replica_Interface {

    private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica_MPI_Master(const replica_MPI_Master & r);
        
        
    public:
        replica_MPI_Master(io::Argument _args, int cont, int globalThreadID, simulation::MpiControl & replica_mpi_control);
        virtual ~replica_MPI_Master();
        
        /**
         * run MD
         * @return 
         */
        void run_MD() override;
        double calculateEnergies() override;

    private:
        void send_coordinates();
        
        /**
         * Attributes
         */


    };
}
#endif /* REPLICA_MPI_MASTER_H */

