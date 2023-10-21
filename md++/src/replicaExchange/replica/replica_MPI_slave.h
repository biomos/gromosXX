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

//
// Created by bschroed on 10/28/19.
//

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

//used in mpi splitting on nodes
#include <math/periodicity.h>
#include <algorithm/constraints/shake.h>
#include <algorithm/constraints/m_shake.h>
#include <algorithm/integration/monte_carlo.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include "_replica_Interface.h"

#ifndef GROMOSXX_REPLICA_MPI_Slave_H
#define GROMOSXX_REPLICA_MPI_Slave_H

namespace re {
    class replica_MPI_Slave: public virtual replica_Interface {
    protected:
        /**
         * SLAVE THINGS
         */
        //slave directives
        bool do_nonbonded;  
        bool do_cmc;
        bool do_shake;
        bool do_m_shake;
        bool do_shake_twice;

        //attributes
        interaction::Forcefield * ff;
        interaction::Interaction * nb;
        algorithm::Shake * shake;
        algorithm::M_Shake * m_shake;
        algorithm::Monte_Carlo * monte_carlo;

        void receive_coords();
        
    public:
        replica_MPI_Slave(io::Argument _args, int cont, int globalThreadID, simulation::MpiControl & replica_mpi_control);
        ~replica_MPI_Slave();
        
        /**
         * simulation
         */
        //simulation::Simulation sim;

        /**
         * run MD
         * @return 
         */
        void run_MD() override;

        /**
         * help updating rep Energies
         * @return
         */
        void calculateEnergiesHelper();

        
        /**
         *  init -  this is an intentaionally empty func.
         */
        void init() override {};
    };
}//namespace re

#endif //GROMOSXX_REPLICA_MPI_H
