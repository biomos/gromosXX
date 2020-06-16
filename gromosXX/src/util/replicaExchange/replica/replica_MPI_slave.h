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
#include <util/replicaExchange/repex_mpi.h>
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

namespace util {
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
         *  init -  this is an intentaionally empty func.
         */
        void init() override {};
    };
}//namespace util

#endif //GROMOSXX_REPLICA_MPI_H
