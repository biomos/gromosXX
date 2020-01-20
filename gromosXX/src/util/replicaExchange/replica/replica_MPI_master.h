/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_MPI_master.h
 * Author: bschroed
 *
 * Created on November 8, 2019, 12:48 PM
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
#include <util/replicaExchange/repex_mpi.h>
#include <util/usage.h>
#include "replica_Interface.h"


#ifndef REPLICA_MPI_MASTER_H
#define REPLICA_MPI_MASTER_H

namespace util{
    class replica_MPI_Master: public virtual replica_Interface {
    protected:
        int simulation_ID;
        int simulation_rank;
        int simulation_num_threads;
                
    private:
    /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
    replica_MPI_Master(const replica_MPI_Master & r);
        
        
    public:
        replica_MPI_Master(io::Argument _args, int cont, int _ID, int _rank, int simulation_rank, int simulation_ID, int simulation_num_threads);
        virtual ~replica_MPI_Master();
        
#ifdef XXMPI
        /**
         * run MD
         * @return 
         */
        void run_MD() override;

#endif

    private:

    };
}
#endif /* REPLICA_MPI_MASTER_H */

