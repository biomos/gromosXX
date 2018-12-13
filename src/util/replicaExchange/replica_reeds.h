/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_reeds.h
 * Author: bschroed
 *
 * Created on April 23, 2018, 3:25 PM
 */

#include <util/replicaExchange/replica.h>

#include <io/argument.h>
#include <util/error.h>
#include <math/volume.h>
#include <simulation/simulation.h>
#include <simulation/parameter.h>

#ifdef XXMPI
#include <mpi.h>
#endif

#ifndef REPLICA_REEDS_H
#define REPLICA_REEDS_H

namespace util{
    class replica_reeds : virtual public replica {
    public:
        /**
         * Constructor
         * Reads in all input parameters and assigns rank and ID of the replica
         * @param _args io::Argument, copy of the one initialized in repex_mpi.cc 
         * @param _ID integer, unique identifier of the replica
         * @param _rank integer, for knowing where to send the data
         */
        replica_reeds(io::Argument _args, int cont, int _ID, int _rank);
        /*
        * energy calculation for statistical purposes of eds_stat() in replica_exchange_base.cc
        * for given configuration with new smoothing parameter s.
        */
        double calc_energy_eds_stat(double s);
        double calculate_energy(const int partner);
        
    private:
     /**
     * copy constructor
     * don't allow copies because there's a bug in the copy constructor of configuration::configurtiaon
     */
     replica_reeds(const replica_reeds& orig);
     /**
     * Sets eds_struct() parameters to original value of replica
     */
    void reset_eds();
    /**
     * Sets  eds_struct() parameters to value of partner (for energy calculations)
     */
    void change_eds(const unsigned int partner);
    
    public:
    /**
    * contains all parameters for RE-EDS
    */
    simulation::Parameter::eds_struct eds_para;
    /**
     * contains original forces which have to be reset after RE-EDS exchange energy calculation
     */
    math::VArray force_orig;
    /**
     * contains original virial which have to be reset after RE-EDS exchange energy calculation
     */
    math::Matrix virial_tensor_orig;

    //configuration::Configuration conf;

    };
}
#endif /* REPLICA_REEDS_H */

