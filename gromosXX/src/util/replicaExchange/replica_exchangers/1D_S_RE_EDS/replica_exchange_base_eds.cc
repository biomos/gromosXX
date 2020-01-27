/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_exchange_base_eds.cc
 * Author: bschroed
 * 
 * Created on April 18, 2018, 3:38 PM
 */
#include "util/replicaExchange/replica_mpi_tools.h"
#include <util/replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_base_eds.h>

//Constructor
#include <util/replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_base.h>
#include <util/replicaExchange/replica/replica.h>
#include "replicaExchange/replica/replica_MPI_master.h"
#include "replicaExchange/replica/replica_MPI_slave.h"
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
#include <math/random.h>
#include <math/volume.h>
#include <string>

#ifdef XXMPI
#include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

/*
 * ADDITIONAL FUNCTIONS for REDS
 */

util::replica_exchange_base_eds::replica_exchange_base_eds(io::Argument _args, 
                                                            unsigned int cont, 
                                                            unsigned int globalThreadID, 
                                                            std::vector<std::vector<unsigned int> >  replica_owned_threads, 
                                                            std::map<ID_t, rank_t> &thread_id_replica_map, 
                                                            simulation::mpi_control_struct replica_mpi_control):  
                            replica_exchange_base(_args, cont, globalThreadID, replica_owned_threads, thread_id_replica_map, replica_mpi_control),
                            reedsParam(replica->sim.param().reeds)
{
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t START");
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t replica Type\t "<< typeid(replica).name());
    createReplicas(cont, globalThreadID, replica_mpi_control);
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t DONE");
}

void util::replica_exchange_base_eds::createReplicas(int cont, int globalThreadID, simulation::mpi_control_struct replica_mpi_control){
  MPI_DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":createReplicas:\t START");
  
  if(!replica){
    delete replica;
  }
  MPI_DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":createReplicas:\t Deleted");

  // create the number of replicas that are assigned to my node
    if(simulationThreads>1){
        if(simulationThreadID == 0){
            replica = new util::replica_MPI_Master(args, cont, globalThreadID, replica_mpi_control);
        }
        else{
            replica = new util::replica_MPI_Slave(args, cont, globalThreadID, replica_mpi_control);
        }
    }
    else{
        replica = new util::replica(args, cont, globalThreadID, replica_mpi_control);
    }
  MPI_DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":createReplicas:\t DONE");
   
  //give the correct replica vals to simulation
  eds_para = replica->sim.param().reeds.eds_para[simulationID];
  replica->sim.param().eds = eds_para;
    //TODO: Still ugly ! l==s bschroed
  l = replica->sim.param().eds.s[0];    //todoAssume only 1s EDS
}


util::replica_exchange_base_eds::~replica_exchange_base_eds() {
    delete replica;
}

void util::replica_exchange_base_eds::init() {
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t init \t START");
  MPI_DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t start init from baseclass \t NEXT");
  replica->init();
  MPI_DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t init_eds_stat \t NEXT");
  this->init_eds_stat();
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t DONE");
}

//initialize output files  
void util::replica_exchange_base_eds::init_eds_stat(){
        DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init_eds_stat:\t START");
        
        ID_t currentID=1000; //error value
        currentID = simulationID;
        replicaStatData[currentID].ID =currentID;
        replicaStatData[currentID].T=T;
        replicaStatData[currentID].s=l; //l==s because of the implementation of hamiltonian replica exchange.
        replicaStatData[currentID].dt=dt;
        replicaStatData[currentID].run=0;
        replicaStatData[currentID].epot_vec.resize(this->numReplicas);
        replicaStatData[currentID].prob_vec.resize(this->numReplicas);
        
        DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init_eds_stat:\t DONE");
}

//RE
////exchange params
void util::replica_exchange_base_eds::reset_eds() {//only reset switched parameters of change_eds() function 
  replica->sim.param().eds = eds_para;
  replica->sim.param().step.dt = dt;
  replica->conf.current().force= force_orig;
  replica->conf.current().virial_tensor= virial_tensor_orig;
}

void util::replica_exchange_base_eds::change_eds(const unsigned int partner){//only change parameters, which are needed for energy calculation i.e. 
  
  int idx;
  if (replica->sim.param().reeds.num_l == 1){
    idx = 0;
  }
  else{
    idx = partner;
  }
  
  replica->sim.param().step.dt = replica->sim.param().reeds.dt[idx];
  replica->sim.param().eds= replica->sim.param().reeds.eds_para[idx];
  force_orig = replica->conf.current().force;
  virial_tensor_orig = replica->conf.current().virial_tensor;
}

////calc exchange Energy
 /*
 * calc_energy_eds_stat() is only used for statistical purposes in eds_stat()
 * In order to avoid any adjustment of the mpi communication and thus reducing the complexity, the 
 * energy_calculation and probability calculations from replica.cc are not adjusted to work 
 * for non-pairwise exchanges. Instead, calc_energy_eds_stat() calculates the potential energy
 * of the current configuration for a new smoothing parameter s.
 * The exchange probabilities can be calculated in a postprocessing step, using these energies
 * given in the energy_stat output files.
 */
double util::replica_exchange_base_eds::calc_energy_eds_stat(double s){
    double old_dt;
    double old_s;
    double old_eds_vr;
    algorithm::Algorithm * ff;   
    if(replica->sim.param().eds.eds){
          //to reset old state
          old_dt=replica->sim.param().step.dt;
          old_s=replica->sim.param().eds.s[0];
          old_eds_vr=replica->conf.current().energies.eds_vr;
          force_orig = replica->conf.current().force;
          virial_tensor_orig = replica->conf.current().virial_tensor;
          //only temporary change
          replica->sim.param().eds.s[0]=s;
          
          ff = replica->md.algorithm("EDS");
    }
    else {
          print_info("eds_stat() i.e calc_energy_eds_stat() called for non EDS simulation!");
      #ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
      #endif
    }
    
    //Calculate energies
    if (ff->apply(replica->topo, replica->conf, replica->sim)) {
      print_info("Error in Forcefield energy calculation!");
     #ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
     #endif
      return 1;
    }
    
    double energy=replica->conf.current().energies.eds_vr; 
    
    // reset old EDS state
    replica->conf.current().energies.eds_vr=old_eds_vr;
    replica->sim.param().eds.s[0] = old_s;
    replica->sim.param().step.dt = old_dt;
    replica->conf.current().force=force_orig;
    replica->conf.current().virial_tensor=virial_tensor_orig;
    
    return energy;
}

double util::replica_exchange_base_eds::calculate_energy_core() {
    
    double energy = 0.0;
    algorithm::Algorithm * ff;  
    
     ff = replica->md.algorithm("EDS");

    //Calculate energies    
    DEBUG(5, "replica_reeds "<< globalThreadID <<":calculate_energy:\t calc energies"); 
    if (ff->apply(replica->topo, replica->conf, replica->sim)) {
      print_info("Error in Forcefield energy calculation!");
     #ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
    #endif
      return 1;
    }
    
    //return energies
    DEBUG(5, "replica_reeds "<< globalThreadID <<":calculate_energy"
            ":\t return energies"); 
    energy=replica->conf.current().energies.eds_vr; 
    return energy;
}


double util::replica_exchange_base_eds::calculate_energy(const unsigned int selectedReplicaID) {
    DEBUG(4, "replica_reeds "<< globalThreadID <<":calculate_energy:\t START"); 
 

    DEBUG(5, "replica_reeds "<< globalThreadID <<":calculate_energy:\t get Partner settings"); 
    if(selectedReplicaID!=simulationID){ 
        change_eds(selectedReplicaID);
    }
    
    double energy =  calculate_energy_core();
    
    if(selectedReplicaID!=simulationID){
        reset_eds();
    }
    DEBUG(4, "replica_reeds "<< globalThreadID <<":calculate_energy:\t DONE"); 
    return energy;
}


