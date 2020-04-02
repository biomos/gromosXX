/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   replica_exchange_base_2d_s_eoff_eds.cc
 * Author: theosm
 *
 * Created on March 29, 2020, 11:08 AM
 */
#include "util/replicaExchange/replica_mpi_tools.h"
#include <util/replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_base_2d_s_eoff_eds.h>

//Constructor
#include <util/replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
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

util::replica_exchange_base_2d_s_eoff_eds::replica_exchange_base_2d_s_eoff_eds(io::Argument _args,
                                                            unsigned int cont,
                                                            unsigned int globalThreadID,
                                                            replica_graph_mpi_control replicaGraphMPIControl,
                                                            simulation::mpi_control_struct replica_mpi_control):
                            replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
                            reedsParam(replica->sim.param().reeds)
{
    MPI_DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t START" );
    DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t simID "<<simulationID);

    //RE-Vars
    MPI_DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t setParams" );
    setParams();

    MPI_DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t DONE");
}

util::replica_exchange_base_2d_s_eoff_eds::~replica_exchange_base_2d_s_eoff_eds() {
    delete replica;
}

void util::replica_exchange_base_2d_s_eoff_eds::setParams(){
    // set some variables
    stepsPerRun = replica->sim.param().step.number_of_steps;
    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t NUMBER OF STEPS "<<stepsPerRun);

    run = 0;
    total_runs = replica->sim.param().replica.trials + replica->sim.param().replica.equilibrate;
    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t NUMBER OF total_runs "<<total_runs);

    partnerReplicaID = simulationID;
    time = replica->sim.time();
    steps = 0;
    switched = 0;
    replica->curentStepNumber=0;
    replica->totalStepNumber = total_runs*stepsPerRun;
    replica->stepsPerRun= stepsPerRun;

    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t PARAM START");

    T = replica->sim.param().reeds.temperature;
    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t got  T " << T);
    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t got simulationID: "<< simulationID);

    //set position info
    replica->sim.param().reeds.eds_para[simulationID].pos_info = std::make_pair(simulationID, simulationID);
    DEBUG(1, "BASE Constructor with simulationID, replica->pos_info= " << simulationID << ", "
    << replica->sim.param().reeds.eds_para[simulationID].pos_info.first << ", "
    << replica->sim.param().reeds.eds_para[simulationID].pos_info.second << "\n");

    pos_info = replica->sim.param().reeds.eds_para[simulationID].pos_info;
    DEBUG(1, "BASE Constructor with simulationID, pos_info= " << simulationID << ", "
    << pos_info.first << ", " << pos_info.second << "\n");

    //just to check -- theosm
    std::pair<int, int> a = reedsParam.eds_para[simulationID].pos_info;
    DEBUG(1, "JUST TO CHECK: BASE Constructor with simulationID, reedsParam->pos_info= " << simulationID << ", "
    << a.first << ", " << a.second << "\n");

    set_s();
    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t got s" << l);

    dt = replica->sim.param().step.dt;
    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t dt " <<dt);

    MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t PARAM DONE ");
}

void util::replica_exchange_base_2d_s_eoff_eds::set_s() {
  MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":set_s:\t START ");

  eds_para = replica->sim.param().reeds.eds_para[simulationID];
  replica->sim.param().eds = eds_para;
  MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":set_s:\t eds_para s size: " << replica->sim.param().eds.s.size());

  l = replica->sim.param().eds.s[0];    //todoAssume only 1s EDS
  MPI_DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":set_s:\t DONE " );
}

void util::replica_exchange_base_2d_s_eoff_eds::init() {
  MPI_DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t init \t START");
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t start init from baseclass \t NEXT");
  //replica->init();
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t init_eds_stat \t NEXT");
  init_eds_stat();
  MPI_DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t DONE");
}

//initialize output files
void util::replica_exchange_base_2d_s_eoff_eds::init_eds_stat(){
        DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init_eds_stat:\t START");

        ID_t currentID=1000; //error value
        currentID = simulationID;
        replicaStatData[currentID].ID=currentID;
        //assignment of position info of each replica
        replicaStatData[currentID].pos_info.first = pos_info.first;
        replicaStatData[currentID].pos_info.second = pos_info.second;
        DEBUG(3, "init_eds_stat(), replicaStatData[currentID].pos_info.first= " << replicaStatData[currentID].pos_info.first
        << " with currentID= " << replicaStatData[currentID].ID << "\n");
        DEBUG(3, "init_eds_stat(), replicaStatData[currentID].pos_info.second= " << replicaStatData[currentID].pos_info.second
        << " with currentID= " << replicaStatData[currentID].ID << "\n");
        replicaStatData[currentID].T=T;
        replicaStatData[currentID].s=l; //l==s because of the implementation of hamiltonian replica exchange.
        replicaStatData[currentID].dt=dt;
        replicaStatData[currentID].run=0;
        replicaStatData[currentID].epot_vec.resize(replicaGraphMPIControl.numberOfReplicas);
        replicaStatData[currentID].prob_vec.resize(replicaGraphMPIControl.numberOfReplicas);

        DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init_eds_stat:\t DONE");
}

//RE
////exchange params
void util::replica_exchange_base_2d_s_eoff_eds::reset_eds() {//only reset switched parameters of change_eds() function
  replica->sim.param().eds = eds_para;
  replica->sim.param().step.dt = dt;
  replica->conf.current().force= force_orig;
  replica->conf.current().virial_tensor= virial_tensor_orig;
}

void util::replica_exchange_base_2d_s_eoff_eds::change_eds(const unsigned int partner){//only change parameters, which are needed for energy calculation i.e.

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
double util::replica_exchange_base_2d_s_eoff_eds::calc_energy_eds_stat(double s){
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

double util::replica_exchange_base_2d_s_eoff_eds::calculate_energy_core() {

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


double util::replica_exchange_base_2d_s_eoff_eds::calculate_energy(const unsigned int selectedReplicaID) {
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
