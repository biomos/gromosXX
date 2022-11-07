/*
 * File:   replica_exchange_base_eds.cc
 * Author: bschroed
 *
 * Created on April 18, 2018, 3:38 PM
 * Modified June 18, 2021 - bschroed, srieder
 */
#include <replicaExchange/replica_exchangers/1D_S_RE_EDS/replica_exchange_base_eds.h>

//Constructor
#include <replicaExchange/replica_exchangers/replica_exchange_base_interface.h>
#include <replicaExchange/replica/replica.h>
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
#define MODULE re
#define SUBMODULE replica_exchanger

/*
 * ADDITIONAL FUNCTIONS for REDS
 */

re::replica_exchange_base_eds::replica_exchange_base_eds(io::Argument _args,
                                                            unsigned int cont, 
                                                            unsigned int globalThreadID, 
                                                            replica_graph_control & replicaGraphMPIControl, 
                                                            simulation::MpiControl & replica_mpi_control):  
                            replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
                            reedsParam(replica->sim.param().reeds)
{
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t START" );
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t simID "<<simulationID);

    //RE-Vars
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t setParams" );
    setParams();

    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t DONE");
}

re::replica_exchange_base_eds::~replica_exchange_base_eds() {
    delete replica;
}

//SWAPPING Functions
void re::replica_exchange_base_eds::determine_switch_probabilities(){
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":determineSwitchPos:\t START");
    replica->sim.param().reeds.eds_para[partnerReplicaID].pos_info.second = simulationID;
    swap_replicas_1D(partnerReplicaID);
    DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":determineSwitchPos:\t DONE");
}

void re::replica_exchange_base_eds::swap_replicas_1D(const unsigned int partnerReplicaID) {
    DEBUG(4, "replica_exchange_base_eds:replica "<<  globalThreadID <<":swap1D - S:\t  START");
    DEBUG(4, "replica_exchange_base_eds:replica "<<  globalThreadID <<":swap1D - S:\t  sim: "<< simulationID << " \t\t "<< partnerReplicaID);
    unsigned int num_s = replica->sim.param().reeds.num_s;

    // does partner exist?
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaID) {
        DEBUG(4, "replica_exchange_base_eds:replica "<<  globalThreadID <<":swap1D:\t  swap prop calculation");

        // posts a MPI_Recv(...) matching the MPI_Send below
        DEBUG(7, "\nSWAP REPLICAs: " << globalThreadID << "\n");
        DEBUG(7, "replica "<< globalThreadID <<": Partner0 in swap Current ER  = "<< replica->conf.current().energies.eds_vr <<" \n");
        DEBUG(7, "replica "<< globalThreadID <<": Partner0 in swap OLD ER   = "<< replica->conf.old().energies.eds_vr <<" \n");
        DEBUG(7, "replica "<< globalThreadID <<": Partner0 in swap epot = "<< epot <<" \n");

        probability = calc_probability(partnerReplicaID);
        const double randNum = rng.get();

        std::vector<double> prob(2);
        prob[0] = probability;
        prob[1] = randNum;
      #ifdef XXMPI
        MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaID, SENDCOORDS, replicaGraphMPIControl().comm);
      #endif

        if (randNum < probability) {
          switched = true;
        } else
          switched = false;
      } else {    //The Partner sends his data to The calculating Thread
        DEBUG(7,  "replica "<< globalThreadID <<": Partner1 in swap epot   = "<< epot <<" \n");

        //special case if lambda also needs to be exchanged
        bool sameLambda = (l == replica->sim.param().reeds.eds_para[partnerReplicaID].s[0]);
        if(!sameLambda){      //exchange LAMBDA
          // E21: Energy with configuration 2 and lambda 1(of partner)
          const double E21 = calculate_energy(partnerReplicaID);
          // this we can store as the partner energy of the current replica
          epot_partner = E21;
          // E22: Energy with configuration 2 and lambda 2(own one)

      #ifdef XXMPI
          const double E22 = epot;
          // send E21 and E22
          double energies[2] = {E22, E21};
          //this send operation is matched in calc_probability()
          DEBUG(1,"\n\nreplica_exchange_base_eds:replica "<<  globalThreadID <<"replica_exchange_base_interface: SWAP_1D before Send\n");
          MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
          DEBUG(1,"\n\nreplica_exchange_base_eds:replica "<<  globalThreadID <<"replica_exchange_base_interface: SWAP_1D after Send\n");
      #endif
        } else { // sameLambda
      #ifdef XXMPI

          double energies[2] = {epot, 0.0};
          MPI_Send(&energies[0],2,MPI_DOUBLE, partnerReplicaID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
      #endif
       }
        if (replica->sim.param().pcouple.scale != math::pcouple_off) {
      #ifdef XXMPI
          math::Box box_replica = replica->conf.current().box;    //exchange box
          MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerReplicaID, BOX,  replicaGraphMPIControl().comm);
      #endif
        }

      #ifdef XXMPI
        MPI_Status status;
      #endif
        std::vector<double> prob;
        prob.resize(2);
      #ifdef XXMPI
        MPI_Recv(&prob[0], 2, MPI_DOUBLE, partnerReplicaID, SENDCOORDS,  replicaGraphMPIControl().comm, &status);
      #endif
        //Have we been exchanged little partner?
        probability = prob[0];
        double randNum = prob[1];

        if (randNum < probability) {
          switched = true;
        } else {
          switched = false;
        }
    }
    DEBUG(1, "replica_exchange_base_eds:replica_exchange_base_eds:replica "<< globalThreadID <<":swap1D:\t  DONE");
}

int re::replica_exchange_base_eds::find_partner() const{
    DEBUG(1,"replica_exchange_base_eds: FIND_PARTNER\tSTART\n\n");
    unsigned int numS = replica->sim.param().reeds.num_s;
    unsigned int ID = simulationID;
    unsigned int partner = ID;
    bool even = ID % 2 == 0;

    if (run % 2 == 1) {
        if (even){
            partner = ID + 1;
        }
        else{
            partner = ID - 1;
        }
    } else {
        if (even){
            partner = ID - 1;
        }
        else{
            partner = ID + 1;
        }
    }
    if (partner > numS - 1 || partner < 0) {
        partner = ID;
    }

    DEBUG(1,"replica_exchange_base_eds: FIND_PARTNER\tDONE \tpartner: "<< partner <<" \n\n");
    return partner;
}

void re::replica_exchange_base_eds::setParams(){
    // set some variables
    stepsPerRun = replica->sim.param().step.number_of_steps;
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t NUMBER OF STEPS "<<stepsPerRun);

    run = 0;
    total_runs = replica->sim.param().reeds.trials + replica->sim.param().reeds.equilibrate;
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t NUMBER OF total_runs "<<total_runs);

    partnerReplicaID = simulationID;
    time = replica->sim.time();
    steps = 0;
    switched = 0;
    replica->curentStepNumber=0;
    replica->totalStepNumber = total_runs*stepsPerRun;
    replica->stepsPerRun= stepsPerRun;

    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t PARAM START");

    T = replica->sim.param().reeds.temperature;
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t got  T " << T);
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t got simulationID: "<< simulationID);

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
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t got s" << l);

    dt = replica->sim.param().step.dt;
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t dt " <<dt);
    DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":setParams:\t PARAM DONE ");    
}

void re::replica_exchange_base_eds::set_s() {
  DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":set_s:\t START ");

  eds_para = replica->sim.param().reeds.eds_para[simulationID];
  replica->sim.param().eds = eds_para;
  DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":set_s:\t eds_para s size: " << replica->sim.param().eds.s.size());
  l = replica->sim.param().eds.s[0];    //todoAssume only 1s EDS
  DEBUG(4,"replica_exchange_base_eds "<< globalThreadID <<":set_s:\t DONE " );
}

void re::replica_exchange_base_eds::init() {
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t init \t START");
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t start init from baseclass \t NEXT");
  //replica->init();
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t init_eds_stat \t NEXT");
  init_eds_stat();
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init:\t DONE");
}

//initialize output files
void re::replica_exchange_base_eds::init_eds_stat(){
        DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init_eds_stat:\t START");

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
        replicaStatData[currentID].epot_vec.resize(replicaGraphMPIControl().numberOfReplicas);
        replicaStatData[currentID].prob_vec.resize(replicaGraphMPIControl().numberOfReplicas);

        DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":init_eds_stat:\t DONE");
}

//RE
////exchange params
void re::replica_exchange_base_eds::reset_eds() {//only reset switched parameters of change_eds() function
  replica->sim.param().eds = eds_para;
  replica->sim.param().step.dt = dt;
}

void re::replica_exchange_base_eds::change_eds(const unsigned int partner){//only change parameters, which are needed for energy calculation i.e.

  int idx = 0;
  if (replica->sim.param().reeds.num_s == 1){
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


double re::replica_exchange_base_eds::calculate_energy_core() {

    force_orig = replica->conf.current().force;
    virial_tensor_orig = replica->conf.current().virial_tensor;

    double energy = 0.0;
    algorithm::Algorithm * ff = nullptr;

    ff = replica->md.algorithm("EDS");

    //Calculate energies
    DEBUG(5, "replica_reeds_base_eds "<< globalThreadID <<":calculate_energy:\t calc energies");
    if (ff->apply(replica->topo, replica->conf, replica->sim)) {
        print_info("Error in Forcefield energy calculation!");
     #ifdef XXMPI
        MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
    #endif
        return 1;
    }

    //return energies
    DEBUG(5, "replica_reeds_base_edsreplica_reeds "<< globalThreadID <<":calculate_energy"
                                                                       ":\t return energies");
    energy=replica->conf.current().energies.eds_vr;

    replica->conf.current().force= force_orig;
    replica->conf.current().virial_tensor= virial_tensor_orig;
    return energy;
}


double re::replica_exchange_base_eds::calculate_energy(const unsigned int partnerThreadID) {
    DEBUG(4, "replica_exchange_base_eds "<< globalThreadID <<":calculate_energy:\t  START");

    change_eds(partnerThreadID);

    double energy = calculate_energy_core();

    reset_eds();
    DEBUG(4, "replica_exchange_base_eds "<< globalThreadID <<":calculate_energy:\t  DONE");

    return energy;
}

void re::replica_exchange_base_eds::calculate_energy_helper(const unsigned int partnerThreadID){
    DEBUG(4, "replica_exchange_base_eds "<< globalThreadID <<":calculate_energy_helper:\t  START");

    DEBUG(4, "replica_exchange_base_eds "<< globalThreadID <<":calculate_energy_helper:\t  DONE");
}