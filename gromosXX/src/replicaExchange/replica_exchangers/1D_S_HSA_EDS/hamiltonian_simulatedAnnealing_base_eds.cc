/*
 * File:   hamiltonian_simulatedAnnealing_base_eds.cc
 * Author: bschroed
 *
 * Created on April 18, 2018, 3:38 PM
 * Modified June 18, 2021 - bschroed, srieder
 */
#include <replicaExchange/replica_exchangers/1D_S_HSA_EDS/hamiltonian_simulatedAnnealing_base_eds.h>

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

re::hamiltonian_simulatedAnnealing_base_eds::hamiltonian_simulatedAnnealing_base_eds(io::Argument _args,
                                                            unsigned int cont, 
                                                            unsigned int globalThreadID, 
                                                            replica_graph_control & replicaGraphMPIControl, 
                                                            simulation::MpiControl & replica_mpi_control):  
                            replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
                            reedsParam(replica->sim.param().reeds)
{
    DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":Constructor:\t START" );
    DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":Constructor:\t simID "<<simulationID);

    //RE-Vars
    DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":Constructor:\t setParams" );
    setParams();

    DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":Constructor:\t DONE");
}

re::hamiltonian_simulatedAnnealing_base_eds::~hamiltonian_simulatedAnnealing_base_eds() {
    delete replica;
}

//SWAPPING Functions
void re::hamiltonian_simulatedAnnealing_base_eds::determine_switch_probabilities(){
    DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":determineSwitchPos:\t START");
    swap_replicas_1D(partnerReplicaID);
    DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":determineSwitchPos:\t DONE");
}

void re::hamiltonian_simulatedAnnealing_base_eds::swap_replicas_1D(const unsigned int partnerReplicaID) {
    DEBUG(4, "hamiltonian_simulatedAnnealing_base_eds:replica "<<  globalThreadID <<":swap1D - S:\t  START");
    DEBUG(4, "hamiltonian_simulatedAnnealing_base_eds:replica "<<  globalThreadID <<":swap1D - S:\t  sim: "<< simulationID << " \t\t "<< partnerReplicaID);
    unsigned int partnerReplicaMasterThreadID = partnerReplicaID;
    unsigned int num_s = replica->sim.param().reeds.num_l;

    // does partner exist?
    // the one with the higher ID does the probability calculation

    if (simulationID > partnerReplicaID) {
            DEBUG(4, "hamiltonian_simulatedAnnealing_base_eds:replica "<<  globalThreadID <<":swap1D:\t  swap prop calculation :)");

            // posts a MPI_Recv(...) matching the MPI_Send below
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
    DEBUG(1, "hamiltonian_simulatedAnnealing_base_eds:hamiltonian_simulatedAnnealing_base_eds:replica "<< globalThreadID <<":swap1D:\t  DONE");
}

double re::hamiltonian_simulatedAnnealing_base_eds::calc_probability(const unsigned int partnerReplicaID) {
    /*
     * This function is calculating the probability of the the coordinates to get pushed to the next higher s-Value.
     */
    DEBUG(1,"\n\nreplica "<<globalThreadID<<":replica_exchange_base_interface: CALC_PROBABILITY by ID:" << simulationID << "\n");

    unsigned int partnerReplicaMasterThreadID = partnerReplicaID;

    double delta;
    const double b1 = 1.0 / (math::k_Boltzmann * T);
    const double b2 = 1.0 / (math::k_Boltzmann * replica->sim.param().replica.temperature[partnerReplicaID % replica->sim.param().replica.num_T]);

    bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaID / replica->sim.param().replica.num_T]);// horizontal switch with same lambda?

    if (sameLambda) {
        // use simple formula
        // get energy from the other partner
        delta = (b1 - b2)*(0); //*  (E21 - E11 = 0 if same lambda)

    } else {
        // 2D formula
        /*
         * E12: Energy with lambda from 1 and configuration from 2
         * delta = b1 * ( E22 - E11) - b2*(E21  - E12);
         * E22 and E12 needed from partner
         */
        const double E11 = epot;
        const double E21 = calculate_energy(partnerReplicaID);

        // store this as the partner energy
        delta = b1 * E11 - b2 * E21;
    }

    // NPT? add PV term
    if (replica->sim.param().pcouple.scale != math::pcouple_off) {
        math::Box box_partner = replica->conf.current().box;
#ifdef XXMPI
        MPI_Status status;
        MPI_Recv(&box_partner(0)[0], 1, MPI_BOX, partnerReplicaMasterThreadID, BOX,  replicaGraphMPIControl().comm, &status);
#endif
        double V1 = math::volume(replica->conf.current().box, replica->conf.boundary_type);
        double V2 = math::volume(box_partner, replica->conf.boundary_type);
        //std::cerr << "volumes: " << V1 << " " << V2 << std::endl;

        // isotropic!
        double pressure = (replica->sim.param().pcouple.pres0(0,0)
                           + replica->sim.param().pcouple.pres0(1,1)
                           + replica->sim.param().pcouple.pres0(2,2)) / 3.0;
        delta += pressure * (b1 - b2) * (V2 - V1);
    }

    if (delta < 0.0)
        return 1.0;
    else
        return exp(-delta);
}

void re::hamiltonian_simulatedAnnealing_base_eds::setParams(){
    // set some variables
    stepsPerRun = replica->sim.param().step.number_of_steps;
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t NUMBER OF STEPS "<<stepsPerRun);

    run = 0;
    total_runs = replica->sim.param().reeds.trials + replica->sim.param().reeds.equilibrate;
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t NUMBER OF total_runs "<<total_runs);

    partnerReplicaID = simulationID;
    time = replica->sim.time();
    steps = 0;
    switched = 0;
    replica->curentStepNumber=0;
    replica->totalStepNumber = total_runs*stepsPerRun;
    replica->stepsPerRun= stepsPerRun;

    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t PARAM START");

    T = replica->sim.param().reeds.temperature;
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t got  T " << T);
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t got simulationID: "<< simulationID);

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
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t got s" << l);

    dt = replica->sim.param().step.dt;
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t dt " <<dt);
    DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":setParams:\t PARAM DONE ");    
}

void re::hamiltonian_simulatedAnnealing_base_eds::set_s() {
  DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":set_s:\t START ");

  eds_para = replica->sim.param().reeds.eds_para[simulationID];
  replica->sim.param().eds = eds_para;
  DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":set_s:\t eds_para s size: " << replica->sim.param().eds.s.size());
  l = replica->sim.param().eds.s[0];    //todoAssume only 1s EDS
  DEBUG(4,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":set_s:\t DONE " );
}

void re::hamiltonian_simulatedAnnealing_base_eds::init() {
  DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":init:\t init \t START");
  DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":init:\t start init from baseclass \t NEXT");
  //replica->init();
  DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":init:\t init_eds_stat \t NEXT");
  init_eds_stat();
  DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":init:\t DONE");
}

//initialize output files
void re::hamiltonian_simulatedAnnealing_base_eds::init_eds_stat(){
        DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":init_eds_stat:\t START");

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

        DEBUG(3,"hamiltonian_simulatedAnnealing_base_eds "<< globalThreadID <<":init_eds_stat:\t DONE");
}

//RE
////exchange params
void re::hamiltonian_simulatedAnnealing_base_eds::reset_eds() {//only reset switched parameters of change_eds() function
  replica->sim.param().eds = eds_para;
  replica->sim.param().step.dt = dt;
}

void re::hamiltonian_simulatedAnnealing_base_eds::change_eds(const unsigned int partner){//only change parameters, which are needed for energy calculation i.e.

  int idx;
  if (replica->sim.param().reeds.num_l == 1){
    idx = 0;
  }
  else{
    idx = partner;
  }

  replica->sim.param().step.dt = replica->sim.param().reeds.dt[idx];
  replica->sim.param().eds= replica->sim.param().reeds.eds_para[idx];
}

//Execute Swapping
/**OVERRIDE Possibly THIS NICE FUNCTION*/
void re::hamiltonian_simulatedAnnealing_base_eds::execute_swap(const unsigned int partnerReplicaID) {
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":executeSwap:\t START");
    if (simulationID > partnerReplicaID) {
        send_coord(partnerReplicaID);
        replica->conf.exchange_state();
    } else {
        receive_new_coord(partnerReplicaID);
    }
    // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
    exchange_averages();
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":executeSwap:\t DONE");
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
double re::hamiltonian_simulatedAnnealing_base_eds::calc_energy_eds_stat(double s){
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

double re::hamiltonian_simulatedAnnealing_base_eds::calculate_energy_core() {

    force_orig = replica->conf.current().force;
    virial_tensor_orig = replica->conf.current().virial_tensor;

    double energy = 0.0;
    algorithm::Algorithm * ff;

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

double re::hamiltonian_simulatedAnnealing_base_eds::calculate_energy(const unsigned int selectedReplicaID) {
    DEBUG(4, "replica_reeds_base_edsreplica_reeds "<< globalThreadID <<":calculate_energy:\t START"); 
 

    DEBUG(5, "replica_reeds_base_edsreplica_reeds "<< globalThreadID <<":calculate_energy:\t get Partner settings"); 
    if(selectedReplicaID!=simulationID){ 
        change_eds(selectedReplicaID);
    }

    double energy =  calculate_energy_core();

    if(selectedReplicaID!=simulationID){
        reset_eds();
    }


    DEBUG(4, "replica_reeds_base_edsreplica_reeds "<< globalThreadID <<":calculate_energy:\t DONE"); 
    return energy;
}

