/* 
 * File:   replica_exchange_base_2d_l_T_HREMD.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include "replicaExchange/replica/replica.h"
#include "replicaExchange/replica/replica_MPI_master.h"
#include "replicaExchange/replica/replica_MPI_slave.h"
#include "replicaExchange/replica_exchangers/2D_T_lambda_REPEX/replica_exchange_base_2d_l_T_HREMD.h"

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

re::replica_exchange_base_2d_l_T_HREMD::replica_exchange_base_2d_l_T_HREMD(io::Argument _args,
                                                                           unsigned int cont,
                                                                           unsigned int globalThreadID,
                                                                           replica_graph_control & replicaGraphMPIControl,
                                                                           simulation::MpiControl & replica_mpi_control) :
        replica_exchange_base_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control){
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t START ");
  DEBUG(3,"replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t SIMULATIONID:  "<< simulationID);

  //construct replica obj
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t  createReplica START");
  createReplicas(cont, globalThreadID, replica_mpi_control);
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t createdReplica DONE" << replica);
  T=replica->sim.param().replica.temperature[replicaGraphMPIControl.threadID % replica->sim.param().replica.num_T];
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t setTemperature: T= " << T);

  //RE-Vars
  DEBUG(3,"replica_exchange_base_eds "<< globalThreadID <<":Constructor:\t setParams" );
  setParams();  
  
  DEBUG(3,"replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":Constructor:\t Constructor \t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}


re::replica_exchange_base_2d_l_T_HREMD::~replica_exchange_base_2d_l_T_HREMD() {
    delete replica;
}

void re::replica_exchange_base_2d_l_T_HREMD::setParams(){
    // set some variables
    stepsPerRun = replica->sim.param().step.number_of_steps;
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t NUMBER OF STEPS "<<stepsPerRun);

    run = 0;
    total_runs = replica->sim.param().replica.trials + replica->sim.param().replica.equilibrate;
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t NUMBER OF total_runs "<<total_runs);

    partnerReplicaID = simulationID;
    time = replica->sim.time();
    steps = 0;
    switched = 0;
    replica->curentStepNumber=0;
    replica->totalStepNumber = total_runs*stepsPerRun;
    replica->stepsPerRun= stepsPerRun;

    const int numT = replica->sim.param().replica.num_T;
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t PARAM START");

    T = replica->sim.param().replica.temperature[simulationID % numT];
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t T START");

    l = replica->sim.param().replica.lambda[simulationID / numT];
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t Lamba");

    dt = replica->sim.param().replica.dt[simulationID / numT];
    DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t PARAM DONE ");

    set_lambda();
    set_temp();
}

//SWAPPING Functions
void re::replica_exchange_base_2d_l_T_HREMD::determine_switch_probabilities(){
    DEBUG(3,"replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":determineSwitchPos:\t START");
    DEBUG(3,"replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":determineSwitchPos:\t partner: "<<partnerReplicaID);
    swap_replicas_2D(partnerReplicaID);
    DEBUG(3,"replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":determineSwitchPos:\t DONE");
}

void re::replica_exchange_base_2d_l_T_HREMD::swap_replicas_2D(const unsigned int partnerReplicaMasterThreadID) {
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<  globalThreadID <<":swap_replicas_2D:\t  START");
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD " << globalThreadID << ":swap_replicas_2D:\t  sim: " << simulationID << " \t\t " << partnerReplicaMasterThreadID);

  unsigned int numT = replica->sim.param().replica.num_T;
  unsigned int numL = replica->sim.param().replica.num_l;

  if (partnerReplicaMasterThreadID < numT * numL && partnerReplicaMasterThreadID != simulationID) {  // does partner exist? Does it really make sense to check this here?!
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaMasterThreadID) {
        DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<  globalThreadID <<":swap_replicas_2D:\t  swap because smaller :)");

      // posts a MPI_Recv(...) matching the MPI_Send below
      probability = calc_probability(partnerReplicaMasterThreadID);
      const double randNum = rng.get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;
      DEBUG(7, "replica_exchange_base_2d_l_T_HREMD "<<  globalThreadID <<":swap_replicas_2D:\t  probabilities:" <<probability << "  rng: "<< randNum);


#ifdef XXMPI
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS, replicaGraphMPIControl().comm);
#endif

      if (randNum < probability) {
        switched = true;
      } else{
        switched = false;
      }
    } else {    //The Partner sends his data to The calculating Thread
      //special case if lambda also needs to be exchanged
      bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaMasterThreadID / replica->sim.param().replica.num_T]);
      if(!sameLambda){      //exchange LAMBDA
        // E21: Energy with configuration 2 and lambda 1(of partner)
        const double E21 = calculate_energy(partnerReplicaMasterThreadID);
        // this we can store as the partner energy of the current replica
        epot_partner = E21;
        // E22: Energy with configuration 2 and lambda 2(own one)
#ifdef XXMPI
        const double E22 = epot;
        // send E21 and E22
        double energies[2] = {E22, E21};
        //this send operation is matched in calc_probability()
        DEBUG(1,"\n\nreplica "<<  globalThreadID <<"replica_exchange_base_interface: SWAP_2D before Send\n");
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
        DEBUG(1,"\n\nreplica "<<  globalThreadID <<"replica_exchange_base_interface: SWAP_2D after Send\n");
#endif
      } else { // sameLambda
#ifdef XXMPI
        double energies[2] = {epot, 0.0};
        MPI_Send(&energies[0],2,MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
#endif
     }
      if (replica->sim.param().pcouple.scale != math::pcouple_off) {
#ifdef XXMPI
        math::Box box_replica = replica->conf.current().box;    //exchange box
        MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerReplicaMasterThreadID, BOX,  replicaGraphMPIControl().comm);
#endif
      }

      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
      MPI_Status status;
      MPI_Recv(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS,  replicaGraphMPIControl().comm, &status);
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

  } else {//This should be an error!
    //throw "Partner does not exist!";
    DEBUG(1, "replica_exchange_base_2d_l_T_HREMD "<<  globalThreadID <<":swap2D:\t  No swap because edgy");
    switched = false;
    probability = 0.0;
  }
  DEBUG(1, "replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":swap2D:\t resultValues: exchange? "<<switched << " p: "<< probability);
  DEBUG(1, "replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":swap2D:\t  DONE");
}

////Energy Determination
int re::replica_exchange_base_2d_l_T_HREMD::find_partner() const {
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<globalThreadID<<":find_partner \tSTART\n");

    unsigned int numT = replica->sim.param().replica.num_T;
    unsigned int numL = replica->sim.param().replica.num_l;

    unsigned int ID = simulationID;
    unsigned int partner = -1;
    bool even = ID % 2 == 0;

    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<globalThreadID<<":find_partner \tnumL "<<numL<<"\n");
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<globalThreadID<<":find_partner \tnumT "<<numT<<"\n");

    if (numT == 1 || numL == 1) { // 1D-RE with T or lam
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
        if ((numT == 1 && partner > numL - 1) || (numL == 1 && partner > numT - 1) || partner < 0) {    // solves edge cases
            partner = ID;
        }
    } else { // 2D-RE with T and Lam!
        bool firstElement = (ID % numT == 0);
        bool lastElement = (ID % numT == numT - 1);
        bool numTeven = (numT % 2 == 0);
        bool evenRow = (ID / numT) % 2 == 0;

        // determine switch direction
        DEBUG(5, "replica_exchange_base_2d_l_T_HREMD "<<globalThreadID<<":find_partner \t2D-ExchangeCase: "<<run % 4<<"\teven: "<< numTeven << "\n");
        switch (run % 4) {
            case 0:{ // Temperature direction
                if (numTeven) {
                    if (even){
                        partner = ID - 1;
                    }
                    else{
                        partner = ID + 1;
                    }
                    if (firstElement || lastElement){
                        partner = ID;
                    }
                } else {
                    if (evenRow) {
                        if (even){
                            partner = ID - 1;
                        }
                        else{
                            partner = ID + 1;
                        }
                    } else {
                        if (even){
                            partner = ID + 1;
                        }
                        else{
                            partner = ID - 1;
                        }
                    }
                    if (firstElement){
                        partner = ID;
                    }
                }
                break;
            }
            case 1:{ // lam-direction
                if (evenRow){
                    partner = ID + numT;
                }
                else{
                    partner = ID - numT;
                }
                break;
            }
            case 2:{ // temperature direction
                if (numTeven) {// Todo: @bschroed - I do not understand Why this needs to be destinguished / or is this the temperature dir?
                    if (even){
                        partner = ID + 1;
                    }
                    else{
                        partner = ID - 1;
                    }
                } else {
                    if (evenRow) {
                        if (even){
                            partner = ID + 1;
                        }
                        else{
                            partner = ID - 1;
                        }
                        if (lastElement){
                            partner = ID;
                        }
                    } else {
                        if (even){
                            partner = ID - 1;
                        }
                        else{
                            partner = ID + 1;
                        }
                    }
                    if (lastElement){
                        partner = ID;
                    }
                }
                break;
            }
            case 3:{// lam-direction WHY? //generates error
                if (evenRow) {
                    partner = ID - numT;
                }
                else {
                    partner = ID + numT;
                }
                break;
            }
        }

        // lambda-exchange edge cases: This is required for edge cases in lambda exchange
        if (run % 4 == 3 || run % 4 == 1 ){
            if (partner < 0 || partner > numT * numL - 1){
                partner = ID;
            }
        }
    }

    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<globalThreadID<<":find_partner \tpartner: "<<partner<<"\n");
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<<globalThreadID<<":find_partner \tDONE\n\n");

    return partner;
}


//Lambda Exchange
void re::replica_exchange_base_2d_l_T_HREMD::set_lambda() {
    DEBUG(5,"\n\nreplica_exchange_base_interface: SET_LAMBDA\n\n");
    // change Lambda in simulation
    replica->sim.param().perturbation.lambda = l;
    replica->topo.lambda(l);
    // twice, to set old_lambda... STILL NEEDED
    replica->topo.lambda(l);
    replica->topo.update_for_lambda();
    // set corresponding timestep
    replica->sim.param().step.dt = dt;
}

void re::replica_exchange_base_2d_l_T_HREMD::set_temp() {
    DEBUG(5,"\n\nreplica_exchange_base_interface: SET_TEMP\n\n");
    // change T in simulation
    replica->sim.param().stochastic.temp = T;

    for (unsigned int i = 0; i < replica->sim.multibath().size(); ++i) {
        assert(&replica->sim.multibath()[i].temperature != 0);
        replica->sim.multibath()[i].temperature = T;
    }
}

void re::replica_exchange_base_2d_l_T_HREMD::change_lambda(const unsigned int partnerReplicaID) {
    DEBUG(5,"\n\nreplica_exchange_base_interface: CHANGE_LAMBDA\n\n");
    int idx = 0;
    if (replica->sim.param().replica.num_l == 1){
        idx = 0;
    }
    else{
        idx = partnerReplicaID / replica->sim.param().replica.num_T;
    }

    const double lambda = replica->sim.param().replica.lambda[idx];
    const double dt = replica->sim.param().replica.dt[idx];
    // change Lambda in simulation
    replica->sim.param().perturbation.lambda = lambda;
    replica->topo.lambda(lambda);
    // twice, to set old_lambda... STILL NEEDED
    replica->topo.lambda(lambda);
    replica->topo.update_for_lambda();
    // set corresponding timestep
    replica->sim.param().step.dt = dt;
}

void re::replica_exchange_base_2d_l_T_HREMD::change_temp(const unsigned int partnerReplicaID) {
    DEBUG(5,"\n\nreplica_exchange_base_interface: CHANGE_TEMP\n\n");
    int idx = 0;
    if (replica->sim.param().replica.num_T == 1){
        idx = 0;
    }
    else{
        idx = partnerReplicaID % replica->sim.param().replica.num_T;
    }
    const double temp = replica->sim.param().replica.temperature[idx];
    // change T in simulation
    replica->sim.param().stochastic.temp = temp;

    for (unsigned int i = 0; i < replica->sim.multibath().size(); ++i) {
        assert(&replica->sim.multibath()[i].temperature != 0);
        replica->sim.multibath()[i].temperature = temp;
    }
}

double re::replica_exchange_base_2d_l_T_HREMD::calculate_energy(const unsigned int partnerThreadID) {
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":calculate_energy:\t  START");

    change_lambda(partnerThreadID);
    double energy = replica->calculateEnergies();
    set_lambda();
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":calculate_energy:\t  DONE");

    return energy;
}

void re::replica_exchange_base_2d_l_T_HREMD::calculate_energy_helper(const unsigned int partnerThreadID) {
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":calculate_energy_helper:\t  START");

    change_lambda(partnerThreadID);
    //change_temp(partnerThreadID); //Todo: Missing???
    re::replica_MPI_Slave * rs = dynamic_cast<re::replica_MPI_Slave* >(replica);
    rs->calculateEnergiesHelper();
    set_lambda();
    //set_temp(); //Todo: Missing???
    DEBUG(4, "replica_exchange_base_2d_l_T_HREMD "<< globalThreadID <<":calculate_energy_helper:\t  DONE");

}