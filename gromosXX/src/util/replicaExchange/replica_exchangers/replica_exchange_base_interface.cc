/* 
 * File:   replica_exchange_base.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include "replicaExchange/replica/replica.h"
#include "replicaExchange/replica/replica_MPI_master.h"
#include "replicaExchange/replica/replica_MPI_slave.h"
#include "replicaExchange/replica_exchangers/replica_exchange_base_interface.h"

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

util::replica_exchange_base_interface::replica_exchange_base_interface(io::Argument _args,
                                                   unsigned int cont, 
                                                   unsigned int globalThreadID,
                                                   std::vector<std::vector<unsigned int> >  replica_owned_threads, 
                                                   std::map<ID_t, rank_t> &thread_id_replica_map) : 
        args(_args),  rng(-1),
        cont(cont), 
        globalThreadID(globalThreadID), globalMasterThreadID(0), numReplicas(replica_owned_threads.size()),
        simulationID(thread_id_replica_map[globalThreadID]), simulationMasterID(replica_owned_threads[simulationID][0]), simulationThreads(replica_owned_threads[simulationID].size()), 
        simulationThreadID(std::find(replica_owned_threads[simulationID].begin(), replica_owned_threads[simulationID].end(), globalThreadID) - replica_owned_threads[simulationID].begin()), 
        replica_owned_threads(replica_owned_threads), thread_id_replica_map(thread_id_replica_map){
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t START ");
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t SIMULATIONID:  "<< simulationID);

  //construct replica obj
  DEBUG(4,"replica_exchange_base "<< globalThreadID <<":Constructor:\t  createReplica");
  createReplicas(cont, globalThreadID);
  DEBUG(4,"replica_exchange_base "<< globalThreadID <<":Constructor:\t createdReplica T");

  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":Constructor:\t Constructor \t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}

util::replica_exchange_base_interface::~replica_exchange_base_interface() {
    delete replica;
}

void util::replica_exchange_base_interface::createReplicas(int cont, int globalThreadID){
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":createReplicas:\t START");
  // create the number of replicas that are assigned to my node
    if(simulationThreads>1){
        if(simulationThreadID == 0){
            replica = new util::replica_MPI_Master(args, cont, simulationID, globalThreadID, simulationThreadID, simulationID, simulationThreads);
        }
        else{
            replica = new util::replica_MPI_Slave(args, cont, simulationID, globalThreadID, simulationThreadID, simulationID, simulationThreads);
        }
    }
    else{
        replica = new util::replica(args, cont, simulationID, globalThreadID);
        //DEBUG(2,"replica_exchange_base "<< globalThreadID <<":create:\t rep" << replica->sim.param().replica.num_l);
        //DEBUG(2,"replica_exchange_base "<< globalThreadID <<":create:\t rep" << replica->sim.param().replica.lambda[0]);
    }
   
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":createReplicas:\t DONE");
}

void util::replica_exchange_base_interface::init() {
  DEBUG(3, "replica_exchange_base "<< globalThreadID <<":init:\t START");
  // do init for all replica assigned to this node
  DEBUG(4,"replica_exchange_base "<< globalThreadID <<":init:\t initReplicas");
  replica->init();
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":init:\t DONE");
}

void util::replica_exchange_base_interface::run_MD() {
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":run_MD:\t START");
  
    replica->sim.steps() = steps;
    replica->sim.time() = time;

    replica->run_MD();
    
    updateReplica_params();
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":run_MD:\t END");
}

void util::replica_exchange_base_interface::updateReplica_params(){
  // update replica information
  time = replica->sim.time();
  steps = replica->sim.steps();
  ++run;
  epot = calculate_energy_core();   
}
/*
 * REplica Exchanges
 */

void util::replica_exchange_base_interface::swap(){
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":swap:\t START");
  
    partnerReplicaID = find_partner();
    
    if (partnerReplicaID != simulationID) // different replica?
      {
        swap_replicas_priv(partnerReplicaID);
        if (switched) {
          if (globalThreadID < partnerReplicaID) {
            send_coord(partnerReplicaID);
            receive_new_coord(partnerReplicaID);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
            exchange_averages();
          } else {
            receive_new_coord(partnerReplicaID);
            send_coord(partnerReplicaID);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
            exchange_averages();
          }
        }
      }
    else {  // no exchange with replica itself
      probability = 0.0;
      switched = 0;
    }
    if(switched && replica->sim.param().replica.scale) {
      velscale(partnerReplicaID);
    }
     
  DEBUG(3,"replica_exchange_base "<< globalThreadID <<":swap:\t DONE");
}

void util::replica_exchange_base_interface::write_final_conf() {
  // write coordinates to cnf for all replica assigned to this node
   replica->write_final_conf();
}


//RE
//SWAPPING Functions

void util::replica_exchange_base_interface::exchange_averages() {
  // after a swap the averages of current and old are exchanged and have to be switched back
  configuration::Average  dummy = replica->conf.current().averages;
  replica->conf.current().averages=replica->conf.old().averages;
  replica->conf.old().averages=dummy; // <- why not removing?
}

//sending stuff
void util::replica_exchange_base_interface::send_coord(const unsigned int receiverReplicaID) {
#ifdef XXMPI

  unsigned int receiverReplicaMasterThreadID = replica_owned_threads[receiverReplicaID][0];

  configuration::Configuration  conf = replica->conf;
  
  MPI_Send(&conf.current().pos[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, POS, MPI_COMM_WORLD);
  MPI_Send(&conf.current().posV[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, POSV, MPI_COMM_WORLD);
  MPI_Send(&conf.current().vel[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, VEL, MPI_COMM_WORLD);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if (simulationID < receiverReplicaID) {
    MPI_Send(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, LATTSHIFTS, MPI_COMM_WORLD);
  } else {
    MPI_Send(&((*replica->latticeTMP)[0][0]), 1, MPI_VARRAY, receiverReplicaMasterThreadID, LATTSHIFTS, MPI_COMM_WORLD);
    delete replica->latticeTMP;
    replica->latticeTMP = NULL;
  }

  MPI_Send(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, STOCHINT, MPI_COMM_WORLD);
  MPI_Send(&conf.current().box(0)[0], 1, MPI_BOX, receiverReplicaMasterThreadID, BOX, MPI_COMM_WORLD);

  std::vector<double> angles;
  angles.resize(3);
  angles[0] = conf.current().phi;
  angles[1] = conf.current().psi;
  angles[2] = conf.current().theta;

  MPI_Send(&angles[0], angles.size(), MPI_DOUBLE, receiverReplicaMasterThreadID, ANGLES, MPI_COMM_WORLD);
  MPI_Send(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, receiverReplicaMasterThreadID, DF, MPI_COMM_WORLD);
#endif
}

void util::replica_exchange_base_interface::receive_new_coord(const unsigned int senderReplicaID) {
#ifdef XXMPI
   
  unsigned int senderReplicaMasterThreadID = replica_owned_threads[senderReplicaID][0];

  MPI_Status status;
  configuration::Configuration  conf = replica->conf;

  conf.exchange_state();
  MPI_Recv(&conf.current().pos[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, POS, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().posV[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, POSV, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().vel[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, VEL, MPI_COMM_WORLD, &status);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if ( globalThreadID > senderReplicaMasterThreadID){
    replica->latticeTMP = new math::VArray(conf.special().lattice_shifts);
  }

  MPI_Recv(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, LATTSHIFTS, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, STOCHINT, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().box(0)[0], 1, MPI_BOX, senderReplicaMasterThreadID, BOX, MPI_COMM_WORLD, &status);

  std::vector<double> angles;
  angles.resize(3);
  MPI_Recv(&angles[0], angles.size(), MPI_DOUBLE, senderReplicaMasterThreadID, ANGLES, MPI_COMM_WORLD, &status);

  conf.current().phi = angles[0];
  conf.current().psi = angles[1];
  conf.current().theta = angles[2];

  MPI_Recv(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, senderReplicaMasterThreadID, DF, MPI_COMM_WORLD, &status);

  if ( simulationID > senderReplicaID){
    conf.exchange_state();
  }
#endif
}

// IO
void util::replica_exchange_base_interface::print_coords(std::string name) {
    
    io::Out_Configuration received_traj(GROMOSXX " Configuration", std::cout);
    io::Argument args2(args);
    args2.erase("fin");
    std::stringstream tmp;
    tmp << name << "_replica_" << replica->ID << "_run_" << run;
    std::string fin = tmp.str() + ".cnf";
    args2.insert(std::pair<std::string, std::string > ("fin", fin));

    received_traj.init(args2, replica->sim.param());

    received_traj.title("bla");
    received_traj.write(replica->conf, replica->topo, replica->sim, io::final);
}
