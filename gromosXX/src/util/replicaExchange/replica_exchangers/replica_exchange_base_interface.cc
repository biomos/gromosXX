/* 
 * File:   replica_exchange_base.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include "replicaExchange/replica_mpi_tools.h"
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
                                                   replica_graph_mpi_control replicaGraphMPIControl,
                                                   simulation::mpi_control_struct replica_mpi_control) :
        args(_args), rng(-1), cont(cont), 
        globalThreadID(globalThreadID), simulationID(replica_mpi_control.simulationID),
        replicaGraphMPIControl(replicaGraphMPIControl)
{
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t START ");
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t SIMULATIONID:  "<< simulationID);

  //construct replica obj
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t  createReplica");
  createReplicas(cont, globalThreadID, replica_mpi_control);
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t createdReplica T");

  //RE-Vars
  // set some variables
  stepsPerRun = replica->sim.param().step.number_of_steps;
  run = 0;
  total_runs = replica->sim.param().replica.trials + replica->sim.param().replica.equilibrate;
  partnerReplicaID = simulationID;
  time = replica->sim.time();
  steps = 0;
  switched = 0;
  replica->curentStepNumber=0;
  replica->totalStepNumber = total_runs*stepsPerRun;
  replica->stepsPerRun= stepsPerRun;

  const int numT = replica->sim.param().replica.num_T;

  T = replica->sim.param().replica.temperature[simulationID % numT];
  l = replica->sim.param().replica.lambda[simulationID / numT];
  dt = replica->sim.param().replica.dt[simulationID / numT];

  set_lambda();
  set_temp();
  
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t Constructor \t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}


util::replica_exchange_base_interface::~replica_exchange_base_interface() {
    delete replica;
}

//TODO: REMOVE
void util::replica_exchange_base_interface::createReplicas(int cont, int globalThreadID, simulation::mpi_control_struct replica_mpi_control){
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":createReplicas:\t START \t THREADS "<<replica_mpi_control.simulationNumberOfThreads);
  // create the number of replicas that are assigned to my node
    if(replica_mpi_control.simulationNumberOfThreads>1){
        if(replica_mpi_control.simulationThisThreadID == replica_mpi_control.simulationMasterThreadID ){
            replica = new util::replica_MPI_Master(args, cont, globalThreadID, replica_mpi_control);
        }
        else{
            replica = new util::replica_MPI_Slave(args, cont, globalThreadID, replica_mpi_control);
        }
    }
    else{
        replica = new util::replica(args, cont, globalThreadID, replica_mpi_control);
        //DEBUG(2,"replica_exchange_base_interface "<< globalThreadID <<":create:\t rep" << replica->sim.param().replica.num_l);
        //DEBUG(2,"replica_exchange_base_interface "<< globalThreadID <<":create:\t rep" << replica->sim.param().replica.lambda[0]);
    }
   
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":createReplicas:\t DONE");
}

//TODO: REMOVE
void util::replica_exchange_base_interface::init() {
  DEBUG(3, "replica_exchange_base_interface "<< globalThreadID <<":init:\t START");
  // do init for all replica assigned to this node
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":init:\t initReplicas");
  replica->init();
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":init:\t DONE");
}

//TODO: REMOVE
void util::replica_exchange_base_interface::run_MD() {
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":run_MD:\t START");
  
    replica->sim.steps() = steps;
    replica->sim.time() = time;

    replica->run_MD();
    
    updateReplica_params();
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":run_MD:\t END");
}

//TODO: REMOVE
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

//TODO: REMOVE
void util::replica_exchange_base_interface::swap(){
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":swap:\t START");
  
    partnerReplicaID = find_partner();
    
    if (partnerReplicaID != simulationID) // different replica?
    {
      //TODO: RENAME 
      swap_replicas_2D(partnerReplicaID);
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
     
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":swap:\t DONE");
}

//TODO: REMOVE
void util::replica_exchange_base_interface::write_final_conf() {
  // write coordinates to cnf for all replica assigned to this node
   replica->write_final_conf();
}


//RE
//SWAPPING Functions
void util::replica_exchange_base_interface::swap_replicas_2D(const unsigned int partnerReplicaID) {
  DEBUG(4, "replica "<<  globalThreadID <<":swap:\t  START");

  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;
  unsigned int numT = replica->sim.param().replica.num_T;
  unsigned int numL = replica->sim.param().replica.num_l;

  // does partner exist?
  if (partnerReplicaID < numT * numL && partnerReplicaID != simulationID) {
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaID) {
    
      // posts a MPI_Recv(...) matching the MPI_Send below 
      probability = calc_probability(partnerReplicaID);
      const double randNum = rng.get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;

#ifdef XXMPI
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS, replicaGraphMPIControl.replicaGraphCOMM);
#endif

      if (randNum < probability) {
        switched = true;
      } else
        switched = false;
    } else {    //The Partner sends his data to The calculating Thread
      //special case if lambda also needs to be exchanged
      bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaID / replica->sim.param().replica.num_T]);
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
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
      } else { // sameLambda
#ifdef XXMPI
        double energies[2] = {epot, 0.0}; 
        MPI_Send(&energies[0],2,MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
     }
      if (replica->sim.param().pcouple.scale != math::pcouple_off) {
#ifdef XXMPI
        math::Box box_replica = replica->conf.current().box;    //exchange box
        MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerReplicaMasterThreadID, BOX,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
      }
      
#ifdef XXMPI
      MPI_Status status;
#endif
      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
      MPI_Recv(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS,  replicaGraphMPIControl.replicaGraphCOMM, &status);
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
      throw "Partner does not exist!";
    /*
      partner = ID;
    switched = false;
    probability = 0.0;
      
    */
  }
    DEBUG(4, "replica "<< globalThreadID <<":swap:\t  DONE");
}

void util::replica_exchange_base_interface::swap_replicas_1D(const unsigned int partnerReplicaID) {
  DEBUG(4, "replica "<<  globalThreadID <<":swap:\t  START");

  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;
  unsigned int numP = replica->sim.param().replica.num_T;

  // does partner exist?
  if (partnerReplicaID < numP && partnerReplicaID != simulationID) {
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaID) {
    
      // posts a MPI_Recv(...) matching the MPI_Send below 
      probability = calc_probability(partnerReplicaID);
      const double randNum = rng.get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;

#ifdef XXMPI
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS, replicaGraphMPIControl.replicaGraphCOMM);
#endif

      if (randNum < probability) {
        switched = true;
      } else
        switched = false;
    } else {    //The Partner sends his data to The calculating Thread
      
        //special case if lambda also needs to be exchanged
      bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaID / replica->sim.param().replica.num_T]);
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
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
      } else { // sameLambda
#ifdef XXMPI
        double energies[2] = {epot, 0.0}; 
        MPI_Send(&energies[0],2,MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
     }
      if (replica->sim.param().pcouple.scale != math::pcouple_off) {
#ifdef XXMPI
        math::Box box_replica = replica->conf.current().box;    //exchange box
        MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerReplicaMasterThreadID, BOX,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
      }
      
#ifdef XXMPI
      MPI_Status status;
#endif
      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
      MPI_Recv(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS,  replicaGraphMPIControl.replicaGraphCOMM, &status);
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
      throw "Partner does not exist!";
    /*
      partner = ID;
    switched = false;
    probability = 0.0;
      
    */
  }
    DEBUG(4, "replica "<< globalThreadID <<":swap:\t  DONE");
}

int util::replica_exchange_base_interface::find_partner() const {
    //TODO: REWRITE To get Replica ID BSCHROED
  unsigned int numT = replica->sim.param().replica.num_T;
  unsigned int numL = replica->sim.param().replica.num_l;
  
  unsigned int ID = simulationID;

  unsigned int partner = ID; //makes sense why?
  bool even = ID % 2 == 0;
  bool evenRow = (ID / numT) % 2 == 0;
  bool firstElement = (ID % numT == 0);
  bool lastElement = (ID % numT == numT - 1);
  bool numTeven = (numT % 2 == 0);

  // only 1D-RE ?
  if (numT == 1 || numL == 1) {
    if (run % 2 == 1) {
      if (even)
        partner = ID + 1;
      else
        partner = ID - 1;
    } else {
      if (even)
        partner = ID - 1;
      else
        partner = ID + 1;
    }
  } else { // 2D-RE
    // determine switch direction
    switch (run % 4) {
      case 0: // lambda direction 
        if (numTeven) {
          if (even)
            partner = ID - 1;
          else
            partner = ID + 1;
          if (firstElement || lastElement)
            partner = ID;
        } else {
          if (evenRow) {
            if (even)
              partner = ID - 1;
            else
              partner = ID + 1;
          } else {
            if (even)
              partner = ID + 1;
            else
              partner = ID - 1;
          }
          if (firstElement)
            partner = ID;
        }
        break;

      case 1: // temp-direction
        if (evenRow)
          partner = ID + numT;
        else
          partner = ID - numT;
        break;

      case 2: // lambda direction
        if (numTeven) {
          if (even)
            partner = ID + 1;
          else
            partner = ID - 1;
        } else {
          if (evenRow) {
            if (even)
              partner = ID + 1;
            else
              partner = ID - 1;
            if (lastElement)
              partner = ID;
          } else {
            if (even)
              partner = ID - 1;
            else
              partner = ID + 1;
          }
          if (lastElement)
            partner = ID;
        }
        break;
      case 3:// temp-direction WHY? 
        if (evenRow)
          partner = ID - numT;
        else
          partner = ID + numT;
        break;
    }
  }
  // partner out of range ? - Do we really need this or is it more a hack hiding bugs?
  if (partner > numT * numL - 1)
    partner = ID;

  return partner;
}

double util::replica_exchange_base_interface::calc_probability(const unsigned int partnerReplicaID) {
    
  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;

  double delta;
  const double b1 = 1.0 / (math::k_Boltzmann * T);
  const double b2 = 1.0 / (math::k_Boltzmann * replica->sim.param().replica.temperature[partnerReplicaID % replica->sim.param().replica.num_T]);
  
  bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaID / replica->sim.param().replica.num_T]);// horizontal switch with same lambda?

  if (sameLambda) {
    // use simple formula
    // get energy from the other partner
    double energies[2] = {0.0, 0.0};
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl.replicaGraphCOMM, &status);
#endif

    epot_partner = energies[0];
    delta = (b1 - b2)*(epot_partner - epot); //*  (E21 - E11=
    
  } else {
    // 2D formula
    /*
     * E12: Energy with lambda from 1 and configuration from 2
     * delta = b1 * ( E22 - E11) - b2*(E21  - E12);
     * E22 and E12 needed from partner
     */
    double energies[2] = {0.0, 0.0};
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl.replicaGraphCOMM, &status);
#endif
    const double E22 = energies[0];
    const double E12 = energies[1];

    const double E11 = epot;
    const double E21 = calculate_energy(partnerReplicaID);
    
    // store this as the partner energy 
    epot_partner = E21;

    // Chris: I think this is wrong
    // delta = b1 * (E22 - E11) - b2 * (E21 - E12);
    //std::cerr << "b1: " << b1 << " b2: " << b2 << std::endl;
    //std::cerr << "E11: " << E11 << " E22: " << E22 << std::endl;
    //std::cerr << "E21: " << E21 << " E12: " << E12 << std::endl;

    delta = b1 * (E12 - E11) - b2 * (E22 - E21);
  }
  
  // NPT? add PV term
  if (replica->sim.param().pcouple.scale != math::pcouple_off) {
    math::Box box_partner = replica->conf.current().box;
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&box_partner(0)[0], 1, MPI_BOX, partnerReplicaMasterThreadID, BOX,  replicaGraphMPIControl.replicaGraphCOMM, &status);
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

//THIS COULD GO TO REPLICA and into the func above!
double util::replica_exchange_base_interface::calculate_energy_core() {
  double energy = 0.0;
  //chris: you do need to re-evaluate the energy, otherwise you get the energy of before the previous step
  algorithm::Algorithm * ff = replica->md.algorithm("Forcefield");

  if (ff->apply(replica->topo, replica->conf, replica->sim)) {
    print_info("Error in energy calculation!");
 #ifdef XXMPI
    MPI_Abort( replicaGraphMPIControl.replicaGraphCOMM, E_UNSPECIFIED);
#endif
    return 1;
  }

  replica->conf.current().energies.calculate_totals();
  switch (replica->sim.param().xrayrest.replica_exchange_parameters.energy_switcher) {
    case simulation::energy_tot:
      energy = replica->conf.current().energies.potential_total + replica->conf.current().energies.special_total;
      break;

    case simulation::energy_phys:
      energy = replica->conf.current().energies.potential_total;
      break;
    case simulation::energy_special:
      energy = replica->conf.current().energies.special_total;
      break;
    default:
      print_info("Error in energy switching!");
#ifdef XXMPI
      MPI_Abort(replicaGraphMPIControl.replicaGraphCOMM, E_UNSPECIFIED);
#endif
  }
  return energy;
}


double util::replica_exchange_base_interface::calculate_energy(const unsigned int partnerThreadID) {
  DEBUG(4, "replica "<< globalThreadID <<":calculate_energy:\t  START");

  change_lambda(partnerThreadID);

  double energy = calculate_energy_core();
  /*
   //THIS COULD GO be removed!
  algorithm::Algorithm* ff = replica->md.algorithm("Forcefield");

  if (ff->apply(replica->topo, replica->conf, replica->sim)) {
    print_info("Error in energy calculation!");
 #ifdef XXMPI
    MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
#endif
    return 1;
  }

  replica->conf.current().energies.calculate_totals();
  switch (replica->sim.param().xrayrest.replica_exchange_parameters.energy_switcher) {
    case simulation::energy_tot:
      energy = replica->conf.current().energies.potential_total + replica->conf.current().energies.special_total;
      break;
    case simulation::energy_phys:
      energy = replica->conf.current().energies.potential_total;
      break;
    case simulation::energy_special:
      energy = replica->conf.current().energies.special_total;
      break;
    default:
      std::cerr << "Something is wrong in energy calculation" << std::endl;
      print_info("Error in energy calculation!");
#ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
#endif
      return 1;
  }
   * */

  set_lambda();
  DEBUG(4, "replica "<< globalThreadID <<":calculate_energy:\t  DONE");

  return energy;
}


//TODO: REMOVE
//THIS COULD GO TO REPLICA and into the func above!
void util::replica_exchange_base_interface::exchange_averages() {
  // after a swap the averages of current and old are exchanged and have to be switched back
  configuration::Average  dummy = replica->conf.current().averages;
  replica->conf.current().averages=replica->conf.old().averages;
  replica->conf.old().averages=dummy; // <- why not removing?
}

//sending stuff
//TODO: REMOVE
void util::replica_exchange_base_interface::send_coord(const unsigned int receiverReplicaID) {
#ifdef XXMPI

  unsigned int receiverReplicaMasterThreadID = receiverReplicaID;

  configuration::Configuration  conf = replica->conf;
  
  MPI_Send(&conf.current().pos[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, POS,  replicaGraphMPIControl.replicaGraphCOMM);
  MPI_Send(&conf.current().posV[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, POSV,  replicaGraphMPIControl.replicaGraphCOMM);
  MPI_Send(&conf.current().vel[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, VEL,  replicaGraphMPIControl.replicaGraphCOMM);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if (simulationID < receiverReplicaID) {
    MPI_Send(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, LATTSHIFTS,  replicaGraphMPIControl.replicaGraphCOMM);
  } else {
    MPI_Send(&((*replica->latticeTMP)[0][0]), 1, MPI_VARRAY, receiverReplicaMasterThreadID, LATTSHIFTS,  replicaGraphMPIControl.replicaGraphCOMM);
    delete replica->latticeTMP;
    replica->latticeTMP = NULL;
  }

  MPI_Send(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, STOCHINT,  replicaGraphMPIControl.replicaGraphCOMM);
  MPI_Send(&conf.current().box(0)[0], 1, MPI_BOX, receiverReplicaMasterThreadID, BOX,  replicaGraphMPIControl.replicaGraphCOMM);

  std::vector<double> angles;
  angles.resize(3);
  angles[0] = conf.current().phi;
  angles[1] = conf.current().psi;
  angles[2] = conf.current().theta;

  MPI_Send(&angles[0], angles.size(), MPI_DOUBLE, receiverReplicaMasterThreadID, ANGLES,  replicaGraphMPIControl.replicaGraphCOMM);
  MPI_Send(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, receiverReplicaMasterThreadID, DF,  replicaGraphMPIControl.replicaGraphCOMM);
#endif
}

//TODO: REMOVE
void util::replica_exchange_base_interface::receive_new_coord(const unsigned int senderReplicaID) {
#ifdef XXMPI
   
  unsigned int senderReplicaMasterThreadID = senderReplicaID;

  MPI_Status status;
  configuration::Configuration  conf = replica->conf;

  conf.exchange_state();
  MPI_Recv(&conf.current().pos[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, POS,  replicaGraphMPIControl.replicaGraphCOMM, &status);
  MPI_Recv(&conf.current().posV[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, POSV,  replicaGraphMPIControl.replicaGraphCOMM, &status);
  MPI_Recv(&conf.current().vel[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, VEL,  replicaGraphMPIControl.replicaGraphCOMM, &status);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if ( globalThreadID > senderReplicaMasterThreadID){
    replica->latticeTMP = new math::VArray(conf.special().lattice_shifts);
  }

  MPI_Recv(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, LATTSHIFTS,  replicaGraphMPIControl.replicaGraphCOMM, &status);
  MPI_Recv(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, STOCHINT,  replicaGraphMPIControl.replicaGraphCOMM, &status);
  MPI_Recv(&conf.current().box(0)[0], 1, MPI_BOX, senderReplicaMasterThreadID, BOX,  replicaGraphMPIControl.replicaGraphCOMM, &status);

  std::vector<double> angles;
  angles.resize(3);
  MPI_Recv(&angles[0], angles.size(), MPI_DOUBLE, senderReplicaMasterThreadID, ANGLES,  replicaGraphMPIControl.replicaGraphCOMM, &status);

  conf.current().phi = angles[0];
  conf.current().psi = angles[1];
  conf.current().theta = angles[2];

  MPI_Recv(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, senderReplicaMasterThreadID, DF,  replicaGraphMPIControl.replicaGraphCOMM, &status);

  if ( simulationID > senderReplicaID){
    conf.exchange_state();
  }
#endif
}

//TRE
void util::replica_exchange_base_interface::velscale(int unsigned partnerReplica){ 
  double T1 = replica->sim.param().replica.temperature[simulationID]; 
  double T2 = replica->sim.param().replica.temperature[partnerReplica];
  if (T1 != T2) {
    double factor = sqrt(T1/T2);
    for (unsigned int k = 0; k < replica->topo.num_atoms(); ++k) {
     replica->conf.current().vel(k) *= factor;
    }
  } 
}

//Lambda Exchange (Kinda Hamiltonian Exchange)
void util::replica_exchange_base_interface::set_lambda() {
  // change Lambda in simulation
  replica->sim.param().perturbation.lambda = l;
  replica->topo.lambda(l);
  // twice, to set old_lambda... STILL NEEDED
  replica->topo.lambda(l);
  replica->topo.update_for_lambda();
  // set corresponding timestep
  replica->sim.param().step.dt = dt;
}

void util::replica_exchange_base_interface::set_temp() {
  // change T in simulation
  replica->sim.param().stochastic.temp = T;

  for (unsigned int i = 0; i < replica->sim.multibath().size(); ++i) {
    assert(&replica->sim.multibath()[i].temperature != 0);
    replica->sim.multibath()[i].temperature = T;
  }
}

void util::replica_exchange_base_interface::change_lambda(const unsigned int partnerReplicaID) {
  int idx;
  if (replica->sim.param().replica.num_l == 1)
    idx = 0;
  else
    idx = partnerReplicaID / replica->sim.param().replica.num_T;
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

void util::replica_exchange_base_interface::change_temp(const unsigned int partnerReplicaID) {
  int idx;
  if (replica->sim.param().replica.num_T == 1)
    idx = 0;
  else
    idx = partnerReplicaID % replica->sim.param().replica.num_T;
  const double temp = replica->sim.param().replica.temperature[idx];
  // change T in simulation
  replica->sim.param().stochastic.temp = temp;

  for (unsigned int i = 0; i < replica->sim.multibath().size(); ++i) {
    assert(&replica->sim.multibath()[i].temperature != 0);
    replica->sim.multibath()[i].temperature = temp;
  }
}

//TODO: REMOVE
// IO
void util::replica_exchange_base_interface::print_coords(std::string name) {
    
    io::Out_Configuration received_traj(GROMOSXX " Configuration", std::cout);
    io::Argument args2(args);
    args2.erase("fin");
    std::stringstream tmp;
    tmp << name << "_replica_" << simulationID << "_run_" << run;
    std::string fin = tmp.str() + ".cnf";
    args2.insert(std::pair<std::string, std::string > ("fin", fin));

    received_traj.init(args2, replica->sim.param());

    received_traj.title("bla");
    received_traj.write(replica->conf, replica->topo, replica->sim, io::final);
}

unsigned int util::replica_exchange_base_interface::get_stepsPerRun(){
    return stepsPerRun;
};

//TODO: REMOVE
void util::replica_exchange_base_interface::print_info(std::string bla) const {
  std::cout << "\n" << bla << std::endl;
  std::cout << "#"
          << std::setw(5) << "ID"
          << " "
          << std::setw(7) << "partner"
          << std::setw(7) << "run"

          << std::setw(13) << "li"
          << std::setw(13) << "Ti"
          << std::setw(14) << "Epoti"
          << std::setw(13) << "lj"
          << std::setw(13) << "Tj"
          << std::setw(14) << "Epotj"
          << std::setw(13) << "p"
          << std::setw(4) << "s"
          << "\n";
  //And this makes sense why again??? todo: bschroed
  std::cout << std::setw(6) << (simulationID+1)
          << " "
          << std::setw(6) << partnerReplicaID
          << std::setw(6) << run
          << std::setw(13) << l        //TODO into fitting Exchange class
          << std::setw(13) << T        //TODO into fitting Exchange class
          << " "
          << std::setw(18) << epot
          << std::setw(13) << l        //TODO into fitting Exchange class
          << std::setw(13) << T        //TODO into fitting Exchange class
          << " "
          << std::setw(18) << epot_partner
          << std::setw(13) << probability
          << std::setw(4) << switched
          << std::endl;
}
