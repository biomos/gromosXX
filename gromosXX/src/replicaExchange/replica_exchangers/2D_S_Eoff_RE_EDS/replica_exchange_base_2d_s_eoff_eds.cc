/*
 * File:   replica_exchange_base_2d_s_eoff_eds.cc
 * Author: theosm
 *
 * Created on March 29, 2020, 11:08 AM
 */
#include <replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_base_2d_s_eoff_eds.h>

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

re::replica_exchange_base_2d_s_eoff_eds::replica_exchange_base_2d_s_eoff_eds(io::Argument _args,
                                                            unsigned int cont,
                                                            unsigned int globalThreadID,
                                                            replica_graph_control &replicaGraphMPIControl,
                                                            simulation::MpiControl &replica_mpi_control):
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

re::replica_exchange_base_2d_s_eoff_eds::~replica_exchange_base_2d_s_eoff_eds() {
    delete replica;
}

void re::replica_exchange_base_2d_s_eoff_eds::setParams(){
    // set some variables
    stepsPerRun = replica->sim.param().step.number_of_steps;
    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t NUMBER OF STEPS "<<stepsPerRun);

    run = 0;
    total_runs = replica->sim.param().replica.trials + replica->sim.param().replica.equilibrate;
    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t NUMBER OF total_runs "<<total_runs);

    partnerReplicaID = simulationID;
    time = replica->sim.time();
    steps = 0;
    switched = 0;
    replica->curentStepNumber=0;
    replica->totalStepNumber = total_runs*stepsPerRun;
    replica->stepsPerRun= stepsPerRun;

    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t PARAM START");

    T = replica->sim.param().reeds.temperature;
    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t got  T " << T);
    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t got simulationID: "<< simulationID);

    //set position info
    replica->sim.param().reeds.eds_para[simulationID].pos_info = std::make_pair(simulationID, simulationID);
    DEBUG(4, "BASE Constructor with simulationID, replica->pos_info= " << simulationID << ", "
    << replica->sim.param().reeds.eds_para[simulationID].pos_info.first << ", "
    << replica->sim.param().reeds.eds_para[simulationID].pos_info.second << "\n");

    pos_info = replica->sim.param().reeds.eds_para[simulationID].pos_info;
    DEBUG(4, "BASE Constructor with simulationID, pos_info= " << simulationID << ", "
    << pos_info.first << ", " << pos_info.second << "\n");

    //just to check -- theosm
    std::pair<int, int> a = reedsParam.eds_para[simulationID].pos_info;
    DEBUG(4, "JUST TO CHECK: BASE Constructor with simulationID, reedsParam->pos_info= " << simulationID << ", "
    << a.first << ", " << a.second << "\n");

    set_s();
    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t got s" << l);

    dt = replica->sim.param().step.dt;
    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t dt " <<dt);

    DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":setParams:\t PARAM DONE ");
}

void re::replica_exchange_base_2d_s_eoff_eds::set_s() {
  DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":set_s:\t START ");

  eds_para = replica->sim.param().reeds.eds_para[simulationID];
  replica->sim.param().eds = eds_para;
  DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":set_s:\t eds_para s size: " << replica->sim.param().eds.s.size());

  l = replica->sim.param().eds.s[0];    //todoAssume only 1s EDS
  DEBUG(4,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":set_s:\t DONE " );
}

void re::replica_exchange_base_2d_s_eoff_eds::init() {
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t init \t START");
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t start init from baseclass \t NEXT");
  //replica->init();
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t init_eds_stat \t NEXT");
  init_eds_stat();
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init:\t DONE");
}

//initialize output files
void re::replica_exchange_base_2d_s_eoff_eds::init_eds_stat(){
        DEBUG(3,"\n\nreplica_exchange_base_2d_s_eoff_eds: INIT_EDS_STAT\n\n");
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
        replicaStatData[currentID].epot_vec.resize(replicaGraphMPIControl().numberOfReplicas);
        replicaStatData[currentID].prob_vec.resize(replicaGraphMPIControl().numberOfReplicas);

        DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":init_eds_stat:\t DONE");
}

//RE

/**
 * Override Exchange Functions
 */
// top layer 
void re::replica_exchange_base_2d_s_eoff_eds::determine_switch_probabilities(){
    DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":determine_switch_probabilities:\t DONE");
    replica->sim.param().reeds.eds_para[partnerReplicaID].pos_info.second = simulationID;

    if(run % 2 == 1){
        swap_s(partnerReplicaID);
    }
    else{
        swap_eoff(partnerReplicaID);
    }
    
    DEBUG(3,"replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":determine_switch_probabilities:\t DONE");
}

int re::replica_exchange_base_2d_s_eoff_eds::find_partner() const {
  unsigned int num_eoff = replica->sim.param().reeds.num_eoff;
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds:find_partner: num_eoff= " << num_eoff << "\n");
  unsigned int num_l = replica->sim.param().reeds.num_s;
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds:find_partner: num_l= " << num_l << "\n");
  unsigned int numT = replica->sim.param().replica.num_T;
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds:find_partner: numT= " << numT << "\n");
  unsigned int numReps = num_l * num_eoff;
  DEBUG(3,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds:find_partner: numReps= " << numReps << "\n");
  unsigned int ID = simulationID;

  unsigned int partner = ID;
  bool evenRow = (ID % num_l) % 2 == 0;//1st row is here the 0th row and therefore even!
  bool evenCol = (ID / num_l) % 2 == 0;//1st col is here the 0th col and therefore even!
  bool numEoffeven = num_eoff % 2 == 0;//used for periodic boundary


  //edge cases for s dimension
  bool upper = ID % num_l == 0;
  bool lower = ID % num_l == num_l - 1;

  //current s coord == j € [0, num_l -1)
  unsigned int j = ID % num_l;

  //edge cases for eoff dimension
  bool left_edge = ID == j;
  bool right_edge = ID == (numReps - num_l + j);
  DEBUG(3,"ID, j, upper, lower, left_edge, right_edge= " << ID << ", " << j << ", " << upper << ", " << lower
  << ", " << left_edge << ", " << right_edge << "\n");


    // determine switch direction
    switch ((run % 4) - 1) {
      case 0: //s dimension
      DEBUG(5,"find_partner: FIRST case\n");
          if (evenRow) {
            partner = ID + 1;
            //edge case
            if(lower)
              partner = ID;
          }
          else {
            partner = ID - 1;
            //edge case
            if(upper)
              partner = ID;
          }
      DEBUG(1,"find_partner(first case): partner of ID=" << ID << " is " << partner << "\n");
        break;

      case 1: //eoff dimension
      DEBUG(5,"find_partner: SECOND case\n");
        if(numEoffeven){
          DEBUG(5,"find_partner(second case): numEoffeven_firstCase \n");
          partner = partner_eoffDim_numEoffeven_firstCase();
        }

        else{
          DEBUG(5,"find_partner(second case): numEoffodd_firstCase \n");
          partner = partner_eoffDim_numEoffodd_firstCase();
        }
      DEBUG(1,"find_partner(second case): partner of ID=" << ID << " is " << partner << "\n");
        break;

      case 2: //s dimension
      DEBUG(5,"find_partner: THIRD case\n");
        if (evenRow) {
          partner = ID - 1;
          //edge case
          if(upper)
            partner = ID;
        }
        else {
          partner = ID + 1;
          //edge case
          if(lower)
            partner = ID;
        }
      DEBUG(1,"find_partner(third case): partner of ID=" << ID << " is " << partner << "\n");
        break;

      case -1: //eoff dimension
      DEBUG(5,"find_partner: FOURTH case\n");
        if(numEoffeven){
          DEBUG(5,"find_partner(fourth case): numEoffeven_secondCase \n");
          partner = partner_eoffDim_numEoffeven_secondCase();
        }

        else{
          DEBUG(5,"find_partner(fourth case): numEoffodd_secondCase \n");
          partner = partner_eoffDim_numEoffodd_secondCase();
        }
      DEBUG(1,"find_partner(fourth case): partner of ID=" << ID << " is " << partner << "\n");
        break;
    }

  return partner;
}

//find partner functions
int re::replica_exchange_base_2d_s_eoff_eds::partner_eoffDim_numEoffeven_firstCase() const {
  unsigned int ID = simulationID;
  unsigned int partner = ID;

  unsigned int num_eoff = replica->sim.param().reeds.num_eoff;
  unsigned int num_l = replica->sim.param().reeds.num_s;
  unsigned int numReps = num_l * num_eoff;
  bool periodic = replica->sim.param().reeds.periodic;

  bool evenCol = (ID / num_l) % 2 == 0;//1st col is here the 0th col and therefore even!
  //current s coord == j € [0, num_l -1)
  unsigned int j = ID % num_l;

  //edge cases for eoff dimension
  bool left_edge = ID == j;
  bool right_edge = ID == (numReps - num_l + j);

  if (evenCol) {
    partner = ID + num_l;
    //edge case
    if(right_edge && !periodic) partner = ID;
    if(right_edge && periodic) {partner = (ID + num_l) % numReps; DEBUG(5,"\nPERIODIC\n");}
  }
  else {
    partner = ID - num_l;
    //edge case
    if(left_edge && !periodic) partner = ID;
    if(left_edge && periodic) {partner = ID + (numReps - num_l); DEBUG(5,"\nPERIODIC\n");}
  }
  return partner;
}

int re::replica_exchange_base_2d_s_eoff_eds::partner_eoffDim_numEoffodd_firstCase() const {
  unsigned int ID = simulationID;
  unsigned int partner = ID;

  unsigned int num_eoff = replica->sim.param().reeds.num_eoff;
  unsigned int num_l = replica->sim.param().reeds.num_s;
  unsigned int numReps = num_l * num_eoff;
  bool periodic = replica->sim.param().reeds.periodic;

  bool evenCol = (ID / num_l) % 2 == 0;//1st col is here the 0th col and therefore even!
  //current s coord == j € [0, num_l -1)
  unsigned int j = ID % num_l;

  //edge cases for eoff dimension
  bool left_edge = ID == j;
  bool right_edge = ID == (numReps - num_l + j);

  if(periodic){
    partner = partner_eoffDim_numEoffodd_cyclic();
  }

  else{
    if (evenCol) {
      partner = ID + num_l;
      //edge case
      if(right_edge) partner = ID;
    }
    else {
      partner = ID - num_l;
      //edge case
      if(left_edge) partner = ID;
    }
  }
  return partner;
}

int re::replica_exchange_base_2d_s_eoff_eds::partner_eoffDim_numEoffeven_secondCase() const {
  unsigned int ID = simulationID;
  unsigned int partner = ID;

  unsigned int num_eoff = replica->sim.param().reeds.num_eoff;
  unsigned int num_l = replica->sim.param().reeds.num_s;
  unsigned int numReps = num_l * num_eoff;
  bool periodic = replica->sim.param().reeds.periodic;

  bool evenCol = (ID / num_l) % 2 == 0;//1st col is here the 0th col and therefore even!
  //current s coord == j € [0, num_l -1)
  unsigned int j = ID % num_l;

  //edge cases for eoff dimension
  bool left_edge = ID == j;
  bool right_edge = ID == (numReps - num_l + j);

  if (evenCol) {
    partner = ID - num_l;
    //edge case
    if(left_edge && !periodic) partner = ID;
    if(left_edge && periodic) {partner = ID + (numReps - num_l); DEBUG(5,"replica "<<globalThreadID<<":PERIODIC\n");}
  }
  else {
    partner = ID + num_l;
    //edge case
    if(right_edge && !periodic) partner = ID;
    if(right_edge && periodic) {partner = (ID + num_l) % numReps; DEBUG(5,"\nPERIODIC\n");}
  }
  return partner;
}

int re::replica_exchange_base_2d_s_eoff_eds::partner_eoffDim_numEoffodd_secondCase() const {
  unsigned int ID = simulationID;
  unsigned int partner = ID;

  unsigned int num_eoff = replica->sim.param().reeds.num_eoff;
  unsigned int num_l = replica->sim.param().reeds.num_s;
  unsigned int numReps = num_l * num_eoff;
  bool periodic = replica->sim.param().reeds.periodic;

  bool evenCol = (ID / num_l) % 2 == 0;//1st col is here the 0th col and therefore even!
  //current s coord == j € [0, num_l -1)
  unsigned int j = ID % num_l;

  //edge cases for eoff dimension
  bool left_edge = ID == j;
  bool right_edge = ID == (numReps - num_l + j);

  if(periodic){
    partner = partner_eoffDim_numEoffodd_cyclic();
  }

  else{
    if (evenCol) {
      partner = ID - num_l;
      //edge case
      if(left_edge) partner = ID;
    }
    else {
      partner = ID + num_l;
      //edge case
      if(right_edge) partner = ID;
    }
  }
  return partner;
}

int re::replica_exchange_base_2d_s_eoff_eds::partner_eoffDim_numEoffodd_cyclic() const {
  //very important to keep every previous unsigned int as int bc partner might be negative during the computation
  int ID = simulationID;
  int partner = ID;

  int num_eoff = replica->sim.param().reeds.num_eoff;
  int num_l = replica->sim.param().reeds.num_s;
  int numReps = num_l * num_eoff;

  bool evenCol = (ID / num_l) % 2 == 0;//1st col is here the 0th col and therefore even!
  int blocked_col = (run/2 - 1) % num_eoff;//to identify the blocked column
  bool even_blocked_col = blocked_col % 2 == 0;//check whether blocked column is an even column
  int border = blocked_col * num_l + (num_l - 1);//to identify if ID is in column left from blocked_col or right from it
  bool exchanged = false;//for efficiency reasons not to compute irrelevant conditions for IDs in blocked column

  DEBUG(5,"\nblocked_col, border: " << blocked_col << ", " << border << "\n");

  //blocked column does not exchange
  for(int i=0; i<num_l; ++i){
	if(ID == (blocked_col * num_l + i)){
		partner = ID;
		exchanged = true;
		break;
	}
  }

  if(!exchanged){
	if(ID<border){
		if(evenCol && even_blocked_col){
			partner = (ID + num_l) % numReps;
			DEBUG(7,"\ncase a with ID: " << ID << "\n");
		}
		if(!evenCol && even_blocked_col) {
			partner = (numReps + ((ID - num_l) % numReps)) % numReps;//to handle possible negative partner right
			DEBUG(7,"\ncase b with ID: " << ID << "\n");
		}
		if(evenCol && !even_blocked_col){
			partner = (numReps + ((ID - num_l) % numReps)) % numReps;//to handle possible negative partner right
			DEBUG(7,"\ncase e with ID: " << ID << "\n");
		}
		if(!evenCol && !even_blocked_col){
			partner = (ID + num_l) % numReps;
			DEBUG(7,"\ncase f with ID: " << ID << "\n");
		}
	}
	if(ID>border){
		if(evenCol && even_blocked_col){
			partner = (numReps + ((ID - num_l) % numReps)) % numReps;//to handle possible negative partner right
			DEBUG(7,"\ncase c with ID: " << ID << "\n");
		}
		if(!evenCol && even_blocked_col){
			partner = (ID + num_l) % numReps;
			DEBUG(7,"\ncase d with ID: " << ID << "\n");
		}
		if(evenCol && !even_blocked_col){
			partner = (ID + num_l) % numReps;
			DEBUG(7,"\ncase g with ID: " << ID << "\n");
		}
		if(!evenCol && !even_blocked_col){
			partner = (numReps + ((ID - num_l) % numReps)) % numReps;//to handle possible negative partner right
			DEBUG(7,"\ncase h with ID: " << ID << "\n");
		}
	}
  }

  return partner;
}


//determine swaps
void re::replica_exchange_base_2d_s_eoff_eds::swap_s(const unsigned int partnerReplicaID) {
  DEBUG(4, "replica "<<  globalThreadID <<":swap:\t  START");

  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;
  unsigned int numReps = replica->sim.param().reeds.num_s * replica->sim.param().reeds.num_eoff;

  // does partner exist?
  if (partnerReplicaID < numReps && partnerReplicaID != simulationID) {
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaID) {

      // posts a MPI_Recv(...) matching the MPI_Send below
      probability = calc_probability(partnerReplicaID);
      const double randNum = rng.get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;

#ifdef XXMPI
      DEBUG(4,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds: SWAP_S before Send\n");
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS, replicaGraphMPIControl().comm);
      DEBUG(4,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds: SWAP_S after Send\n");
#endif

      if (randNum < probability) {
        switched = true;
      } else
        switched = false;
    } else {    //The Partner sends his data to The calculating Thread
      //special case if lambda also needs to be exchanged
      bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaID]);
      DEBUG(3, "swap_s: simID, s value= " << simulationID << ", " << l << "\n");
      DEBUG(3,"swap_s: simID, bool sameLambda= " << simulationID << ", " << sameLambda << "\n");
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
        DEBUG(4,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds: SWAP_S before Send\n");
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
        DEBUG(4,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds: SWAP_S after Send\n");
#endif
      } else { // sameLambda
#ifdef XXMPI
        double energies[2] = {epot, 0.0};
        DEBUG(4,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds: SWAP_S before Send\n");
        MPI_Send(&energies[0],2,MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
        DEBUG(4,"replica "<<globalThreadID<<":replica_exchange_base_2d_s_eoff_eds: SWAP_S after Send\n");
#endif
     }
      if (replica->sim.param().pcouple.scale != math::pcouple_off) {
#ifdef XXMPI
        math::Box box_replica = replica->conf.current().box;    //exchange box
        MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerReplicaMasterThreadID, BOX,  replicaGraphMPIControl().comm);
#endif
      }

#ifdef XXMPI
      MPI_Status status;
#endif
      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
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
      throw "Partner does not exist!";
    /*
      partner = ID;
    switched = false;
    probability = 0.0;

    */
  }
    DEBUG(4, "replica "<< globalThreadID <<":replica_exchange_base_2d_s_eoff_eds:swap:\t  DONE");
}

void re::replica_exchange_base_2d_s_eoff_eds::swap_eoff(const unsigned int partnerReplicaID) {
  DEBUG(4, "replica "<<  globalThreadID <<":swap_eoff:\t  START");


  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;
  unsigned int numReps = replica->sim.param().reeds.num_s * replica->sim.param().reeds.num_eoff;

  // does partner exist?
  if (partnerReplicaID < numReps && partnerReplicaID != simulationID) {
    // the one with lower ID does probability calculation
    if (simulationID < partnerReplicaID) {

      // posts a MPI_Recv(...) matching the MPI_Send below
      probability = calc_probability_for_eoff_exchange(partnerReplicaID);
      DEBUG(5,"SWAP_EOFF: ID, probability = " << simulationID << ", " << probability << "\n");
      const double randNum = rng.get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;

#ifdef XXMPI
      DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: SWAP_EOFF before Send\n");
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SENDCOORDS, replicaGraphMPIControl().comm);
      DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: SWAP_EOFF after Send\n");
#endif

      if (randNum < probability) {
        switched = true;
      } else
        switched = false;
      DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: ID, switched = " << simulationID << ", " << switched << "\n");
    } else {    //The Partner sends his data to The calculating Thread
      //special case if eoff vector also needs to be exchanged
      bool sameEoffvector = true;
      for(int i=0; i<replica->sim.param().reeds.num_states; ++i){
        if(replica->sim.param().reeds.eds_para[simulationID].eir[i] != replica->sim.param().reeds.eds_para[partnerReplicaID].eir[i]){
          sameEoffvector = false;
          DEBUG(1,"\nSWAP_EOFF: ID " << simulationID << " sets sameEoffvector to false because of " << i << "th position of Eoffvector\n");
          break;
        }
      }
      DEBUG(5,"swap_eoff: simID, bool sameEoffvector= " << simulationID << ", " << sameEoffvector << "\n");
      if(!sameEoffvector){      //exchange EOFF_VEC
        // E21: Energy with configuration 2 and lambda 1(of partner)
        const double E21 = calculate_energy(partnerReplicaMasterThreadID);
        // this we can store as the partner energy of the current replica
        epot_partner = E21;
        // E22: Energy with configuration 2 and lambda 2(own one)
#ifdef XXMPI
        const double E22 = epot;
        // send E21 and E22
        double energies[2] = {E22, E21};
        //this send operation is matched in calc_probability_for_eoff_exchange()
        DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: SWAP_EOFF before Send\n");
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm);
        DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: SWAP_EOFF after Send\n");
#endif
      } else { // sameEoffvector
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

#ifdef XXMPI
      MPI_Status status;
#endif
      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
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
      DEBUG(5,"\nreplica_exchange_base_2d_s_eoff_eds: 2nd time ID, switched = " << simulationID << ", " << switched << "\n");
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

double re::replica_exchange_base_2d_s_eoff_eds::calc_probability_for_eoff_exchange(const unsigned int partnerReplicaID) {
  DEBUG(4,"replica_exchange_base_2d_s_eoff_eds: CALC_PROBABILITY by ID: " << simulationID << "\n");

  unsigned int partnerReplicaMasterThreadID = partnerReplicaID;

  double delta = 0.0;
  const double b1 = 1.0 / (math::k_Boltzmann * T);
  const double b2 = 1.0 / (math::k_Boltzmann * replica->sim.param().replica.temperature[partnerReplicaID % replica->sim.param().reeds.num_eoff]);

  bool sameEoffvector = true;
  for(int i=0; i<replica->sim.param().reeds.num_states; ++i){
    if(replica->sim.param().reeds.eds_para[simulationID].eir[i] != replica->sim.param().reeds.eds_para[partnerReplicaID].eir[i]){
      sameEoffvector = false;
      DEBUG(1,"\nSWAP_EOFF in calc_prob(): ID " << simulationID << " sets sameEoffvector to false because of " << i << "th position of Eoffvector\n");
      break;
    }
  }
  DEBUG(5,"swap_eoff in calc_prob(): simID, bool sameEoffvector= " << simulationID << ", " << sameEoffvector << "\n");

  if (sameEoffvector) {
    // use simple formula
    // get energy from the other partner
    double energies[2] = {0.0, 0.0};
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm, &status);
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
    DEBUG(1,"\n\nreplica_exchange_base_2d_s_eoff_eds: CALC_PROBABILITY before Recv\n");
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm, &status);
    DEBUG(1,"\n\nreplica_exchange_base_2d_s_eoff_eds: CALC_PROBABILITY after Recv\n");
#endif
    const double E22 = energies[0];
    const double E12 = energies[1];

    const double E11 = epot;
    const double E21 = calculate_energy(partnerReplicaID);

    // store this as the partner energy
    epot_partner = E21;

    delta = b1 * (E12 - E11) - b2 * (E22 - E21);
  }

  DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: CALC_PROBABILITY before NPT\n");
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

  DEBUG(4,"replica_exchange_base_2d_s_eoff_eds: CALC_PROBABILITY by ID: " << simulationID << " done\n");
}


////exchange params
void re::replica_exchange_base_2d_s_eoff_eds::reset_eds() {//only reset switched parameters of change_eds() function
  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds: RESET_EDS\n\n");
  replica->sim.param().eds = eds_para;
  replica->sim.param().step.dt = dt;
  replica->conf.current().force= force_orig;
  replica->conf.current().virial_tensor= virial_tensor_orig;
}

void re::replica_exchange_base_2d_s_eoff_eds::change_eds(const unsigned int partner){//only change parameters, which are needed for energy calculation i.e.

  DEBUG(3,"replica_exchange_base_2d_s_eoff_eds: CHANGE_EDS\n\n");
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
double re::replica_exchange_base_2d_s_eoff_eds::calc_energy_eds_stat(double s){
    double old_dt = 0.0;
    double old_s = 0.0;
    double old_eds_vr = 0.0;
    algorithm::Algorithm * ff = nullptr;
    DEBUG(5,"replica_exchange_base_2d_s_eoff_eds: CALC_ENERGY_EDS_STAT\n\n");
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


double re::replica_exchange_base_2d_s_eoff_eds::calculate_energy(const unsigned int selectedReplicaID) {
    DEBUG(4, "replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":calculate_energy:\t START");

    DEBUG(5, "replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":calculate_energy:\t get Partner settings");
    if(selectedReplicaID!=simulationID){
        change_eds(selectedReplicaID);
    }

    double energy = calculate_energy_core();

    if(selectedReplicaID!=simulationID){
        reset_eds();
    }
    DEBUG(4, "replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":calculate_energy:\t DONE");
    return energy;
}


double re::replica_exchange_base_2d_s_eoff_eds::calculate_energy_core() {

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


void re::replica_exchange_base_2d_s_eoff_eds::calculate_energy_helper(const unsigned int selectedReplicaID) {
    DEBUG(4, "replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":calculate_energy:\t START");

    DEBUG(4, "replica_exchange_base_2d_s_eoff_eds "<< globalThreadID <<":calculate_energy:\t DONE");
}
