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


#include <util/replicaExchange/replica_exchangers/replica_exchange_base_eds.h>
#include <util/replicaExchange/replica/replica_reeds.h>

//Constructor
#include <util/replicaExchange/replica_exchangers/replica_exchange_base.h>
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

util::replica_exchange_base_eds::replica_exchange_base_eds(io::Argument _args, int cont, int rank, int simulationRank, int simulationID, int simulationThreads,
        std::vector<int> repIDs, std::map<ID_t, rank_t> &_repMap):  
        replica_exchange_base(_args, cont, rank, simulationRank, simulationID, simulationThreads, repIDs, _repMap){
    DEBUG(3,"replica_exchange_base_eds "<< rank <<":Constructor:\t START");
    DEBUG(4,"replica_exchange_base_eds "<< rank <<":Constructor:\t replica Type\t "<< typeid(replicas).name());
    createReplicas(cont, repIDs, rank);
    DEBUG(3,"replica_exchange_base_eds "<< rank <<":Constructor:\t DONE");
}

void util::replica_exchange_base_eds::createReplicas(int cont, std::vector<int>  repIDs, int rank){
  MPI_DEBUG(3,"replica_exchange_base_eds "<< rank <<":createReplicas:\t START");
  
  if(!replica){
    delete replica;
  }
  MPI_DEBUG(3,"replica_exchange_base_eds "<< rank <<":createReplicas:\t Deleted");

  // create the number of replicas that are assigned to my node
    if(simulationThreads>1){
        if(simulationRank == 0){
            replica = new util::replica_MPI_Master(args, cont, simulationID, rank, simulationRank, simulationID, simulationThreads);
        }
        else{
            replica = new util::replica_MPI_Slave(args, cont, simulationID, rank, simulationRank, simulationID, simulationThreads);
        }
    }
    else{
        replica = new util::replica_reeds(args, cont, simulationID, rank);
    }
  MPI_DEBUG(3,"replica_exchange_base_eds "<< rank <<":createReplicas:\t DONE");
    /*
  replicas.resize(numReplicas);
  // create the number of replicas that are assigned to my node
   int i = 0;
   for (repIterator it = replicas.begin(); it < replicas.end(); ++it, ++i) {
    *it = new util::replica_reeds(args, cont, repIDs[i], rank);
    DEBUG(4, "replica_exchange_base_eds "<< rank <<":initReplicas:\tConstructor \ttopo\t" << (*it)->topo.check_state())
    DEBUG(4, "replica_exchange_base_eds "<< rank <<":initReplicas:\tConstructor \tconf:\t" << (*it)->conf.check((*it)->topo, (*it)->sim))
   }  
   */
}


util::replica_exchange_base_eds::~replica_exchange_base_eds() {
    delete replica;
  /*
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    delete *it;
  }
   * */
}

void util::replica_exchange_base_eds::swap() {
  DEBUG(3,"replica_exchange_base_eds "<< rank <<":swap:\t Start");
  
      unsigned int partner = replica->find_partner();
      unsigned int partnerRank = repMap.find(partner)->second;
      if (partner != (replica)->ID) {
      // are the replicas on the same node?
      DEBUG(4,"replica_exchange_base_eds "<< rank <<":swap:\t sending "<< replica->ID <<" \t Start");
      if (partnerRank != replica->rank) // replicas on different nodes
      {
        replica->swap(partner, partnerRank);
        if (replica->switched) {
          if (replica->ID < partner) {
            replica->send_coord(partner, partnerRank);
            replica->receive_new_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back bschroed:that is new!
            //(*it)->exchange_averages();
          } else {
            replica->receive_new_coord(partner, partnerRank);
            replica->send_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back bschroed:that is new!
            //(*it)->exchange_averages();
          }
        }
      }
    } else {
      replica->partner = replica->ID;
      replica->probability = 0.0;
      replica->switched = 0;
    }
  
  
  
  /*
  // do this for all replicas; if two of them on this node, do special swap
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    unsigned int partner = (*it)->find_partner();
    unsigned int partnerRank = repMap.find(partner)->second;
    bool partnerOnThisNode = (partnerRank == (*it)->rank && (*it)->ID < partner);

    // attempt switch in this run?
    if (partner != (*it)->ID) {
      // are the replicas on the same node?
      DEBUG(4,"replica_exchange_base_eds "<< rank <<":swap:\t sending "<< (*it)->ID <<" \t Start");
      if (partnerOnThisNode) {
        swap_on_node(it, partner);
        if ((*it)->switched)
          switch_coords_on_node(it, partner);
      } else if (partnerRank != (*it)->rank) // replicas on different nodes
      {
        (*it)->swap(partner, partnerRank);
        if ((*it)->switched) {
          if ((*it)->ID < partner) {
            (*it)->send_coord(partner, partnerRank);
            (*it)->receive_new_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back bschroed:that is new!
            //(*it)->exchange_averages();
          } else {
            (*it)->receive_new_coord(partner, partnerRank);
            (*it)->send_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back bschroed:that is new!
            //(*it)->exchange_averages();
          }
        }
      }
    } else {
      (*it)->partner = (*it)->ID;
      (*it)->probability = 0.0;
      (*it)->switched = 0;
    }
  }
   * */
    DEBUG(3,"replica_exchange_base_eds "<< rank <<":swap:\t Done");
}

void util::replica_exchange_base_eds::write_final_conf() {
  // write coordinates to cnf for all replica assigned to this node
  for (std::vector< util::replica_reeds * >::iterator it = replicas.begin(); it < replicas.end(); ++it) {
    (*it)->write_final_conf();
  }
}



void util::replica_exchange_base_eds::run_MD() {
    DEBUG(3, "replica_exchange_base_eds "<< rank <<":run_MD:\t START")

     // do a md run for all replica assigned to this node
    for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
        DEBUG(4, "replica_exchange_base_eds "<< rank <<":run_MD:\t replica "<< (*it)->ID << ":")
        DEBUG(4, "replica_exchange_base_eds "<< rank <<":run_MD:\t replica\t Start")
        (*it)->run_MD();
        DEBUG(4, "replica_exchange_base_eds "<< rank <<":run_MD:\t replica\t Done")
    }
    DEBUG(3, "replica_exchange_base_eds "<< rank <<":run_MD:\t DONE")
}


/*SID: eds_stat is used to calculate and output additional acceptance
 probability information in order to develop an optimization procedure
 for the s distribution of replica exchange eds simulation.
 * This function increases the time for a replica exchange trial and should 
 * not be called if the output is not needed.
 */
void util::replica_exchange_base_eds::eds_stat(){
    //Calculate exchange probabilities for each replica
        /*
         * loop over all replica on current node
         */
        for (util::replica_reeds* x : replicas){
            ID_t currentID = x->ID;
            for(unsigned int rep=0; rep< this->numReplicas;rep++){
               double rep_s=x->sim.param().reeds.lambda[rep];
               replicaStatData[currentID].epot_vec[rep]= x->calc_energy_eds_stat(rep_s);
            }
            replicaStatData[currentID].run=replicaStatData[currentID].run +1;

          //Write data to output file
          *(eds_stat_out[currentID]) << std::setw(6) << (replicaStatData[currentID].ID + 1)
                  << " "
                  << std::setw(1) << replicaStatData[currentID].run
                  << std::setw(1) << replicaStatData[currentID].s
                  << std::setw(1) << replicaStatData[currentID].T
                  << " ";
            for(unsigned int rep=0; rep<this->numReplicas; rep++){
                    *(eds_stat_out[currentID]) << std::setw(10) << replicaStatData[currentID].epot_vec[rep]
                    << " ";
            }
          *(eds_stat_out[currentID])   << std::endl;
            
        }
    }

void util::replica_exchange_base_eds::init() {
  DEBUG(3,"replica_exchange_base_eds "<< rank <<":init:\t init \t START");
  MPI_DEBUG(3,"replica_exchange_base_eds "<< rank <<":init:\t start init from baseclass \t NEXT");
  replica->init();
  /*
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    (*it)->init();
   }
   * */
  //replica_exchange_base::init();
  MPI_DEBUG(3,"replica_exchange_base_eds "<< rank <<":init:\t init_eds_stat \t NEXT");
  this->init_eds_stat();
  DEBUG(3,"replica_exchange_base_eds "<< rank <<":init:\t DONE");
}

//initialize output files  
void util::replica_exchange_base_eds::init_eds_stat(){
        DEBUG(3,"replica_exchange_base_eds "<< rank <<":init_eds_stat:\t START");
        
        ID_t currentID=1000; //error value
        currentID = replica->ID;
        replicaStatData[currentID].ID =currentID;
        replicaStatData[currentID].T=replica->T;
        replicaStatData[currentID].s=replica->l; //l==s because of the implementation of hamiltonian replica exchange.
        replicaStatData[currentID].dt=replica->dt;
        replicaStatData[currentID].run=0;
        replicaStatData[currentID].epot_vec.resize(this->numReplicas);
        replicaStatData[currentID].prob_vec.resize(this->numReplicas);
        
        /*
        for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
          currentID = (*it)->ID;
          replicaStatData[currentID].ID =currentID;
          replicaStatData[currentID].T=(*it)->T;
          replicaStatData[currentID].s=(*it)->l; //l==s because of the implementation of hamiltonian replica exchange.
          replicaStatData[currentID].dt=(*it)->dt;
          replicaStatData[currentID].run=0;
          replicaStatData[currentID].epot_vec.resize(this->numReplicas);
          replicaStatData[currentID].prob_vec.resize(this->numReplicas);
        }
         * */
        DEBUG(3,"replica_exchange_base_eds "<< rank <<":init_eds_stat:\t DONE");
}
