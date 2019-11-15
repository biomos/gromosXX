/* 
 * File:   replica_exchange_base.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */

#include "replicaExchange/replica/replica.h"
#include "replicaExchange/replica/replica_MPI_master.h"
#include "replicaExchange/replica/replica_MPI_slave.h"
#include "replicaExchange/replica_exchangers/replica_exchange_base.h"

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

util::replica_exchange_base::replica_exchange_base(io::Argument _args, int cont, int rank, int simulationRank, int simulationID, int simulationThreads,
        std::vector<int> repIDs, std::map<ID_t, rank_t> &_repMap) : args(_args), repMap(_repMap), numReplicas(repIDs.size()),
        rank(rank), simulationRank(simulationRank), simulationID(simulationID), simulationThreads(simulationThreads), cont(cont), repIDs(repIDs), rng(-1){
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base "<< rank <<":Constructor:\t START ");
  DEBUG(3,"replica_exchange_base "<< rank <<":Constructor:\t SIMULATIONID:  "<< simulationID);

  //construct replica obj
  DEBUG(4,"replica_exchange_base "<< rank <<":Constructor:\t  createReplica");
  createReplicas(cont, repIDs, rank);
  DEBUG(4,"replica_exchange_base "<< rank <<":Constructor:\t createdReplica T");

  DEBUG(3,"replica_exchange_base "<< rank <<":Constructor:\t Constructor \t DONE");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!"; 
  #endif
}

util::replica_exchange_base::~replica_exchange_base() {
    delete replica;
    
  //TODO: Think about if this implementation is practical at any point? or needs to be removed. bschroed
    /*
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    delete *it;
  }
   */
}

void util::replica_exchange_base::run_MD() {
  DEBUG(3,"replica_exchange_base "<< rank <<":run_MD:\t START");
    replica->run_MD();
    replica->conf.current().energies.eds_vr;
    //TODO: Think about if this implementation is practical at any point? or needs to be removed. bschroed
   /*
  // do a md run for all replica assigned to this node
  for (std::vector< util::replica * >::iterator it = replicas.begin(); it < replicas.end(); ++it) {
    (*it)->run_MD();
    (*it)->conf.current().energies.eds_vr;
  }
    * */
  DEBUG(3,"replica_exchange_base "<< rank <<":run_MD:\t END");

}

void util::replica_exchange_base::write_final_conf() {
  // write coordinates to cnf for all replica assigned to this node
   replica->write_final_conf();
  //TODO: Think about if this implementation is practical at any point? or needs to be removed. bschroed
  //for (std::vector< util::replica * >::iterator it = replicas.begin(); it < replicas.end(); ++it) {
  //  (*it)->write_final_conf();
  //}
}

void util::replica_exchange_base::init() {
  DEBUG(3, "replica_exchange_base "<< rank <<":init:\t START");

  // do init for all replica assigned to this node
  DEBUG(4,"replica_exchange_base "<< rank <<":init:\t initReplicas");
  replica->init();

    //TODO: Think about if this implementation is practical at any point? or needs to be removed. bschroed
  //for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
  //  (*it)->init();
  //}
  DEBUG(3,"replica_exchange_base "<< rank <<":init:\t DONE");

}

void util::replica_exchange_base::createReplicas(int cont, std::vector<int> repIDs, int rank){
  DEBUG(3,"replica_exchange_base "<< rank <<":createReplicas:\t START");
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
        replica = new util::replica(args, cont, simulationID, rank);
        //DEBUG(2,"replica_exchange_base "<< rank <<":create:\t rep" << replica->sim.param().replica.num_l);
        //DEBUG(2,"replica_exchange_base "<< rank <<":create:\t rep" << replica->sim.param().replica.lambda[0]);
    }
  
  //TODO: Think about if this implementation is practical at any point? or needs to be removed. bschroed
  /*
   *   util::replica* replica;

  int i = 0;
  replicas.resize(numReplicas);
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it, ++i) {
        
      //select correct replica implementation - Polymorphism.
    if(_subSimulationThreadedMPI>1){
        replica = new util::replica_MPI(args, cont, repIDs[i++], rank, _subSimulationThreadedMPI);
    }
    else{
        replica = new util::replica(args, cont, repIDs[i++], rank);
    }
    // one could use a vector<util::replica> and theref ore we needed a copy constructor
    // for the replica class but there seems to be something wrong with those of conf/topo/...
    *it = replica;

  
    DEBUG(4, "replica_exchange_base "<< rank <<":createReplicas:\t topo\t" << (*it)->topo.check_state())
    DEBUG(4, "replica_exchange_base "<< rank <<":createReplicas:\t conf:\t" << (*it)->conf.check((*it)->topo, (*it)->sim))
  }
  */
  DEBUG(3,"replica_exchange_base "<< rank <<":createReplicas:\t DONE");
}

void util::replica_exchange_base::swap(){
  DEBUG(3,"replica_exchange_base "<< rank <<":swap:\t START");
  
    unsigned int partner = replica->find_partner();
    unsigned int partnerRank = repMap.find(partner)->second;
    
    if (partnerRank != replica->rank) // replicas on different nodes
      {
        replica->swap(partner, partnerRank);
        if (replica->switched) {
          if (replica->ID < partner) {
            replica->send_coord(partner, partnerRank);
            replica->receive_new_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
            replica->exchange_averages();
          } else {
            replica->receive_new_coord(partner, partnerRank);
            replica->send_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
            replica->exchange_averages();
          }
        }
      }
    else {
      replica->partner = replica->ID;
      replica->probability = 0.0;
      replica->switched = 0;
    }
    if(replica->switched && replica->sim.param().replica.scale) {
      replica->velscale(replica->partner);
    }
     
  DEBUG(3,"replica_exchange_base "<< rank <<":swap:\t DONE");
    //TODO: Think about if this implementation is practical at any point? or needs to be removed. bschroed
  // do this for all replicas; if two of them on this node, do special swap
    /*
     * 
   

 // for (std::vector< util::replica* >::iterator it = replicas.begin(); it < replicas.end(); ++it) {
    //unsigned int partner = (*it)->find_partner();
    unsigned int partner = (*replica)->find_partner();
    unsigned int partnerRank = repMap.find(partner)->second;
    bool partnerOnThisNode = (partnerRank == (*it)->rank && (*it)->ID < partner);

    // attempt switch in this run?
    \/\* 
    if (partner != (*it)->ID) {
      // are the replicas on the same node?
  
    if (partnerOnThisNode) {
        swap_on_node(it, partner);
        if ((*it)->switched)
          switch_coords_on_node(it, partner);
      } else if (partnerRank != (*it)->rank) // replicas on different nodes
      {
        
      if (partnerRank != (*it)->rank) // replicas on different nodes
      {
        (*it)->swap(partner, partnerRank);
        if ((*it)->switched) {
          if ((*it)->ID < partner) {
            (*it)->send_coord(partner, partnerRank);
            (*it)->receive_new_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
            (*it)->exchange_averages();
          } else {
            (*it)->receive_new_coord(partner, partnerRank);
            (*it)->send_coord(partner, partnerRank);
            // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
            (*it)->exchange_averages();
          }
        }
      }
    } else {
      (*it)->partner = (*it)->ID;
      (*it)->probability = 0.0;
      (*it)->switched = 0;
    }

  }

  // scale the velocities?
  for (std::vector< util::replica* >::iterator it = replicas.begin();it < replicas.end(); ++it) {
    if((*it)->switched && (*it)->sim.param().replica.scale) {
      (*it)->velscale((*it)->partner);
    }
  }
  DEBUG(3,"replica_exchange_base "<< rank <<":swap:\t DONE");
}
*/
}
/*
void util::replica_exchange_base::swap_on_node(repIterator it1, const unsigned int partner) {
  replica* rep1 = *it1;
  replica* rep2;

  // find partner of rep1 on my node
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    if ((*it)->ID == partner) {
      rep2 = (*it);
      break;
    }
  }

  rep1->probability = calc_probability(rep1, rep2);
  rep2->probability = rep1->probability;

  const double randNum = rng.get();

  if (randNum < rep1->probability) {
    rep1->switched = true;
    rep2->switched = true;
  } else {
    rep1->switched = false;
    rep2->switched = false;
  }
}
*/

double util::replica_exchange_base::calc_probability(util::replica_Interface * rep1, util::replica_Interface * rep2) {
  double delta;
  bool sameLambda = (rep2->l == rep1->l); // horizontal switch with same lambda?
  const double b1 = 1.0 / (math::k_Boltzmann * rep1->T);
  const double b2 = 1.0 / (math::k_Boltzmann * rep2->T);

  if (sameLambda) {
    // use simple formula
    rep1->epot_partner = rep2->epot;
    rep2->epot_partner = rep1->epot;
    delta = (b1 - b2)*(rep2->epot - rep1->epot); //*  (E21 - E11=
  } else {
    // 2D formula
    /*
     * E12: Energy with lambda from 1 and configuration from 2
     * delta = b1 * ( E22 - E11) - b2*(E21  - E12);
     * E22 and E12 needed from partner
     */
    const double E22 = rep2->epot;
    const double E12 = rep2->calculate_energy(rep1->ID);

    const double E11 = rep1->epot;
    const double E21 = rep1->calculate_energy(rep2->ID);
    // Chris: I think this was wrong
    // delta = b1 * (E22 - E11) - b2 * (E21 - E12);
    //std::cerr << "b1: " << b1 << " b2: " << b2 << std::endl;
    //std::cerr << "E11: " << E11 << " E22: " << E22 << std::endl;
    //std::cerr << "E21: " << E21 << " E12: " << E12 << std::endl;

    delta = b1 * (E12 - E11) - b2 * (E22 - E21);
  }
  
  // NPT? add PV term
  if (rep1->sim.param().pcouple.scale != math::pcouple_off) {
    // isotropic! p0 is the same for rep1 and rep2
    double pressure = (rep1->sim.param().pcouple.pres0(0,0)
              + rep1->sim.param().pcouple.pres0(1,1)
              + rep1->sim.param().pcouple.pres0(2,2)) / 3.0;
    // get the volume
    double V1 = math::volume(rep1->conf.current().box, rep1->conf.boundary_type);
    double V2 = math::volume(rep2->conf.current().box, rep2->conf.boundary_type);
    // add the PV term to delta
    delta += pressure * (b1 - b2) * (V2 - V1);
  }

  if (delta < 0.0)
    return 1.0;
  else {
    return exp(-delta);
  }
}

/*
void util::replica_exchange_base::switch_coords_on_node(repIterator it1, const unsigned int partner) {
  replica* rep1 = *it1;
  replica* rep2;

  // find my partner
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    if ((*it)->ID == partner) {
      rep2 = (*it);
      break;
    }
  }

  rep1->conf.exchange_state();

  rep1->conf.current().pos = rep2->conf.current().pos;
  rep2->conf.old().pos = rep1->conf.old().pos;

  rep1->conf.current().posV = rep2->conf.current().posV;
  rep2->conf.old().posV = rep1->conf.old().posV;

  rep1->conf.current().vel = rep2->conf.current().vel;
  rep2->conf.old().vel = rep1->conf.old().vel;

  rep1->conf.current().stochastic_integral = rep2->conf.current().stochastic_integral;
  rep2->conf.old().stochastic_integral = rep1->conf.old().stochastic_integral;

  rep1->conf.current().box = rep2->conf.current().box;
  rep2->conf.old().box = rep1->conf.old().box;

  // temp array for lattice shifts
  math::VArray latticeTMP(rep1->conf.special().lattice_shifts);

  rep1->conf.special().lattice_shifts = rep2->conf.special().lattice_shifts; // lattice shifts?
  rep2->conf.special().lattice_shifts = latticeTMP;

  rep1->conf.current().phi = rep2->conf.current().phi;
  rep2->conf.old().phi = rep1->conf.old().phi;

  rep1->conf.current().psi = rep2->conf.current().psi;
  rep2->conf.old().psi = rep1->conf.old().psi;

  rep1->conf.current().theta = rep2->conf.current().theta;
  rep2->conf.old().theta = rep1->conf.old().theta;

  rep2->conf.exchange_state();
  // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
  rep1->exchange_averages();
  rep2->exchange_averages();
}
 * */

void util::replica_exchange_base::print_coords(std::string name) {
    
    io::Out_Configuration received_traj(GROMOSXX " Configuration", std::cout);
    io::Argument args2(args);
    args2.erase("fin");
    std::stringstream tmp;
    tmp << name << "_replica_" << replica->ID << "_run_" << replica->run;
    std::string fin = tmp.str() + ".cnf";
    args2.insert(std::pair<std::string, std::string > ("fin", fin));

    received_traj.init(args2, replica->sim.param());

    received_traj.title("bla");
    received_traj.write(replica->conf, replica->topo, replica->sim, io::final);

  /*
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    io::Out_Configuration received_traj(GROMOSXX " Configuration", std::cout);

    io::Argument args2(args);
    args2.erase("fin");
    std::stringstream tmp;
    tmp << name << "_replica_" << (*it)->ID << "_run_" << (*it)->run;
    std::string fin = tmp.str() + ".cnf";
    args2.insert(std::pair<std::string, std::string > ("fin", fin));

    received_traj.init(args2, (*it)->sim.param());

    received_traj.title("bla");
    received_traj.write((*it)->conf, (*it)->topo, (*it)->sim, io::final);
  }
  */
}
