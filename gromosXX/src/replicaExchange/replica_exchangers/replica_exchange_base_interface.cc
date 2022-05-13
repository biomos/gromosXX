/*
 * File:   replica_exchange_base_2d_l_T_HREMD.cc
 * Author: wissphil, sriniker
 *
 * Created on April 29, 2011, 2:06 PM
 * Modified June 18, 2021 - bschroed, srieder
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
#include "util/debug.h"
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


re::replica_exchange_base_interface::replica_exchange_base_interface(io::Argument _args,
                                                   unsigned int cont,
                                                   unsigned int globalThreadID,
                                                   replica_graph_control &replicaGraphMPIControl,
                                                   simulation::MpiControl &replica_mpi_control) :
        args(_args), rng(-1), cont(cont),
        globalThreadID(globalThreadID), simulationID(replica_mpi_control.simulationID),
        m_replicaGraphMPIControl(replicaGraphMPIControl)
{
#ifdef XXMPI
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t START ");\
  createReplicas(cont, globalThreadID, replica_mpi_control);

  // random generator
  std::stringstream seed;
  seed << replica->sim.param().start.ig*(simulationID+1);
  rng.seed(seed.str());

  
  //Copy by hand....
  m_replicaGraphMPIControl.graphID = replicaGraphMPIControl.graphID;
  m_replicaGraphMPIControl.masterID = replicaGraphMPIControl.masterID;
  m_replicaGraphMPIControl.threadID = replicaGraphMPIControl.threadID;
  m_replicaGraphMPIControl.numberOfThreads = replicaGraphMPIControl.numberOfThreads;
  m_replicaGraphMPIControl.numberOfReplicas = replicaGraphMPIControl.numberOfReplicas;
  m_replicaGraphMPIControl.mpiColor = replicaGraphMPIControl.mpiColor;
  m_replicaGraphMPIControl.replicaMasterIDs = replicaGraphMPIControl.replicaMasterIDs;
  m_replicaGraphMPIControl.replicaThreads = replicaGraphMPIControl.replicaThreads;
  m_replicaGraphMPIControl.threadReplicaMap = replicaGraphMPIControl.threadReplicaMap;
  m_replicaGraphMPIControl.comm = replicaGraphMPIControl.comm;
  

  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":Constructor:\t Constructor \t DONE ");
  #else
    throw "Cannot use send_to_master from replica_exchange_slave_eds without MPI!";
  #endif
}


re::replica_exchange_base_interface::~replica_exchange_base_interface() {
    delete replica;
}

void re::replica_exchange_base_interface::createReplicas(int cont, int globalThreadID, simulation::MpiControl &replica_mpi_control){
  MPI_DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":createReplicas:\t START \t THREADS "<<replica_mpi_control.numberOfThreads);

  DEBUG(3,"\n\nreplica_exchange_base_interface: CREATEREPLICAS\n\n");
  // create the number of replicas that are assigned to my node
    if(replica_mpi_control.numberOfThreads>1){

        if(replica_mpi_control.threadID == replica_mpi_control.masterID ){
            replica = new re::replica_MPI_Master(args, cont, globalThreadID, replica_mpi_control);
        }
        else{
            replica = new re::replica_MPI_Slave(args, cont, globalThreadID, replica_mpi_control);
        }
    }
    else{
        replica = new re::replica(args, cont, globalThreadID, replica_mpi_control);
    }
  DEBUG(4,"replica_exchange_base_interface "<< globalThreadID <<":create:\t replica numS " << replica->sim.param().replica.num_l);
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":createReplicas:\t DONE");
}

void re::replica_exchange_base_interface::init() {
  DEBUG(3, "replica_exchange_base_interface "<< globalThreadID <<":init:\t START EMPTY");
  DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":init:\t DONE");
}

void re::replica_exchange_base_interface::run_MD() {
  MPI_DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":run_MD:\t START");
    replica->sim.steps() = steps;
    replica->sim.time() = time;
    replica->run_MD();

    if(replicaInfoSender){
        updateReplica_params();
    }else{
        re::replica_MPI_Slave * rs = dynamic_cast<re::replica_MPI_Slave* >(replica);
        rs->calculateEnergiesHelper();
    }
    ++run;

  MPI_DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":run_MD:\t DONE");
}


//TODO: Maybe REMOVE in future?
void re::replica_exchange_base_interface::updateReplica_params(){
    DEBUG(3,"\n\nreplica_exchange_base_interface "<< globalThreadID <<": UPDATEREPLICA_PARAMS\n\n");
  // update replica information
  time = replica->sim.time();
  steps = replica->sim.steps();
  epot = replica->calculateEnergies();
}
/*
 * REplica Exchanges
 */

void re::replica_exchange_base_interface::swap(){
        DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":swap:    START");
        #ifdef XXMPI
            // find partner
            DEBUG(3,"replica_exchange_base_interface"<< globalThreadID <<":determine_swaps:\t mySwapBuddy: "<< partnerReplicaID)
            partnerReplicaID = find_partner();

            if (replicaInfoSender) // different replica?
            {
                //determine swaps
                determine_swaps();

                //final execution of swap
                if (switched) {
                    execute_swap(partnerReplicaID);
                }
                if(switched && replica->sim.param().replica.scale) {
                    velscale(partnerReplicaID);
                }
            }
            else {  // no exchange for helper replicas (they are only used in md_mpi run)
                swapHelper();
            }
        #endif
        DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":swap:    DONE");
}

void re::replica_exchange_base_interface::swapHelper(){
    DEBUG(5,"replica_exchange_base_interface "<< globalThreadID <<":swap:\t I'm not a swapper!");
    #ifdef XXMPI
        //do we potentially need to do FF calcs?
        if(replica->sim.param().replica.num_l>1){
            //first update
            DEBUG(5, "replica_exchange_base_interface"<<globalThreadID<<":  Slave\texchangeSTART " << exchange)
            MPI_Bcast(&exchange , 1, MPI_INT, replica->replica_mpi_control.masterID, replica->replica_mpi_control.comm);
            DEBUG(5, "replica_exchange_base_interface"<<globalThreadID<<":  Slave\texchangeEND " << exchange)

            if(exchange){  // TH-RE-MD
                DEBUG(5, "replica"<<globalThreadID<<":replica_exchange_base_interface:  \t Slave\tff help_energy")
                if((replica->sim.param().replica.num_l>1 && replica->sim.param().replica.num_T == 1) ||               //case of 1D H-RE MD
                (replica->sim.param().replica.num_l>1 && replica->sim.param().replica.num_T>1 && this->run%4%2==1))   //case of 2D TH-RE MD only Lam ex.
                {
                    calculate_energy_helper(partnerReplicaID);
                }
            }
        }
        if(replica->sim.param().reeds.reeds>0){   //1D/2D RE-EDS
            replica->sim.param().reeds.eds_para[partnerReplicaID].pos_info.second = simulationID;
        }
    #endif

    probability = 0.0;
    switched = 0;

}

void re::replica_exchange_base_interface::write_final_conf() {
  // write coordinates to cnf for all replica assigned to this node
    if(replicaInfoSender){
       replica->write_final_conf();
    }
}

//RE
//SWAPPING Functions
void re::replica_exchange_base_interface::determine_swaps() {
    DEBUG(3,"replica_exchange_base_interface"<< globalThreadID <<":determine_swap:\t START");
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":partnerID:partnerRep: "<<partnerReplicaID<<"\tsimulationID"<<simulationID);
    #ifdef XXMPI
        // get swap likelihood and discretize to bool
        exchange = int(partnerReplicaID != simulationID);
        if(replica->sim.param().replica.num_l>1 && replica->replica_mpi_control.numberOfThreads>1){ // if lambdaSims on, we need to do a forceField eval, to calc lam change
            DEBUG(5, "replica_exchange_base_interface"<<globalThreadID<<":  Master\texchangeSTART " << exchange);
            MPI_Bcast(&exchange, 1, MPI_INT, replica->replica_mpi_control.masterID, replica->replica_mpi_control.comm);
        }

        if(exchange){
            DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":swap:\t DO swap!!:)");
            determine_switch_probabilities();
        } else {
            DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":swap:\t No swap here!!:)");
            probability = 0.0;
            switched = 0;
        }
    #endif
    DEBUG(3,"replica_exchange_base_interface"<< globalThreadID <<":determine_swap:\t DONE");
}

/**OVERRIDE THIS NICE FUNCTION*/
void re::replica_exchange_base_interface::determine_switch_probabilities(){
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":determineSwitchPos:\t START");
    throw "Implement the exchange probability Function!";
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":determineSwitchPos:\t DONE");
}

/**OVERRIDE THIS NICE FUNCTION*/
int re::replica_exchange_base_interface::find_partner() const {
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":find_partner:\t START");
    throw "Implement the find_partner Function!";
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":find_partner:\t DONE");
}


//Execute Swapping
void re::replica_exchange_base_interface::execute_swap(const unsigned int partnerReplicaID) {
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":executeSwap:\t START");
    if (simulationID < partnerReplicaID) {
      send_coord(partnerReplicaID);
      receive_new_coord(partnerReplicaID);
    } else {
      receive_new_coord(partnerReplicaID);
      send_coord(partnerReplicaID);
    }
    // the averages of current and old are interchanged after calling exchange_state() and have to be switched back
    exchange_averages();
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":executeSwap:\t DONE");
}

/**OVERRIDE Possibly THIS NICE FUNCTION*/
double re::replica_exchange_base_interface::calc_probability(const unsigned int partnerReplicaMasterThreadID) {

  DEBUG(1,"\n\nreplica "<<globalThreadID<<":replica_exchange_base_interface: CALC_PROBABILITY by ID:" << simulationID << "\n");
  double delta = 0.0;
  const double b1 = 1.0 / (math::k_Boltzmann * T);
  const double b2 = 1.0 / (math::k_Boltzmann * replica->sim.param().replica.temperature[partnerReplicaMasterThreadID % replica->sim.param().replica.num_T]);

  bool sameLambda = (l == replica->sim.param().replica.lambda[partnerReplicaMasterThreadID / replica->sim.param().replica.num_T]);// horizontal switch with same lambda?

  if (sameLambda) {
    // use simple formula
    // get energy from the other partner
    double energies[2] = {0.0, 0.0};
#ifdef XXMPI
    MPI_Status status;
    DEBUG(1,"\n\nreplica_exchange_base_interface: CALC_PROBABILITY before Recv\n");
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES,  replicaGraphMPIControl().comm, &status);
    DEBUG(1,"\n\nreplica_exchange_base_interface: CALC_PROBABILITY after Recv\n");
#endif
    epot_partner = energies[0];

    //DEBUG
    DEBUG(7, "replica "<<globalThreadID<<" Exchange: E1" << epot<< "\n");
    DEBUG(7, "replica "<<globalThreadID<<" Exchange: E2" << epot_partner<< "\n");
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
    DEBUG(1,"\n\nreplica_exchange_base_interface: CALC_PROBABILITY before Recv\n");
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerReplicaMasterThreadID, SWITCHENERGIES, replicaGraphMPIControl().comm, &status);
    DEBUG(1,"\n\nreplica_exchange_base_interface: CALC_PROBABILITY after Recv\n");
#endif
    const double E22 = energies[0];
    const double E12 = energies[1];

    const double E11 = epot;
    const double E21 = calculate_energy(partnerReplicaMasterThreadID);

    // store this as the partner energy
    epot_partner = E21;

    //DEBUG
    DEBUG(7, "\n\nb1: " << b1 << " b2: " << b2 << std::endl);
    DEBUG(7, "replica "<<globalThreadID<<" Exchange:E11: " << E11 << " E22: " << E22 << std::endl);
    DEBUG(7, "replica "<<globalThreadID<<" Exchange: E21: " << E21 << " E12: " << E12 << std::endl<<std::endl);


    delta = b1 * (E12 - E11) - b2 * (E22 - E21);
  }
  DEBUG(5,"replica_exchange_base_interface: CALC_PROBABILITY - before pressure delta: "<<delta << "\n");

  // NPT? add PV term
  if (replica->sim.param().pcouple.scale != math::pcouple_off)
  {
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
  DEBUG(5,"replica_exchange_base_interface: CALC_PROBABILITY - delta: "<<delta << "\n");

  if (delta < 0.0){
    return 1.0;
  }
  else {
    return exp(-delta);
  }
}


void re::replica_exchange_base_interface::setParams(){
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t START");
    throw "Implement the setParams Function!";
    DEBUG(3,"replica_exchange_base_interface "<< globalThreadID <<":setParams:\t DONE");
}

double re::replica_exchange_base_interface::calculate_energy(const unsigned int partnerThreadID) {
    throw "Implement the calculate_energy Function!";
}

void re::replica_exchange_base_interface::calculate_energy_helper(const unsigned int partnerThreadID) {
    throw "Implement the calculate_energy_helper Function!";
}

void re::replica_exchange_base_interface::exchange_averages() {
  DEBUG(3,"replica"<<globalThreadID<<":replica_exchange_base_interface: EXCHANGE_AVERAGES\n\n");
  // after a swap the averages of current and old are exchanged and have to be switched back
  configuration::Average  dummy = replica->conf.current().averages;
  replica->conf.current().averages=replica->conf.old().averages;
  replica->conf.old().averages=dummy; // <- why not removing?
}

//sending stuff
void re::replica_exchange_base_interface::send_coord(const unsigned int receiverReplicaID) {
  DEBUG(5,"replica"<<globalThreadID<<":replica_exchange_base_interface: SEND_COORD\n\n");
#ifdef XXMPI
  unsigned int receiverReplicaMasterThreadID = receiverReplicaID;

  configuration::Configuration& conf = replica->conf;
  DEBUG(5,"replica"<<globalThreadID<<":replica_exchange_base_interface: get conf\n\n");

  MPI_Send(&conf.current().pos[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, POS,  replicaGraphMPIControl().comm);
  MPI_Send(&conf.current().posV[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, POSV, replicaGraphMPIControl().comm);
  MPI_Send(&conf.current().vel[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, VEL,  replicaGraphMPIControl().comm);
  DEBUG(5,"replica"<<globalThreadID<<":replica_exchange_base_interface: SEND_COORD POS,POSV, VEL\n\n");

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if (simulationID < receiverReplicaID) {
    MPI_Send(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, LATTSHIFTS,  replicaGraphMPIControl().comm);
  } else {
    MPI_Send(&((*replica->latticeTMP)[0][0]), 1, MPI_VARRAY, receiverReplicaMasterThreadID, LATTSHIFTS,  replicaGraphMPIControl().comm);
    delete replica->latticeTMP;
    replica->latticeTMP = NULL;
  }
  DEBUG(5,"replica"<<globalThreadID<<":replica_exchange_base_interface: SEND_COORD latticeShifts\n\n");

  MPI_Send(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, receiverReplicaMasterThreadID, STOCHINT,  replicaGraphMPIControl().comm);
  MPI_Send(&conf.current().box(0)[0], 1, MPI_BOX, receiverReplicaMasterThreadID, BOX,  replicaGraphMPIControl().comm);

  std::vector<double> angles;
  angles.resize(3);
  angles[0] = conf.current().phi;
  angles[1] = conf.current().psi;
  angles[2] = conf.current().theta;

  DEBUG(5,"replica"<<globalThreadID<<":replica_exchange_base_interface: SEND_COORD angles\n\n");
  MPI_Send(&angles[0], angles.size(), MPI_DOUBLE, receiverReplicaMasterThreadID, ANGLES,  replicaGraphMPIControl().comm);
  MPI_Send(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, receiverReplicaMasterThreadID, DF,  replicaGraphMPIControl().comm);
  
  if ( simulationID > receiverReplicaID){
    conf.exchange_state();
  }
  DEBUG(5,"replica"<<globalThreadID<<":replica_exchange_base_interface: DONE angles\n\n");

#endif
}

void re::replica_exchange_base_interface::receive_new_coord(const unsigned int senderReplicaID) {
  DEBUG(3,"replica"<<globalThreadID<<":replica_exchange_base_interface: RECEIVE_NEW_COORD\n\n");
#ifdef XXMPI
  unsigned int senderReplicaMasterThreadID = senderReplicaID;

  MPI_Status status;
  configuration::Configuration&  conf = replica->conf;

  conf.exchange_state();
  MPI_Recv(&conf.current().pos[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, POS,  replicaGraphMPIControl().comm, &status);
  MPI_Recv(&conf.current().posV[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, POSV,  replicaGraphMPIControl().comm, &status);
  MPI_Recv(&conf.current().vel[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, VEL,  replicaGraphMPIControl().comm, &status);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if ( globalThreadID > senderReplicaMasterThreadID){
    replica->latticeTMP = new math::VArray(conf.special().lattice_shifts);
  }

  MPI_Recv(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, LATTSHIFTS,  replicaGraphMPIControl().comm, &status);
  MPI_Recv(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, senderReplicaMasterThreadID, STOCHINT,  replicaGraphMPIControl().comm, &status);
  MPI_Recv(&conf.current().box(0)[0], 1, MPI_BOX, senderReplicaMasterThreadID, BOX,  replicaGraphMPIControl().comm, &status);

  std::vector<double> angles;
  angles.resize(3);
  MPI_Recv(&angles[0], angles.size(), MPI_DOUBLE, senderReplicaMasterThreadID, ANGLES,  replicaGraphMPIControl().comm, &status);

  conf.current().phi = angles[0];
  conf.current().psi = angles[1];
  conf.current().theta = angles[2];

  MPI_Recv(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, senderReplicaMasterThreadID, DF,  replicaGraphMPIControl().comm, &status);

  if ( simulationID > senderReplicaID){
    conf.exchange_state();
  }
#endif
}

//TRE
void re::replica_exchange_base_interface::velscale(int unsigned partnerReplica){
  DEBUG(5,"\n\nreplica_exchange_base_interface: VELSCALE\n\n");
  double T1 = replica->sim.param().replica.temperature[simulationID];
  double T2 = replica->sim.param().replica.temperature[partnerReplica];
  if (T1 != T2) {
    double factor = sqrt(T1/T2);
    for (unsigned int k = 0; k < replica->topo.num_atoms(); ++k) {
     replica->conf.current().vel(k) *= factor;
    }
  }
}

// IO
void re::replica_exchange_base_interface::print_coords(std::string name) {

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

// Getter&Setter
unsigned int re::replica_exchange_base_interface::get_stepsPerRun(){
    return stepsPerRun;
};


void re::replica_exchange_base_interface::print_info(std::string bla) const {
  std::cout << "\n" << bla << std::endl;
  std::cout << "#"
          << std::setw(5) << "ID"
          << " "
          << std::setw(7) << "partner"
          << std::setw(7) << "run"
          << std::setw(14) << "Epoti"
          << std::setw(14) << "Epotj"
          << std::setw(13) << "p"
          << std::setw(4) << "s"
          << "\n";

  std::cout << std::setw(6) << (simulationID+1)
          << " "
          << std::setw(6) << partnerReplicaID
          << std::setw(6) << run
          << " "
          << std::setw(18) << epot
          << " "
          << std::setw(18) << epot_partner
          << std::setw(13) << probability
          << std::setw(4) << switched
          << std::endl;
}
