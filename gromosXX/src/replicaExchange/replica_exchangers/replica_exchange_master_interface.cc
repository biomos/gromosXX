/*
 * File:   replica_exchange_master.cc
 * Author: wissphil, sriniker
 *
 * Created on April 29, 2011, 2:18 PM
 * Modified June 18, 2021 - bschroed, srieder
 */

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
#include <util/replicaExchange/replica/replica.h>
#include <util/replicaExchange/replica_exchangers/replica_exchange_master_interface.h>

#include <string>

#ifdef XXMPI
    #include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

re::replica_exchange_master_interface::replica_exchange_master_interface(io::Argument & args,
        unsigned int cont,
        unsigned int globalThreadID,
        replica_graph_control &replicaGraphMPIControl,
        simulation::MpiControl &replica_mpi_control) :
        replica_exchange_base_interface(args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control),
        repParams(replica->sim.param().replica),
        repdatName(args["repdat"])
{
#ifdef XXMPI
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t START");

  assert(replicaGraphMPIControl.masterID == replicaGraphMPIControl.threadID);    //TODO: This can be removed in future! bscrhoed
  assert(replicaGraphMPIControl.numberOfReplicas > 0);
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t rep_params THERE?");
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t replica Num " << replica->sim.param().replica.num_l);
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t replica S" << replica->sim.param().replica.lambda[0]);

  assert(repParams.num_l > 0);

  DEBUG(4,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t Init Replicas \t Next");
  replicaData.resize(replicaGraphMPIControl.numberOfReplicas);
  DEBUG(4,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t Replica_data type \t " << typeid(replicaData).name());

  //initialize data of replicas
  int ID = 0;
  for (int i = 0; i < repParams.num_l; ++i) {
    for (int j = 0; j < repParams.num_T; ++j) {
      replicaData[ID].ID = ID;
      replicaData[ID].T = repParams.temperature[j];
      DEBUG(4,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t Init Replicas \t "<< repParams.temperature[j]);
      replicaData[ID].l = repParams.lambda[i];
      replicaData[ID].dt = repParams.dt[i];
      ++ID;
    }
  }

  // set output file
 DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":Constructor:\t DONE");
#else
   throw "Cannot initialize replica_exchange_master without MPI!";
#endif
}


re::replica_exchange_master_interface::~replica_exchange_master_interface() {
   repOut.close();
}

void re::replica_exchange_master_interface::receive_from_all_slaves() {
    #ifdef XXMPI
    DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":receive_from_all_slaves:\t START\n");
    double start = MPI_Wtime();

    MPI_Status status;
    re::repInfo info;

    // receive all information from slaves
    for (unsigned int slaveReplicaID = 0; slaveReplicaID < replicaGraphMPIControl().numberOfReplicas; ++slaveReplicaID) {
      if (slaveReplicaID != replicaGraphMPIControl().masterID) {
        unsigned int replicaMasterThreadID = slaveReplicaID;
        DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":receive_from_all_slaves:\t get_MPI replicaMasterThread: "<< replicaMasterThreadID << "\n");
        MPI_Recv(&info, 1, MPI_REPINFO, replicaMasterThreadID, REPINFO,  replicaGraphMPIControl().comm, &status);
        replicaData[slaveReplicaID].run = info.run;
        replicaData[slaveReplicaID].epot = info.epot;
        replicaData[slaveReplicaID].epot_partner = info.epot_partner;
        replicaData[slaveReplicaID].probability = info.probability;
        replicaData[slaveReplicaID].switched = info.switched;
        replicaData[slaveReplicaID].partner = info.partner;
      }
      DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":receive_from_all_slaves:\t got all MPI reps\n");

      // write all information from master node to data structure
      replicaData[simulationID].run = run;
      replicaData[simulationID].partner = partnerReplicaID;
      replicaData[simulationID].epot = epot;
      replicaData[simulationID].epot_partner = epot_partner;
      replicaData[simulationID].probability = probability;
      replicaData[simulationID].switched = switched;
    }
    DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":receive_from_all_slaves:\t " << "time used for receiving all messages: " << MPI_Wtime() - start << " seconds\n");
    DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":receive_from_all_slaves:\t DONE: \n");
    #else
     throw "Cannot use replica_exchange_master without MPI!";
    #endif
}


void re::replica_exchange_master_interface::init_repOut_stat_file() {
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t START");
  repOut.open(repdatName.c_str());
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t  repdat file open ");

  repOut << "Number of temperatures:\t" << repParams.num_T << "\n"
         << "Number of lambda values:\t" << repParams.num_l << "\n";

  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t set precision ");
  repOut.precision(4);
  repOut.setf(std::ios::fixed, std::ios::floatfield);

  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t write Temperatures ");
  repOut << "T    \t";
  for (int t = 0; t < repParams.num_T; ++t){
    DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t it: "<<  t);
    DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t T: "<<  repParams.temperature[t]);
    repOut << std::setw(12) << repParams.temperature[t];
  }

  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":init_repOut_stat_file:\t write lambdas ");
  repOut << "\nlambda    \t";
  for (int l = 0; l < repParams.num_l; ++l){
    repOut << std::setw(12) << repParams.lambda[l];
  }

  repOut << "\n\n";

  repOut << "#"
          << std::setw(6) << "ID"
          << " "
          << std::setw(6) << "partner"
          << std::setw(6) << "run"
          << " "
          << std::setw(13)  << "li"
          << std::setw(13)  << "Ti"
          << std::setw(18)  << "Epoti"
          << std::setw(13)  << "lj"
          << std::setw(13)  << "Tj"
          << std::setw(18)  << "Epotj"
          << std::setw(13)  << "p"
          << std::setw(6) << "exch";
  repOut << "\n";
}


void re::replica_exchange_master_interface::write() {
   DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":write:\t START");

  for (unsigned int treplicaID = 0; treplicaID < replicaGraphMPIControl().numberOfReplicas; ++treplicaID) {
    repOut << std::setw(6) << (replicaData[treplicaID].ID + 1)
            << " "
            << std::setw(6) << (replicaData[treplicaID].partner + 1)
            << std::setw(6) << replicaData[treplicaID].run
            << std::setw(13) << replicaData[treplicaID].l
            << std::setw(13) << replicaData[treplicaID].T
            << " "
            << std::setw(18) << replicaData[treplicaID].epot
            << std::setw(13) << replicaData[replicaData[treplicaID].partner].l
            << std::setw(13) << replicaData[replicaData[treplicaID].partner].T
            << " ";
    if(replicaData[treplicaID].l == replicaData[replicaData[treplicaID].partner].l)
	repOut << std::setw(18) << replicaData[replicaData[treplicaID].partner].epot;
    else
        repOut << std::setw(18) << replicaData[treplicaID].epot_partner;
    repOut  << std::setw(13) << replicaData[treplicaID].probability
            << std::setw(6) << replicaData[treplicaID].switched
            << std::endl;
  }
  DEBUG(2,"replica_exchange_master_interface "<< globalThreadID <<":write:\t DONE");

}
