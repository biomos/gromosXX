/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_exchange_master_eds.cc
 * Author: bschroed
 * 
 * Created on April 18, 2018, 3:20 PM
 */

#include <util/replicaExchange/replica_exchangers/replica_exchange_master_eds.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

//using namespace util;
util::replica_exchange_master_eds::replica_exchange_master_eds(io::Argument _args,
        int cont,
        int rank,
        int simulationRank,
        int simulationID,
        int simulationThreads,
        int _size,
        int _numReplicas,
        std::vector<int> repIDs,
        std::map<ID_t, rank_t> & repMap) :
        replica_exchange_base(_args, cont, rank, simulationRank, simulationID, simulationThreads, repIDs, repMap),
        replica_exchange_base_eds(_args, cont, rank, simulationRank, simulationID, simulationThreads, repIDs, repMap),
        replica_exchange_master(_args, cont, rank, simulationRank, simulationID, simulationThreads, _size, _numReplicas, repIDs, repMap),
        reedsParam(replica->sim.param().reeds) {
    #ifdef XXMPI
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":Constructor:\t START");
    
    DEBUG(4,"replica_exchange_master_eds "<< rank <<":Constructor:\t Next");
    //initialize data of replicas    
    replicaData.resize(repParams.num_l);      
    DEBUG(4,"replica_exchange_master_eds "<< rank <<":Constructor:\t replicaDatasize\t "<< replicaData.size());
    //DEBUG(4,"replica_exchange_master_eds "<< rank <<":Constructor:\t reeds- lambda\t "<< replica->sim.param().reeds.num_l);
    //DEBUG(4,"replica_exchange_master_eds "<< rank <<":Constructor:\t eds \t "<< replica->sim.param().eds.s.size());

        
    int ID = 0;
    for (int i = 0; i < repParams.num_l; ++i) {
        replicaData[ID].ID = ID;
        replicaData[ID].T = repParams.temperature[i];
        DEBUG(4,"replica_exchange_master_eds "<< rank <<":Constructor:\t Init Replicas ID"<<ID<<"\t "<< repParams.temperature[i]);
        replicaData[ID].l = repParams.lambda[i];
        replicaData[ID].dt = repParams.dt[i];
        DEBUG(4,"replica_exchange_master_eds "<< rank <<":Constructor:\t Init Replicas eds_param"<<ID<<"\t "<< reedsParam.eds_para[0].numstates);
        replicaData[ID].Vi.assign(reedsParam.eds_para[0].numstates,0);
        ++ID;
    }
    
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":Constructor:\t DONE\n");
    #else
   throw "Cannot initialize replica_exchange_master_eds without MPI!"; 
    #endif
}

void util::replica_exchange_master_eds::receive_from_all_slaves() {
  #ifdef XXMPI
  DEBUG(2,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t START\n");

  double start = MPI_Wtime();

  MPI_Status status;
  MPI_Status status_eds;

  util::repInfo info;

  // receive all information from slaves
  DEBUG(4,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t receive from slaves");
  DEBUG(4, "replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t numReps: "<< numReplicas)
  for (unsigned int rep = 0; rep < numReplicas; ++rep) {
    unsigned int rank = repMap.find(rep)->second;
    if (rank != 0) {
        DEBUG(2,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t \t out of rank ID: " << rank);
        MPI_Recv(&info, 1, MPI_REPINFO, rank, REPINFO, MPI_COMM_WORLD, &status);
        DEBUG(2,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t \t ID: " << rank);
        replicaData[rep].run = info.run;
        replicaData[rep].epot = info.epot;
        replicaData[rep].epot_partner = info.epot_partner;
        replicaData[rep].probability = info.probability;
        replicaData[rep].switched = info.switched;
        replicaData[rep].partner = info.partner;
        DEBUG(4,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t REP:" <<rep<< " EpotTot: "<< replicaData[rep].epot);

        MPI_Recv(&replicaData[rep].Vi[0],1, MPI_EDSINFO, rank, EDSINFO, MPI_COMM_WORLD, &status_eds);
        for(unsigned int s=0;s< replicaData[rep].Vi.size(); s++){
            DEBUG(4,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t "<< s << " En: "<< replicaData[rep].Vi[s]);
        }        
    }    
  }
  
  // write all information from master node to data structure
  DEBUG(2,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t  t\twrite_own data");
  for (repIterator it = replicas.begin(); it < replicas.end(); it++) {
    int ID = (*it)->ID;
    replicaData[ID].run = (*it)->run;
    replicaData[ID].partner = (*it)->partner;
    replicaData[ID].epot = (*it)->epot;
    replicaData[ID].epot_partner = (*it)->epot_partner;
    replicaData[ID].probability = (*it)->probability;
    replicaData[ID].switched = (*it)->switched;
    replicaData[ID].Vi = (*it)->conf.current().energies.eds_vi;
  }
  DEBUG(4,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t Master:\n" << "time used for receiving all messages: " << MPI_Wtime() - start 
            << " seconds");
  DEBUG(2,"replica_exchange_master_eds "<< rank <<":receive_from_all_slaves:\t DONE")
  #else
    throw "Cannot use receive_from_all_slaves from replica_exchange_master_eds without MPI!"; 
  #endif
}

int util::replica_exchange_master_eds::getSValPrecision(double minLambda){
  // this method basically counts the zeros between the dot and the first significant digit.
  // 0.0002560 => 4 digits=2
  int numDigits = 0;
  int minLambdaInt= int(minLambda);
  double minLambdaDecimals = minLambda - minLambdaInt;

  int num =100; //kill infloop!
  while (minLambdaDecimals > 0 && num > 0){
      num--;
      minLambdaInt = int(minLambdaDecimals);
      if(minLambdaInt == 0){
          minLambdaDecimals = minLambdaDecimals*10;
          numDigits++;
          continue;
      }
      else{
          break;
      }
  }
  return numDigits;
}

void util::replica_exchange_master_eds::write() {
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":write:\t START");
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":write:\t svalPrecision "<<svalPrecision);
    for (unsigned int r = 0; r < numReplicas; ++r) {
      repOut << std::setw(6) << (replicaData[r].ID + 1)
              << " "
              << std::setw(6) << (replicaData[r].partner + 1)
              << std::setw(6) << replicaData[r].run;
      repOut.precision(svalPrecision);
      DEBUG(2,"replica_exchange_master_eds "<< rank <<":write:\t s_val raus! "<<replicaData[r].l);

      repOut  << std::setw(13) << replicaData[r].l;
      repOut.precision(generalPrecision);
      repOut  << std::setw(13) << replicaData[r].T
              << " "
              << std::setw(18) << replicaData[r].epot;
      repOut.precision(svalPrecision);
      repOut  << std::setw(13) << replicaData[replicaData[r].partner].l;
      
      repOut.precision(generalPrecision);
      repOut  << std::setw(13) << replicaData[replicaData[r].partner].T;

      if(replicaData[r].l == replicaData[replicaData[r].partner].l){
          repOut << std::setw(18) << replicaData[replicaData[r].partner].epot;
      }
      else{
          repOut << std::setw(18) << replicaData[r].epot_partner;
      }
      repOut  << std::setw(13) << replicaData[r].probability
              << std::setw(6) << replicaData[r].switched;

    for(unsigned int s=0;s<reedsParam.eds_para[0].numstates; s++){
        DEBUG(2, "replica_exchange_master_eds "<< rank <<":write:\t WRITE out POTS " <<s<< "\t" << replicaData[r].Vi[s]);
        repOut   << std::setw(18) << std::min(replicaData[r].Vi[s], 10000000.0);//Output potential energies for each state 
    }
      repOut << std::endl;
    }
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":write:\t DONE");
}

void util::replica_exchange_master_eds::init_repOut_stat_file() {
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":init_repOut_stat_file:\t START");

    repOut.open(repdatName.c_str());
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":init_repOut_stat_file:\t repdat file open ");

    repOut << "#=========================================================================================================\n"
           << "#\tREPLICA Exchange - Enveloping Distribution Sampling - Exchange Output\n"
           << "#=========================================================================================================\n#\n"
           << "#Replica Exchange Parameters:\n"
           << "#Number of temperatures:\t" << repParams.num_T << "\n"
           << "#Number of s values:\t" << repParams.num_l << "\n";

    DEBUG(4,"replica_exchange_master_eds "<< rank <<":init_repOut_stat_file:\t set precision ");
    
    //double minS = *(std::min_element(reedsParam.eds_para[0].s.begin(), reedsParam.eds_para[0].s.end()));
    //svalPrecision = getSValPrecision(minS);
    repOut.setf(std::ios::fixed, std::ios::floatfield);

    repOut.precision(generalPrecision);
    repOut << "#T    \t";
    for (int t = 0; t < repParams.num_T; ++t){
        repOut << std::setw(12) << repParams.temperature[t];
    }
    repOut << "\n";
    repOut << "#s\t";
    repOut.precision(svalPrecision);
    int num_l = reedsParam.num_l;    
    for (int i = 0; i < num_l; ++i){
            repOut << std::setw(12) << reedsParam.eds_para[i].s[0];
    }

    repOut.precision(generalPrecision);
    for (int j = 0; j < reedsParam.num_states; ++j) {
        repOut << "\n# E"<< (j+1)<<"R(s)\t";
        for (int i = 0; i < num_l; ++i)
            repOut << "\t" << std::setw(4) << reedsParam.eds_para[i].eir[j];
    }  

    repOut << "\n#\n";
    repOut << ""
          << std::setw(6) << "ID"
          << " "
          << std::setw(6) << "partner"
          << std::setw(6) << "run"

          << std::setw(13) << "si"
          << std::setw(13) << "Ti"
          << std::setw(18) << "Epoti"

          << std::setw(13) << "sj"
          << std::setw(13) << "Tj"
          << std::setw(18) << "Epotj"
          << std::setw(13) << "p"
          << std::setw(6) << "exch";
    for (int x=1; x<= reedsParam.num_states; x++){
        repOut<<std::setw(17)<<"Vr"<<x;
    }

    repOut << "\n";
    repOut.flush();    
    DEBUG(2,"replica_exchange_master_eds "<< rank <<":init_repOut_stat_file:\t DONE");

}
