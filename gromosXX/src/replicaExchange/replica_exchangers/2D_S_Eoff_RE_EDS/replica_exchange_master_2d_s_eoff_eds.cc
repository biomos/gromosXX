/*
 * File:   replica_exchange_master_2d_s_eoff_eds.cc
 * Author: theosm
 *
 * Created on March 29, 2020, 11:00 AM
 */
#include <replicaExchange/replica_exchangers/2D_S_Eoff_RE_EDS/replica_exchange_master_2d_s_eoff_eds.h>

#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica_exchanger

re::replica_exchange_master_2d_s_eoff_eds::replica_exchange_master_2d_s_eoff_eds(io::Argument _args,
                                                                unsigned int cont,
                                                                unsigned int globalThreadID,
                                                                replica_graph_control &replicaGraphMPIControl,
                                                                simulation::MpiControl &replica_mpi_control) :
        replica_exchange_base_interface(_args, cont, globalThreadID,  replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_base_2d_s_eoff_eds(_args, cont, globalThreadID,  replicaGraphMPIControl, replica_mpi_control),
        replica_exchange_master_interface(_args, cont, globalThreadID, replicaGraphMPIControl, replica_mpi_control)
{
    #ifdef XXMPI
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t START");
    DEBUG(3,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t Replicas: "<<replicaGraphMPIControl.numberOfReplicas);
    DEBUG(3,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t ReplicasOLD: "<<repParams.num_l);
    DEBUG(3,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t ReplicasMASTER: "<< replicaGraphMPIControl.masterID);


    //initialize data of replicas
    replicaData.resize(replicaGraphMPIControl.numberOfReplicas);
    DEBUG(3,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t replicaDatasize\t "<< replicaData.size());
    //DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< rank <<":Constructor:\t reeds- lambda\t "<< replica->sim.param().reeds.num_s);
    //DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< rank <<":Constructor:\t eds \t "<< replica->sim.param().eds.s.size());

    //initialize bookkeeping data structure
    coordIDPositionsVector.resize(replicaGraphMPIControl.numberOfReplicas);
    DEBUG(3,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t coordIDPositionsVectorsize\t "<< coordIDPositionsVector.size());


    for (int replicaID = 0; replicaID<  replicaData.size(); ++replicaID) {
        replicaData[replicaID].ID = replicaID;
        //set position info
        replicaData[replicaID].pos_info = std::make_pair(replicaID, replicaID);
        DEBUG(1, "MASTER Constructor with replicaID, pos_info= " << replicaID << ", "
        << replicaData[replicaID].pos_info.first << ", " << replicaData[replicaID].pos_info.second << "\n");
        reedsParam.eds_para[replicaID].pos_info = replicaData[replicaID].pos_info;

        //just to check
        std::pair<int, int> a = replica->sim.param().reeds.eds_para[replicaID].pos_info;
        DEBUG(1, "JUST TO CHECK: MASTER Constructor with replicaID, replica->pos_info= " << replicaID << ", "
        << a.first << ", " << a.second << "\n");

        replicaData[replicaID].T = repParams.temperature[replicaID];
        DEBUG(5,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t Init Replicas ID"<<replicaID<<"\t "<< repParams.temperature[replicaID]);
        replicaData[replicaID].l = repParams.lambda[replicaID];
        replicaData[replicaID].dt = repParams.dt[replicaID];
        DEBUG(5,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t Init Replicas eds_param"<<replicaID<<"\t "<< reedsParam.eds_para[0].numstates);
        replicaData[replicaID].Vi.assign(reedsParam.eds_para[0].numstates,0);
    }
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t SIMID "<< simulationID <<"\n");
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":Constructor:\t DONE\n");
    #else
        throw "Cannot initialize replica_exchange_master_2d_s_eoff_eds without MPI!";
    #endif
}

void re::replica_exchange_master_2d_s_eoff_eds::receive_from_all_slaves() {
  #ifdef XXMPI
  DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t START\n");

  double start = MPI_Wtime();

  MPI_Status status;
  MPI_Status status_eds;

  re::repInfo info;
  //theosm -- new
  bool begin = 0;
  int real_pos;

  // receive all information from slaves
  DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t receive from slaves");
  DEBUG(4, "replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t numReps: "<< replicaGraphMPIControl().numberOfReplicas)
  for (unsigned int slaveReplicaID = 0; slaveReplicaID < replicaGraphMPIControl().numberOfReplicas; ++slaveReplicaID) {
    if (slaveReplicaID != replicaGraphMPIControl().masterID) {
        DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t \t out of rank ID: " << replicaGraphMPIControl().masterID);
        MPI_Recv(&info, 1, MPI_REPINFO, slaveReplicaID, REPINFO, replicaGraphMPIControl().comm, &status);
        DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t \t ID: " << replicaGraphMPIControl().masterID);
        replicaData[slaveReplicaID].run = info.run;
        replicaData[slaveReplicaID].epot = info.epot;
        replicaData[slaveReplicaID].epot_partner = info.epot_partner;
        replicaData[slaveReplicaID].probability = info.probability;
        replicaData[slaveReplicaID].switched = info.switched;
        replicaData[slaveReplicaID].partner = info.partner;

        if(replicaData[slaveReplicaID].partner == 0 && info.run == 1){
          begin = 1;
          real_pos = slaveReplicaID;
          DEBUG(3, "begin, real_pos & slaveReplicaID= " << begin << ", " << real_pos << ", " << slaveReplicaID << "\n");
        }

        //check whether switched or not for the first trial
        //assuming switched is treated as a boolean
        if(info.switched && info.run == 1){
          replicaData[slaveReplicaID].pos_info.second = reedsParam.eds_para[replicaData[slaveReplicaID].partner].pos_info.second;
        }


        //normal case
        if(info.run > 1 && info.switched){
          begin = 0;
          replicaData[slaveReplicaID].pos_info.second = coordIDPositionsVector[info.partner];
        }


        DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t REP:" <<slaveReplicaID<< " EpotTot: "<< replicaData[slaveReplicaID].epot);

        MPI_Recv(&replicaData[slaveReplicaID].Vi[0],1, MPI_EDSINFO, slaveReplicaID, EDSINFO, replicaGraphMPIControl().comm, &status_eds);
        for(unsigned int s=0;s< replicaData[slaveReplicaID].Vi.size(); s++){
            DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< replicaGraphMPIControl().masterID <<":receive_from_all_slaves:\t "<< s << " En: "<< replicaData[slaveReplicaID].Vi[s]);
        }
    }

  }

  // write all information from master node to data structure
  DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t  t\twrite_own data");
  replicaData[simulationID].run = run;
  replicaData[simulationID].partner = partnerReplicaID;
  replicaData[simulationID].epot = epot;
  replicaData[simulationID].epot_partner = epot_partner;
  replicaData[simulationID].probability = probability;
  replicaData[simulationID].switched = switched;
  //theosm
  if(switched){
    replicaData[simulationID].pos_info.second = coordIDPositionsVector[partnerReplicaID];
  }

  //check whether switched or not for the first trial
  if(switched && begin){
    replicaData[simulationID].pos_info.second = real_pos;
  }

  if (switched){
    replicaData[simulationID].Vi = replica->conf.old().energies.eds_vi;
  } else {
    replicaData[simulationID].Vi = replica->conf.current().energies.eds_vi;
  }


  DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t Master:\n" << "time used for receiving all messages: " << MPI_Wtime() - start
            << " seconds");
  DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":receive_from_all_slaves:\t DONE")
  #else
    throw "Cannot use receive_from_all_slaves from replica_exchange_master_2d_s_eoff_eds without MPI!";
  #endif
}

int re::replica_exchange_master_2d_s_eoff_eds::getSValPrecision(double minLambda){
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

void re::replica_exchange_master_2d_s_eoff_eds::write() {
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":write:\t START");
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":write:\t svalPrecision "<<svalPrecision);

    for (unsigned int treplicaID = 0; treplicaID < replicaGraphMPIControl().numberOfReplicas; ++treplicaID) {
      repOut << std::setw(6) << (replicaData[treplicaID].ID)//removed  + 1 for consistency reasons
              << " "
              << std::setw(6) << replicaData[treplicaID].pos_info.first
              << " "
              << std::setw(6) << replicaData[treplicaID].pos_info.second
              << " "
              << std::setw(6) << (replicaData[treplicaID].partner) //removed  + 1 for consistency reasons
              << "   "
              << std::setw(6) << replicaData[replicaData[treplicaID].partner].pos_info.first
              << "\t\t"
              << std::setw(6) << replicaData[replicaData[treplicaID].partner].pos_info.second
              << "\t\t"
              << std::setw(6) << replicaData[treplicaID].run << "  ";

      coordIDPositionsVector[treplicaID] = replicaData[treplicaID].pos_info.second;
      DEBUG(1,"write: ID, coordIDPositionsVector "<< treplicaID << ", " << coordIDPositionsVector[treplicaID] << "\n");
      /*
      repOut.precision(svalPrecision);
      DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":write:\t s_val raus! "<<replicaData[treplicaID].l);

      //do not use it anymore since the info is in header
      //repOut  << std::setw(13) << replicaData[treplicaID].l;
      */

      repOut.precision(generalPrecision);
      //do not use temperatures anymore
      //repOut  << std::setw(13) << replicaData[treplicaID].T
      //be careful when commenting in -- extra repOut
      repOut  << " "
              << std::setw(18) << replicaData[treplicaID].epot;

      /*
      //do not use it anymore since the info is in header
      repOut.precision(svalPrecision);
      repOut  << std::setw(13) << replicaData[replicaData[treplicaID].partner].l;
      */

      repOut.precision(generalPrecision);
      //do not use temperatures anymore
      //repOut  << std::setw(13) << replicaData[replicaData[treplicaID].partner].T;

      if(replicaData[treplicaID].l == replicaData[replicaData[treplicaID].partner].l){
          repOut << std::setw(18) << replicaData[replicaData[treplicaID].partner].epot;
      }
      else{
          repOut << std::setw(18) << replicaData[treplicaID].epot_partner;
      }
      repOut  << std::setw(13) << replicaData[treplicaID].probability
              << std::setw(6) << replicaData[treplicaID].switched;

    for(unsigned int s=0;s<reedsParam.eds_para[0].numstates; s++){
        DEBUG(2, "replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":write:\t WRITE out POTS " <<s<< "\t" << replicaData[treplicaID].Vi[s]);
        repOut   << std::setw(18) << std::min(replicaData[treplicaID].Vi[s], 10000000.0);//Output potential energies for each state
    }
      repOut << std::endl;
    }
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":write:\t DONE");
}

void re::replica_exchange_master_2d_s_eoff_eds::init_repOut_stat_file() {
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":init_repOut_stat_file:\t START");

    repOut.open(repdatName.c_str());
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":init_repOut_stat_file:\t repdat file open ");

    repOut << "#=========================================================================================================\n"
           << "#\tREPLICA Exchange - Enveloping Distribution Sampling - Exchange Output\n"
           << "#=========================================================================================================\n#\n"
           << "#Replica Exchange Parameters:\n"
           << "#Number of temperatures:\t" << repParams.num_T << "\n"
           << "#Number of s values:\t" << repParams.num_l << "\n";

    DEBUG(4,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":init_repOut_stat_file:\t set precision ");

    //double minS = *(std::min_element(reedsParam.eds_para[0].s.begin(), reedsParam.eds_para[0].s.end()));
    //svalPrecision = getSValPrecision(minS);
    repOut.setf(std::ios::fixed, std::ios::floatfield);

    repOut.precision(generalPrecision);
    repOut << "#T    \t";
    for (int t = 0; t < repParams.num_T; ++t){
        repOut << std::setw(12) << repParams.temperature[t];
    }
    repOut << "\n";
    repOut << "#Start coordinate == Position coordinate\n";
    repOut << "#Position\t";
    for(int i=0; i < replicaGraphMPIControl().numberOfReplicas; ++i){
      repOut << std::setw(12) << i;
    }
    repOut << "\n";
    repOut << "#s\t\t";
    repOut.precision(svalPrecision);
    int num_l = reedsParam.num_l;
    for (int i = 0; i < replicaGraphMPIControl().numberOfReplicas; ++i){
            repOut << std::setw(12) << reedsParam.eds_para[i].s[0];
    }

    repOut.precision(generalPrecision);
    for (int j = 0; j < reedsParam.num_states; ++j) {
        repOut << "\n# E"<< (j+1)<<"R(s)\t";
        for (int i = 0; i < replicaGraphMPIControl().numberOfReplicas; ++i)
            repOut << std::setw(12) << reedsParam.eds_para[i].eir[j];
    }

    repOut << "\n#\n";
    repOut << ""
          << std::setw(6) << "pos"
          << " "
          << std::setw(6) << "start"
          << " "
          << std::setw(6) << "coord_ID"
          << "   "
          << std::setw(6) << "partner"
          << " "
          << std::setw(6) << "partner_start"
          << " "
          << std::setw(6) << "partner_coord_ID"
          << "  "
          << std::setw(6) << "run"

          //do not use it anymore since the info is in header
          //<< std::setw(13) << "si"
          //do not use temperatures anymore
          //<< std::setw(13) << "Ti"
          << std::setw(18) << "Epoti"
          //do not use it anymore since the info is in header
          //<< std::setw(13) << "sj"
          //do not use temperatures anymore
          //<< std::setw(13) << "Tj"

          << std::setw(18) << "Epotj"
          << std::setw(13) << "p"
          << "  "
          << std::setw(6) << "exch";
    for (int x=1; x<= reedsParam.num_states; x++){
        repOut<<std::setw(17)<<"Vr"<<x;
    }

    repOut << "\n";
    repOut.flush();
    DEBUG(2,"replica_exchange_master_2d_s_eoff_eds "<< globalThreadID <<":init_repOut_stat_file:\t DONE");
}
