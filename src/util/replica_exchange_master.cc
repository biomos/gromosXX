/* 
 * File:   replica_exchange_master.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:18 PM
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
#include <util/replica.h>
#include <util/replica_exchange_master.h>
#include <util/replica_exchange_base.h>
#include <string>

#ifdef XXMPI
#include <mpi.h>
#endif

util::replica_exchange_master::replica_exchange_master(io::Argument & args,
        int cont,
        int rank,
        int _size,
        int _numReplicas,
        std::vector<int> repIDs,
        std::map<ID_t, rank_t> & repMap)
:
replica_exchange_base(args, cont, rank, repIDs, repMap),
repParams(replicas[0]->sim.param().replica),
size(_size),
numReplicas(_numReplicas) {
  assert(rank == 0);
  assert(numReplicas > 0);
  replicaData.resize(numReplicas);

  //initialize data of replicas
  int ID = 0;
  for (int i = 0; i < repParams.num_l; ++i) {
    for (int j = 0; j < repParams.num_T; ++j) {
      replicaData[ID].ID = ID;
      replicaData[ID].T = repParams.temperature[j];
      replicaData[ID].l = repParams.lambda[i];
      replicaData[ID].dt = repParams.dt[i];
      ++ID;
    }
  }

  // set output file
  std::string repdatName = args["repdat"];
  repOut.open(repdatName.c_str());
  repOut << "Number of temperatures:\t" << repParams.num_T << "\n"
          << "Number of lambda values:\t" << repParams.num_l << "\n";

  repOut.precision(4);
  repOut.setf(std::ios::fixed, std::ios::floatfield);

  repOut << "T    \t";
  for (int t = 0; t < repParams.num_T; ++t)
    repOut << std::setw(12) << repParams.temperature[t];

  repOut << "\nlambda    \t";
  for (int l = 0; l < repParams.num_l; ++l)
    repOut << std::setw(12) << repParams.lambda[l];

  repOut << "\n\n";

  repOut << "#"
          << std::setw(5) << "ID"
          << " "
          << std::setw(7) << "partner"
          << std::setw(7) << "run"

          << std::setw(13) << "li"
          << std::setw(13) << "Ti"
          << std::setw(14) << "Epoti"
          << std::setw(13) << "lj"
          << std::setw(13) << "Tj"
          << std::setw(15) << "Epotj"
          << std::setw(15) << "p"
          << std::setw(8) << "s"
          << "\n";
}

util::replica_exchange_master::~replica_exchange_master() {
  repOut.close();
}

void util::replica_exchange_master::receive_from_all_slaves() {
#ifdef XXMPI
  double start = MPI_Wtime();

  MPI_Status status;
  util::repInfo info;

  // receive all information from slaves
  for (unsigned int rep = 0; rep < numReplicas; ++rep) {
    unsigned int rank = repMap.find(rep)->second;
    if (rank != 0) {
      MPI_Recv(&info, 1, MPI_REPINFO, rank, REPINFO, MPI_COMM_WORLD, &status);
      replicaData[rep].run = info.run;
      replicaData[rep].epot = info.epot;
      replicaData[rep].epot_partner = info.epot_partner;
      replicaData[rep].probability = info.probability;
      replicaData[rep].switched = info.switched;
      replicaData[rep].partner = info.partner;
    }
  }

  // write all information from master node to data structure
  for (repIterator it = replicas.begin(); it < replicas.end(); ++it) {
    int ID = (*it)->ID;
    replicaData[ID].run = (*it)->run;
    replicaData[ID].partner = (*it)->partner;
    replicaData[ID].epot = (*it)->epot;
    replicaData[ID].epot_partner = (*it)->epot_partner;
    replicaData[ID].probability = (*it)->probability;
    replicaData[ID].switched = (*it)->switched;
  }

  std::cout << "Master:\n" << "time used for receiving all messages: " << MPI_Wtime() - start 
            << " seconds" << std::endl;
#endif
}

void util::replica_exchange_master::write() {
  for (unsigned int r = 0; r < numReplicas; ++r) {

    repOut << std::setw(6) << (replicaData[r].ID + 1)
            << " "
            << std::setw(6) << (replicaData[r].partner + 1)
            << std::setw(6) << replicaData[r].run
            << std::setw(13) << replicaData[r].l
            << std::setw(13) << replicaData[r].T
            << " "
            << std::setw(18) << replicaData[r].epot
            << std::setw(13) << replicaData[replicaData[r].partner].l
            << std::setw(13) << replicaData[replicaData[r].partner].T
            << " ";
    // the following is a little bit clumsy. Because for temperature REMD, we did not have to recalculate
    // any energies using different potentials, we could not properly store the epot_partner value (only
    // for the first of the switching pair. Therefor we just print the epot of the partner (which is exactly
    // what epot_partner would have been if we had taken the effort to send it back to the partner.
    if(replicaData[r].l == replicaData[replicaData[r].partner].l)
	repOut << std::setw(18) << replicaData[replicaData[r].partner].epot;
    else
        repOut << std::setw(18) << replicaData[r].epot_partner;
    repOut  << std::setw(13) << replicaData[r].probability
            << std::setw(4) << replicaData[r].switched
            << std::endl;
  }
}
