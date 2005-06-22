/**
 * @file replica_exchange_control.cc
 * replica exchange control client
 */

#include <stdheader.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>
#include <unistd.h>

#include <io/configuration/out_configuration.h>

#include "replica_exchange.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica

#ifdef XXMPI


////////////////////////////////////////////////////////////////////////////////
// replica control /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Control::Replica_Exchange_Control()
{
}

int util::Replica_Exchange_Control::run
(
 io::Argument & args
 )
{
  char port_name[MPI::MAX_PORT_NAME];

  if (args.count("control") < 1){
    io::messages.add("control: connection name required",
		     "replica exchange",
		     io::message::error);
    MPI_Finalize();
    return 1;
  }
  
  MPI::Lookup_name(args["control"].c_str(), MPI_INFO_NULL, port_name);

  DEBUG(8, "control: connecting..");
  MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
  
  int i = 1, nr;
  MPI_Status status;
  
  DEBUG(8, "control: requesting replica information");
  MPI_Send(&i, 1, MPI_INT, 0, 4, master);
  
  // number of replicas
  MPI_Recv(&nr, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, master, &status);
  
  replica_data.resize(nr);
  
  MPI_Recv(&replica_data[0], sizeof(Replica_Data) * nr,
	   MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
	   master, &status);
  
  MPI_Comm_disconnect(&master);

  if (args.count("control") >= 2){
    std::map<std::string, std::string>::const_iterator
      it = args.lower_bound("control"),
      to = args.upper_bound("control");

    ++it;
    if (it->second == "change" && args.count("control") == 11){
      ++it;
      std::istringstream is(it->second);
      unsigned int nr;
      is >> nr;
      
      --nr;
      if (nr >= replica_data.size()){
	std::cout << "wrong replica ID selected!" << std::endl;
	MPI_Finalize();
	return 1;
      }
      
      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].run;
      }
      
      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].Ti;
	--replica_data[nr].Ti;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].li;
	--replica_data[nr].li;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].epot_i;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].Tj;
	--replica_data[nr].Tj;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].lj;
	--replica_data[nr].lj;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].epot_j;
      }

      is.clear();
      ++it;
      if (it->second == "wait") replica_data[nr].state = waiting;
      if (it->second == "rdy") replica_data[nr].state = ready;
      if (it->second == "run") replica_data[nr].state = running;
      if (it->second == "err") replica_data[nr].state = st_error;
      if (it->second == "term") replica_data[nr].state = terminate;

      DEBUG(8, "control: connecting..");
      MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
  
      int i = 2;

      DEBUG(8, "control: requesting replica change");
      MPI_Send(&i, 1, MPI_INT, 0, 4, master);
  
      DEBUG(8, "control: sending change ID");
      MPI_Send(&nr, 1, MPI_INT, 0, 4, master);
  
      MPI_Send(&replica_data[nr], sizeof(Replica_Data),
	       MPI_CHAR, 0, 4, master);
  
      MPI_Comm_disconnect(&master);
      
    }
    else if (it->second == "stop"){
      std::cout << "send stop request to master" << std::endl;

      DEBUG(8, "control: connecting..");
      MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
  
      int i = 3;

      DEBUG(8, "control: requesting server quit");
      MPI_Send(&i, 1, MPI_INT, 0, 4, master);
      MPI_Comm_disconnect(&master);

    }
    else{
      std::cout << "invalid command\n" << std::endl;
    }
  }
  
  std::cout.precision(4);
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  
  std::cout << "\n"
	    << std::setw(6) << "ID"
	    << std::setw(6) << "run"
	    << std::setw(10) << "T"
	    << std::setw(10) << "l"
	    << std::setw(18) << "Epot"
	    << std::setw(10) << "sT"
	    << std::setw(10) << "sl"
	    << std::setw(18) << "sEpot"
	    << std::setw(10) << "p"
	    << std::setw(5) << "s"
	    << std::setw(10) << "ste"
	    << "\n";
  
  for(int r=0; r<nr; ++r){
    std::cout << std::setw(6) << r + 1
	      << std::setw(6) << replica_data[r].run
	      << std::setw(10) << replica_data[r].Ti + 1
	      << std::setw(10) << replica_data[r].li + 1
	      << std::setw(18) << replica_data[r].epot_i
	      << std::setw(10) << replica_data[r].Tj + 1
	      << std::setw(10) << replica_data[r].lj + 1
	      << std::setw(18) << replica_data[r].epot_j
	      << std::setw(10) << replica_data[r].probability
	      << std::setw(5) << replica_data[r].switched;

    switch(replica_data[r].state){
      case waiting:
	std::cout << std::setw(10) << "wait";
	break;
      case ready:
	std::cout << std::setw(10) << "rdy";
	break;
      case running:
	std::cout << std::setw(10) << "run";
	break;
      case st_error:
	std::cout << std::setw(10) << "err";
	break;
      case terminate:
	std::cout << std::setw(10) << "term";
	break;
      default:
	std::cout << std::setw(10) << "???";
    }
    std::cout << std::endl;
  }
  
  std::cout << "\nuse\n\t@control change ID run T l Epot sT sl sEpot state" 
	    << "\n\tto change replicas or stop to shutdown master\n" << std::endl;
  
  MPI_Finalize();
  return 0;
}

#endif
