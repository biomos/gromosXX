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

#include <io/configuration/out_configuration.h>

#include "replica_exchange.h"

#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica

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
  std::string server_name;
  {
    char buffer[256];
    gethostname(buffer, 255);
    server_name = buffer;
  }
  std::cerr << "running on host " << server_name << std::endl;
  
  int server_port = 29375;
  
  bool cmd_change = false;
  bool cmd_quit = false;

  if (args.count("control") > 0){
    io::Argument::const_iterator iter=args.lower_bound("control");
    if(iter!=args.upper_bound("control")){

      if (iter->second == "change")
	cmd_change = true;
      else if (iter->second == "quit")
	cmd_quit = true;
      else{

	std::cerr << "setting server name to: " << iter->second << std::endl;
	server_name = iter->second;
	++iter;
	if(iter!=args.upper_bound("control")){

	  if (iter->second == "change")
	    cmd_change = true;
	  else if (iter->second == "quit")
	    cmd_quit = true;
	  else{
	    std::istringstream is(iter->second);
	    if (!(is >> server_port)){
	      io::messages.add("control [server [port number]]",
			       "replica_exchange",
			       io::message::error);
	      return 1;
	    }

	    ++iter;
	    if(iter!=args.upper_bound("control")){
	      
	      if (iter->second == "change")
		cmd_change = true;
	      else if (iter->second == "quit")
		cmd_quit = true;
	    }
	  }
	}
      }
    }
  }

  struct addrinfo *addrinfo_p;
  struct addrinfo hints;
  hints.ai_flags = 0;
  hints.ai_family = PF_INET;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_protocol = 0;
  hints.ai_addrlen = 0;
  hints.ai_addr = NULL;
  hints.ai_canonname = NULL;
  hints.ai_next = NULL;

  int error = getaddrinfo(server_name.c_str(), NULL,
			  &hints, &addrinfo_p);
  
  if (error){
    std::cerr << "getaddrinfo error!\n"
	      << gai_strerror(error)
	      << std::endl;
    return 1;
  }

  ((sockaddr_in *)addrinfo_p->ai_addr)->sin_port = htons(server_port);

  sockaddr * s_addr_p = addrinfo_p->ai_addr;
  int len = addrinfo_p->ai_addrlen;

  cl_socket = socket(addrinfo_p->ai_family, addrinfo_p->ai_socktype,
		     addrinfo_p->ai_protocol);

  DEBUG(8, "control: connecting..");

  int result;
  
  result = connect(cl_socket, s_addr_p, len);
  freeaddrinfo(addrinfo_p);
  
  if (result == -1){
    io::messages.add("could not connect to server",
		     "replica_exchange",
		     io::message::error);
    close(cl_socket);
    return 1;
  }

  std::cerr << "connected!" << std::endl;

  double magic[4] = { 3.1415927, 29375, 243, 8.3116 };
  double magic_buff[4];
  
  // magic cookie exchange
  if (read(cl_socket, magic_buff, 4 * sizeof(double)) != 4 * sizeof(double)){
    std::cerr << "could not read magic cookie" << std::endl;
    close(cl_socket);
    return 1;
  }
  if (write(cl_socket, &magic, 4 * sizeof(double)) != 4 * sizeof(double)){
    std::cerr << "could not write magic cookie" << std::endl;
    close(cl_socket);
    return 1;
  }
  if (magic[0] != magic_buff[0] || magic[1] != magic_buff[1] ||
      magic[2] != magic_buff[2] || magic[3] != magic_buff[3]){
    
    std::cerr << "magic cookie exchange failed" << std::endl;
    close(cl_socket);
    return 1;
  }
  
  // interactive
  char ch = 4;
  write(cl_socket, &ch, 1);
  read(cl_socket, &ch, 1);
  
  if(ch != 0){
    io::messages.add("master reported error",
		     "replica_exchange",
		     io::message::error);
    close(cl_socket);
    return 1;
  }
  
  if (cmd_quit){
    ch = 3;
    write(cl_socket, &ch, 1);
    close(cl_socket);
    return 0;
  }

  // ask for replicas
  std::cerr << "requesting replicas" << std::endl;
  if (cmd_change) ch = 2;
  else ch = 1;
  
  write(cl_socket, &ch, 1);
  
  DEBUG(8, "control: requesting replica information");
  int nr = 0;
  read(cl_socket, (char *) &nr, sizeof(int));
  std::cerr << nr << " replicas." << std::endl;

  replica_data.resize(nr);

  read(cl_socket,  (char *) &replica_data[0],
       nr * sizeof(Replica_Data));

  //////////////////////////////////////////////////////////////////////
  // change a replica
  //////////////////////////////////////////////////////////////////////
  if (cmd_change){
    std::map<std::string, std::string>::const_iterator
      it = args.lower_bound("control"),
      to = args.upper_bound("control");

    //////////////////////////////////////////////////
    // parse change arguments
    //////////////////////////////////////////////////
    for( ; it != to && it->second != "change"; ++it){
    }

    nr = -1;
    while (it != to){
      ++it;
      if (it == to) break;

      std::istringstream is(it->second);
      is >> nr;
      
      --nr;
      if (nr >= int(replica_data.size())){
	std::cout << "wrong replica ID selected!" << std::endl;
	return 1;
      }
      
      std::cerr << "trying to change replica " << nr << std::endl;

      is.clear();
      ++it;
      if (it == to) break;
	  
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].run;
      }
      
      is.clear();
      ++it;
      
      if (it == to) break;
      
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].Ti;
	--replica_data[nr].Ti;
      }
      
      is.clear();
      ++it;
      
      if (it == to) break;

      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].li;
	--replica_data[nr].li;
      }
      
      is.clear();
      ++it;
      
      if (it == to) break;

      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].epot_i;
      }
      
      is.clear();
      ++it;
      
      if (it == to) break;
      
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].Tj;
	--replica_data[nr].Tj;
      }
      
      is.clear();
      ++it;
      
      if (it == to) break;
      
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].lj;
	--replica_data[nr].lj;
      }
      
      is.clear();
      ++it;
      if (it == to) break;
		      
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].epot_j;
      }
      
      is.clear();
      ++it;
      if (it == to) break;
			
      if (it->second == "wait") replica_data[nr].state = waiting;
      if (it->second == "rdy") replica_data[nr].state = ready;
      if (it->second == "run") replica_data[nr].state = running;
      if (it->second == "err") replica_data[nr].state = st_error;
      if (it->second == "term") replica_data[nr].state = terminate;
    }
    //////////////////////////////////////////////////

    std::cerr << "writing to master! (nr = " << nr << ")" << std::endl;

    if (nr >= 0 && nr < int(replica_data.size())){

      DEBUG(8, "control: sending change ID");
      if (write(cl_socket, (char *) &nr, sizeof(int)) != sizeof(int)){
	std::cerr << "could not write ID" << std::endl;
      }
      
      std::cerr << "and now replica data" << std::endl;
      
      if (write(cl_socket, (char *) &replica_data[nr], sizeof(Replica_Data)) != sizeof(Replica_Data)){
	std::cerr << "could not write replica data" << std::endl;
      }
    }
  }
  
  std::cerr << "control: disconnecting..." << std::endl;
  close(cl_socket);
 
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
  
  for(unsigned int r=0; r<replica_data.size(); ++r){
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
	    << "\n\tto change replicas or"
	    << "\n\t@control stop to shutdown master\n" << std::endl;
  
  return 0;
}
