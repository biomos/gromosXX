/**
 * @file replica_exchange.cc
 * replica exchange common routines
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
#include <interaction/special/external_interaction.h>

#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>
#include <util/virtual_grain.h>

#include <io/read_input.h>
#include <io/print_block.h>

#include <time.h>

#include <io/configuration/out_configuration.h>

#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
// linux includes?
#include <sys/types.h>
#include <sys/socket.h>
// end linux includes
#include <netdb.h>

#include "replica_exchange.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica



////////////////////////////////////////////////////////////////////////////////
// get / put replica data
////////////////////////////////////////////////////////////////////////////////
int util::Replica_Exchange::put_replica_data(Replica_Data & r)
{
  if (write(cl_socket, (char *) &r.ID, sizeof(int)) != sizeof(int)){
    std::cerr << "could not write to socket" << std::endl;
    return 1;
  }
  if (write(cl_socket, (char *) &r, sizeof(Replica_Data))
      != sizeof(Replica_Data)){

    std::cerr << "could not write Replica_Data to socket" << std::endl;
    return 1;
  }   
  return 0;
}

int util::Replica_Exchange::get_replica_data(Replica_Data & r, int & i)
{
  if (read(cl_socket, (char *) &i, sizeof(int)) != sizeof(int)){
    std::cerr << "could not read from socket" << std::endl;
    i = -1;
    throw std::runtime_error("could not read replica index");
  }

  // it's a matter of trust...
  if (i < 0){
    std::cerr << "received invalid replica index! (" << i << ")" <<  std::endl;
    i = -1;
    throw std::runtime_error("could not read replica index");
  }
  
  if (read(cl_socket, (char *) &r, sizeof(Replica_Data)) != sizeof(Replica_Data)){
    std::cerr << "could not read replica data " << i << std::endl;
    throw std::runtime_error("could not read replica index");
  }
  return 0;
}

int util::Replica_Exchange::get_replica_data(std::vector<Replica_Data> & r, int & i)
{
  if (read(cl_socket, (char *) &i, sizeof(int)) != sizeof(int)){
    std::cerr << "could not read from socket" << std::endl;
    i = -1;
    throw std::runtime_error("could not read replica index");
  }

  if (i < 0 || unsigned(i) >= r.size()){
    std::cerr << "received invalid replica index! (" << i << ")" <<  std::endl;
    i = -1;
    throw std::runtime_error("could not read replica data");
  }

  if (read(cl_socket, (char *) &r[i], sizeof(Replica_Data)) != sizeof(Replica_Data)){
    std::cerr << "could not read replica data " << i  << std::endl;
    throw std::runtime_error("could not read replica data");
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// get / put configuration
////////////////////////////////////////////////////////////////////////////////

int util::Replica_Exchange::get_configuration
(
 configuration::Configuration & conf
 )
{
  const ssize_t num = 3 * conf.current().pos.size() * sizeof(double);
  DEBUG(10, "receiving " << num / sizeof(double) << " coords");

  if (num >= SSIZE_MAX){
    std::cerr << "chunk size not large enough to exchange configuration" << std::endl;
    close(cl_socket);
    return 1;
  }

  try{
    readblock((char *) &conf.current().pos(0)(0), num);
    readblock((char *) &conf.current().vel(0)(0), num);   

    if (read(cl_socket, (char *) &conf.current().box(0)(0),
		      9 * sizeof(double))
	!= 9 * sizeof(double))

      throw std::runtime_error("could not read box");      
  }
  catch (std::runtime_error e){
    std::cerr << "Exception: " << e.what() << std::endl;
    std::cout << "Exception: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}

int util::Replica_Exchange::put_configuration
(
 configuration::Configuration & conf
 )
{
  // positions
  DEBUG(9, "sending " << 3 * conf.current().pos.size() << " coords");
  const ssize_t num = 3 * conf.current().pos.size() * sizeof(double);
  
  if (num >= SSIZE_MAX){
    std::cerr << "chunk size not large enough to exchange configuration" << std::endl;
    return 1;
  }

  try{
    writeblock((char *) &conf.current().pos(0)(0), num);
    writeblock((char *) &conf.current().vel(0)(0), num);  
    if (write(cl_socket, (char *) &conf.current().box(0)(0),
	      9 * sizeof(double)) != 9 * sizeof(double))
      throw std::runtime_error("could not write box");
  }
  catch(std::runtime_error e){
    std::cout << "Exception: " << e.what() << std::endl;
    std::cerr << "Exception: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// socket read / write
////////////////////////////////////////////////////////////////////////////////

ssize_t util::Replica_Exchange::readblock(char * source, ssize_t size)
{
  ssize_t current;
  ssize_t window = 1024;
  while(size > 0){
    if (size > window)
      current = window;
    else current = size;
      
    ssize_t count;
    if ((count = read(cl_socket, source, current)) == 0){
      std::cerr << "received zero bytes instead of "
		<< current << " !!!" << std::endl;
      throw std::runtime_error("could not read data block");
    }
      
    if (current != count){
      std::cerr << "received only " << count << " bytes..." << std::endl;
      // make window smaller...
      window = count;
    }

    char c = 0;
    if (write(cl_socket, &c, 1) != 1){
      std::cerr << "sending ACK failed" << std::endl;
      throw std::runtime_error("could not send ACK");
    }

    source += count;
    size -= count;
  }
  return 0;
}
  
ssize_t util::Replica_Exchange::writeblock(char * dest, ssize_t size)
{
  ssize_t current;
  ssize_t window = 1024;
  while(size > 0){
    if (size > window)
      current = window;
    else current = size;
      
    ssize_t count;
    if ((count = write(cl_socket, dest, current)) == 0){
      std::cerr << "could not write a single byte!\n"
		<< "tried to send " << current << std::endl;
      throw std::runtime_error("could not write data block");
    }

    if (current != count){
      std::cerr << "sent only " << count << " bytes..." << std::endl;
      // make window smaller...
      window = count;
    }
      
    char c = 0;
    if (read(cl_socket, &c, 1) != 1){
      std::cerr << "getting ACK failed" << std::endl;
      throw std::runtime_error("could not read ACK");
    }
    if (c != 0){
      std::cerr << "wrong ACK received" << std::endl;
      throw std::runtime_error("wrong ACK received");
    }
      
    dest += count;
    size -= count;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// magic cookie
////////////////////////////////////////////////////////////////////////////////

bool util::Replica_Exchange::magic_cookie(bool master)
{
  double magic[4] = { 3.1415927, 29375, 243, 8.3116 };
  double magic_buff[4];
    
  // magic cookie exchange

  if (master){
    if (write(cl_socket, &magic, 4 * sizeof(double)) != 4 * sizeof(double)){
      std::cerr << "could not write magic cookie" << std::endl;
      return false;
    }
    if (read(cl_socket, magic_buff, 4 * sizeof(double)) != 4 * sizeof(double)){
      std::cerr << "could not read magic cookie" << std::endl;
      return false;
    }
  }
  else{
    if (read(cl_socket, magic_buff, 4 * sizeof(double)) != 4 * sizeof(double)){
      std::cerr << "could not read magic cookie" << std::endl;
      return false;
    }
    if (write(cl_socket, &magic, 4 * sizeof(double)) != 4 * sizeof(double)){
      std::cerr << "could not write magic cookie" << std::endl;
      return false;
    }
  }
  
  if (magic[0] != magic_buff[0] || magic[1] != magic_buff[1] ||
      magic[2] != magic_buff[2] || magic[3] != magic_buff[3]){
    
    std::cerr << "magic cookie exchange failed" << std::endl;
    return false;
  }
  else{
    std::cout << "magic cookie test succeeded!" << std::endl;
    return true;
  }
}

addrinfo * util::Replica_Exchange::get_server(io::Argument & args, addrinfo &hints, std::string server_name)
{
  {
    char buffer[256];
    gethostname(buffer, 255);
    server_name = buffer;
  }

  int server_port = 29375;
  
  if (args.count("slave") > 0){
    io::Argument::const_iterator iter=args.lower_bound("slave");
    if(iter!=args.upper_bound("slave")){

      server_name = iter->second;
      ++iter;
      if(iter!=args.upper_bound("slave")){
	
	std::istringstream is(iter->second);
	if (!(is >> server_port)){
	  io::messages.add("control [server [port number]]",
			   "replica_exchange",
			   io::message::error);

	  throw std::runtime_error("could not get server port");
	}
      }
    }
  }

  std::cout << "\trunning on host " << server_name << " : " << server_port << std::endl;

  struct addrinfo *addrinfo_p;
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
    io::messages.add("could not get server address info",
		     "replica_exchange",
		     io::message::error);
    std::cerr << "getaddrinfo error!\n"
	      << gai_strerror(error)
	      << std::endl;

    throw std::runtime_error("could not get server address info");
    
  }

  ((sockaddr_in *)addrinfo_p->ai_addr)->sin_port = htons(server_port);
  

  return addrinfo_p;
  
}
