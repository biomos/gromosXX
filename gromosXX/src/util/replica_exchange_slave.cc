/**
 * @file replica_exchange_slave.cc
 * replica exchange
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
// replica slave    ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Slave::Replica_Exchange_Slave()
{
}

int util::Replica_Exchange_Slave::run
(
 io::Argument & args)
{
  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  // read the files (could also be copied from the master)
  io::read_input(args, topo, conf, sim,  md, std::cout);

  // initialises everything
  md.init(topo, conf, sim, std::cout);

  // aliases
  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;

  io::Out_Configuration traj("GromosXX\n", std::cout);
  traj.title("GromosXX\n" + sim.param().title);
  traj.init(args["fin"], args["trj"], args["trv"], args["trf"], 
	    args["tre"], args["trg"], args["bae"], args["bag"],
	    sim.param());

  // should be the same as from master...
  std::cout << "\nMESSAGES FROM (SLAVE) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }
  io::messages.clear();
  
  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

  std::cout << "slave initialised" << std::endl;
  
  std::string server_name;
  int server_port;

  if (args.count("slave") > 0){
    io::Argument::const_iterator iter=args.lower_bound("slave");
    if(iter!=args.upper_bound("slave")){

      std::cerr << "setting server name to: " << iter->second << std::endl;
      server_name = iter->second;
      ++iter;
      if(iter!=args.upper_bound("slave")){

	std::istringstream is(iter->second);
	if (!(is >> server_port)){
	  io::messages.add("control [server [port number]]",
			   "replica_exchange",
			   io::message::error);
	  return 1;
	}
      }
    }
  }

  struct hostent *hostinfo;
  int error;
  hostinfo = getipnodebyname(server_name.c_str(), AF_INET, AI_DEFAULT, &error);
  if (hostinfo == NULL){
    io::messages.add("could not get hostinfo on server",
		     "replica_exchange",
		     io::message::error);
    std::cerr << "could not get hostinfo: error = " << error << std::endl;
    return 1;
  }
  if (hostinfo->h_addrtype != AF_INET){
    io::messages.add("host is not an IP host",
		     "replica_exchange",
		     io::message::error);
    return 1;
  }
  
  cl_socket = socket(AF_INET, SOCK_STREAM, 0);
  struct sockaddr_in address;

  address.sin_family = AF_INET;
  address.sin_port = htons(server_port);
  address.sin_addr = *(struct in_addr *) *hostinfo->h_addr_list;
  socklen_t len = sizeof(address);

  freehostent(hostinfo);

  for(int run=0; run < sim.param().replica.slave_runs; ++run){
    
    DEBUG(8, "slave: connecting..");
    int result = connect(cl_socket, (sockaddr *) &address, len);

    if (result == -1){
      std::cout << "could not connect to master. master finished?"
		<< std::endl;
      return 1;
    }

    DEBUG(9, "slave: connected");

    // request a job
    int server_response = get_replica_data();

    if (server_response == 0){
      // accept job
      std::cout << "slave: got a job! (replica "
		<< replica_data.ID << ")" << std::endl;
      
      std::cout << "slave:  running replica " << replica_data.ID
		<< " at T=" << T[replica_data.Ti]
		<< " and l=" << l[replica_data.li]
		<< " (run = " << replica_data.run << ")"
		<< std::endl;

      // init replica parameters (t, conf, T, lambda)
      if (init_replica(topo, conf, sim)){
	std::cerr << "slave: disconnecting..." << std::endl;
	close(cl_socket);
	return 1;
      }

      // close connection
      close(cl_socket);

      // run it
      int error = run_md(topo, conf, sim, md, traj);

      if (!error){
	// do we need to reevaluate the potential energy ?
	// yes! 'cause otherwise it's just the energy for
	// the previous configuration...
	algorithm::Algorithm * ff = md.algorithm("Forcefield");
	
	if (ff == NULL){
	  std::cerr << "forcefield not found in MD algorithm sequence"
		    << std::endl;
	  return 1;
	}
	
	ff->apply(topo, conf, sim);
	conf.current().energies.calculate_totals();
	replica_data.epot_i = conf.current().energies.potential_total;
	
	if (replica_data.Ti != replica_data.Tj)
	  std::cout << "replica_energy " << replica_data.epot_i 
		    << " @ " << T[replica_data.Ti] << "K" << std::endl;
	else if (replica_data.li != replica_data.lj)
	  std::cout << "replica_energy " << replica_data.epot_i 
		    << " @ " << l[replica_data.li] << std::endl;
	
	if (replica_data.li != replica_data.lj){
	  
	  // change the lambda value
	  sim.param().perturbation.lambda = l[replica_data.lj];
	  topo.lambda(l[replica_data.lj]);
	  // twice: set old lambda as well
	  topo.lambda(l[replica_data.lj]);
	  topo.update_for_lambda();
	  std::cout << "\tlambda = " << topo.lambda() << "\n";
	  
	  // recalc energy
	  ff->apply(topo, conf, sim);
	  conf.current().energies.calculate_totals();
	  replica_data.epot_j = conf.current().energies.potential_total;
	} 
	else{
	  replica_data.epot_j = replica_data.epot_i;
	}
	
	++replica_data.run;
	replica_data.state = waiting;
      }
      else{
	replica_data.state = st_error;
      }
	
      DEBUG(8, "slave: connecting (after run)...");
      // std::cerr << "salve: connecting..." << std::endl;
      cl_socket = socket(AF_INET, SOCK_STREAM, 0);
      result = connect(cl_socket, (sockaddr *) &address, len);
      if (result == -1){
	std::cout << "could not (re-)connect to master. master finished?"
		  << std::endl;
	return 1;
      }

      DEBUG(9, "slave: connected");
      DEBUG(8, "slave: finished job " << replica_data.ID);

      char ch = 2;
      write(cl_socket, &ch, 1);
      read(cl_socket, &ch, 1);
      if (ch != 0){
	io::messages.add("server reported error",
			 "replica_exchange",
			 io::message::error);
	close(cl_socket);
	return 1;
      }
      
      update_replica_data();
      update_configuration(topo, conf);
      
      // and disconnect
      close(cl_socket);
    }
    else if (server_response == 9){
      std::cout << "server has finished!\n"
		<< "exiting...\n"
		<< std::endl;

      std::cerr << "slave: disconnecting..." << std::endl;
      close(cl_socket);

      return 0;
    }
    else{
      close(cl_socket);
      std::cout << "slave waiting..." << std::endl;
      // wait apropriate time for master to prepare
      sleep(timeout);
    }
  } // for slave_runs
  
  std::cerr << "slave: finished runs. terminating..." << std::endl;
  return 0;
}

int util::Replica_Exchange_Slave::run_md
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 algorithm::Algorithm_Sequence & md,
 io::Out_Configuration & traj
 )
{
  // do we scale the initial temperatures?
  if (sim.param().replica.scale){
    algorithm::Temperature_Calculation tcalc;
    tcalc.apply(topo, conf, sim);
    
    algorithm::Berendsen_Thermostat tcoup;
    tcoup.calc_scaling(topo, conf, sim, true);
    tcoup.scale(topo, conf, sim);
  }

  double end_time = sim.time() + 
    sim.time_step_size() * (sim.param().step.number_of_steps - 1);
    
  int error;

  while(sim.time() < end_time + math::epsilon){
      
    traj.write_replica_step(sim, replica_data);
    traj.write(conf, topo, sim, io::reduced);

    // run a step
    if ((error = md.run(topo, conf, sim))){
      
      std::cout << "\nError during MD run!\n" << std::endl;
      break;
    }

    traj.print(topo, conf, sim);

    sim.time() += sim.time_step_size();
    ++sim.steps();

  }
    
  return error;
}

int util::Replica_Exchange_Slave::get_replica_data()
{
  DEBUG(8, "slave: requesting job");
  char ch = 1;
  write(cl_socket, &ch, 1);

  read(cl_socket, &ch, 1);
  if (ch != 0){
    std::cout << "no job received from server (" << ch << ")" << std::endl;
  }
  else{
    
    DEBUG(8, "slave: waiting for replica data");
    read(cl_socket, (char *) &replica_data, sizeof(Replica_Data));

    DEBUG(9, "slave: got replica " << replica_data.ID
	  << " temperature=" << replica_data.Ti
	  << " lambda=" << replica_data.li);
  }

  return int(ch);
}

int util::Replica_Exchange_Slave::update_replica_data()
{
  DEBUG(8, "slave: updating replica data");

  write(cl_socket, (char *) &replica_data.ID, sizeof(int));
  write(cl_socket, (char *) &replica_data, sizeof(Replica_Data));

  return 0;
}

int util::Replica_Exchange_Slave::get_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  DEBUG(10, "receiving " << 3 * conf.current().pos.size() << " coords");
  int error = 0;

  read(cl_socket, (char *) &conf.current().pos(0)(0),
       conf.current().pos.size() * 3 * sizeof(double));

  read(cl_socket, (char *) &conf.current().vel(0)(0),
       conf.current().vel.size() * 3 * sizeof(double));

  read(cl_socket, (char *) &conf.current().box(0)(0),
       9 * sizeof(double));

  return error;
}

int util::Replica_Exchange_Slave::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  // positions
  DEBUG(9, "sending " << 3 * conf.current().pos.size() << " coords");

  write(cl_socket, (char *) &conf.current().pos(0)(0),
	conf.current().pos.size() * 3 * sizeof(double));

  write(cl_socket, (char *) &conf.current().vel(0)(0),
	conf.current().vel.size() * 3 * sizeof(double));

  write(cl_socket, (char *) &conf.current().box(0)(0),
	9 * sizeof(double));

  return 0;
}

int util::Replica_Exchange_Slave::init_replica
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim
 )
{
  // get configuration from master
  if (get_configuration(topo, conf))
    return 1;
  
  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;

  // change all the temperature coupling temperatures
  for(unsigned int i=0; i<sim.multibath().size(); ++i){
    sim.multibath()[i].temperature = T[replica_data.Ti];
  }
  
  // change the lambda value
  sim.param().perturbation.lambda = l[replica_data.li];
  topo.lambda(l[replica_data.li]);
  // twice, to set old_lambda...
  topo.lambda(l[replica_data.li]);
  topo.update_for_lambda();
  std::cout << "\tlambda = " << topo.lambda() << "\n";
  
  // change simulation time
  sim.time() = replica_data.run * 
    sim.param().step.number_of_steps * sim.param().step.dt +
    sim.param().step.t0;
  
  sim.steps() = replica_data.run * sim.param().step.number_of_steps;

  return 0;
}
