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
// replica slave    ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Slave::Replica_Exchange_Slave()
{
}

int util::Replica_Exchange_Slave::run
(
 io::Argument & args)
{
  multigraining = false;
  if (args.count("cg_topo") >= 0){
    multigraining = true;
  }

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  // and the coarse-grained stuff
  topology::Topology cg_topo;
  configuration::Configuration cg_conf;
  algorithm::Algorithm_Sequence cg_md;
  simulation::Simulation cg_sim;

  // read the files (could also be copied from the master)

  // add an external interaction
  if (multigraining)
    sim.param().force.external_interaction = 1;

  std::cout << "reading (standard) input" << std::endl;
  io::read_input(args, topo, conf, sim,  md, std::cout);
  md.init(topo, conf, sim, std::cout);

  interaction::Forcefield * cg_ff;
  if (multigraining){
    interaction::Forcefield * ff = 
      dynamic_cast<interaction::Forcefield *>(md.algorithm("Forcefield"));
    if (ff == NULL){
      std::cout << "Error: no forcefield in MD" << std::endl;
      return 1;
    }
    interaction::External_Interaction * ei = 
      dynamic_cast<interaction::External_Interaction *>(ff->interaction("External"));
    if (ei == NULL){
      std::cout << "Error: no external interaction in forcefield" << std::endl;
      return 1;
    }
    
    ei->set_coarsegraining(cg_topo, cg_conf);
    
    io::argname_conf = "cg_conf";
    io::argname_topo = "cg_topo";
    io::argname_pttopo = "cg_pttopo";
    io::argname_input = "cg_input";
    
    io::read_input(args, cg_topo, cg_conf, cg_sim, cg_md);
    cg_ff = dynamic_cast<interaction::Forcefield *>(cg_md.algorithm("Forcefield"));
    if (cg_ff == NULL){
      std::cout << "Error: no forcefield in cg_MD" << std::endl;
      return 1;
    }

    cg_md.init(cg_topo, cg_conf, cg_sim);
  }

  // aliases
  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;
  std::vector<double> const & dt = sim.param().replica.dt;

  io::Out_Configuration traj("GromosXX\n", std::cout);
  traj.title("GromosXX\n" + sim.param().title);

  traj.init(args, sim.param());
  /*
  traj.init(args["fin"], args["trj"], args["trv"], args["trf"], 
	    args["tre"], args["trg"], args["bae"], args["bag"],
	    sim.param());
  */

  // should be the same as from master... (more or less, multigraining)
  std::cout << "\nMESSAGES FROM (SLAVE) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }
  io::messages.clear();
  
  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

  std::cout << "\n\nslave initialised\n" << std::endl;
  if (multigraining){
    std::cout <<"\tmultigrained simulation" << std::endl;
  }
  
  std::string server_name;
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
	  return 1;
	}
      }
    }
  }

  std::cout << "\trunning on host " << server_name << " : " << server_port << std::endl;

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
    io::messages.add("could not get server address info",
		     "replica_exchange",
		     io::message::error);
    std::cerr << "getaddrinfo error!\n"
	      << gai_strerror(error)
	      << std::endl;
    return 1;
  }

  ((sockaddr_in *)addrinfo_p->ai_addr)->sin_port = htons(server_port);
  sockaddr * s_addr_p = addrinfo_p->ai_addr;
  int len = addrinfo_p->ai_addrlen;

  for(int run=0; run < sim.param().replica.slave_runs; ++run){

    cl_socket = socket(addrinfo_p->ai_family, addrinfo_p->ai_socktype,
		       addrinfo_p->ai_protocol);
    if (cl_socket < 0){
      std::cerr << "could not create client socket" << std::endl;
      return 1;
    }
    
    DEBUG(8, "slave: connecting..");
    int result = connect(cl_socket, s_addr_p, len);

    if (result == -1){
      std::cout << "could not connect to master. master finished?"
		<< std::endl;
      return 1;
    }

    double magic[4] = { 3.1415927, 29375, 243, 8.3116 };
    double magic_buff[4];
    
    // magic cookie exchange
    if (read(cl_socket, magic_buff, 4 * sizeof(double)) != 4 * sizeof(double)){
      std::cerr << "could not read magic cookie" << std::endl;
      close(cl_socket);
      continue;
    }
    if (write(cl_socket, &magic, 4 * sizeof(double)) != 4 * sizeof(double)){
      std::cerr << "could not write magic cookie" << std::endl;
      close(cl_socket);
      continue;
    }
    if (magic[0] != magic_buff[0] || magic[1] != magic_buff[1] ||
	magic[2] != magic_buff[2] || magic[3] != magic_buff[3]){

      std::cerr << "magic cookie exchange failed" << std::endl;
      close(cl_socket);
      continue;
    }
    else{
      std::cout << "magic cookie test succeeded!" << std::endl;
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
		<< " (run = " << replica_data.run 
		<< " and dt = " << dt[replica_data.li] << ")"
		<< std::endl;

      // init replica parameters (t, conf, T, lambda)
      if (init_replica(topo, conf, sim, cg_topo, cg_sim)){
	std::cerr << "init_replica returned error!" << std::endl;
	std::cerr << "slave: disconnecting..." << std::endl;
	close(cl_socket);
	return 1;
      }

      // close connection
      close(cl_socket);

      // run it
      std::cerr << "running md" << std::endl;
      int error = run_md(topo, conf, sim, md, cg_topo, cg_conf, cg_sim, cg_ff, traj);
      std::cerr << "run finished" << std::endl;
      
      std::cout << "\nslave: run finished!\n\n";
      std::cout << "\tcalculating potential energies of last configuration\n\n";
      
      if (!error){
	// do we need to reevaluate the potential energy ?
	// yes! 'cause otherwise it's just the energy for
	// the previous configuration...

	if (multigraining){
	  // coarse grained atom positions are based upon
	  // real atom positions

	  // std::cerr << "update virtual pos" << std::endl;
	  util::update_virtual_pos(cg_topo, cg_conf, topo, conf);
	  
	  // calculate the cg forces first!
	  // std::cerr << "cg step" << std::endl;
	  if ((error = cg_ff->apply(cg_topo, cg_conf, cg_sim))){
	    cg_conf.current().energies.calculate_totals();
	    io::print_ENERGY(traj.output(), cg_conf.current().energies,
			     cg_topo.energy_groups(),
			     "CGOLDERROR", "CGOLDERR_");
	    
	    io::print_ENERGY(traj.output(), cg_conf.old().energies, 
			     cg_topo.energy_groups(),
			     "CGERROR", "CGERR_");
	    
	    std::cout << "\nError during CG final energy calc!\n" << std::endl;
	    return 1;
	  }

	  cg_conf.current().energies.calculate_totals();
	  io::print_ENERGY(traj.output(), cg_conf.current().energies,
			   cg_topo.energy_groups(),
			   "CGENERGY_Li", "CGELi_");
	} // multigraining
	
	algorithm::Algorithm * ff = md.algorithm("Forcefield");
	
	if (ff == NULL){
	  std::cerr << "forcefield not found in MD algorithm sequence"
		    << std::endl;
	  return 1;
	}
	
	if (ff->apply(topo, conf, sim)){
	  io::print_ENERGY(traj.output(), conf.current().energies,
			   topo.energy_groups(),
			   "ERROR", "ERR_");
	  
	  io::print_ENERGY(traj.output(), conf.old().energies, 
			   topo.energy_groups(),
			   "OLDERROR", "OLDERR_");
	  
	  std::cout << "\nError during final energy calc!\n" << std::endl;
	  return 1;
	}

	conf.current().energies.calculate_totals();
	replica_data.epot_i = conf.current().energies.potential_total +
	  conf.current().energies.special_total;
	
	io::print_ENERGY(traj.output(), conf.current().energies,
			 topo.energy_groups(),
			 "ENERGY_Li", "ELi_");
	
	traj.write_replica_energy(replica_data, sim,
				  conf.current().energies,
				  0);

	if (replica_data.Ti != replica_data.Tj)
	  std::cout << "replica_energy " << replica_data.epot_i 
		    << " @ " << T[replica_data.Ti] << "K" << std::endl;
	else if (replica_data.li != replica_data.lj)
	  std::cout << "replica_energy " << replica_data.epot_i 
		    << " @ " << l[replica_data.li] << std::endl;
	
	if (replica_data.li != replica_data.lj){
	  
	  std::cout << "energies at switched lambda:\n\n";

	  // change the lambda value
	  sim.param().perturbation.lambda = l[replica_data.lj];
	  topo.lambda(l[replica_data.lj]);
	  // twice: set old lambda as well
	  topo.lambda(l[replica_data.lj]);
	  topo.update_for_lambda();
	  std::cout << "\tlambda = " << topo.lambda() << "\n";

	  if (multigraining){
	    cg_sim.param().perturbation.lambda = l[replica_data.lj];
	    cg_topo.lambda(l[replica_data.lj]);
	    cg_topo.lambda(l[replica_data.lj]);
	    cg_topo.update_for_lambda();

	    if ((error = cg_ff->apply(cg_topo, cg_conf, cg_sim))){
	      cg_conf.current().energies.calculate_totals();
	      io::print_ENERGY(traj.output(), cg_conf.current().energies,
			       cg_topo.energy_groups(),
			       "CGOLDERROR", "CGOLDERR_");
	      
	      io::print_ENERGY(traj.output(), cg_conf.old().energies, 
			       cg_topo.energy_groups(),
			       "CGERROR", "CGERR_");
	      
	      std::cout << "\nError during CG force recalc at different lambda!\n"
			<< std::endl;
	      return 1;
	    }

	    cg_conf.current().energies.calculate_totals();
	    io::print_ENERGY(traj.output(), cg_conf.current().energies,
			     cg_topo.energy_groups(),
			     "CGENERGY_Lj", "CGELj_");

	  }

	  // recalc energy
	  if (ff->apply(topo, conf, sim)){
	    io::print_ENERGY(traj.output(), conf.current().energies,
			     topo.energy_groups(),
			     "OLDERROR", "OLDERR_");
	    
	    io::print_ENERGY(traj.output(), conf.old().energies, 
			     topo.energy_groups(),
			     "ERROR", "ERR_");
	    
	    std::cout << "\nError during force recalc at different lambda!\n" << std::endl;
	  }
	  
	  conf.current().energies.calculate_totals();
	  replica_data.epot_j = conf.current().energies.potential_total +
	    conf.current().energies.special_total;

	  io::print_ENERGY(traj.output(), conf.current().energies,
			   topo.energy_groups(),
			   "ENERGY_Lj", "ELj_");

	  traj.write_replica_energy(replica_data, sim,
				    conf.current().energies,
				    1);

	} 
	else{
	  replica_data.epot_j = replica_data.epot_i;
	}
	
	++replica_data.run;
	replica_data.time = sim.time();
	replica_data.state = waiting;
      }
      else{
	replica_data.state = st_error;
      }
	
      DEBUG(8, "slave: connecting (after run)...");
      cl_socket = socket(addrinfo_p->ai_family, addrinfo_p->ai_socktype,
			 addrinfo_p->ai_protocol);

      result = connect(cl_socket, s_addr_p, len);
      if (result == -1){
	std::cout << "could not (re-)connect to master. master finished?"
		  << std::endl;
	return 1;
      }

      DEBUG(9, "slave: connected");

      double magic[4] = { 3.1415927, 29375, 243, 8.3116 };
      double magic_buff[4];
      
      // magic cookie exchange
      if (read(cl_socket, magic_buff, 4 * sizeof(double)) != 4 * sizeof(double)){
	std::cerr << "could not read magic cookie" << std::endl;
	close(cl_socket);
	continue;
      }
      if (write(cl_socket, &magic, 4 * sizeof(double)) != 4 * sizeof(double)){
	std::cerr << "could not write magic cookie" << std::endl;
	close(cl_socket);
	continue;
      }
      if (magic[0] != magic_buff[0] || magic[1] != magic_buff[1] ||
	  magic[2] != magic_buff[2] || magic[3] != magic_buff[3]){
	
	std::cerr << "magic cookie exchange failed" << std::endl;
	close(cl_socket);
	continue;
      }
      else{
	std::cout << "magic cookie test succeeded" << std::endl;
      }

      DEBUG(8, "slave: finished job " << replica_data.ID);

      char ch = 2;
      if (write(cl_socket, &ch, 1) != 1){
	std::cerr << "could not write to socket" << std::endl;
	close(cl_socket);
	return 1;
      }
      
      if (read(cl_socket, &ch, 1) != 1){
	std::cerr << "could not read" << std::endl;
	close(cl_socket);
	return 1;
      }
      
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
      std::cout << "slave sleeping ..." << std::endl;
      // wait apropriate time for master to prepare
      sleep(timeout);
    }
  } // for slave_runs
  
  std::cerr << "slave: finished runs. terminating..." << std::endl;

  freeaddrinfo(addrinfo_p);

  return 0;
}

int util::Replica_Exchange_Slave::run_md
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 algorithm::Algorithm_Sequence & md,
 topology::Topology & cg_topo,
 configuration::Configuration & cg_conf,
 simulation::Simulation & cg_sim,
 interaction::Forcefield *cg_ff,
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

    if (multigraining){
      // coarse grained atom positions are based upon
      // real atom positions
      util::update_virtual_pos(cg_topo, cg_conf, topo, conf);

      // calculate the cg forces first!
      if ((error = cg_ff->apply(cg_topo, cg_conf, cg_sim))){
	io::print_ENERGY(traj.output(), cg_conf.current().energies,
			 cg_topo.energy_groups(),
			 "CGOLDERROR", "CGOLDERR_");
	
	io::print_ENERGY(traj.output(), cg_conf.old().energies, 
			 cg_topo.energy_groups(),
			 "CGERROR", "CGERR_");
	
	std::cout << "\nError during CG MD run!\n" << std::endl;
	break;
      }
    }

    // run a step
    if ((error = md.run(topo, conf, sim))){

      io::print_ENERGY(traj.output(), conf.current().energies,
		       topo.energy_groups(),
		       "OLDERROR", "OLDERR_");
      
      io::print_ENERGY(traj.output(), conf.old().energies, 
		       topo.energy_groups(),
		       "ERROR", "ERR_");
      
      std::cout << "\nError during MD run!\n" << std::endl;
      break;
    }
    
    traj.print(topo, conf, sim);

    sim.time() += sim.time_step_size();
    ++sim.steps();

    if (multigraining){
      cg_sim.time() += cg_sim.time_step_size();
      ++cg_sim.steps();
    }
    
  }
    
  return error;
}

int util::Replica_Exchange_Slave::get_replica_data()
{
  DEBUG(8, "slave: requesting job");
  char ch = 1;
  if (write(cl_socket, &ch, 1) != 1){
    std::cerr << "could not write to socket" << std::endl;
    close(cl_socket);
    return 1;
  }  

  if (read(cl_socket, &ch, 1) != 1){
    std::cerr << "could not read" << std::endl;
    close(cl_socket);
    return 1;
  }
    
  if (ch != 0){
    std::cout << "no job received from server (" << int(ch) << ")"
	      << std::endl;
  }
  else{
    
    DEBUG(8, "slave: waiting for replica data");
    if (read(cl_socket, (char *) &replica_data, sizeof(Replica_Data)) != sizeof(Replica_Data)){
      std::cerr << "could not read" << std::endl;
      close(cl_socket);
      return 1;
    }

    DEBUG(9, "slave: got replica " << replica_data.ID
	  << " temperature=" << replica_data.Ti
	  << " lambda=" << replica_data.li);
  }

  return int(ch);
}

int util::Replica_Exchange_Slave::update_replica_data()
{
  DEBUG(8, "slave: updating replica data");

  if (write(cl_socket, (char *) &replica_data.ID, sizeof(int)) != sizeof(int)){
    std::cerr << "could not write to socket" << std::endl;
    close(cl_socket);
    return 1;
  }

  if (write(cl_socket, (char *) &replica_data, sizeof(Replica_Data))
      != sizeof(Replica_Data)){

    std::cerr << "could not write to socket" << std::endl;
    close(cl_socket);
    return 1;
  }    

  return 0;
}

int util::Replica_Exchange_Slave::get_configuration
(
 topology::Topology const & topo,
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

  ssize_t n_rec = 0;

  readblock((char *) &conf.current().pos(0)(0), num);

  readblock((char *) &conf.current().vel(0)(0), num);   

  if ((n_rec = read(cl_socket, (char *) &conf.current().box(0)(0),
	   9 * sizeof(double))) != 9 * sizeof(double)){
    std::cerr << "could not read box" << std::endl;
    std::cerr << "got: " << n_rec << "\texpected: " << num << std::endl;
    close(cl_socket);
    return 1;
  }

  return 0;
}

int util::Replica_Exchange_Slave::update_configuration
(
 topology::Topology const & topo,
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

  writeblock((char *) &conf.current().pos(0)(0), num);
    
  writeblock((char *) &conf.current().vel(0)(0), num);  

  if (write(cl_socket, (char *) &conf.current().box(0)(0),
	    9 * sizeof(double)) != 9 * sizeof(double)){
    std::cerr << "could not write to socket" << std::endl;
    close(cl_socket);
    return 1;
  }

  return 0;
}

int util::Replica_Exchange_Slave::init_replica
(
 topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim,
 topology::Topology & cg_topo,
 simulation::Simulation & cg_sim
 )
{
  // get configuration from master
  if (get_configuration(topo, conf)){
    std::cerr << "get configuration failed" << std::endl;
    return 1;
  }
  
  std::vector<double> const & T = sim.param().replica.temperature;
  std::vector<double> const & l = sim.param().replica.lambda;

  // change all the temperature coupling temperatures
  for(unsigned int i=0; i<sim.multibath().size(); ++i){
    sim.multibath()[i].temperature = T[replica_data.Ti];
    if (multigraining){
      cg_sim.multibath()[i].temperature = T[replica_data.Ti];
      // temperature constraining
      sim.multibath()[i].tau = sim.param().replica.dt[replica_data.li];
      cg_sim.multibath()[i].tau = sim.param().replica.dt[replica_data.li];
    }
  }
  
  // change the lambda value
  sim.param().perturbation.lambda = l[replica_data.li];
  topo.lambda(l[replica_data.li]);
  // twice, to set old_lambda...
  topo.lambda(l[replica_data.li]);
  topo.update_for_lambda();
  std::cout << "\tlambda = " << topo.lambda() << "\n";

  if (multigraining){
    cg_sim.param().perturbation.lambda = l[replica_data.li];
    cg_topo.lambda(l[replica_data.li]);
    cg_topo.lambda(l[replica_data.li]);
    cg_topo.update_for_lambda();
  }
  
  // change simulation time
  sim.time() = replica_data.time;
  sim.time_step_size() = sim.param().replica.dt[replica_data.li];
  sim.steps() = replica_data.run * sim.param().step.number_of_steps;

  return 0;
}
