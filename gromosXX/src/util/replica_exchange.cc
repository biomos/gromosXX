/**
 * @file replica_exchange.cc
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
// replica exchange ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange::Replica_Exchange()
  : m_ID(-1)
{
}

////////////////////////////////////////////////////////////////////////////////
// replica master   ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Master::Replica_Exchange_Master()
{
  DEBUG(7, "creating replica master");
  
  // enable control via environment variables
  gsl_rng_env_setup();
  const gsl_rng_type * rng_type = gsl_rng_default;

  // get the random number generator
  m_rng = gsl_rng_alloc(rng_type);

  // master
  m_ID = 0;
}

int util::Replica_Exchange_Master::run
(
 io::Argument & args)
{
  // master writes to std::cout
  // MPI::Intracomm client;
  MPI_Comm client;
  MPI_Status status;
  
  char port_name[MPI::MAX_PORT_NAME];

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  if(mpi_size != 1) 
    std::cerr << "Server consists of multiple MPI processes" << std::endl;
  
  MPI::Open_port(MPI_INFO_NULL, port_name);
  printf("server available at %s or by name: gromosXX\n", port_name);

  if (args.count("master") != 1){
    io::messages.add("master: connection name required",
		     "replica exchange",
		     io::message::error);
    MPI_Finalize();
    return 1;
  }
  MPI::Publish_name(args["master"].c_str(), MPI_INFO_NULL, port_name);

  DEBUG(8, "replica master registered: " << port_name);

  // create the simulation classes (necessary to store the configurations)
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  // read the files
  io::read_input(args, topo, conf, sim,  md, std::cout);

  // initialises everything
  md.init(topo, conf, sim, std::cout);

  // check input
  if (sim.param().replica.number <= 1){
    io::messages.add("replica exchange with less than 2 replicas?!",
		     "Replica_Exchange",
		     io::message::error);
  }

  if (sim.param().replica.trials < 1){
    io::messages.add("replica exchange with less than 1 trials?!",
		     "Replica_Exchange",
		     io::message::error);
  }

  // setup replica information
  replica_data.resize(sim.param().replica.number);
  neighbour.push_back(-1);

  for(int i=0; i<sim.param().replica.number; ++i){
    replica_data[i].temperature = sim.param().replica.temperature[i];
    replica_data[i].lambda = sim.param().replica.lambda[i];
    replica_data[i].energy = 0.0;
    replica_data[i].switch_temperature = 0.0;
    replica_data[i].switch_lambda = 0.0;
    replica_data[i].switch_energy = 0.0;
    replica_data[i].ID = i;
    replica_data[i].run = 0;
    replica_data[i].state = waiting;
    replica_data[i].probability = 0.0;
    replica_data[i].switched = false;

    neighbour.push_back(i);
    neighbour_pos.push_back(i+1);

    // store the positions of the all replicas
    // Change : maybe start from different initial positions!
    m_conf.push_back(conf);
    
  }

  neighbour.push_back(-1);

  std::vector<int> snapshot = neighbour;
  std::vector<double> temp_snapshot(neighbour.size(), 0.0);
  std::vector<double> lambda_snapshot(neighbour.size(), 0.0);

  std::cout << "\nMESSAGES FROM (MASTER) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    std::cout << "\nErrors during initialization!\n" << std::endl;
    return 1;
  }
  
  io::messages.clear();

  // set seed if not set by environment variable
  // (get from sim.param()...)
  if (gsl_rng_default_seed == 0)
    gsl_rng_set (m_rng, 123456789);

  std::cout << "master thread initialised" << std::endl;

  std::cout << "Replica Exchange\n"
	    << "\treplicas:\t" << sim.param().replica.number << "\n"
	    << "\tthreads:\t" << "unspecified" << "\n\n"
	    << "\t" << std::setw(10) << "ID"
	    << std::setw(20) << "Temp"
	    << std::setw(20) << "lambda\n";
    
  for(int i=0; i<sim.param().replica.number; ++i){
    std::cout << "\t" << std::setw(10) << i+1
	      << std::setw(20) << replica_data[i].temperature
	      << std::setw(20) << replica_data[i].lambda
	      << "\n";
  }
  std::cout << "\n\ttrials:\t" << sim.param().replica.trials << "\n"
	    << "\nEND\n" << std::endl;

  int trials = 0;
  int runs = 0;
    
  std::ofstream out_neighbour("neighbour.dat");
  out_neighbour << "#" << std::setw(4) << "time"
		<< std::setw(5) << "rep1"
		<< std::setw(5) << "..."
		<< std::endl;

  std::ofstream out_temp("temperature.dat");
  out_temp << "#" << std::setw(4) << "time"
	   << std::setw(5) << "T1"
	   << std::setw(5) << "..."
	   << std::endl;
  
  std::ofstream out_lambda("lambda.dat");
  out_lambda << "#" << std::setw(4) << "time"
	     << std::setw(5) << "l1"
	     << std::setw(5) << "..."
	     << std::endl;

  std::ofstream rep_out("replica.dat");
  rep_out << "#"
	  << std::setw(5) << "ID"
	  << std::setw(5) << "run"
	  << std::setw(6) << "pos"
	  << std::setw(13) << "T"
	  << std::setw(13) << "l"
	  << std::setw(13) << "Epot"
	  << std::setw(13) << "sT"
	  << std::setw(13) << "sl"
	  << std::setw(13) << "sEpot"
	  << std::setw(13) << "p"
	  << std::setw(4) << "s"
	  << "\n";

  while(true){
    
    if (runs == sim.param().replica.number){

      DEBUG(9, "master: finished trial period " << trials);

      out_neighbour << std::setw(8) << trials;
      for(unsigned int i=1; i<snapshot.size()-1; ++i){
	out_neighbour << std::setw(5) << snapshot[i];
      }
      out_neighbour << std::endl;

      out_temp << std::setw(8) << trials;
      for(unsigned int i=1; i<temp_snapshot.size()-1; ++i){
	out_temp << std::setw(5) << temp_snapshot[i];
      }
      out_temp << std::endl;

      out_lambda << std::setw(8) << trials;
      for(unsigned int i=1; i<lambda_snapshot.size()-1; ++i){
	out_lambda << std::setw(5) << lambda_snapshot[i];
      }
      out_lambda << std::endl;

      ++trials;
      runs = 0;
    }
      
    if(trials > sim.param().replica.trials){
      DEBUG(8, "master: finished all trials...");
      break;
    }
    
    // wait for a thread to connect
    DEBUG(9, "master: accepting connection...");
    MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD,
		    &client);

    DEBUG(9, "client connected!");

    int i;
    MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE,
	     MPI_ANY_TAG, client, &status);

    switch(status.MPI_TAG){
      case 0:
	std::cout << "process " << i << " says hello\n";
	break;
      case 1:
	std::cout << "process " << i << " requests job\n";

	// select a replica to run
	int r;
	for(r=0; r<sim.param().replica.number; ++r){
	  
	  if(replica_data[r].state == waiting){
	    // try a switch
	    switch_replica(r);
	  }
	    
	  if(replica_data[r].state == ready && replica_data[r].run == trials){
	    std::cout << "\tmaster: replica " << r << " is ready for run "
		      << trials << "..." << std::endl;
	    
	    // print state
	    std::cout << std::setw(6) << r
		      << std::setw(6) << replica_data[r].run
		      << std::setw(6) << neighbour_pos[r]
		      << std::setw(13) << replica_data[r].temperature
		      << std::setw(13) << replica_data[r].lambda
		      << std::setw(13) << replica_data[r].energy
		      << std::setw(13) << replica_data[r].switch_temperature
		      << std::setw(13) << replica_data[r].switch_lambda
		      << std::setw(13) << replica_data[r].switch_energy
		      << std::setw(13) << replica_data[r].probability
		      << std::setw(4) << replica_data[r].switched
		      << std::endl;

	    rep_out << std::setw(6) << r
		    << std::setw(6) << replica_data[r].run
		    << std::setw(6) << neighbour_pos[r]
		    << std::setw(13) << replica_data[r].temperature
		    << std::setw(13) << replica_data[r].lambda
		    << std::setw(13) << replica_data[r].energy
		    << std::setw(13) << replica_data[r].switch_temperature
		    << std::setw(13) << replica_data[r].switch_lambda
		    << std::setw(13) << replica_data[r].switch_energy
		    << std::setw(13) << replica_data[r].probability
		    << std::setw(4) << replica_data[r].switched
		    << "\n";
	    
	    // assign it!
	    replica_data[r].state = running;

	    snapshot[neighbour_pos[r]] = r;
	    temp_snapshot[r+1] = replica_data[r].temperature;
	    lambda_snapshot[r+1] = replica_data[r].lambda;

	    // send parameters
	    MPI_Send(&replica_data[r], sizeof(Replica_Data), MPI_CHAR,
		     0, 0, client);

	    // positions
	    DEBUG(9, "sending " << 3 * m_conf[r].current().pos.size() << " coords");

	    MPI_Send(&m_conf[r].current().pos(0)(0), m_conf[r].current().pos.size()*3,
		     MPI::DOUBLE, 0, 0, client);
	    
	    // and velocities
	    MPI_Send(&m_conf[r].current().vel(0)(0), m_conf[r].current().vel.size()*3,
		     MPI::DOUBLE, 0, 0, client);

	    DEBUG(10, "master: pos = " << math::v2s(m_conf[r].current().pos(0)));
	    DEBUG(10, "master: vel = " << math::v2s(m_conf[r].current().vel(0)));

	    break;
	  }
	} // replica selected

	if (r==sim.param().replica.number){
	  // error!
	  std::cout << "could not select replica!!!" << std::endl;
	  MPI_Send(&replica_data[0], sizeof(Replica_Data), MPI_CHAR,
		   0, 1, client);
	}

	break;

      case 2:
	std::cout << "process " << i << " has finished run\n";
	
	MPI_Status status;
  
	DEBUG(8, "master: waiting for replica data " << i);
	MPI_Recv(&replica_data[i], sizeof(Replica_Data), MPI_CHAR,
		 MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
  
	if (replica_data[i].state != st_error){
	  // get configuration
	  MPI_Recv(&m_conf[i].current().pos(0)(0),
		   m_conf[i].current().pos.size() * 3,
		   MPI::DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
	  
	  if (status.MPI_ERROR != MPI_SUCCESS){
	    std::cout << "MPI ERROR!!! " << status.MPI_ERROR << std::endl;
	  }
	  
	  MPI_Recv(&m_conf[i].current().vel(0)(0),
		   m_conf[i].current().vel.size()*3,
		   MPI::DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client, &status);
	  
	  if (status.MPI_ERROR != MPI_SUCCESS){
	    std::cout << "MPI ERROR!!! " << status.MPI_ERROR << std::endl;
	  }
	}
	
	DEBUG(9, "master: got replica " << replica_data[i].ID
	      << " temperature=" << replica_data[i].temperature
	      << " lambda=" << replica_data[i].lambda);
	++runs;
	break;

      case 3:
	std::cout << "process " << i << " has aborted run\n";
	break;

      case 4: // interactive session
	switch(i){
	  case 1: // replica information
	    {
	      int i = replica_data.size();
	      MPI_Send(&i, 1, MPI_INT, 0, 4, client);
	      
	      MPI_Send(&replica_data[0],
		       sizeof(Replica_Data)*replica_data.size(),
		       MPI_CHAR, 0, 4, client);
	      break;
	    }
	  case 2: // replica change
	    {
	      int i;
	      MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE,
		       MPI_ANY_TAG, client, &status);
	      
	      MPI_Recv(&replica_data[i], sizeof(Replica_Data),
		       MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
		       client, &status);
	      break;
	    }
	  case 3: // quit
	    {
	      std::cout << "master: stopping" << std::endl;
	      trials = sim.param().replica.trials + 1;
	      break;
	    }
	}
	
	break;

      default:
	std::cout << "message not understood\n";
    }

    // std::cout << "disconnecting..." << std::endl;
    MPI_Comm_disconnect(&client);

  } // while
    
  // simulation done
  DEBUG(9, "master: done");

  MPI_Unpublish_name((char *) args["master"].c_str(), MPI_INFO_NULL, port_name);
  MPI_Close_port(port_name);

  out_neighbour.flush();
  out_neighbour.close();

  out_temp.flush();
  out_temp.close();

  out_lambda.flush();
  out_lambda.close();

  rep_out.flush();
  rep_out.close();
  
  return 0;
}

int util::Replica_Exchange_Master::switch_replica(int i)
{
  assert(i>=0 && unsigned(i)<replica_data.size());
  DEBUG(8, "switch replica: " << i);

  if (replica_data[i].state != waiting){
    // why here???
    assert(false);
    return 0;
  }
  
  const int j = ((replica_data[i].run % 2) ==  (neighbour_pos[i] % 2) ) ?
    neighbour[neighbour_pos[i] - 1] : neighbour[neighbour_pos[i] + 1];
  
  if (replica_data[i].run == 0 || j == -1){
    
    const int offset = (((replica_data[i].run + 1) % 2) ==  (neighbour_pos[i] % 2) ) ?
      -1 : +1;
    
    const int neighbour_i = neighbour[neighbour_pos[i] + offset];
    
    if (neighbour_i != -1){
      replica_data[i].switch_temperature = replica_data[neighbour_i].temperature;
      replica_data[i].switch_lambda = replica_data[neighbour_i].lambda;
    }
    else{
      replica_data[i].switch_temperature = replica_data[i].temperature;
      replica_data[i].switch_lambda = replica_data[i].lambda;
    }

    replica_data[i].probability = 0.0;
    replica_data[j].probability = 0.0;
    replica_data[i].switched = false;
    replica_data[j].switched = false;
    
    replica_data[i].state = ready;
    return 0;
  }

  // j running or not yet at same run
  if (replica_data[j].state != waiting ||
      replica_data[i].run != replica_data[j].run){

    return 0;
  }

  // try switch...
  // calculate probability
  double delta = 0;
  if (replica_data[i].lambda != replica_data[j].lambda){
    // 2D formula
    delta = 1.0 / (math::k_Boltzmann * replica_data[i].temperature) *
      (replica_data[j].switch_energy - replica_data[i].energy) -
      1.0 / (math::k_Boltzmann * replica_data[j].temperature) *
      (replica_data[j].energy - replica_data[i].switch_energy);
  }
  else{
    // standard formula
    delta = (1.0 / (math::k_Boltzmann * replica_data[i].temperature) -
	     1.0 / (math::k_Boltzmann * replica_data[j].temperature)) *
      (replica_data[j].energy - replica_data[i].energy);
  }
  
  double probability = 1.0;
  if (delta > 0.0)
    probability = exp(-delta);
  
  replica_data[i].probability = probability;
  replica_data[j].probability = probability;
  replica_data[i].switched = false;
  replica_data[j].switched = false;
  
  const double r = gsl_rng_uniform(m_rng);

  if (r < probability){
    
    std::cout << "\t\tswitching i: " << replica_data[i].temperature
	      << " to j: " << replica_data[j].temperature << std::endl;
    
    replica_data[i].switched = true;
    replica_data[j].switched = true;
    
    const double l = replica_data[i].lambda;
    const double t = replica_data[i].temperature;
    replica_data[i].lambda = replica_data[j].lambda;
    replica_data[i].temperature = replica_data[j].temperature;
    replica_data[j].lambda = l;
    replica_data[j].temperature = t;
      
    const int p = neighbour_pos[i];
    neighbour_pos[i] = neighbour_pos[j];
    neighbour_pos[j] = p;
    
    neighbour[neighbour_pos[i]] = i;
    neighbour[neighbour_pos[j]] = j;
  }
    
  // get temperature / lambda for next switch
  const int offset = (((replica_data[i].run + 1) % 2) ==  (neighbour_pos[i] % 2) ) ?
    -1 : +1;

  const int neighbour_i = neighbour[neighbour_pos[i] + offset];
  
  if (neighbour_i != -1){
    replica_data[i].switch_temperature = replica_data[neighbour_i].temperature;
    replica_data[i].switch_lambda = replica_data[neighbour_i].lambda;
  }
  else{
    replica_data[i].switch_temperature = replica_data[i].temperature;
    replica_data[i].switch_lambda = replica_data[i].lambda;
  }

  const int neighbour_j = neighbour[neighbour_pos[j] - offset];
  
  if (neighbour_j != -1){
    replica_data[j].switch_temperature = replica_data[neighbour_j].temperature;
    replica_data[j].switch_lambda = replica_data[neighbour_j].lambda;
  }
  else{
    replica_data[j].switch_temperature = replica_data[j].temperature;
    replica_data[j].switch_lambda = replica_data[j].lambda;
  }
  
  replica_data[i].state = ready;
  replica_data[j].state = ready;
  
  return 0;
}

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

  io::Out_Configuration traj("GromosXX\n", std::cout);
  traj.title("GromosXX\n" + sim.param().title);
  traj.init(args["fin"], args["trj"], args["trv"], args["trf"], 
	    args["tre"], args["trg"], args["bae"], args["bag"],
	    sim.param());

  // should be the same as from master...
  std::cout << "\nMESSAGES FROM (SLAVE) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    std::cout << "\nErrors during initialization!\n" << std::endl;
    MPI_Finalize();
    return 1;
  }
  io::messages.clear();
  
  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

  std::cout << "slave " << m_ID << " initialised" << std::endl;

  char port_name[MPI::MAX_PORT_NAME];

  if (args.count("slave") != 1){
    io::messages.add("slave: connection name required",
		     "replica exchange",
		     io::message::error);
    MPI_Finalize();
    return 1;
  }
  
  MPI::Lookup_name(args["slave"].c_str(), MPI_INFO_NULL, port_name);

  while(true){
    
    DEBUG(8, "slave: connecting..");
    //////////////////////////////////////////////////////////////////////////////////
    // wait apropriate time for master to prepare
    sleep(3);
    //////////////////////////////////////////////////////////////////////////////////

    MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);

    DEBUG(9, "slave: connected");

    // request a job
    int server_response = get_replica_data();

    if (server_response == 0){
      // accept job
      std::cout << "slave: got a job! (replica "
		<< replica_data.ID << ")" << std::endl;
	
      std::cout << "slave:  running replica " << replica_data.ID
		<< " at T=" << replica_data.temperature
		<< " and l=" << replica_data.lambda
		<< " (run = " << replica_data.run << ")"
		<< std::endl;

      // init replica parameters (t, conf, T, lambda)
      init_replica(topo, conf, sim);

      // close connection
      MPI_Comm_disconnect(&master);

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
	  MPI_Finalize();
	  MPI_Comm_free(&master);
	  return 1;
	}
	
	ff->apply(topo, conf, sim);
	conf.current().energies.calculate_totals();
	replica_data.energy = conf.current().energies.potential_total;
	
	std::cout << "replica_energy " << replica_data.energy 
		  << " @ " << replica_data.temperature << "K" << std::endl;
	
	if (replica_data.lambda != replica_data.switch_lambda){
	  
	  // change the lambda value
	  sim.param().perturbation.lambda = replica_data.switch_lambda;
	  topo.lambda(replica_data.switch_lambda);
	  topo.update_for_lambda();
	  
	  // recalc energy
	  ff->apply(topo, conf, sim);
	  conf.current().energies.calculate_totals();
	  replica_data.switch_energy = conf.current().energies.potential_total;
	} 
	else{
	  replica_data.switch_energy = replica_data.energy;
	}
	
	++replica_data.run;
	replica_data.state = waiting;
      }
      else{
	replica_data.state = st_error;
      }
	
      DEBUG(8, "slave: connecting (after run)...");
      MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
      
      DEBUG(9, "slave: connected");

      DEBUG(8, "slave: finished job " << replica_data.ID);
      MPI_Send(&replica_data.ID, 1, MPI_INT, 0, 2, master);

      if (!error){
	// store configuration on master
	update_replica_data();
	update_configuration(topo, conf);
      }
      
      // and disconnect
      DEBUG(9, "disconnecting...");
      MPI_Comm_disconnect(&master);

    }
    else{
      MPI_Comm_disconnect(&master);
      std::cout << "slave " << m_ID << " waiting..." << std::endl;
    }
    
  } // while
    
  // MPI_Comm_free(&master);
  
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
  double end_time = sim.time() + 
    sim.time_step_size() * (sim.param().step.number_of_steps - 1);
    
  int error;

  while(sim.time() < end_time + math::epsilon){
      
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
  int i = m_ID;
  MPI_Status status;
  
  DEBUG(8, "slave: requesting job");
  MPI_Send(&i, 1, MPI_INT, 0, 1, master);
  
  DEBUG(8, "slave: waiting for replica data");
  MPI_Recv(&replica_data, sizeof(Replica_Data), MPI_CHAR,
	   MPI_ANY_SOURCE, MPI_ANY_TAG, master, &status);
  
  if (status.MPI_ERROR != MPI_SUCCESS){
    std::cout << "MPI ERROR!!! " << status.MPI_ERROR << std::endl;
  }

  DEBUG(9, "tag = " << status.MPI_TAG);
  DEBUG(9, "slave: got replica " << replica_data.ID
	<< " temperature=" << replica_data.temperature
	<< " lambda=" << replica_data.lambda);
  
  return status.MPI_TAG;
}

int util::Replica_Exchange_Slave::update_replica_data()
{
  DEBUG(8, "slave: updating replica data");
  MPI_Send(&replica_data, sizeof(Replica_Data), MPI_CHAR,
	   0, 0, master);

  return 0;
}

int util::Replica_Exchange_Slave::get_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  MPI_Status status;

  DEBUG(10, "receiving " << 3 * conf.current().pos.size() << " coords");
  MPI_Recv(&conf.current().pos(0)(0), conf.current().pos.size() * 3,
	   MPI::DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, master, &status);

  if (status.MPI_ERROR != MPI_SUCCESS){
    std::cout << "MPI ERROR!!! " << status.MPI_ERROR << std::endl;
  }
  
  MPI_Recv(&conf.current().vel(0)(0), conf.current().vel.size()*3,
	   MPI::DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, master, &status);

  if (status.MPI_ERROR != MPI_SUCCESS){
    std::cout << "MPI ERROR!!! " << status.MPI_ERROR << std::endl;
  }

  DEBUG(10, "slave: pos = " << math::v2s(conf.current().pos(0)));
  DEBUG(10, "slave: vel = " << math::v2s(conf.current().vel(0)));
  
  return status.MPI_TAG;
}

int util::Replica_Exchange_Slave::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  // positions
  DEBUG(9, "sending " << 3 * conf.current().pos.size() << " coords");
  
  MPI_Send(&conf.current().pos(0)(0), conf.current().pos.size()*3,
	   MPI::DOUBLE, 0, 0, master);
  
  // and velocities
  MPI_Send(&conf.current().vel(0)(0), conf.current().vel.size()*3,
	   MPI::DOUBLE, 0, 0, master);
  
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
  get_configuration(topo, conf);
  
  // change all the temperature coupling temperatures
  for(unsigned int i=0; i<sim.multibath().size(); ++i){
    sim.multibath()[i].temperature = replica_data.temperature;
  }
  
  // change the lambda value
  sim.param().perturbation.lambda = replica_data.lambda;
  topo.lambda(replica_data.lambda);
  topo.update_for_lambda();
  
  // change simulation time
  sim.time() = replica_data.run * 
    sim.param().step.number_of_steps * sim.param().step.dt +
    sim.param().step.t0;
  
  sim.steps() = replica_data.run * sim.param().step.number_of_steps;

  // ready to run ???

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// replica interactive    //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

util::Replica_Exchange_Interactive::Replica_Exchange_Interactive()
{
}

int util::Replica_Exchange_Interactive::run
(
 io::Argument & args
 )
{
  char port_name[MPI::MAX_PORT_NAME];

  if (args.count("control") < 1){
    io::messages.add("interactive: connection name required",
		     "replica exchange",
		     io::message::error);
    MPI_Finalize();
    return 1;
  }
  
  MPI::Lookup_name(args["control"].c_str(), MPI_INFO_NULL, port_name);

  DEBUG(8, "interactive: connecting..");
  MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
  
  int i = 1, nr;
  MPI_Status status;
  
  DEBUG(8, "interactive: requesting replica information");
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
      if (nr >= replica_data.size() || nr <0){
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
	is >> replica_data[nr].temperature;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].lambda;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].energy;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].switch_temperature;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].switch_lambda;
      }

      is.clear();
      ++it;
      if (it->second != "."){
	is.str(it->second);
	is >> replica_data[nr].switch_energy;
      }

      is.clear();
      ++it;
      if (it->second == "wait") replica_data[nr].state = waiting;
      if (it->second == "rdy") replica_data[nr].state = ready;
      if (it->second == "run") replica_data[nr].state = running;
      if (it->second == "err") replica_data[nr].state = st_error;
      if (it->second == "term") replica_data[nr].state = terminate;

      DEBUG(8, "interactive: connecting..");
      MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
  
      int i = 2;

      DEBUG(8, "interactive: requesting replica change");
      MPI_Send(&i, 1, MPI_INT, 0, 4, master);
  
      DEBUG(8, "interactive: sending change ID");
      MPI_Send(&nr, 1, MPI_INT, 0, 4, master);
  
      MPI_Send(&replica_data[nr], sizeof(Replica_Data),
	       MPI_CHAR, 0, 4, master);
  
      MPI_Comm_disconnect(&master);
      
    }
    else if (it->second == "stop"){
      std::cout << "send stop request to master" << std::endl;

      DEBUG(8, "interactive: connecting..");
      MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
  
      int i = 3;

      DEBUG(8, "interactive: requesting server quit");
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
	    << std::setw(10) << "Epot"
	    << std::setw(10) << "sT"
	    << std::setw(10) << "sl"
	    << std::setw(10) << "sEpot"
	    << std::setw(10) << "p"
	    << std::setw(5) << "s"
	    << std::setw(10) << "ste"
	    << "\n";
  
  for(int r=0; r<nr; ++r){
    std::cout << std::setw(6) << r + 1
	      << std::setw(6) << replica_data[r].run
	      << std::setw(10) << replica_data[r].temperature
	      << std::setw(10) << replica_data[r].lambda
	      << std::setw(10) << replica_data[r].energy
	      << std::setw(10) << replica_data[r].switch_temperature
	      << std::setw(10) << replica_data[r].switch_lambda
	      << std::setw(10) << replica_data[r].switch_energy
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
