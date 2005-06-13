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

  // replica output files
  std::vector<std::ofstream *> rep_out;
  
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
    
    /*
      should be written after accepting a job
    // print header
    (*rep_out[i]) << "#"
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
    */
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
	    // std::cout << "\tmaster: replica " << r 
	    // << " is waiting (trying switch)" << std::endl;
	    switch_replica(r);
	  }
	    
	  if(replica_data[r].state == ready && replica_data[r].run == trials){
	    std::cout << "\tmaster: replica " << r << " is ready for run "
		      << trials << "..." << std::endl;
	    
	    // print state
	    std::cout << std::setw(5) << r
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

	    // assign it!
	    replica_data[r].state = running;

	    snapshot[neighbour_pos[r]] = r;
	    temp_snapshot[r+1] = replica_data[r].temperature;
	    lambda_snapshot[r+1] = replica_data[r].lambda;

	    MPI_Send(&replica_data[r], sizeof(Replica_Data), MPI_CHAR,
		     0, 0, client);
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
	++runs;
	break;
      case 3:
	std::cout << "process " << i << " has aborted run\n";
	break;
      default:
	std::cout << "message not understood\n";
    }

    std::cout << "disconnecting..." << std::endl;
    MPI_Comm_disconnect(&client);

    // see what the thread wants...

    /*
      if (slave_data[i].state == waiting){
	// std::cout << "\tmaster: thread " << i << " is waiting (run = "
	// << runs << ")..." << std::endl;
	
      } // thread waiting
    */

  } // while
    
  // simulation done
  DEBUG(9, "master: done");

  MPI_Comm_free(&client);
  MPI_Unpublish_name((char *) args["master"].c_str(), MPI_INFO_NULL, port_name);
  MPI_Close_port(port_name);

  out_neighbour.flush();
  out_neighbour.close();

  out_temp.flush();
  out_temp.close();

  out_lambda.flush();
  out_lambda.close();

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
  // m_ID = tid;
  
  // prepare output file
  /**
  std::ostringstream oss;
  oss << "repex_" << m_ID << ".out";
  std::ofstream os(oss.str().c_str());
  */
  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  // read the files (could also be copied from the master)
  io::read_input(args, topo, conf, sim,  md, std::cout);

  // initialises everything
  md.init(topo, conf, sim, std::cout);

  // trajectory (per thread!)
  // create output files...
  // rename them appropriately
  // (not necessary anymore...)
  std::ostringstream suff;
  // suff << "." << m_ID;
  
  std::string fin = "";
  if (args.count("fin")) fin = args["fin"] + suff.str();
  std::string trj = "";
  if (args.count("trj")) trj = args["trj"] + suff.str();
  std::string tre = "";
  if (args.count("tre")) tre = args["tre"] + suff.str();
  std::string trg = "";
  if (args.count("trg")) trg = args["trg"] + suff.str();
  std::string bae = "";
  if (args.count("bae")) bae = args["bae"] + suff.str();
  std::string bag = "";
  if (args.count("bag")) bag = args["bag"] + suff.str();
  std::string trv = "";
  if (args.count("trv")) trv = args["trv"] + suff.str();
  std::string trf = "";
  if (args.count("trf")) trf = args["trf"] + suff.str();

  io::Out_Configuration traj("GromosXX\n", std::cout);
  traj.title("GromosXX\n" + sim.param().title);
  traj.init(fin, trj, trv, trf, tre, trg, bae, bag, sim.param());

  // should be the same as from master...
  std::cout << "\nMESSAGES FROM (SLAVE) INITIALIZATION\n";
  if (io::messages.display(std::cout) >= io::message::error){
    // exit
    std::cout << "\nErrors during initialization!\n" << std::endl;
    // os.flush();
    // os.close();
    return 1;
  }
  io::messages.clear();
  
  std::cout.precision(5);
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

  std::cout << "slave " << m_ID << " initialised" << std::endl;
  // std::cerr << "slave " << m_ID << " initialised" << std::endl;

  char port_name[MPI::MAX_PORT_NAME];

  if (args.count("slave") != 1){
    io::messages.add("slave: connection name required",
		     "replica exchange",
		     io::message::error);
    MPI_Finalize();
    return 1;
  }
  
  MPI::Lookup_name(args["slave"].c_str(), MPI_INFO_NULL, port_name);

  while(slave_data.state != terminate){
    
    DEBUG(8, "slave: connecting..");
    sleep(1);
    MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);

    DEBUG(9, "slave: connected");

    // request a job
    int server_response = get_replica_data();

    // get state
    // get_state();

    // close connection
    MPI_Comm_disconnect(&master);

    if (server_response == 0){
      // accept job
      std::cout << "slave " << m_ID << " got a job! (replica "
		<< slave_data.replica << ")" << std::endl;
	
      std::cout << "slave " << m_ID << " running replica " << slave_data.replica
		<< " at T=" << replica_data.temperature
		<< " and l=" << replica_data.lambda
		<< " (run = " << replica_data.run << ")"
		<< std::endl;

      // init replica parameters (t, conf, T, lambda)
      init_replica(topo, conf, sim);
      // run it
      run_md(topo, conf, sim, md, traj);
      std::cout << "replica_energy final "
		<< conf.old().energies.potential_total << std::endl;
      
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
      replica_data.energy = conf.current().energies.potential_total;
      
      std::cout << "replica_energy " << replica_data.energy 
		<< " @ " << replica_data.temperature << "K" << std::endl;
      std::cout << "pos(0) " << math::v2s(conf.current().pos(0)) << std::endl;
      
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

      ++replica_data.run;
      replica_data.state = waiting;
      slave_data.state = waiting;
	
      DEBUG(8, "slave: connecting (after run)...");
      MPI_Comm_connect(port_name, MPI::INFO_NULL, 0, MPI::COMM_WORLD, &master);
      
      DEBUG(9, "slave: connected");

      // store configuration on master
      update_configuration(topo, conf);
      update_replica_data();
      // update_slave_data();

      // and disconnect
      std::cout << "disconnecting..." << std::endl;
      MPI_Comm_disconnect(&master);

    }
    else{
      std::cout << "slave " << m_ID << " waiting..." << std::endl;
    }
    
  } // while
    
  MPI_Comm_free(&master);
  
  // os.flush();
  // os.close();

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
      // try to save the final structures...
      break;
    }

    traj.print(topo, conf, sim);

    sim.time() += sim.time_step_size();
    ++sim.steps();

  }
    
  // std::cout << "writing final configuration" << std::endl;
    
  // traj.write(conf, topo, sim, io::final);
  // traj.print_final(topo, conf, sim);
    
  // std::cout << "\nMESSAGES FROM SIMULATION\n";
  // io::message::severity_enum err_msg = io::messages.display(std::cout);

  // std::cout << "\n\n";
    
  return error;
}

int util::Replica_Exchange_Slave::get_slave_data()
{
  // slave_data = replica_master->slave_data[m_ID];

  return 0;
}

int util::Replica_Exchange_Slave::update_slave_data()
{
  // replica_master->slave_data[m_ID] = slave_data;

  return 0;
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
  
  DEBUG(9, "tag = " << status.MPI_TAG);
  DEBUG(9, "slave: got replica " << replica_data.ID
	<< " temperature=" << replica_data.temperature
	<< " lambda=" << replica_data.lambda);
  
  return status.MPI_TAG;
}

int util::Replica_Exchange_Slave::update_replica_data()
{
  int i = m_ID;
  MPI_Status status;
  
  DEBUG(8, "slave: finished job");
  MPI_Send(&i, 1, MPI_INT, 0, 2, master);

  return 0;
}

int util::Replica_Exchange_Slave::get_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  // conf = replica_master->conf(slave_data.replica);
  
  return 0;
}

int util::Replica_Exchange_Slave::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  // replica_master->conf(slave_data.replica) = conf;

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

#endif
