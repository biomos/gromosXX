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

#ifdef OMP
#include <omp.h>
#endif

#include "replica_exchange.h"

namespace util
{
  Replica_Exchange_Master * replica_master = NULL;
}

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
 io::Argument & args,
 int tid, 
 int num_threads)
{
  if (tid != 0){
    std::cerr << "master thread should run with thread ID 0" << std::endl;
    return 1;
  }
  
  // master writes to std::cout

  // create the simulation classes (necessary to store the configurations)
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  // read the files
  io::read_input(args, topo, conf, sim,  md, std::cout);

  // initialises everything
  md.init(topo, conf, sim, std::cout);

  // initialise thread state, and assigned replicas
  slave_data.resize(num_threads, Slave_Data(waiting, -1));
  
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
    replica_data[i].TID = -1;
    replica_data[i].run = 0;
    replica_data[i].state = waiting;
    replica_data[i].probability = 0.0;
    replica_data[i].switched = false;

    neighbour.push_back(i);
    neighbour_pos.push_back(i+1);

    // store the positions of the all replicas
    // Change : maybe start from different initial positions!
    m_conf.push_back(conf);
    
    std::ostringstream oss;
    oss << "replica_" << i << ".dat";
    rep_out.push_back(new std::ofstream(oss.str().c_str()));

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

    rep_out[i]->precision(4);
    rep_out[i]->setf(std::ios::scientific, std::ios::floatfield);

  }

  neighbour.push_back(-1);
  std::vector<int> snapshot = neighbour;
  std::vector<double> temp_snapshot(neighbour.size(), 0.0);
  std::vector<double> lambda_snapshot(neighbour.size(), 0.0);

  std::cout << "\nMESSAGES FROM INITIALIZATION\n";
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
	    << "\tthreads:\t" << num_threads << "\n\n"
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

      out_neighbour << std::setw(8) << trials;
      for(int i=1; i<snapshot.size()-1; ++i){
	out_neighbour << std::setw(5) << snapshot[i];
      }
      out_neighbour << std::endl;

      out_temp << std::setw(8) << trials;
      for(int i=1; i<temp_snapshot.size()-1; ++i){
	out_temp << std::setw(5) << temp_snapshot[i];
      }
      out_temp << std::endl;

      out_lambda << std::setw(8) << trials;
      for(int i=1; i<lambda_snapshot.size()-1; ++i){
	out_lambda << std::setw(5) << lambda_snapshot[i];
      }
      out_lambda << std::endl;

      ++trials;
      runs = 0;
    }
      
    if(trials > sim.param().replica.trials) break;
    
    // is a thread waiting?
    // std::cout << "master: selecting thread..." << std::endl;

    for(int i=1; i<num_threads; ++i){

      if (slave_data[i].state == waiting){
	// std::cout << "\tmaster: thread " << i << " is waiting (run = "
	// << runs << ")..." << std::endl;
	
	// select a replica to run
	for(int r=0; r<sim.param().replica.number; ++r){
	  
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
	    (*rep_out[r]) << std::setw(6) << replica_data[r].run
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

	    slave_data[i].replica = r;
	    slave_data[i].state = ready;
	    ++runs;
	    break;
	  }
	} // replica selected
      } // thread waiting
    } // threads

    sleep(3);
  } // while
    
  // simulation done
  // stop all threads
  for(int i=1; i<num_threads; ++i){
    std::cerr << "terminating " << i << std::endl;
    slave_data[i].state = terminate;
  }

  // close output files
  for(int i=0; i<rep_out.size(); ++i){
    rep_out[i]->flush();
    rep_out[i]->close();
    delete rep_out[i];
  }

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
  assert(i>=0 && i<replica_data.size());

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

    /*
    std::cout << "\tswitch: no switch: state i = " << replica_data[i].state
	      << " state j = " << replica_data[j].state
	      << " run i = " << replica_data[i].run
	      << " run j = " << replica_data[j].run
	      << std::endl;
    */

    return 0;
  }

  // try switch...
  /*
  std::cout << "\tswitch: trying switch " << i << " - " << j 
	    << "(" << replica_data[i].temperature << " - "
	    << replica_data[j].temperature << ")" << std::endl;
  */

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

  // std::cout << "\t\tprob = " << probability << "\tr = " << r << std::endl;
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
 io::Argument & args,
 int tid, 
 int num_threads)
{
  m_ID = tid;
  
  // prepare output file
  std::ostringstream oss;
  oss << "repex_" << tid << ".out";
  std::ofstream os(oss.str().c_str());

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  // read the files (could also be copied from the master)
  io::read_input(args, topo, conf, sim,  md, os);

  // initialises everything
  md.init(topo, conf, sim, os);

  // trajectory (per thread!)
  // create output files...
  // rename them appropriately
  std::ostringstream suff;
  suff << "." << m_ID;
  
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

  io::Out_Configuration traj("GromosXX\n", os);
  traj.title("GromosXX\n" + sim.param().title);
  traj.init(fin, trj, trv, trf, tre, trg, bae, bag, sim.param());

  // should be the same as from master...
  os << "\nMESSAGES FROM INITIALIZATION\n";
  if (io::messages.display(os) >= io::message::error){
    // exit
    os << "\nErrors during initialization!\n" << std::endl;
    os.flush();
    os.close();
    return 1;
  }
  io::messages.clear();
  
  os.precision(5);
  os.setf(std::ios_base::fixed, std::ios_base::floatfield);

  os << "slave " << m_ID << " initialised" << std::endl;
  // std::cerr << "slave " << m_ID << " initialised" << std::endl;

  while(slave_data.state != terminate){

    get_slave_data();

    if (slave_data.state == ready){
      // accept job
      os << "slave " << m_ID << " got a job! (replica "
	 << slave_data.replica << ")" << std::endl;
	
      get_replica_data();
	
      os << "slave " << m_ID << " running replica " << slave_data.replica
	 << " at T=" << replica_data.temperature
	 << " and l=" << replica_data.lambda
	 << " (run = " << replica_data.run << ")"
	 << std::endl;
      /*
      std::cerr << "slave " << m_ID << " running replica " << slave_data.replica
		<< " at T=" << replica_data.temperature
		<< " and l=" << replica_data.lambda
		<< " (run = " << replica_data.run << ")"
		<< std::endl;
      */

      // init replica parameters (t, conf, T, lambda)
      init_replica(topo, conf, sim);
      // run it
      run_md(topo, conf, sim, md, traj);
      os << "replica_energy final " << conf.old().energies.potential_total << std::endl;
      
      // store configuration on master
      // (only necessary if more replicas than threads...)
      update_configuration(topo, conf);

      // do we need to reevaluate the potential energy ?
      // yes! 'cause otherwise it's just the energy for the previous configuration...
      algorithm::Algorithm * ff = md.algorithm("Forcefield");
      
      if (ff == NULL){
	std::cerr << "forcefield not found in MD algorithm sequence" << std::endl;
	return 1;
      }
      
      ff->apply(topo, conf, sim);
      conf.current().energies.calculate_totals();
      replica_data.energy = conf.current().energies.potential_total;
      
      os << "replica_energy " << replica_data.energy 
	 << " @ " << replica_data.temperature << "K" << std::endl;
      os << "pos(0) " << math::v2s(conf.current().pos(0)) << std::endl;
      
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
	
      update_replica_data();
      update_slave_data();

    }
    else{
      os << "slave " << m_ID << " waiting..." << std::endl;
      sleep(1);
    }
    
  } // while
    
  os.flush();
  os.close();

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
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  assert(replica_master->slave_data.size() > m_ID);
  
  slave_data = replica_master->slave_data[m_ID];

  // std::cerr << "slave " << m_ID << " state after synch = "
  // << slave_data.state << std::endl;
  return 0;
}

int util::Replica_Exchange_Slave::update_slave_data()
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  assert(replica_master->slave_data.size() > m_ID);
  
  replica_master->slave_data[m_ID] = slave_data;

  // std::cerr << "slave " << m_ID << " slave_data updated" << std::endl;
  return 0;
}

int util::Replica_Exchange_Slave::get_replica_data()
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that

  assert(replica_master != NULL);
  assert(replica_master->replica_data.size() > slave_data.replica);
  
  replica_data = replica_master->replica_data[slave_data.replica];
  // std::cerr << "slave " << m_ID << " got replica data" << std::endl;
  return 0;
}

int util::Replica_Exchange_Slave::update_replica_data()
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that

  assert(replica_master != NULL);
  assert(replica_master->replica_data.size() > slave_data.replica);
  
  replica_master->replica_data[slave_data.replica] = replica_data;
  
  // std::cerr << "slave " << m_ID << " replica data updated" << std::endl;
  return 0;
}

int util::Replica_Exchange_Slave::get_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that

  assert(replica_master != NULL);
  
  conf = replica_master->conf(slave_data.replica);
  
  // std::cerr << "slave " << m_ID 
  //           << " got replica configuration" << std::endl;
  return 0;
  
}

int util::Replica_Exchange_Slave::update_configuration
(
 topology::Topology const & topo,
 configuration::Configuration & conf
 )
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that

  assert(replica_master != NULL);
  
  replica_master->conf(slave_data.replica) = conf;
  
  // std::cerr << "slave " << m_ID 
  //           << " replica configuration updated" << std::endl;
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
  for(int i=0; i<sim.multibath().size(); ++i){
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
