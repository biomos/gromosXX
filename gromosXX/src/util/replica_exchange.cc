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
  Replica_Exchange * replica_master = NULL;
}

util::Replica_Exchange::Replica_Exchange()
{
  // enable control via environment variables
  gsl_rng_env_setup();
  const gsl_rng_type * rng_type = gsl_rng_default;

  // get the rundom number generator
  m_rng = gsl_rng_alloc(rng_type);

}

int util::Replica_Exchange::run(io::Argument & args, int tid, int num_threads)
{
  std::ostringstream oss;
  oss << "repex_" << tid << ".out";
  
  std::ofstream os(oss.str().c_str());

  // create the simulation classes
  topology::Topology topo;
  configuration::Configuration conf;
  algorithm::Algorithm_Sequence md;
  simulation::Simulation sim;

  io::Out_Configuration traj("GromosXX\n");

  io::read_input(args, topo, conf, sim,  md, os);

  traj.title("GromosXX\n" + sim.param().title);

  // create output files...
  traj.init(args, sim.param());

  // initialises all algorithms (and therefore also the forcefield)
  md.init(topo, conf, sim);

  // initialise thread state, and assigned replicas
  thread_state.resize(num_threads, waiting);
  thread_replica.resize(num_threads, -1);
  
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
    replica_data[i].switch_replica = -1;
    replica_data[i].switch_temperature = 0.0;
    replica_data[i].switch_lambda = 0.0;
    replica_data[i].switch_energy = 0.0;
    replica_data[i].TID = -1;
    replica_data[i].run = 0;
    replica_data[i].state = waiting;

    neighbour.push_back(i);
    neighbour_pos.push_back(i+1);
    
  }

  neighbour.push_back(-1);

  os << "\nMESSAGES FROM INITIALIZATION\n";
  if (io::messages.display(os) >= io::message::error){
    // exit
    os << "\nErrors during initialization!\n" << std::endl;
    os.flush();
    os.close();
    return 1;
  }
  
  os.precision(5);
  os.setf(std::ios_base::fixed, std::ios_base::floatfield);

  if (tid == 0){

    io::messages.clear();

    // store the positions of the all replicas
    m_conf.push_back(conf);

    // set seed if not set by environment variable
    // (get from sim.param()...)
    if (gsl_rng_default_seed == 0)
      gsl_rng_set (m_rng, 123456789);

    os << "starting master thread" << std::endl;

    os << "Replica Exchange\n"
       << "\treplicas:\t" << sim.param().replica.number << "\n"
       << "\tthreads:\t" << num_threads << "\n\n"
       << "\t" << std::setw(10) << "ID"
       << std::setw(20) << "Temp"
       << std::setw(20) << "lambda\n";
    
    for(int i=0; i<sim.param().replica.number; ++i){
      os << "\t" << std::setw(10) << i+1
	 << std::setw(20) << replica_data[i].temperature
	 << std::setw(20) << replica_data[i].lambda
	 << "\n";
    }
    os << "\n\ttrials:\t" << sim.param().replica.trials << "\n"
       << "\nEND\n\n";

    os.flush();
    
    int runs = 0;
    
    while(true){
      
      if(runs >= sim.param().replica.trials * sim.param().replica.number) break;
    
      // is a thread waiting?
      std::cerr << "master: selecting thread..." << std::endl;
      // master does not accept jobs...
      for(int i=1; i<num_threads; ++i){
	if (thread_state[i] == waiting){
	
	  std::cerr << "thread " << i << " is waiting..." << std::endl;
	  
	  // select a replica to run
	  for(int r=0; r<sim.param().replica.number; ++r){
	    
	    if(replica_data[r].state == waiting){
	      // try a switch
	      std::cerr << "replica " << r << " is waiting..." << std::endl;
	      switch_replica(r, replica_data[r].switch_replica);
	    }
	    
	    if(replica_data[r].state == ready){
	      std::cerr << "replica " << r << " is ready..." << std::endl;
	      
	      // assign it!
	      thread_replica[i] = r;
	      thread_state[i] = ready;
	      ++runs;
	      break;
	    }
	  } // replica selected
	} // thread waiting
      } // threads

      sleep(10);
    } // while
    
    // stop all threads
    for (int i=0; i<thread_state.size(); ++i){
      std::cerr << "terminating " << i << std::endl;
      thread_state[i] = terminate;
    }

  }
  else{

    os << "client " << tid << " waiting..." << std::endl;
    std::cerr << "client " << tid << " waiting..." << std::endl;

    assert(thread_state.size() > tid);
    while(thread_state[tid] != terminate){

      synchronise_thread_state(tid);

      if (thread_state[tid] == ready){
	// accept job
	std::cerr << "thread " << tid << " got a job!" << std::endl;
	
	synchronise_thread_replica(tid);
	synchronise_replica(thread_replica[tid]);
	
	assert(thread_replica.size() > tid);
	assert(replica_data.size() > thread_replica[tid]);
	
	std::cerr << "thread " << tid << " running replica " << thread_replica[tid]
		  << " at T=" << replica_data[thread_replica[tid]].temperature
		  << " and l=" << replica_data[thread_replica[tid]].lambda
		  << std::endl;
	os << "thread " << tid << " running replica " << thread_replica[tid]
	   << " at T=" << replica_data[thread_replica[tid]].temperature
	   << " and l=" << replica_data[thread_replica[tid]].lambda
	   << std::endl;
	
	sleep(10);
	
	++replica_data[thread_replica[tid]].run;
	replica_data[thread_replica[tid]].state = waiting;
	thread_state[tid] = waiting;
	
	update_replica(thread_replica[tid]);
	update_thread_state(tid);

      }
      else{
	std::cerr << "thread " << tid << " waiting (state " << thread_state[tid] << ")" << std::endl;
	sleep(5);
      }

    } // while
    

  }

  os.flush();
  os.close();

  return 0;
}

int util::Replica_Exchange::switch_replica(int i, int j)
{
  if (i < 0){
    assert(j>=0 && j<replica_data.size());
    if (replica_data[j].state != waiting) return 0;

    replica_data[j].switch_replica = ((replica_data[j].run % 2) ==  (neighbour_pos[j] % 2) ) ?
      neighbour[neighbour_pos[j] - 1] : neighbour[neighbour_pos[j] + 1];
    std::cerr << "run " << replica_data[j].run << " switch " 
	      << j << " - " << replica_data[j].switch_replica << std::endl;
    
    replica_data[j].state = ready;
  }
  else if(j < 0){
    assert(i>=0 && i<replica_data.size());
    if (replica_data[i].state != waiting) return 0;

    replica_data[i].switch_replica = ((replica_data[i].run % 2) ==  (neighbour_pos[i] % 2) ) ?
      neighbour[neighbour_pos[i] - 1] : neighbour[neighbour_pos[i] + 1];
    std::cerr << "run " << replica_data[i].run << " switch " 
	      << i << " - " << replica_data[i].switch_replica << std::endl;

    replica_data[i].state = ready;
  }
  else{
    // both finished? same number of runs?
    if (replica_data[i].state != waiting ||
	replica_data[j].state != waiting ||
	replica_data[i].run != replica_data[j].run) return 0;

    // try switch...

    replica_data[j].switch_replica = ((replica_data[j].run % 2) ==  (neighbour_pos[j] % 2) ) ?
      neighbour[neighbour_pos[j] - 1] : neighbour[neighbour_pos[j] + 1];
    replica_data[i].switch_replica = ((replica_data[i].run % 2) ==  (neighbour_pos[i] % 2) ) ?
      neighbour[neighbour_pos[i] - 1] : neighbour[neighbour_pos[i] + 1];

    std::cerr << "run " << replica_data[i].run << " switch " 
	      << i << " - " << replica_data[i].switch_replica << std::endl;
    std::cerr << "run " << replica_data[j].run << " switch " 
	      << j << " - " << replica_data[j].switch_replica << std::endl;

    replica_data[i].state = ready;
    replica_data[j].state = ready;
    
  }
  
  return 0;
}

int util::Replica_Exchange::synchronise_thread_state(int tid)
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  assert(thread_state.size() > tid);
  assert(replica_master->thread_state.size() > tid);
  
  thread_state[tid] = replica_master->thread_state[tid];
  std::cerr << "thread " << tid << " state after synch = " << thread_state[tid] << std::endl;
  
  return 0;
}

int util::Replica_Exchange::synchronise_thread_replica(int tid)
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  assert(thread_replica.size() > tid);
  assert(replica_master->thread_replica.size() > tid);
  
  thread_replica[tid] = replica_master->thread_replica[tid];
  return 0;
}

int util::Replica_Exchange::synchronise_replica(int r)
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  assert(replica_data.size() > r);
  assert(replica_master->replica_data.size() > r);
  
  replica_data[r].temperature = replica_master->replica_data[r].temperature;
  replica_data[r].lambda = replica_master->replica_data[r].lambda;
  replica_data[r].energy = replica_master->replica_data[r].energy;
  replica_data[r].switch_replica = replica_master->replica_data[r].switch_replica;
  replica_data[r].switch_temperature = replica_master->replica_data[r].switch_temperature;
  replica_data[r].switch_lambda = replica_master->replica_data[r].switch_lambda;
  replica_data[r].switch_energy = replica_master->replica_data[r].switch_energy;
  replica_data[r].TID = replica_master->replica_data[r].TID;
  replica_data[r].run = replica_master->replica_data[r].run;
  replica_data[r].state = replica_master->replica_data[r].state;
  
  return 0;
}

int util::Replica_Exchange::update_thread_state(int tid)
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  replica_master->thread_state[tid] = thread_state[tid];
  return 0;
}

int util::Replica_Exchange::update_replica(int r)
{
  // for OMP the master is set in a variable
  // MPI would use messages to 0 for that
  
  assert(replica_master != NULL);
  replica_master->replica_data[r].temperature = replica_data[r].temperature;
  replica_master->replica_data[r].lambda = replica_data[r].lambda;
  replica_master->replica_data[r].energy = replica_data[r].energy;
  replica_master->replica_data[r].switch_replica = replica_data[r].switch_replica;
  replica_master->replica_data[r].switch_temperature = replica_data[r].switch_temperature;
  replica_master->replica_data[r].switch_lambda = replica_data[r].switch_lambda;
  replica_master->replica_data[r].switch_energy = replica_data[r].switch_energy;
  replica_master->replica_data[r].TID = replica_data[r].TID;
  replica_master->replica_data[r].run = replica_data[r].run;
  replica_master->replica_data[r].state = replica_data[r].state;
  
  return 0;
}

