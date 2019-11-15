/* 
 * File:   replica.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 */



#include <util/replicaExchange/replica/replica.h>
#include <io/argument.h>
#include <util/error.h>
#include <math/volume.h>

#ifdef XXMPI
#include <mpi.h>
#endif


#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica::replica(io::Argument _args, int cont, int _ID, int _rank) : replica_Interface(_ID, _rank, _args){
  // read input again. If copy constructors for topo, conf, sim, md work, one could
  // also pass them down from repex_mpi.cc ...
  
  DEBUG(4, "replica "<< rank <<":Constructor:\t  "<< rank <<":\t START");
  
  // do continuation run?
  // change name of input coordinates
  if(cont == 1){
    std::multimap< std::string, std::string >::iterator it = args.lower_bound(("conf"));
    size_t pos = (*it).second.find_last_of(".");
    std::stringstream tmp;
    tmp << "_" << (ID+1);
    (*it).second.insert(pos, tmp.str());
  }

  // set output file
  std::stringstream tmp;
  tmp << "_" << (ID+1);
  std::string out;
  std::multimap< std::string, std::string >::iterator it = args.lower_bound(("repout"));
  size_t pos = (*it).second.find_last_of(".");
  (*it).second.insert(pos, tmp.str());
  os = new std::ofstream((*it).second.c_str());
  
  util::print_title(true, *os, true); // printing read in.

  // set trajectory
  std::stringstream trajstr;
  trajstr << GROMOSXX << "\n\tReplica Exchange with Replica ID " << (ID+1) << std::endl;
  std::string trajname = trajstr.str();

  traj = new io::Out_Configuration(trajname, *os);
  
  if (io::read_input(args, topo, conf, sim, md, *os, true)) { 
    io::messages.display(*os);
    std::cerr << "\nErrors during initialization!\n" << std::endl;
#ifdef XXMPI
    MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
#endif
  }
  
  // set some variables
  maxSteps = sim.param().step.number_of_steps;
  run = 0;
  total_runs = sim.param().replica.trials + sim.param().replica.equilibrate;
  partner = ID;
  time = sim.time();
  steps = 0;
  switched = 0;

  const int numT = sim.param().replica.num_T;

  T = sim.param().replica.temperature[ID % numT];
  l = sim.param().replica.lambda[ID / numT];
  dt = sim.param().replica.dt[ID / numT];

  set_lambda();
  set_temp();

  // manipulate trajectory files
  // just inserting ID: NAME_ID.cnf
  std::string fin;
  it = args.lower_bound(("fin"));
  pos = (*it).second.find_last_of(".");
  (*it).second.insert(pos, tmp.str());
  
  if (sim.param().write.position && args.count("trc") > 0) {
    it = args.lower_bound("trc");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if (sim.param().write.energy && args.count("tre") > 0) {
    it = args.lower_bound("tre");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if (sim.param().write.free_energy && args.count("trg") > 0) {
    it = args.lower_bound("trg");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if (sim.param().write.velocity && args.count("trv") > 0) {
    it = args.lower_bound("trv");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if ((sim.param().polarise.write || sim.param().jvalue.write || sim.param().xrayrest.write
          || sim.param().distanceres.write || sim.param().distancefield.write
          || sim.param().dihrest.write || sim.param().localelev.write
          || sim.param().electric.dip_write || sim.param().electric.cur_write
          || sim.param().addecouple.write || sim.param().nemd.write
          || sim.param().orderparamrest.write || sim.param().print.monitor_dihedrals ) && args.count("trs") > 0) {
    it = args.lower_bound("trs");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if (sim.param().write.force && args.count("trf") > 0) {
    it = args.lower_bound("trf");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if (sim.param().write.block_average && sim.param().write.energy && args.count("bae") > 0) {
    it = args.lower_bound("bae");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }
  if (sim.param().write.block_average && sim.param().write.free_energy && args.count("bag") > 0) {
    it = args.lower_bound("bag");
    pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());
  }

  // Chris: setting the title after init does not make much sense. The init function already prints it
  std::stringstream trajtitle;
  trajtitle << GROMOSXX << "\n" << sim.param().title << "\n\tReplica " << (ID+1) << "on Node " << rank;
  traj->title(trajtitle.str());

  traj->init(args, sim.param());
  
  *os << "\nMESSAGES FROM INITIALISATION\n";
  if (io::messages.display(*os) >= io::message::error) {
    *os << "\nErrors during initialization!\n" << std::endl;
#ifdef XXMPI
    MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
#endif
  }

  // random generator
  std::stringstream seed;
  seed << sim.param().start.ig*ID;
  rng = new math::RandomGeneratorGSL(seed.str(), -1);

  *os << "==================================================\n"
      << " MAIN MD LOOP\n"
      << "==================================================\n\n";

    DEBUG(4, "replica "<< rank <<":Constructor:\t Temp of replica  "<< rank <<": " << ID << " \t" << sim.param().multibath.multibath.bath(0).temperature);
    DEBUG(4, "replica "<< rank <<":Constructor:\t replica Constructor  "<< rank <<": \t DONE");
}

void util::replica::init() {
  // init MD simulation
  md.init(topo, conf, sim, *os, true);
}

void util::replica::run_MD() {
  // run MD simulation
  int error;
  sim.steps() = steps;
  sim.time() = time;
  while ((unsigned int)(sim.steps()) < maxSteps + steps) {
    DEBUG(5, "replica "<< rank <<":run_MD:\t Start");      
    traj->write(conf, topo, sim, io::reduced);
    // run a step
    DEBUG(5, "replica "<< rank <<":run_MD:\t simulation!:");
    if ((error = md.run(topo, conf, sim))) {
      switch (error) {
        case E_SHAKE_FAILURE:
          std::cerr << "SHAKE FAILURE in Replica " << (ID+1) << " on node " << rank << std::endl;
          io::messages.display();
          print_info("Info:");
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        case E_SHAKE_FAILURE_SOLUTE:
          std::cerr << "SHAKE FAILURE SOLUTE in Replica " << (ID+1) << " on node " << rank << std::endl;
          io::messages.display();
          print_info("Info:");
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        case E_SHAKE_FAILURE_SOLVENT:
          std::cerr << "SHAKE FAILURE SOLVENT in Replica " << (ID+1) << " on node " << rank << std::endl;
          io::messages.display();
          print_info("Info:");
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        case E_NAN:
          std::cerr << "NAN error in Replica " << (ID+1) << " on node " << rank << std::endl;
          io::messages.display();
          print_info("Info:");
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        default:
          std::cerr << "Unknown error in Replica " << (ID+1) << " on node " << rank << std::endl;
          io::messages.display();
          print_info("Info:");
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
      }
      error = 0; // clear error condition
      break;
    }
    
    
    DEBUG(5, "replica "<< rank <<":run_MD:\t clean up:");      
    traj->print(topo, conf, sim);

    ++sim.steps();
    sim.time() = sim.param().step.t0 + sim.steps() * sim.time_step_size();

  } // main md loop
  DEBUG(4, "replica "<< rank <<":run_MD:\t  DONE:");      
  // update replica information
  time = sim.time();
  steps = sim.steps();
  ++run;
  epot = calculate_energy();

  // print final data of run
  if (run ==  total_runs) {
    traj->print_final(topo, conf, sim);
  }
}
