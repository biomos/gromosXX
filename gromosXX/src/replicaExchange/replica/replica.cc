/* 
 * File:   replica.cc
 * Author: wissphil, sriniker
 * 
 * Created on April 29, 2011, 2:06 PM
 * Modified June 18, 2021 - bschroed, srieder
 */


#ifdef XXMPI
    #include <mpi.h>
#endif


#include <io/argument.h>
#include <util/error.h>
#include <math/volume.h>
#include "replica.h"


#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica

re::replica::replica(io::Argument _args, int cont, int globalThreadID, simulation::MpiControl & replica_mpi_control) : replica_Interface(globalThreadID, replica_mpi_control, _args){
  // read input again. If copy constructors for topo, conf, sim, md work, one could
  // also pass them down from repex_mpi.cc ...
  
  DEBUG(4, "replica "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t START");

  // do continuation run?
  // change name of input coordinates
  if(cont == 1){
    std::multimap< std::string, std::string >::iterator it = args.lower_bound(("conf"));
    size_t pos = (*it).second.find_last_of(".");
    std::stringstream tmp;
    tmp << "_" << (simulationID+1);
    (*it).second.insert(pos, tmp.str());
  }

  // set output file
  std::stringstream tmp;
  tmp << "_" << (simulationID+1);
  std::string out;
  std::multimap< std::string, std::string >::iterator it = args.lower_bound(("repout"));
  size_t pos = (*it).second.find_last_of(".");
  (*it).second.insert(pos, tmp.str());
  os = new std::ofstream((*it).second.c_str());
  
  util::print_title(true, *os, true); // printing read in.

  
  // set trajectory
  std::stringstream trajstr;
  trajstr << GROMOSXX << "\n\tReplica Exchange with Replica ID " << (simulationID+1) << std::endl;
  std::string trajname = trajstr.str();

  traj = new io::Out_Configuration(trajname, *os);
    
    //sim.mpiControl() = ;
  if (io::read_input(args, topo, conf, sim, md, *os, true)) { 
    io::messages.display(*os);
    std::cerr << "\nErrors during initialization!\n" << std::endl;
#ifdef XXMPI
    MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
#endif
  }
  
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

  std::stringstream trajtitle;
  trajtitle << GROMOSXX << "\n" << sim.param().title << "\n\tReplica " << (simulationID+1) << "on Node " << globalThreadID;
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
  seed << sim.param().start.ig*(simulationID+1);
  rng = new math::RandomGeneratorGSL(seed.str(), -1);

  //init REPLICA ATTRIBUTES
  curentStepNumber = 0;
  stepsPerRun = sim.param().step.number_of_steps;
  totalStepNumber = stepsPerRun*(sim.param().replica.trials + sim.param().replica.equilibrate);
  
  *os << "==================================================\n"
      << " MAIN MD LOOP\n"
      << "==================================================\n\n";

    DEBUG(4, "replica "<< globalThreadID <<":Constructor:\t Temp of replica  "<< globalThreadID <<": " << simulationID << " \t" << sim.param().multibath.multibath.bath(0).temperature);
    DEBUG(4, "replica "<< globalThreadID <<":Constructor:\t replica Constructor  "<< globalThreadID <<": \t DONE");
    
     
    //TODO: HERE?
    // init MD simulation
    md.init(topo, conf, sim, *os, true);

}

void re::replica::init() {
  // init MD simulation
  md.init(topo, conf, sim, *os, true);
}

void re::replica::run_MD() {
  // run MD simulation
  int error = 0;
  DEBUG(4,  "replica "<< globalThreadID <<": run_MD:\t Start");      
  DEBUG(5, "replica "<< globalThreadID <<":run_MD:\t doing steps: "<<stepsPerRun<< " till: "<< stepsPerRun + curentStepNumber << " starts at: " << curentStepNumber << "TOTAL RUNS: "<< totalStepNumber );
  DEBUG(7, "replica "<<globalThreadID<<" BeforeSimulation "<<sim.time()<<": Conf Current Epot" << conf.current().energies.potential_total<< "\n");
  DEBUG(7, "replica "<<globalThreadID<<" BeforeSimulation "<<sim.time()<<": Conf OLD Epot" << conf.old().energies.potential_total<< "\n");

  if(sim.steps() == 0){ //initial write out
    traj->write(conf, topo, sim, io::reduced);
  }
  while ((unsigned int)(sim.steps()) < stepsPerRun + curentStepNumber ){
    // run a step
    DEBUG(5, "replica "<< globalThreadID <<":run_MD:\t simulation Step: "<< sim.steps());
    if ((error = md.run(topo, conf, sim))) {
      switch (error) {
        case E_SHAKE_FAILURE:
          std::cerr << "SHAKE FAILURE in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
          io::messages.display();
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        case E_SHAKE_FAILURE_SOLUTE:
          std::cerr << "SHAKE FAILURE SOLUTE in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
          io::messages.display();
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        case E_SHAKE_FAILURE_SOLVENT:
          std::cerr << "SHAKE FAILURE SOLVENT in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
          io::messages.display();
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        case E_NAN:
          std::cerr << "NAN error in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
          io::messages.display();
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
        default:
          std::cerr << "Unknown error in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
          io::messages.display();
#ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, error);
#endif
          break;
      }
      error = 0; // clear error condition
      break;
    }

    DEBUG(6, "replica "<< globalThreadID <<":run_MD:\t clean up:");      
    traj->print(topo, conf, sim);

    ++sim.steps();
    sim.time() = sim.param().step.t0 + sim.steps() * sim.time_step_size();
    //DEBUG
    DEBUG(7, "replica "<<globalThreadID<<" STEP"<< sim.steps() <<": Conf Current Epot" << conf.current().energies.potential_total<< "\n");
    DEBUG(7, "replica "<<globalThreadID<<" STEP"<< sim.steps() <<": Conf OLD Epot" << conf.old().energies.potential_total<< "\n");

    traj->write(conf, topo, sim, io::reduced);
  } // main md loop
  DEBUG(6, "replica "<< globalThreadID <<":run_MD:\t  DONE:");

  //DEBUG
  DEBUG(7,"replica "<<globalThreadID<<" AfterSimulation "<<sim.time()<<": Conf Current Epot" << conf.current().energies.potential_total<< "\n");
  DEBUG(7,"replica "<<globalThreadID<<" AfterSimulation "<<sim.time()<<": Conf OLD Epot" << conf.old().energies.potential_total<< "\n");


  //calculateEnergies(); //ONLY for Debugging
  DEBUG(7,"replica "<<globalThreadID<<" AfterReCAlc: Conf Current Epot" << conf.current().energies.potential_total<< "\n");
  DEBUG(7,"replica "<<globalThreadID<<" AfterReCAlc: Conf OLD Epot" << conf.old().energies.potential_total<< "\n");

  // print final data of run
  curentStepNumber = stepsPerRun + curentStepNumber;
  
  if (curentStepNumber ==  totalStepNumber) {
    traj->print_final(topo, conf, sim);
  }

}

double re::replica::calculateEnergies(){
    double energy = 0.0;
    algorithm::Algorithm * ff = nullptr;

    //Calculate energies
    ff = md.algorithm("Forcefield");
    DEBUG(5, "replica "<< globalThreadID <<":calculate_energy:\t calc energies");
    if (ff->apply(topo, conf, sim)) {
        std::cerr << "Error in Forcefield calculation in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
    #ifdef XXMPI
        MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
    #endif
    }

    //first calc energies for eds
    if(sim.param().reeds.reeds >0)
    {
        force_orig = conf.current().force;
        virial_tensor_orig = conf.current().virial_tensor;

        ff = md.algorithm("EDS");
        //Calculate energies
        DEBUG(5, "replica "<< globalThreadID <<":calculate_energy:\t calc energies");
        if (ff->apply(topo, conf, sim)) {
            std::cerr << "Error in EDS calculation in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
        #ifdef XXMPI
            MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
        #endif
        }

        conf.current().energies.calculate_totals();
        conf.current().force= force_orig;
        conf.current().virial_tensor= virial_tensor_orig;
        energy=conf.current().energies.eds_vr;
        conf.current().energies.calculate_totals();
    }
    else {
        conf.current().energies.calculate_totals();
        switch ( sim.param().xrayrest.replica_exchange_parameters.energy_switcher) {
            case simulation::energy_tot:
                energy =  conf.current().energies.potential_total +  conf.current().energies.special_total;
                break;
            case simulation::energy_phys:
                energy =  conf.current().energies.potential_total;
                break;
            case simulation::energy_special:
                energy =  conf.current().energies.special_total;
                break;
            default:
                std::cerr << ("Error in energy switching!");
                #ifdef XXMPI
                    MPI_Abort(0, E_UNSPECIFIED);
                #endif
        }
        conf.current().energies.calculate_totals();
    }
    DEBUG(5, "replica "<< globalThreadID <<":calculate_energy:\t energy: "<< energy);
    return energy;
}
