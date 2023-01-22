/* 
 * File:   replica_MPI_master.cpp
 * Author: bschroed
 * 
 * Created on November 8, 2019, 12:48 PM
 * Modified June 18, 2021 - bschroed, srieder
 */

#include "replica_MPI_master.h"
#include <io/argument.h>
#include <util/error.h>
#include <util/debug.h>
#include <math/volume.h>

#ifdef XXMPI
    #include <mpi.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE re
#define SUBMODULE replica

re::replica_MPI_Master::replica_MPI_Master(io::Argument _args, int cont,  int globalThreadID,
        simulation::MpiControl & replica_mpi_control) : replica_Interface( globalThreadID, replica_mpi_control, _args){

    /**
     * Build up replica - reads in the input again and build output files.
     * 
     * @param _args
     * @param cont
     * @param _ID
     * @param _globalThreadID
     * @param simulation_globalThreadID
     * @param simulationID
     * @param simulation_num_threads
     */
    MPI_DEBUG(5, "replica_MPI_MASTER "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t START");
#ifdef XXMPI
    /**
     * READ INPUT
     */
    // do continuation run?
    // change name of input coordinates
    if(cont == 1){
      std::multimap< std::string, std::string >::iterator it = args.lower_bound(("conf"));
      size_t pos = (*it).second.find_last_of(".");
      std::stringstream tmp;
      tmp << "_" << (simulationID+1);
      (*it).second.insert(pos, tmp.str());
    }
    
    DEBUG(6, "replica_MPI_MASTER "<< globalThreadID <<":Constructor:\t  "<< simulationID <<":\t start read in");
    //Build structure
    sim.mpi = true;
    sim.mpiControl() = replica_mpi_control;  //build MPI parallelism
    
    if (io::read_input(args, topo, conf, sim, md, *os, true)) { 
      io::messages.display(*os);
      std::cerr << "\nErrors during initialization!\n" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
    }
    DEBUG(6, "replica_MPI_MASTER "<< globalThreadID <<":Constructor:\t  "<< globalThreadID <<":\t REad in input already");

    /**
     * INIT OUTPUT
     */
    // set output file
    std::stringstream tmp;
    tmp << "_" << (simulationID+1);
    std::string out;
    std::multimap< std::string, std::string >::iterator it = args.lower_bound(("repout"));
    size_t pos = (*it).second.find_last_of(".");
    (*it).second.insert(pos, tmp.str());

    os = new std::ofstream((*it).second.c_str());
    MPI_Barrier(MPI_COMM_WORLD); // give signal to child processes that output file was created

    util::print_title(true, *os, true); // printing read in.
    
    // init md
    md.init(topo, conf, sim, *os, true);

    
    /**
     * MASTER SPECIFICS:
     */

    // set trajectory
    std::stringstream trajstr;
    trajstr << GROMOSXX << "\n\tReplica Exchange with Replica ID " << (simulationID+1) << std::endl;
    std::string trajname = trajstr.str();

    traj = new io::Out_Configuration(trajname, *os);
    
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
      MPI_Abort(MPI_COMM_WORLD, E_INPUT_ERROR);
    }

    // random generator
    std::stringstream seed;
    seed << sim.param().start.ig*(simulationID+1);
    rng = new math::RandomGeneratorGSL(seed.str(), -1);

    *os << "==================================================\n"
        << " MAIN MD LOOP\n"
        << "==================================================\n\n";
    

    DEBUG(6, "replica_MPI_MASTER "<< globalThreadID <<":Constructor:\t Temp of replica  "<< globalThreadID <<": " << simulationID << " \t" << sim.param().multibath.multibath.bath(0).temperature);
    MPI_DEBUG(5, "replica_MPI_MASTER "<< globalThreadID <<":Constructor:\t replica Constructor  "<< globalThreadID <<": \t DONE");
#else
    throw "Can not construct Replica_MPI as MPI is not enabled!";
#endif    
}

re::replica_MPI_Master::~replica_MPI_Master() {
  delete rng;
  delete traj;
  delete os;
}

void re::replica_MPI_Master::run_MD(){
    MPI_DEBUG(5, "replica_MPI_Master "<< globalThreadID <<":runMD:\t thread  "<< globalThreadID <<": \t START");
    #ifdef XXMPI

    // run MD simulation
    int error;
    
    //next_stepf for mpi slaves  
    int next_step = 1;  //bool that signalises if next step is fine.
    //after an Replica coordinate exchange, update the coordinates of the slaves
    send_coordinates();
    MPI_DEBUG(5, "replica_MPI_Master "<< globalThreadID <<":runMD:\t\t sent Coords");
    MPI_DEBUG(5, "replica_MPI_Master "<< globalThreadID <<":runMD:\t\t steps: current step: "<<sim.steps()<< "  totalsteps: "<< stepsPerRun << " + " << curentStepNumber << " + 1 = "<< stepsPerRun+curentStepNumber+1);
    DEBUG(7,  "replica "<<globalThreadID<<" BeforeSimulation: Conf Current Epot" << conf.current().energies.potential_total<< "\n");
    DEBUG(7,  "replica "<<globalThreadID<<" BeforeSimulation: Conf OLD Epot" << conf.old().energies.potential_total<< "\n");
    DEBUG(7,  "replica "<<globalThreadID<<" BeforeSimulation: Conf cur CNF" << conf.current().pos[0][0] << "\t" << conf.current().pos[0][1] << "\t"  << conf.current().pos[0][2]<< "\n");

    if(sim.steps() == 0){ //initial write out
        traj->write(conf, topo, sim, io::reduced);
    }
    while ((unsigned int)(sim.steps()) <  stepsPerRun + curentStepNumber) {
      MPI_DEBUG(6, "replica_MPI_SLAVE " << globalThreadID << ":run_MD:\t step: "<< sim.steps() << " \tmaximal \t" << curentStepNumber+stepsPerRun);

      // run a step
      DEBUG(5, "replica_MPI_MASTER "<< globalThreadID <<":run_MD:\t simulation!:");
        if ((error = md.run(topo, conf, sim))) {
            DEBUG(5, "replica_MPI_MASTER "<< globalThreadID <<":run_MD:\t ERROER!:");
            switch (error) {
              case E_SHAKE_FAILURE:
                std::cerr << "SHAKE FAILURE in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
                io::messages.display();
                MPI_Abort(MPI_COMM_WORLD, error);
                break;
              case E_SHAKE_FAILURE_SOLUTE:
                std::cerr << "SHAKE FAILURE SOLUTE in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
                io::messages.display();
                MPI_Abort(MPI_COMM_WORLD, error);
                break;
              case E_SHAKE_FAILURE_SOLVENT:
                std::cerr << "SHAKE FAILURE SOLVENT in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
                io::messages.display();
                MPI_Abort(MPI_COMM_WORLD, error);
                break;
              case E_NAN:
                std::cerr << "NAN error in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
                io::messages.display();
                MPI_Abort(MPI_COMM_WORLD, error);
                break;
              default:
                std::cerr << "Unknown error in Replica " << (simulationID+1) << " on node " << globalThreadID << std::endl;
                io::messages.display();
                MPI_Abort(MPI_COMM_WORLD, error);
                break;
            }
            error = 0; // clear error condition
            std::cout << "\nError during MD run!\n" << std::endl;
            // send error status to slaves
            next_step = 0;
            std::cout << "Telling slaves to quit." << std::endl;
            break;
        }
        
        // tell the slaves to continue
        MPI_Bcast(&next_step, 1, MPI::INT, sim.mpiControl().masterID, sim.mpiControl().comm);
           
        traj->print(topo, conf, sim);
        ++sim.steps();
        sim.time() = sim.param().step.t0 + sim.steps() * sim.time_step_size();
        //DEBUG
        DEBUG(7,  "replica "<<globalThreadID<<" STEP"<< sim.steps() <<": Conf Current Epot" << conf.current().energies.potential_total<< "\n");
        DEBUG(7,  "replica "<<globalThreadID<<" STEP"<< sim.steps() <<": Conf OLD Epot" << conf.old().energies.potential_total<< "\n");

        traj->write(conf, topo, sim, io::reduced);

    } // main md loop
    MPI_DEBUG(5, "replica_MPI_MASTER "<< globalThreadID <<":run_MD:\t after step while");
    //DEBUG
    DEBUG(7, "replica "<<globalThreadID<<" AfterSimulation: Conf Current Epot" << conf.current().energies.potential_total<< "\n");
    DEBUG(7,  "replica "<<globalThreadID<<" AfterSimulation: Conf OLD Epot" << conf.old().energies.potential_total<< "\n");
    DEBUG(7,  "replica "<<globalThreadID<<" AfterSimulation: Conf cur CNF" << conf.current().pos[0][0] << "\t" << conf.current().pos[0][1] << "\t"  << conf.current().pos[0][2]<< "\n");

    //calculateEnergies(); // ONLY for DEBUGGING
    //DEBUG(7, "replica "<<globalThreadID<<" AfterReCAlc: Conf Current Epot" << conf.current().energies.potential_total<< "\n");
    //DEBUG(7,  "replica "<<globalThreadID<<" AfterReCAlc: Conf OLD Epot" << conf.old().energies.potential_total<< "\n");

    curentStepNumber +=  stepsPerRun;

    #endif
    MPI_DEBUG(5, "replica_MPI_MASTER "<< globalThreadID <<":run_MD:\t  DONE: at step= " << curentStepNumber);
}

void re::replica_MPI_Master::send_coordinates(){
  #ifdef XXMPI
  DEBUG(4, "replica_MPI_Master " << globalThreadID << " ::send_coordinates::\t START");
  //EXCHANGE conf parts
  MPI_Bcast(&conf.current().pos[0][0], 1, MPI_VARRAY, replica_mpi_control.masterID , replica_mpi_control.comm);
  MPI_Bcast(&conf.current().posV[0][0], 1, MPI_VARRAY, replica_mpi_control.masterID , replica_mpi_control.comm);
  MPI_Bcast(&conf.current().vel[0][0], 1, MPI_VARRAY, replica_mpi_control.masterID , replica_mpi_control.comm);

  MPI_Bcast(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, replica_mpi_control.masterID, replica_mpi_control.comm);
  MPI_Bcast(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, replica_mpi_control.masterID, replica_mpi_control.comm);
  MPI_Bcast(&conf.current().box(0)[0], 1, MPI_BOX, replica_mpi_control.masterID, replica_mpi_control.comm);

  //EXCHANGE ANGLES
  std::vector<double> angles;
  angles.resize(3);
  angles[0] = conf.current().phi;
  angles[1] = conf.current().psi;
  angles[2] = conf.current().theta;

  MPI_Bcast(&angles[0], angles.size(), MPI_DOUBLE, replica_mpi_control.masterID, replica_mpi_control.comm);
  
  //Exchange STUFF
  MPI_Bcast(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, replica_mpi_control.masterID, replica_mpi_control.comm);
  
  DEBUG(4, "replica_MPI_Master " << globalThreadID << " ::send_coordinates::\t DONE");
  #endif
}

double re::replica_MPI_Master::calculateEnergies(){
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
                            std::cerr << "Error in energy switching!";
                #ifdef XXMPI
                            MPI_Abort(0, E_UNSPECIFIED);
                #endif
        }
    }
    return energy;
}
