/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_reeds.cc
 * Author: bschroed
 * 
 * Created on April 23, 2018, 3:25 PM
 */

#include <util/replicaExchange/replica_reeds.h>



#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_reeds::replica_reeds(io::Argument _args, int cont, int _ID, int _rank): replica(_args, cont, _ID, _rank){
  // read input again. If copy constructors for topo, conf, sim, md work, one could
  // also pass them down from repex_mpi.cc ...  
 DEBUG(4, "replica_reeds "<< rank <<":constructor:\t START");
  eds_para=sim.param().reeds.eds_para[ID];
  sim.param().eds = eds_para;
  
  const int num_s = sim.param().reeds.num_l;
  l == eds_para.s[0];
  DEBUG(5, "replica_reeds "<< rank <<":constructor:\t Temp of replica "<< rank <<" " << _ID << " \t" << sim.param().multibath.multibath.bath(0).temperature)
  DEBUG(5, "replica_reeds "<< rank <<":constructor:\t S of replica "<< rank <<" " << _ID << " \t" << sim.param().eds.s[0])
  DEBUG(5, "replica_reeds "<< rank <<":constructor:\t S size of replica "<< rank <<" " << _ID << " \t" << sim.param().eds.s.size())

  assert(0.0 <= sim.param().multibath.multibath.bath(0).temperature);
  DEBUG(4, "replica_reeds "<< rank <<":constructor:\t DONE"); 

}

void util::replica_reeds::reset_eds() {//only reset switched parameters of change_eds() function 
  sim.param().eds = eds_para;
  sim.param().step.dt = dt;
  conf.current().force=force_orig;
  conf.current().virial_tensor=virial_tensor_orig;
}

void util::replica_reeds::change_eds(const unsigned int partner){//only change parameters, which are needed for energy calculation i.e. 
  int idx;
  if (sim.param().reeds.num_l == 1)
    idx = 0;
  else
    idx = partner;
  sim.param().step.dt = sim.param().reeds.dt[idx];
  sim.param().eds=sim.param().reeds.eds_para[idx];
  force_orig = conf.current().force;
  virial_tensor_orig = conf.current().virial_tensor;
}

/*
 * calc_energy_eds_stat() is only used for statistical purposes in eds_stat()
 * In order to avoid any adjustment of the mpi communication and thus reducing the complexity, the 
 * energy_calculation and probability calculations from replica.cc are not adjusted to work 
 * for non-pairwise exchanges. Instead, calc_energy_eds_stat() calculates the potential energy
 * of the current configuration for a new smoothing parameter s.
 * The exchange probabilities can be calculated in a postprocessing step, using these energies
 * given in the energy_stat output files.
 */
double util::replica_reeds::calc_energy_eds_stat(double s){
    double old_dt;
    double old_s;
    double old_eds_vr;
    algorithm::Algorithm * ff;   
    if(sim.param().eds.eds){
          //to reset old state
          old_dt=sim.param().step.dt;
          old_s=sim.param().eds.s[0];
          old_eds_vr=conf.current().energies.eds_vr;
          force_orig = conf.current().force;
          virial_tensor_orig = conf.current().virial_tensor;
          //only temporary change
          sim.param().eds.s[0]=s;
          
          ff = md.algorithm("EDS");
    }
    else {
          print_info("eds_stat() i.e calc_energy_eds_stat() called for non EDS simulation!");
      #ifdef XXMPI
          MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
      #endif
    }
    
    //Calculate energies
    if (ff->apply(topo, conf, sim)) {
      print_info("Error in Forcefield energy calculation!");
     #ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
     #endif
      return 1;
    }
    
    double energy=conf.current().energies.eds_vr; 
    
    // reset old EDS state
    conf.current().energies.eds_vr=old_eds_vr;
    sim.param().eds.s[0] = old_s;
    sim.param().step.dt = old_dt;
    conf.current().force=force_orig;
    conf.current().virial_tensor=virial_tensor_orig;
    
    return energy;
}

double util::replica_reeds::calculate_energy(const int partner) {
    DEBUG(4, "replica_reeds "<< rank <<":calculate_energy:\t START"); 

    double energy = 0.0;
    algorithm::Algorithm * ff;   

    DEBUG(5, "replica_reeds "<< rank <<":calculate_energy:\t get Partner settings"); 
    if(partner!=ID) change_eds(partner);
    ff = md.algorithm("EDS");

    //Calculate energies    
    DEBUG(5, "replica_reeds "<< rank <<":calculate_energy:\t calc energies"); 
    if (ff->apply(topo, conf, sim)) {
      print_info("Error in Forcefield energy calculation!");
  #ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
  #endif
      return 1;
    }
    
    //return energies
    DEBUG(5, "replica_reeds "<< rank <<":calculate_energy:\t return energies"); 
    energy=conf.current().energies.eds_vr; 
    if(partner!=ID) reset_eds();
    
    DEBUG(4, "replica_reeds "<< rank <<":calculate_energy:\t DONE"); 
    return energy;
}

double util::replica_reeds::calculate_energy() {
    return calculate_energy(ID);
}