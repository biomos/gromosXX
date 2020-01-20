/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   replica_Interface.cpp
 * Author: bschroed
 * 
 * Created on November 8, 2019, 12:53 PM
 */


#include <io/argument.h>
#include <util/error.h>
#include <math/volume.h>

#ifdef XXMPI
#include <mpi.h>
#endif
#include <util/replicaExchange/replica/replica_Interface.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE replica_exchange

util::replica_Interface::replica_Interface(int ID, int rank, io::Argument args): ID(ID), rank(rank), args(args){};

util::replica_Interface::~replica_Interface(){
    delete rng;
    delete traj;
    delete os;
}
void util::replica_Interface::init() {
  // init MD simulation
  md.init(topo, conf, sim, *os, true);
}

void util::replica_Interface::write_final_conf() {
    traj->write(conf, topo, sim, io::final);
    delete  traj;   //bschroed needs to be done, so the flush functions are triggered
}

// this function is deadlock safe, because it's made sure, that the replicas
// are not on the same node

void util::replica_Interface::swap(const unsigned int partnerID, const unsigned int partnerRank) {
  DEBUG(4, "replica "<< rank <<":swap:\t  START");

  partner = partnerID;
  unsigned int numT = sim.param().replica.num_T;
  unsigned int numL = sim.param().replica.num_l;

  // does partner exist?
  if (partner < numT * numL && partner != ID) {
    // the one with lower ID does probability calculation
    if (ID < partner) {
      // posts a MPI_Recv(...) matching the MPI_Send below 
      probability = calc_probability(partnerID, partnerRank);
      const double randNum = rng->get();

      std::vector<double> prob(2);
      prob[0] = probability;
      prob[1] = randNum;

#ifdef XXMPI
      MPI_Send(&prob[0], 2, MPI_DOUBLE, partnerRank, SENDCOORDS, MPI_COMM_WORLD);
#endif

      if (randNum < probability) {
        switched = true;
      } else
        switched = false;
    } else {
      bool sameLambda = (l == sim.param().replica.lambda[partner / sim.param().replica.num_T]);
      if(!sameLambda){
        // E21: Energy with configuration 2 and lambda 1(of partner)
        const double E21 = calculate_energy(partner);
        // this we can store as the partner energy of the current replica
        epot_partner = E21;
        // E22: Energy with configuration 2 and lambda 2(own one)
#ifdef XXMPI
        const double E22 = epot;
        // send E21 and E22
        double energies[2] = {E22, E21};
        //this send operation is matched in calc_probability()
        MPI_Send(&energies[0], 2, MPI_DOUBLE, partnerRank, SWITCHENERGIES, MPI_COMM_WORLD);
#endif
      } else { // sameLambda
#ifdef XXMPI
        double energies[2] = {epot, 0.0};
        MPI_Send(&energies[0],2,MPI_DOUBLE, partnerRank, SWITCHENERGIES, MPI_COMM_WORLD);
#endif
     }
      if (sim.param().pcouple.scale != math::pcouple_off) {
#ifdef XXMPI
        math::Box box_replica = conf.current().box;
        MPI_Send(&box_replica(0)[0], 1, MPI_BOX, partnerRank, BOX, MPI_COMM_WORLD);
#endif
      }
      
#ifdef XXMPI
      MPI_Status status;
#endif
      std::vector<double> prob;
      prob.resize(2);
#ifdef XXMPI
      MPI_Recv(&prob[0], 2, MPI_DOUBLE, partnerRank, SENDCOORDS, MPI_COMM_WORLD, &status);
#endif

      probability = prob[0];
      double randNum = prob[1];

      if (randNum < probability) {
        switched = true;
      } else {
        switched = false;
      }
    }

  } else {
    partner = ID;
    switched = false;
    probability = 0.0;
  }
    DEBUG(4, "replica "<< rank <<":swap:\t  DONE");

}

int util::replica_Interface::find_partner() const {
  unsigned int numT = sim.param().replica.num_T;
  unsigned int numL = sim.param().replica.num_l;
  unsigned int partner = ID;

  bool even = ID % 2 == 0;
  bool evenRow = (ID / numT) % 2 == 0;
  bool firstElement = (ID % numT == 0);
  bool lastElement = (ID % numT == numT - 1);
  bool numTeven = (numT % 2 == 0);

  // only 1D-RE ?
  if (numT == 1 || numL == 1) {
    if (run % 2 == 1) {
      if (even)
        partner = ID + 1;
      else
        partner = ID - 1;
    } else {
      if (even)
        partner = ID - 1;
      else
        partner = ID + 1;
    }
  } else { // 2D-RE
    // determine switch direction
    switch (run % 4) {
      case 0: // lambda direction 
        if (numTeven) {
          if (even)
            partner = ID - 1;
          else
            partner = ID + 1;
          if (firstElement || lastElement)
            partner = ID;
        } else {
          if (evenRow) {
            if (even)
              partner = ID - 1;
            else
              partner = ID + 1;
          } else {
            if (even)
              partner = ID + 1;
            else
              partner = ID - 1;
          }
          if (firstElement)
            partner = ID;
        }
        break;

      case 1: // temp-direction
        if (evenRow)
          partner = ID + numT;
        else
          partner = ID - numT;
        break;

      case 2: // lambda direction
        if (numTeven) {
          if (even)
            partner = ID + 1;
          else
            partner = ID - 1;
        } else {
          if (evenRow) {
            if (even)
              partner = ID + 1;
            else
              partner = ID - 1;
            if (lastElement)
              partner = ID;
          } else {
            if (even)
              partner = ID - 1;
            else
              partner = ID + 1;
          }
          if (lastElement)
            partner = ID;
        }
        break;
      case 3:// temp-direction
        if (evenRow)
          partner = ID - numT;
        else
          partner = ID + numT;
        break;
    }
  }
  // partner out of range ?
  if (partner > numT * numL - 1)
    partner = ID;

  return partner;
}

double util::replica_Interface::calc_probability(const int partner, const int partnerRank) {
  double delta;
//  bool sameLambda = (ID == sim.param().replica.lambda[partner / sim.param().replica.num_T]); // horizontal switch with same lambda?
  bool sameLambda = (l == sim.param().replica.lambda[partner / sim.param().replica.num_T]);
  const double b1 = 1.0 / (math::k_Boltzmann * T);
  const double b2 = 1.0 / (math::k_Boltzmann * sim.param().replica.temperature[partner % sim.param().replica.num_T]);
  
  if (sameLambda) {
    // use simple formula
    // get energy from the other partner
    double energies[2] = {0.0, 0.0};
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerRank, SWITCHENERGIES, MPI_COMM_WORLD, &status);
#endif

    epot_partner = energies[0];
    
    delta = (b1 - b2)*(epot_partner - epot); //*  (E21 - E11=
  } else {
    // 2D formula
    /*
     * E12: Energy with lambda from 1 and configuration from 2
     * delta = b1 * ( E22 - E11) - b2*(E21  - E12);
     * E22 and E12 needed from partner
     */
    double energies[2] = {0.0, 0.0};
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&energies[0], 2, MPI_DOUBLE, partnerRank, SWITCHENERGIES, MPI_COMM_WORLD, &status);
#endif
    const double E22 = energies[0];
    const double E12 = energies[1];

    const double E11 = epot;
    const double E21 = calculate_energy(partner);
    
    // store this as the partner energy 
    epot_partner = E21;

    // Chris: I think this is wrong
    // delta = b1 * (E22 - E11) - b2 * (E21 - E12);
    //std::cerr << "b1: " << b1 << " b2: " << b2 << std::endl;
    //std::cerr << "E11: " << E11 << " E22: " << E22 << std::endl;
    //std::cerr << "E21: " << E21 << " E12: " << E12 << std::endl;

    delta = b1 * (E12 - E11) - b2 * (E22 - E21);
  }
  
  // NPT? add PV term
  if (sim.param().pcouple.scale != math::pcouple_off) {
    math::Box box_partner = conf.current().box;
#ifdef XXMPI
    MPI_Status status;
    MPI_Recv(&box_partner(0)[0], 1, MPI_BOX, partnerRank, BOX, MPI_COMM_WORLD, &status);
#endif
    double V1 = math::volume(conf.current().box, conf.boundary_type);
    double V2 = math::volume(box_partner, conf.boundary_type);
    //std::cerr << "volumes: " << V1 << " " << V2 << std::endl;
    
    // isotropic!
    double pressure = (sim.param().pcouple.pres0(0,0)
              + sim.param().pcouple.pres0(1,1)
              + sim.param().pcouple.pres0(2,2)) / 3.0;
    delta += pressure * (b1 - b2) * (V2 - V1);
  }

  if (delta < 0.0)
    return 1.0;
  else
    return exp(-delta);
}

double util::replica_Interface::calculate_energy(const int partner) {
  DEBUG(4, "replica "<< rank <<":calculate_energy:\t  START");

  change_lambda(partner);

  double energy = 0.0;

  algorithm::Algorithm * ff = md.algorithm("Forcefield");

  if (ff->apply(topo, conf, sim)) {
    print_info("Error in energy calculation!");
 #ifdef XXMPI
    MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
#endif
    return 1;
  }

  conf.current().energies.calculate_totals();
  switch (sim.param().xrayrest.replica_exchange_parameters.energy_switcher) {
    case simulation::energy_tot:
      energy = conf.current().energies.potential_total + conf.current().energies.special_total;
      break;

    case simulation::energy_phys:
      energy = conf.current().energies.potential_total;
      break;
    case simulation::energy_special:
      energy = conf.current().energies.special_total;
      break;
    default:
      std::cerr << "Something is wrong in energy calculation" << std::endl;
      print_info("Error in energy calculation!");
#ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
#endif
      return 1;
  }

  set_lambda();
  DEBUG(4, "replica "<< rank <<":calculate_energy:\t  DONE");

  return energy;
}

double util::replica_Interface::calculate_energy() {
  double energy = 0.0;
  //chris: you do need to re-evaluate the energy, otherwise you get the energy of before the previous step
  algorithm::Algorithm * ff = md.algorithm("Forcefield");

  if (ff->apply(topo, conf, sim)) {
    print_info("Error in energy calculation!");
 #ifdef XXMPI
    MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
#endif
    return 1;
  }


  conf.current().energies.calculate_totals();
  switch (sim.param().xrayrest.replica_exchange_parameters.energy_switcher) {
    case simulation::energy_tot:
      energy = conf.current().energies.potential_total + conf.current().energies.special_total;
      break;

    case simulation::energy_phys:
      energy = conf.current().energies.potential_total;
      break;
    case simulation::energy_special:
      energy = conf.current().energies.special_total;
      break;
    default:
      print_info("Error in energy switching!");
#ifdef XXMPI
      MPI_Abort(MPI_COMM_WORLD, E_UNSPECIFIED);
#endif
  }
  return energy;
}

void util::replica_Interface::exchange_averages() {
  // after a swap the averages of current and old are exchanged and have to be switched back
  configuration::Average  dummy = conf.current().averages;
  conf.current().averages=conf.old().averages;
  conf.old().averages=dummy;
}

void util::replica_Interface::send_coord(const int receiverID, const int receiverRank) {
#ifdef XXMPI
  MPI_Send(&conf.current().pos[0][0], 1, MPI_VARRAY, receiverRank, POS, MPI_COMM_WORLD);
  MPI_Send(&conf.current().posV[0][0], 1, MPI_VARRAY, receiverRank, POSV, MPI_COMM_WORLD);
  MPI_Send(&conf.current().vel[0][0], 1, MPI_VARRAY, receiverRank, VEL, MPI_COMM_WORLD);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if (ID < partner) {
    MPI_Send(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, receiverRank, LATTSHIFTS, MPI_COMM_WORLD);
  } else {
    MPI_Send(&((*latticeTMP)[0][0]), 1, MPI_VARRAY, receiverRank, LATTSHIFTS, MPI_COMM_WORLD);
    delete latticeTMP;
    latticeTMP = NULL;
  }

  MPI_Send(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, receiverRank, STOCHINT, MPI_COMM_WORLD);
  MPI_Send(&conf.current().box(0)[0], 1, MPI_BOX, receiverRank, BOX, MPI_COMM_WORLD);

  std::vector<double> angles;
  angles.resize(3);
  angles[0] = conf.current().phi;
  angles[1] = conf.current().psi;
  angles[2] = conf.current().theta;

  MPI_Send(&angles[0], angles.size(), MPI_DOUBLE, receiverRank, ANGLES, MPI_COMM_WORLD);
  MPI_Send(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, receiverRank, DF, MPI_COMM_WORLD);
  //if ((int) ID > receiverID){
  //  conf.exchange_state();
  //}
#endif
}

void util::replica_Interface::receive_new_coord(const int senderID, const int senderRank) {
#ifdef XXMPI
  MPI_Status status;
  conf.exchange_state();
  MPI_Recv(&conf.current().pos[0][0], 1, MPI_VARRAY, senderRank, POS, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().posV[0][0], 1, MPI_VARRAY, senderRank, POSV, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().vel[0][0], 1, MPI_VARRAY, senderRank, VEL, MPI_COMM_WORLD, &status);

  // we need to store the lattice shifts from the previous configuration and to send them
  // if we are the second one to send the data
  if ((int) ID > senderID)
    latticeTMP = new math::VArray(conf.special().lattice_shifts);

  MPI_Recv(&conf.special().lattice_shifts[0][0], 1, MPI_VARRAY, senderRank, LATTSHIFTS, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().stochastic_integral[0][0], 1, MPI_VARRAY, senderRank, STOCHINT, MPI_COMM_WORLD, &status);
  MPI_Recv(&conf.current().box(0)[0], 1, MPI_BOX, senderRank, BOX, MPI_COMM_WORLD, &status);

  std::vector<double> angles;
  angles.resize(3);
  MPI_Recv(&angles[0], angles.size(), MPI_DOUBLE, senderRank, ANGLES, MPI_COMM_WORLD, &status);

  conf.current().phi = angles[0];
  conf.current().psi = angles[1];
  conf.current().theta = angles[2];

  MPI_Recv(&conf.special().distancefield.distance[0], conf.special().distancefield.distance.size(), MPI_DOUBLE, senderRank, DF, MPI_COMM_WORLD, &status);

  if ((int) ID > senderID)
    conf.exchange_state();
#endif
}

void util::replica_Interface::velscale(int i){ 
  double T1 = sim.param().replica.temperature[ID]; 
  double T2 = sim.param().replica.temperature[i];
  if (T1 != T2) {
    double factor = sqrt(T1/T2);
    for (unsigned int k = 0; k < topo.num_atoms(); ++k) {
      conf.current().vel(k) *= factor;
    }
  } 
}

void util::replica_Interface::set_lambda() {
  // change Lambda in simulation
  sim.param().perturbation.lambda = l;
  topo.lambda(l);
  // twice, to set old_lambda... STILL NEEDED
  topo.lambda(l);
  topo.update_for_lambda();
  // set corresponding timestep
  sim.param().step.dt = dt;
}

void util::replica_Interface::set_temp() {
  // change T in simulation
  sim.param().stochastic.temp = T;

  for (unsigned int i = 0; i < sim.multibath().size(); ++i) {
    assert(&sim.multibath()[i].temperature != 0);
    sim.multibath()[i].temperature = T;
  }
}

void util::replica_Interface::change_lambda(const unsigned int partner) {
  int idx;
  if (sim.param().replica.num_l == 1)
    idx = 0;
  else
    idx = partner / sim.param().replica.num_T;
  const double lambda = sim.param().replica.lambda[idx];
  const double dt = sim.param().replica.dt[idx];
  // change Lambda in simulation
  sim.param().perturbation.lambda = lambda;
  topo.lambda(lambda);
  // twice, to set old_lambda... STILL NEEDED
  topo.lambda(lambda);
  topo.update_for_lambda();
  // set corresponding timestep
  sim.param().step.dt = dt;
}

void util::replica_Interface::change_temp(const unsigned int partner) {
  int idx;
  if (sim.param().replica.num_T == 1)
    idx = 0;
  else
    idx = partner % sim.param().replica.num_T;
  const double temp = sim.param().replica.temperature[idx];
  // change T in simulation
  sim.param().stochastic.temp = temp;

  for (unsigned int i = 0; i < sim.multibath().size(); ++i) {
    assert(&sim.multibath()[i].temperature != 0);
    sim.multibath()[i].temperature = temp;
  }
}

void util::replica_Interface::print_info(std::string bla) const {
  std::cout << "\n" << bla << std::endl;
  std::cout << "#"
          << std::setw(5) << "ID"
          << " "
          << std::setw(7) << "partner"
          << std::setw(7) << "run"

          << std::setw(13) << "li"
          << std::setw(13) << "Ti"
          << std::setw(14) << "Epoti"
          << std::setw(13) << "lj"
          << std::setw(13) << "Tj"
          << std::setw(14) << "Epotj"
          << std::setw(13) << "p"
          << std::setw(4) << "s"
          << "\n";

  std::cout << std::setw(6) << (ID+1)
          << " "
          << std::setw(6) << partner
          << std::setw(6) << run
          << std::setw(13) << l
          << std::setw(13) << T
          << " "
          << std::setw(18) << epot
          << std::setw(13) << l
          << std::setw(13) << T
          << " "
          << std::setw(18) << epot_partner
          << std::setw(13) << probability
          << std::setw(4) << switched
          << std::endl;
}
