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
#include <util/replicaExchange/replica/_replica_Interface.h>

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

