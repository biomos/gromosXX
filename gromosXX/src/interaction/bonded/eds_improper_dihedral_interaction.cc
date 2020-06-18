/**
 * @file eds_improper_dihedral_interaction.cc
 * template methods of eds_Improper_Dihedral_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// interactions
#include "../../interaction/interaction_types.h"
#include "improper_dihedral_interaction.h"
#include "eds_improper_dihedral_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE bonded

/**
 * calculate improper dihedral forces and energies.
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_eds_improper_interactions
( topology::Topology & topo,
  configuration::Configuration & conf,
  simulation::Simulation & sim,
  interaction::Improper_Dihedral_Interaction const & m_interaction)
{
  // this is repeated code from Improper_Dihedral_Interaction !!!

  DEBUG(5, "Multiple perturbed improper dihedral interaction");
  DEBUG(7, "using the improper dihedral interaction: " 
	<< m_interaction.name);
  DEBUG(7, std::setprecision(5));
  
  // loop over the angles
  std::vector<topology::multiple_perturbed_four_body_term_struct>::const_iterator i_it =
    topo.eds_perturbed_solute().improper_dihedrals().begin(),
    i_to = topo.eds_perturbed_solute().improper_dihedrals().end();

  math::VArray &pos   = conf.current().pos;
  math::VArray &force = conf.current().force;
  math::Vec rij, rkj, rlj, rkl, rmj, rnk, fi, fj, fk, fl;
  double dkj2, dkj, dmj2, dmj, dnk2, dnk, ip, q;
  double energy;


  // get beta
  assert(sim.param().multibath.multibath.bath(0).temperature != 0.0);
  const double beta = 1.0 / (sim.param().multibath.multibath.bath(0).temperature * math::k_Boltzmann);

  int numstates = sim.param().eds.numstates;
  std::vector<double> tmpe_vi;
  tmpe_vi.resize(numstates, 0.0);
  unsigned int itcount = 0; //keep track of the actual improper that we are dealing with

  DEBUG(7, "eds_pertubed_solute_size " << topo.eds_perturbed_solute().improper_dihedrals().size());
  DEBUG(7, "eds_improper_num " << sim.param().eds.numimpropers);
  assert(topo.eds_perturbed_solute().improper_dihedrals().size() == sim.param().eds.numimpropers);

  math::Periodicity<B> periodicity(conf.current().box);

  //calculate acceleration if necessary
  if (sim.param().eds.initimpsearch){

    
    for( ; i_it != i_to; ++i_it){

      std::vector<unsigned int> t_types = i_it->T_types;
      std::sort( t_types.begin(), t_types.end() );
      t_types.erase( unique( t_types.begin(), t_types.end() ), t_types.end() );
      double ebarrier = 0.0;

      // search the highest barrier
      for (unsigned int i = 0; i < t_types.size() - 1; i++){
        unsigned int j = i + 1;
        for ( ; j < t_types.size(); j++){
          // state A
          double Ka = topo.impdihedral_types()[t_types[i]].K;
          double q0a = topo.impdihedral_types()[t_types[i]].q0;

          // state B
          double Kb = topo.impdihedral_types()[t_types[j]].K;
          double q0b = topo.impdihedral_types()[t_types[j]].q0;

          // calculate intersection
          double kdiff = Ka - Kb;
          double intersect;
          if (kdiff != 0){
            double maxq0 =  std::max( q0a , q0b );
            double minq0 =  std::min( q0a , q0b );
            double b = -2 * q0a * Ka + 2 * q0b * Kb;
            double c = q0a * q0a * Ka - q0b * q0b * Kb;
            double root = pow( b * b - 4 * kdiff * c, 0.5 );
            double x_1 = ( -b + root )/ ( 2 * kdiff );
		        double x_2 = ( -b - root )/ ( 2 * kdiff );
            intersect = x_1 > minq0 && x_1 < maxq0 ? x_1 : x_2;
          }
          else
          {
            intersect = ( q0b * q0b - q0a * q0a ) / ( 2 * ( q0b - q0a ) );
          }

          // calculate energie barrier height
          double partA = -beta *  0.5 * Ka * (intersect - q0a) * (intersect - q0a);
          double partB = -beta *  0.5 * Kb * (intersect - q0b) * (intersect - q0b);
          DEBUG(7, "partA " << partA);
          DEBUG(7, "partB " << partB);
          double sum_prefactors = std::max(partA, partB) 
          + log(1 + exp(std::min(partA, partB) - std::max(partA, partB)));

          // apply eds hamiltonian
          double eds_e = -1.0 / beta * sum_prefactors;

          // get acceleration limits
          if (std::max(eds_e, ebarrier) != ebarrier){
            sim.param().eds.impaedslimit[itcount][0] = std::min(q0a, q0b);
            sim.param().eds.impaedslimit[itcount][1] = std::max(q0a, q0b);
          }
          ebarrier = std::max(eds_e, ebarrier);
        }
      }
      
      // compute required acceleration
      double target_barrier = sim.param().eds.impbar[itcount];
      if (ebarrier > target_barrier){
        double emin = ( 2 * target_barrier * ebarrier - ebarrier * ebarrier ) / ( 2 * target_barrier );
        sim.param().eds.impemaxs[itcount] += ebarrier;
        sim.param().eds.impemins[itcount] += emin;
      }
      itcount++;
    }
    sim.param().eds.initimpsearch = false;
    itcount = 0;
    i_it = topo.eds_perturbed_solute().improper_dihedrals().begin();
  } //end calculate acceleration if necessary

  // loop EDS dihedrals
  for( ; i_it != i_to; ++i_it){

    periodicity.nearest_image(pos(i_it->k), pos(i_it->j), rkj);
    periodicity.nearest_image(pos(i_it->i), pos(i_it->j), rij);
    periodicity.nearest_image(pos(i_it->k), pos(i_it->l), rkl);

    tmpe_vi.assign(numstates, 0.0);
  
    rmj = cross(rij, rkj);
    rnk = cross(rkj, rkl);
    
    dkj2 = abs2(rkj);
    dmj2 = abs2(rmj);
    dnk2 = abs2(rnk);
    dkj  = sqrt(dkj2);
    dmj  = sqrt(dmj2);
    dnk  = sqrt(dnk2);
    
    DEBUG(15,"dkj="<<dkj<<" dmj="<<dmj<<" dnk="<<dnk);
   

    assert(dmj != 0.0);
    assert(dnk != 0.0);

    ip = dot(rmj, rnk);
   
    double acs = ip / (dmj*dnk);
    if (acs > 1.0) {
      if (acs < 1.0 + math::epsilon) {
        acs = 1.0;
      } else {
        io::messages.add("improper dihedral",
                "acs > 1.0",
                io::message::critical);
      }
    }
    
    if (acs < -1.0) {
      if (acs > -1.0 - math::epsilon) {
        acs = -1.0;
      } else {
        io::messages.add("improper dihedral",
                "acs < -1.0",
                io::message::critical);
      }
    }
    
    q  = acos(acs);

    DEBUG(10, "zeta="<<q);
    
    ip = dot(rij, rnk);
    if(ip < 0) q *= -1.0;

    // first calculate individual energys 
    for (unsigned int state = 0; state < numstates; state++){
      double K = topo.impdihedral_types()[i_it->T_types[state]].K;
      double q0 = topo.impdihedral_types()[i_it->T_types[state]].q0;
      tmpe_vi[state] = 0.5 * K * (q-q0) * (q-q0);
    }


    // calculate EDS hamiltonian
    energy = -beta * tmpe_vi[0];
    for (unsigned int state = 1; state < numstates; state++){
      double part = -beta * tmpe_vi[state];
      energy = std::max(energy, part)
          + log(1 + exp(std::min(energy, part) - std::max(energy, part)));
    }

    // accelerate eds Hamiltonian
    double eds_energy = energy * -1.0/beta;
    double kfac, fkfac;
    if (sim.param().eds.edsimp == 1){
      if (eds_energy <= sim.param().eds.impemins[itcount] || 
       sim.param().eds.impaedslimit[itcount][0] < q ||  sim.param().eds.impaedslimit[itcount][1] > q ) {
        conf.current().energies.eds_it_vr = eds_energy;
      }
      else if (eds_energy >= sim.param().eds.impemaxs[itcount]) {
         conf.current().energies.eds_it_vr = eds_energy - 0.5 * (sim.param().eds.impemaxs[itcount] - sim.param().eds.impemins[itcount]);
      }
      else {
        double demix = eds_energy - sim.param().eds.impemins[itcount];
        kfac = 1.0 / (sim.param().eds.impemaxs[itcount] - sim.param().eds.impemins[itcount]);
        fkfac = 1.0 - kfac * demix;
        conf.current().energies.eds_it_vr = eds_energy - 0.5 * kfac * demix * demix;
      }
    }
    else{
      conf.current().energies.eds_it_vr = eds_energy;
    } // end accelerate
    
    // calculate forces
    for (unsigned int state = 0; state < numstates; state++){
      double K = topo.impdihedral_types()[i_it->T_types[state]].K;
      double q0 = topo.impdihedral_types()[i_it->T_types[state]].q0;
      double ki = -K * (q - q0) * dkj;
      double kl = -ki;

      const double kj1 = dot(rij, rkj) / dkj2 - 1.0;
      const double kj2 = dot(rkl, rkj) / dkj2;
    
      fi = ki * rmj;
      fl = kl * rnk;
      fj = kj1 * fi - kj2 * fl;
      fk = -1.0*(fi + fj + fl);

      if (sim.param().eds.edsimp == 1){
        const long double pi = exp(-beta * tmpe_vi[state] - energy);

        if (eds_energy <= sim.param().eds.impemins[itcount] || 
          eds_energy >= sim.param().eds.impemaxs[itcount] ||
          sim.param().eds.impaedslimit[itcount][0] < q ||  
          sim.param().eds.impaedslimit[itcount][1] > q ) {

            force(i_it->i) += pi * fi;
            force(i_it->j) += pi * fj;
            force(i_it->k) += pi * fk;
            force(i_it->l) += pi * fl;

            // if (V == math::atomic_virial){
            periodicity.nearest_image(pos(i_it->l), pos(i_it->j), rlj);

            for(int a=0; a<3; ++a)
	            for(int bb=0; bb<3; ++bb)
	              conf.current().virial_tensor(a, bb) += 
                  pi * rij(a) * fi(bb) +
	                pi * rkj(a) * fk(bb) +
	                pi * rlj(a) * fl(bb);
            // }
        }
        else {
            force(i_it->i) += pi * fi * fkfac;
            force(i_it->j) += pi * fj * fkfac;
            force(i_it->k) += pi * fk * fkfac;
            force(i_it->l) += pi * fl * fkfac;

            // if (V == math::atomic_virial){
            periodicity.nearest_image(pos(i_it->l), pos(i_it->j), rlj);

            for(int a=0; a<3; ++a)
	            for(int bb=0; bb<3; ++bb)
	              conf.current().virial_tensor(a, bb) += 
                  pi * rij(a) * fi(bb) * fkfac +
	                pi * rkj(a) * fk(bb) * fkfac +
	                pi * rlj(a) * fl(bb) * fkfac;
            // }
        }
      }
      else if (sim.param().eds.edsimp == 2){
        // save forces and virial in eds data holders
        conf.current().energies.eds_vi[state] = tmpe_vi[state];

        conf.special().eds.force_endstates[state](i_it->i) += fi;
        conf.special().eds.force_endstates[state](i_it->j) += fj;
        conf.special().eds.force_endstates[state](i_it->k) += fk;
        conf.special().eds.force_endstates[state](i_it->l) += fl;

        // if (V == math::atomic_virial){
        periodicity.nearest_image(pos(i_it->l), pos(i_it->j), rlj);

        for(int a=0; a<3; ++a)
	        for(int bb=0; bb<3; ++bb)
	          conf.special().eds.virial_tensor_endstates[state](a, bb) += 
              rij(a) * fi(bb) +
	            rkj(a) * fk(bb) +
	            rlj(a) * fl(bb);
        // }
      }
      else {
        io::messages.add("EDS improper dihedral","functional form unknow", io::message::critical);
      }
    } // end calculate forces

  } // end loop EDS dihedrals

  return 0;
  
}

int interaction::EDS_Improper_Dihedral_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  m_timer.start();
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_eds_improper_interactions,
			topo, conf, sim, m_interaction);

  m_timer.stop();

  return 0;
}

// int interaction::EDS_Improper_Dihedral_Interaction
//     ::init(topology::Topology &topo, 
//              configuration::Configuration &conf,
// 		     simulation::Simulation &sim,
// 		     std::ostream &os = std::cout,
// 		     bool quiet = false)
// {
//     //m_timer.start();
  
//     // loop over the angles
//     std::vector<topology::multiple_perturbed_four_body_term_struct>::const_iterator i_it =
//         topo.eds_perturbed_solute().improper_dihedrals().begin(),
//         i_to = topo.eds_perturbed_solute().improper_dihedrals().end();

//     for( ; i_it != i_to; ++i_it){
//         // Iterate over the improper dihedrals
//         DEBUG(7, "improper dihedral " << i_it->i << "-" << i_it->j << "-" << i_it->k << "-" << i_it->l);
//         std::vector<unsigned int>::const_iterator i_it2 = i_it->T_types.begin(), i_to2 = i_it->T_types.end();
//         // iterate over the end states dihedral types

//     }
//     //m_timer.stop();

//     return 0;

// }