/**
 * @file perturbed_distance_field_interaction.cc
 * template methods of Perturbed_Distance_Field_Interaction
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/perturbed_distance_field_interaction.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

template<math::boundary_enum B>
static void _update_grid
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim,
 util::Virtual_Atom &va_i);


int neighbour(int i, int j, std::vector<int> &ngrid);
double assignment_1d(const int &p, const double &xi); 

/**
 * calculate perturbed distance field interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_perturbed_field_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, 
 util::Virtual_Atom &va_i,
 util::Virtual_Atom &va_j)
{
  math::Vec v, f;
  double energy = 0.0, en_term = 0.0;
  std::vector<int> &ngrid = conf.special().distancefield.ngrid;
  math::Box &box = conf.current().box;
  double grid = sim.param().distancefield.grid;
  std::vector<double> &distance = conf.special().distancefield.distance;
  
  // in principle, the grid is correct for the periodicity
  math::Periodicity<B> periodicity(conf.current().box);
 
  math::Vec pos_j = va_j.pos(conf, topo);
  periodicity.put_into_box(pos_j);
  
  DEBUG(6, "DF CALC, position j  in box " << pos_j[0] << " " << pos_j[1] << " " << pos_j[2]);
  
  //pos(sim.param().distancefield.atom_j);
  
  // determine the nearest grid center
  // first we get the position of the particle in terms of grid coordinates (double)
  // we also get from that the coordinates of the lowest gridpoint of the eight surrounding points;
  math::Vec gpos_j;
  math::GenericVec<int> grid_j;
  
  for(int i=0; i<3; i++){
    gpos_j[i] = (pos_j[i] + box(i,i)/2)/grid;
    grid_j[i] = int(gpos_j[i]);
  }

  // fill an array with the indices of the eight neighbouring gridpoints
  std::vector<int> eightpoints(8);
  eightpoints[0] = grid_j[2] * ngrid[0] * ngrid[1] + grid_j[1] * ngrid[0] + grid_j[0];
  eightpoints[1] = neighbour(eightpoints[0], 1, ngrid);
  eightpoints[2] = neighbour(eightpoints[0], 3, ngrid);
  eightpoints[3] = neighbour(eightpoints[1], 3, ngrid);
  eightpoints[4] = neighbour(eightpoints[0], 5, ngrid);
  eightpoints[5] = neighbour(eightpoints[1], 5, ngrid);
  eightpoints[6] = neighbour(eightpoints[2], 5, ngrid);
  eightpoints[7] = neighbour(eightpoints[3], 5, ngrid);

  DEBUG(9, "DF CALC, eightpoints\n\t" << eightpoints[0] << " " << distance[eightpoints[0]] << "\n\t"
        << eightpoints[1] << " " << distance[eightpoints[1]] << "\n\t"
        << eightpoints[2] << " " << distance[eightpoints[2]] << "\n\t"
        << eightpoints[3] << " " << distance[eightpoints[3]] << "\n\t"
        << eightpoints[4] << " " << distance[eightpoints[4]] << "\n\t"
        << eightpoints[5] << " " << distance[eightpoints[5]] << "\n\t"
        << eightpoints[6] << " " << distance[eightpoints[6]] << "\n\t"
        << eightpoints[7] << " " << distance[eightpoints[7]]);

  // and one with their coordinates in grids
  std::vector<math::Vec> eightpos(8);
  eightpos[0][0] = double(grid_j[0]);
  eightpos[0][1] = double(grid_j[1]);
  eightpos[0][2] = double(grid_j[2]);
  eightpos[1] = eightpos[0] + math::Vec(1.0,0.0,0.0);
  eightpos[2] = eightpos[0] + math::Vec(0.0,1.0,0.0);
  eightpos[3] = eightpos[0] + math::Vec(1.0,1.0,0.0);
  eightpos[4] = eightpos[0] + math::Vec(0.0,0.0,1.0);
  eightpos[5] = eightpos[0] + math::Vec(1.0,0.0,1.0);
  eightpos[6] = eightpos[0] + math::Vec(0.0,1.0,1.0);
  eightpos[7] = eightpos[0] + math::Vec(1.0,1.0,1.0);

  // calculate the derivatives in the point. For this, we need to go even beyond
  // the points. Hmmm. Does this carry the risk of getting a force from the protein
  // at very long distances? But it has to be if we want to get a smooth surface.
  std::vector<math::Vec> eightderiv(8);
  for(unsigned int i=0; i< eightderiv.size(); i++){
    eightderiv[i][0] = (distance[neighbour(eightpoints[i], 1, ngrid)] - 
                        distance[neighbour(eightpoints[i], 0, ngrid)]);
    eightderiv[i][1] = (distance[neighbour(eightpoints[i], 3, ngrid)] - 
                        distance[neighbour(eightpoints[i], 2, ngrid)]); 
    eightderiv[i][2] = (distance[neighbour(eightpoints[i], 5, ngrid)] - 
                        distance[neighbour(eightpoints[i], 4, ngrid)]); 
    eightderiv[i] /= (2*grid);
    DEBUG(9, "DF CALC, 8 derivative: " << i << " : " << eightderiv[i][0] << " " << eightderiv[i][1] << " " << eightderiv[i][2]);
 }
  
  // assign the average distance using an assignment function of order 2.
  double dist=0;
  math::Vec deriv(0.0,0.0,0.0);
  
  for(unsigned int i=0; i< eightpoints.size(); i++){
    double P = assignment_1d(2,gpos_j[0] - eightpos[i][0]) * 
      assignment_1d(2,gpos_j[1] - eightpos[i][1]) *
      assignment_1d(2,gpos_j[2] - eightpos[i][2]);
    dist += distance[eightpoints[i]] * P;
    deriv += eightderiv[i] * P;
  }
  DEBUG(9, "DF CALC, assigned distance: " << dist);

  DEBUG(9, "DF CALC, derivative: " << deriv[0] << " " << deriv[1] << " " << deriv[2]);
  

  // the first atom of the first virtual atom 
  // determines the energy group for the output. 
  // we use the same definition for the individual lambdas
  const double l = topo.individual_lambda(simulation::disfield_lambda)
    [topo.atom_energy_group()[va_i.atom(0)]]
    [topo.atom_energy_group()[va_i.atom(0)]];
  const double l_deriv = topo.individual_lambda_derivative
    (simulation::disfield_lambda)
    [topo.atom_energy_group()[va_i.atom(0)]]
    [topo.atom_energy_group()[va_i.atom(0)]];


  const double r0 = (1-l) * topo.perturbed_disfield_restraints().A_r0 + l * topo.perturbed_disfield_restraints().B_r0;
  const double K = (1-l) * topo.perturbed_disfield_restraints().K_A + l * topo.perturbed_disfield_restraints().K_B;
  const double D_r0 = topo.perturbed_disfield_restraints().B_r0 - topo.perturbed_disfield_restraints().A_r0;
  const double D_K = topo.perturbed_disfield_restraints().K_B - topo.perturbed_disfield_restraints().K_A;
  const double r_l = sim.param().distancefield.r_l;
  const double n = topo.perturbed_disfield_restraints().n; 
  const double m = topo.perturbed_disfield_restraints().m; 

  double prefactor = pow(2.0, n + m) * pow(l, n) * pow(1.0-l, m);

  if(fabs(r0 - dist) < r_l){

    DEBUG(9, "PERTDISTANCERES  :  harmonic");
    f = -K * (dist - (r0)) * deriv;
    en_term = 0.5 * K * (dist - r0) * (dist - r0);
  }
  else if(dist < (r0 - r_l)){
    DEBUG(9, "PERTDISTANCERES  : (rep / attr) linear negative");
    f = K * r_l * deriv;
    en_term = -K * (dist - r0 + 0.5 * r_l) * r_l;
  }
  else if(dist > (r0 + r_l)){
    DEBUG(9, "PERTDISTANCERES   : (rep / attr) linear positive");
    f = -K * r_l * deriv;
    en_term = K * (dist - r0 - 0.5 * r_l ) * r_l;
  }
  
  energy = en_term * prefactor;
  f *= prefactor;

  // created an additional special energy, so energies are no longer put to distance restraints
  conf.current().energies.disfieldres_energy[topo.atom_energy_group()
                                             [va_i.atom(0)]] += energy;


  va_i.force(conf, topo, -f);
  va_j.force(conf, topo,  f);

  double dlam_term = 0.0, energy_derivative = 0.0, dprefndl = 0.0, dprefmdl = 0.0;

  // the derivative of the prefactor
  // division by zero precaution
  if (n==0) dprefndl = 0;
  else dprefndl = n * pow(l, n-1) * pow(1.0 - l, m);
  
  if (m == 0) dprefmdl = 0;
  else dprefmdl = m * pow(l, n) * pow(1.0 - l, m-1);

  double dprefdl = pow(2.0, m + n) * 
    (dprefndl - dprefmdl) * en_term;
  
  // calculate energy derivative
  if(fabs(r0 - dist) < r_l){
    DEBUG(9, "PERTDISTANCERES  :  harmonic");
    dlam_term = 0.5 * D_K * (dist - r0) * (dist - r0) - K * (dist - r0) * D_r0;
  }
  else if(dist < (r0 - r_l)){
    DEBUG(9, "PERTDISTANCERES  : (rep / attr) linear negative");
    dlam_term = -D_K * (dist - r0 + 0.5 * r_l) * r_l + K * D_r0 * r_l;
  }
  else if(dist > (r0 + r_l)){
    DEBUG(9, "PERTDISTANCERES   : (rep / attr) linear positive");
    dlam_term = D_K * (dist - r0 - 0.5 * r_l) * r_l - K * D_r0 * r_l;
  }

  double dpotdl = prefactor * dlam_term;

  energy_derivative = l_deriv * (dprefdl + dpotdl);

  // not so nice: we print it based on printout to its own block
  // this went to the energy trajectories and the special trajectory
  /*if(sim.param().print.stepblock && sim.steps() % sim.param().print.stepblock ==0){
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << "PERTDISFIELD dist energy energy_derivative\n"
              << dist << "\t\t" << energy << "\t\t" << energy_derivative << "\n"
              << "END\n";
	      }*/

  // created an additional special energy, so energies are no longer put to distance restraints
  conf.current().perturbed_energy_derivatives.disfieldres_energy[topo.atom_energy_group()
                                             [va_i.atom(0)]] += energy_derivative;

  conf.special().distancefield.dist = dist;
  conf.special().distancefield.energy = energy;
  conf.special().distancefield.energy_deriv = energy_derivative;


  /**
   * ext_TI code - Betty
   */
  
  if (sim.param().precalclam.nr_lambdas &&
      ((sim.steps() % sim.param().write.free_energy) == 0)){
    
    double lambda_step = (sim.param().precalclam.max_lam -
			  sim.param().precalclam.min_lam) /
                          (sim.param().precalclam.nr_lambdas-1);
    
    //loop over nr_lambdas
    for (unsigned int lam_index = 0; lam_index < sim.param().precalclam.nr_lambdas; ++lam_index){
      
      double lam = (lam_index * lambda_step) + sim.param().precalclam.min_lam;
      double prefactorlam = pow(2.0, n + m) * pow(lam, n) * pow(1.0-lam, m);
      
      double r0lam = (1-lam) * topo.perturbed_disfield_restraints().A_r0 + lam * topo.perturbed_disfield_restraints().B_r0;
      double Klam = (1-lam) * topo.perturbed_disfield_restraints().K_A + lam * topo.perturbed_disfield_restraints().K_B;
      double difflam = dist - r0lam;
      double difflam2 = difflam * difflam;

      double en_termlam = 0.0;
      if(fabs(r0lam - dist) < r_l){
	en_termlam = 0.5 * Klam * difflam;
      }
      else if(dist < (r0lam - r_l)){
	en_termlam = -Klam * (difflam + 0.5 * r_l) * r_l;
      }
      else if(dist > (r0lam + r_l)){
	en_termlam = Klam * (difflam - 0.5 * r_l ) * r_l;
      }
      double energylam = prefactorlam * en_termlam;
      
      double dlam_termlam = 0.0, dprefmdlam = 0.0, dprefndlam = 0.0;
      if (n==0) dprefndlam = 0;
      else dprefndlam = n * pow(lam, n-1) * pow(1.0 - lam, m);
      if (m == 0) dprefmdlam = 0;
      else dprefmdlam = m * pow(lam, n) * pow(1.0 - lam, m-1);
      double dprefdlam = pow(2.0, m + n) * (dprefndlam - dprefmdlam) * en_termlam;
      if(fabs(r0lam - dist) < r_l){
	dlam_termlam = 0.5 * D_K * difflam2 - Klam * difflam * D_r0;
      }
      else if(dist < (r0lam - r_l)){
	dlam_termlam = -D_K * (difflam + 0.5 * r_l) * r_l + Klam * D_r0 * r_l;
      }
      else if(dist > (r0lam + r_l)){
	dlam_termlam = D_K * (difflam - 0.5 * r_l) * r_l - Klam * D_r0 * r_l;
      }
      double dpotdlam = prefactorlam * dlam_termlam;
      double energy_derivativlam = dprefdlam + dpotdlam;
      
      conf.current().energies.AB_disfld[lam_index] += energylam;
      conf.current().perturbed_energy_derivatives.AB_disfld[lam_index] += 
	energy_derivativlam;
    }
  // Betty
  }
  return 0;
}

int interaction::Perturbed_Distance_Field_Interaction
::calculate_interactions(topology::Topology &topo,
                         configuration::Configuration &conf,
                         simulation::Simulation &sim)
{

  // it could be that this interaction was created but that we do not need
  // the perturbed version.
  if(!topo.perturbed_disfield_restraints().on){
    return 0;
  }

  // check if the grid needs to be updated
  if(sim.steps()  % sim.param().distancefield.update == 0)
    SPLIT_BOUNDARY(_update_grid, topo, conf, sim, va_i);
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_perturbed_field_restraint_interactions,
                        topo, conf, sim, va_i, va_j);
  
  return 0;
}


/**
 * update distance field grid
 */
template<math::boundary_enum B>
static void _update_grid
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation &sim,
 util::Virtual_Atom &va_i)
{
  math::Periodicity<B> periodicity(conf.current().box);
  std::vector<int> &ngrid = conf.special().distancefield.ngrid;
  std::vector<double> &distance = conf.special().distancefield.distance;
  
  double grid = sim.param().distancefield.grid;
  
  std::set<int> visited;
  std::set<int> todo;

  DEBUG(5, "DF UPDATE");

  // find the gridpoint that is nearest to atom_i. This will be our seed
  math::Box box = conf.current().box;
  int nx = 0, ny = 0, nz = 0;
  math::Vec ppos = va_i.pos(conf, topo);
  periodicity.put_into_box(ppos);

  DEBUG(6, "DF UPDATE, va_i: " << ppos[0] << " " << ppos[1] << " " << ppos[2]);

  nx = int(-(-ppos(0) - box(0,0)/2)/grid+0.5);
  ny = int(-(-ppos(1) - box(1,1)/2)/grid+0.5);
  nz = int(-(-ppos(2) - box(2,2)/2)/grid+0.5);
  math::Vec startpos(nx*grid - box(0,0)/2,ny*grid - box(1,1)/2,nz*grid - box(2,2)/2);
  int start=int(nx + ngrid[0]*ny + ngrid[0]*ngrid[1]*nz);
  double dist_to_i=math::abs(startpos-ppos);

  //std::cerr << "DISTANCEFIELD UPDATE at " << sim.steps() << " with lambda value " << topo.lambda() << std::endl;

  // if protection radius is 0, check if the start pos is in the protein
  if (sim.param().distancefield.protect==0) {
   for(unsigned int a =0; a <= topo.perturbed_disfield_restraints().proteinatoms; a++){
    math::Vec ppos = conf.current().pos(a);
    periodicity.put_into_box(ppos);
    math::Vec d = startpos - ppos;
    if(math::abs(d) < sim.param().distancefield.proteincutoff){
      if(sim.steps()!=0){
        // we're in the protein, skip this update
        std::cerr << "SKIPPING distancefield update at timestep " 
		  << sim.steps() << " with lambda value " << topo.lambda() 
		  << ", because we're in the protein" << std::endl;
        return;
      }
      else{
        // it's the first step, initialize the grid to the reference value
        const double l = topo.individual_lambda(simulation::disfield_lambda)
          [topo.atom_energy_group()[va_i.atom(0)]]
          [topo.atom_energy_group()[va_i.atom(0)]];
        const double r0 = (1-l) * topo.perturbed_disfield_restraints().A_r0 + l * topo.perturbed_disfield_restraints().B_r0;
        for(unsigned int j=0; j < distance.size(); j++){
          distance[j]=r0;
        }
        std::cerr << "SKIPPING distancefield update at first timestep "
		  << "with lambda value " << topo.lambda() 
		  << ", because we're in the protein. "
		  << "Setting all distances to r0 (" << r0 << ")" 
		  << std::endl;
        return;
      }
    }
   }
  }


  DEBUG(10, "DF UPDATE, start: " << start << " at " << dist_to_i);
  // reinitialize all points to a high value (4*(box+box+box))
  // and flag all gridpoints that are within the protein
  // except the ones that are within the protection radius
  for(unsigned int j=0; j < distance.size(); j++)
    distance[j] = 4*(box(0,0)+box(1,1)+box(2,2));

  std::set<int> protein;
  for(unsigned int a =0; a <= topo.perturbed_disfield_restraints().proteinatoms; a++){
    math::Vec ppos = conf.current().pos(a);
    periodicity.put_into_box(ppos);
    int nx_min=int(-( sim.param().distancefield.proteincutoff - ppos(0) - box(0,0)/2)/grid);
    int nx_max=int(-(-sim.param().distancefield.proteincutoff - ppos(0) - box(0,0)/2)/grid)+1;
    int ny_min=int(-( sim.param().distancefield.proteincutoff - ppos(1) - box(1,1)/2)/grid);
    int ny_max=int(-(-sim.param().distancefield.proteincutoff - ppos(1) - box(1,1)/2)/grid)+1;
    int nz_min=int(-( sim.param().distancefield.proteincutoff - ppos(2) - box(2,2)/2)/grid);
    int nz_max=int(-(-sim.param().distancefield.proteincutoff - ppos(2) - box(2,2)/2)/grid)+1;

    for(int ix=nx_min; ix < nx_max; ix++){
      for(int iy=ny_min; iy < ny_max; iy++){
        for(int iz=nz_min; iz < nz_max; iz++){
          math::Vec gpos(ix*grid - box(0,0)/2, iy*grid - box(1,1)/2, iz*grid - box(2,2)/2);
          math::Vec d = gpos - ppos;
	      math::Vec e = startpos - gpos;
          if(math::abs(d) < sim.param().distancefield.proteincutoff && math::abs(e) > sim.param().distancefield.protect){
            int j=int(ix + ngrid[0]*iy + ngrid[0]*ngrid[1]*iz);
            protein.insert(j);
          }
        }
      }
    }
  }

  DEBUG(10, "DF UPDATE, protein gridpoints: " << protein.size());

  int current=start;
  distance[current] = dist_to_i;
  
  if(protein.count(current)) 
    distance[current]+=sim.param().distancefield.proteinoffset;
  
  DEBUG(10, "DF UPDATE, current: " << current );
  
  while(visited.size() != distance.size()) {
    
    visited.insert(current);
    todo.erase(current);
    
    DEBUG(15, "DF UPDATE, visited: " << visited.size() << " current: " << current);

    double currdist = distance[current];
    
    // loop over the six neighbouring gridpoints    
    for(unsigned int i=0; i<6; i++){
      double newdist = currdist + grid;

      int j=neighbour(current, i, ngrid);
      DEBUG(15, "DF UPDATE, neighbour " << i << ": " << j);

      // we don't revisit grids
      if(!visited.count(j)){
          
        // check if it is in the protein
        if(protein.count(j) && !protein.count(current)){ 
          newdist += sim.param().distancefield.proteinoffset;
          DEBUG(10, "hit a protein gridpoint");
        }
        
        if(newdist < distance[j]){
          distance[j]=newdist;
          DEBUG(10, "updated");
        }
        DEBUG(10, "distance " << j << " " << distance[j]);
        // and j is a candidate for the next round
        todo.insert(j);
        
      }
    }
    DEBUG(15, "DF UPDATE, done neighbours")
    // find the next gridpoint that has not been visited yet
    int next=current;
    // ugly initialization of the minimum, distance was initialized to 2*(box+box+box)
    double min=4*(box(0,0)+box(1,1)+box(2,2))+1;
    for(std::set<int>::const_iterator iter=todo.begin(), to=todo.end();
        iter != to; iter++){
      if(distance[*iter] < min){
        next=*iter;
        min=distance[*iter];
      }
    }
    current=next;
  }

  // to avoid the artificialforce due to the protein, do one round of smoothening
  // the first layer of protein grid-points gets a new distance, based on
  // the nearest solvent (you should not be getting deeper 'into' the protein anyway)
  // this is dangerous if your protein is only one grid-point deep
  for(int s = 0; s < sim.param().distancefield.smooth; s++){
    DEBUG(5, "DF UPDATE, smoothening " << s);
    
    std::set<int> remove_from_protein;
    // loop over the protein gridpoints
    for(std::set<int>::const_iterator iter=protein.begin(), to = protein.end();
        iter!=to; ++iter){
      DEBUG(10, "DF smooth: " << *iter);
      
      // find its neighbours
      for(unsigned int i=0; i<6; i++){
        int j=neighbour(*iter, i, ngrid);
        DEBUG(10, "DF smooth neighbour: " << j);
        
        // check that it is in the solvent
        if(!protein.count(j)){
          DEBUG(10, "DF smooth neighbour in solvent " << distance[*iter] << " " << distance[j]);
          
          if(distance[*iter] > (distance[j] + grid)){
            // OK, this can be shortened
            distance[*iter]=distance[j] + grid;
            remove_from_protein.insert(*iter);
            DEBUG(10, "DF shortened");
          }
        }
      }
    }
    DEBUG(10, "number of gridpoints shortened " << remove_from_protein.size());
    
    for(std::set<int>::const_iterator iter=remove_from_protein.begin(), to = remove_from_protein.end();
        iter!=to; ++iter)
      protein.erase(*iter);
 
  }
  
  DEBUG(10, "DF UPDATE done");

#ifndef NDEBUG
  for(unsigned int j=0; j< distance.size(); j++)
    {    
      nz = int(double(j) / double(ngrid[0] * ngrid[1]));
      ny = int(double( j % (ngrid[0] * ngrid[1]) ) / double(ngrid[0]));
      nx = (j % (ngrid[0] * ngrid[1])) % ngrid[0];
      math::Vec gpos(nx*grid - box(0,0)/2, ny*grid - box(1,1)/2, nz*grid - box(2,2)/2);

      DEBUG(10,"dist " << gpos[0] << " " << gpos[1] << " " << gpos[2] << " " << distance[j]);
    }
  
#endif
  
  
}


/**
 * initiate distance field restraint interactions
 */
template<math::boundary_enum B>
static void _init_grid
(configuration::Configuration & conf,
 topology::Topology & topo,
 simulation::Simulation &sim, 
 util::Virtual_Atom &va_i,
 util::Virtual_Atom &va_j)
{
/*
  util::virtual_type ti = util::virtual_type(sim.param().distancefield.vtype_i);
  util::virtual_type tj = util::virtual_type(sim.param().distancefield.vtype_j);
  
  va_i = util::Virtual_Atom(ti, sim.param().distancefield.atom_i, 
                            sim.param().distancefield.dish,
                            sim.param().distancefield.disc);
  va_j = util::Virtual_Atom(tj, sim.param().distancefield.atom_j, 
                            sim.param().distancefield.dish,
                            sim.param().distancefield.disc); 
*/
  va_i = topo.perturbed_disfield_restraints().v1;
  va_j = topo.perturbed_disfield_restraints().v2;

  DEBUG(5, "DF Init, start");
  
  // this will only work for rectangular boxes
  conf.special().distancefield.ngrid.resize(3);
  
  for(int i=0; i<3; i++)
    conf.special().distancefield.ngrid[i] = int(conf.current().box(i,i) / sim.param().distancefield.grid)+1;
  DEBUG(10, "DF Init ngrid: " << conf.special().distancefield.ngrid[0] << " "<< conf.special().distancefield.ngrid[1] << " " << conf.special().distancefield.ngrid[2])
  const int ndim=conf.special().distancefield.ngrid[0] * 
    conf.special().distancefield.ngrid[1] *
    conf.special().distancefield.ngrid[2];

  DEBUG(10, "DF Init ndim:" << ndim);
  
  // we initialize the distances to something large
  conf.special().distancefield.distance.resize(ndim,4*(conf.current().box(0,0) + conf.current().box(1,1) + conf.current().box(2,2)));

  DEBUG(10, "DF Init distance resized");
  
}

int interaction::Perturbed_Distance_Field_Interaction::init(topology::Topology &topo, 
                     configuration::Configuration &conf,
                     simulation::Simulation &sim,
                     std::ostream &os,
                     bool quiet) 
{
  // it could be that this interaction was created but that we do not need
  // the perturbed version.
  if(!topo.perturbed_disfield_restraints().on){
    return 0;
  }

  DEBUG(10, "Perturbed distance field INIT");

  SPLIT_BOUNDARY(_init_grid, conf, topo, sim, va_i, va_j);

  // the first update is called from the interaction function

  if(!quiet) {
    os << "DISTANCEFIELD\n"
       << "\tPerturbed Distance field restraints\n"
       << "\t\tvirtual atom i (type " << topo.perturbed_disfield_restraints().v1.type() << "): ";
    for(int i=0; i< va_i.size(); i++)
      os << va_i.atom(i)+1 << " ";
    os << "\n\t\tvirtual atom j (type " << topo.perturbed_disfield_restraints().v2.type() << "): ";
    for(int i=0; i< va_j.size(); i++)
      os << va_j.atom(i)+1 << " ";
    os << "\n\tForce constant A:\t" << topo.perturbed_disfield_restraints().K_A << "\n"
       << "\tForce constant B:\t" << topo.perturbed_disfield_restraints().K_B << "\n"
       << "\tReference distance A:\t" << topo.perturbed_disfield_restraints().A_r0 << "\n"
       << "\tReference distance B:\t" << topo.perturbed_disfield_restraints().B_r0 << "\n"
       << "\tGrid spacing:\t\t" << sim.param().distancefield.grid << "\n"
       << "\tNumber of grid points:\t" << conf.special().distancefield.distance.size() << "\n"
       << "\tGrid is updated every\t" << sim.param().distancefield.update << " steps\n"
       << "\tGrid points within " << sim.param().distancefield.proteincutoff << " of atoms 1.." 
       << topo.perturbed_disfield_restraints().proteinatoms+1 << " get a distance offset of " 
       << sim.param().distancefield.proteinoffset << "\n"
       << "END\n";
  }
  return 0;
  
}

