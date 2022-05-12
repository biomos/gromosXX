/**
 * @file distance_field_interaction.cc
 * template methods of Distance_Field_Interaction
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

#include "../../interaction/special/distance_field_interaction.h"

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



/**
 * a function to determine the six neighbouring gridpoints of gridpoint i
 */


int neighbour(int i, int j, std::vector<int> &ngrid){
  int n = 0;
  
  switch(j){
    
    case 0: 
      // backward in x
      if(i % ngrid[0] == 0) return i -1 + ngrid[0];
      else return i-1;
    case 1:
      // forward in x
      if((i+1) % ngrid[0] ==0) return i + 1 - ngrid[0];
      else return i+1;
    case 2:
      // backward in y
      n = i - (i % ngrid[0]);
      if(n % (ngrid[0]*ngrid[1]) == 0) return i - ngrid[0] + (ngrid[0]*ngrid[1]);
      else return i - ngrid[0];
    case 3:
      // forward in y
      n = i + ngrid[0] - (i% ngrid[0]);
      if(n % (ngrid[0]*ngrid[1]) == 0) return i + ngrid[0] - (ngrid[0]*ngrid[1]);
      else return i + ngrid[0];
    case 4:
      // backward in z
      n = i - (i % (ngrid[0]*ngrid[1]));
      if(n==0) return i - (ngrid[0]*ngrid[1]) + ngrid[0]*ngrid[1]*ngrid[2];
      else return i - ngrid[0] * ngrid[1];
    case 5:
      // forward in z
      n = i + ngrid[0]*ngrid[1] - (i % (ngrid[0]*ngrid[1]));
      if(n==(ngrid[0]*ngrid[1]*ngrid[2])) return i + ngrid[0]*ngrid[1] - (ngrid[0]*ngrid[1]*ngrid[2]);
      else return i + ngrid[0]*ngrid[1];
    default:
      return 0;
      
  }
  return 0;
}

double assignment_1d(const int &p, const double &xi) {
  double absf = fabs(xi);
  switch (p) {
    case 1:
      return absf < 0.5 ? 1 : 0;
    case 2:
      {
	return absf < 1.0 ? (1.0 - absf) : 0;
      }
    case 3:
      {
	if (absf < 0.5) {
	  return 0.75 - xi * xi;
	} else if (absf >= 0.5 && absf < 1.5) {
	  const double term = 2.0 * absf - 3.0;
	  return 0.125 * term * term;
	} else {
	  return 0.0;
	}
      }
    default:
      io::messages.add("assignment function not implemented", "Lattice Sum",
		       io::message::critical);
      return 0;
  }
}


/**
 * calculate distance field interactions
 */
template<math::boundary_enum B, math::virial_enum V>
static int _calculate_field_restraint_interactions
(topology::Topology & topo,
 configuration::Configuration & conf,
 simulation::Simulation & sim, 
 util::Virtual_Atom &va_i,
 util::Virtual_Atom &va_j)
{
  math::Vec v, f;
  double energy = 0.0;
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

  if(fabs(topo.disfield_restraints().r0 - dist) < sim.param().distancefield.r_l){
    energy = 0.5 * topo.disfield_restraints().K *
      (dist - topo.disfield_restraints().r0) * (dist - topo.disfield_restraints().r0);
    f = -topo.disfield_restraints().K * (dist - topo.disfield_restraints().r0) * deriv;
  }
  else if(dist < (topo.disfield_restraints().r0 - sim.param().distancefield.r_l)){
    energy = -topo.disfield_restraints().K * (dist - topo.disfield_restraints().r0 + 
      0.5 * sim.param().distancefield.r_l) * sim.param().distancefield.r_l;
    f = topo.disfield_restraints().K * sim.param().distancefield.r_l * deriv; 
  }
  else if(dist > (topo.disfield_restraints().r0 + sim.param().distancefield.r_l)){
    energy = topo.disfield_restraints().K * (dist - topo.disfield_restraints().r0 - 
      0.5 * sim.param().distancefield.r_l) * sim.param().distancefield.r_l;
    f = -topo.disfield_restraints().K * sim.param().distancefield.r_l * deriv;
  }

  // not so nice: we print it based on printout to its own block
  // this went to the special trajectory (distance) and energy trajectory
  /*  if(sim.param().print.stepblock && sim.steps() % sim.param().print.stepblock ==0){
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << "DISTANCEFIELD\n"
              << dist << "\t\t" << energy << "\n"
              << "END\n";
  }*/ 
  // created additional special energy, energies are no longer written to distance restraints
  conf.current().energies.disfieldres_energy[topo.atom_energy_group()
					     [va_i.atom(0)]] += energy;

  conf.special().distancefield.dist = dist;
  conf.special().distancefield.energy = energy;


  va_i.force(conf, topo, -f);
  va_j.force(conf, topo,  f);
  
  return 0;
  
}

int interaction::Distance_Field_Interaction
::calculate_interactions(topology::Topology &topo,
			 configuration::Configuration &conf,
			 simulation::Simulation &sim)
{
  // it could be that this interaction was created but that we only need
  // the perturbed version.
  if(!topo.disfield_restraints().on){
    return 0;
  }
  
  // check if the grid needs to be updated
  // slight shift of grid update for use in REMD
  if(sim.steps() % sim.param().distancefield.update == 0)
    SPLIT_BOUNDARY(_update_grid, topo, conf, sim, va_i);
  
  SPLIT_VIRIAL_BOUNDARY(_calculate_field_restraint_interactions,
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

  // if protection radius is 0, check if the start pos is in the protein
  if (sim.param().distancefield.protect==0) {
   for(unsigned int a =0; a <= topo.disfield_restraints().proteinatoms; a++){
    math::Vec ppos = conf.current().pos(a);
    periodicity.put_into_box(ppos);
    math::Vec d = startpos - ppos;
    if(math::abs(d) < sim.param().distancefield.proteincutoff){
      if(sim.steps()!=0){
	// we're in the protein, skip this update
	std::cerr << "SKIPPING distancefield update at timestep "
		  << sim.steps()
		  << ", because we're in the protein" << std::endl;
	return;
      }
      else{
        // it's the first step, initialize the grid to the reference value
        for(unsigned int j=0; j < distance.size(); j++){
          distance[j]=topo.disfield_restraints().r0;
        }
        std::cerr << "SKIPPING distancefield update at first timestep, "
		  << " because we're in "
		  << "the protein. Setting all distances to r0 (" 
		  << topo.disfield_restraints().r0 << ")" 
		  << std::endl;
        return;
      }

    }
   }
  }

  DEBUG(10, "DF UPDATE, start: " << start << " at " << dist_to_i);
  // reinitialize all points to a high value (4*(box+box+box))
  // and flag all gridpoints that are within the protein
  for(unsigned int j=0; j < distance.size(); j++)
    distance[j] = 4*(box(0,0)+box(1,1)+box(2,2));

  std::set<int> protein;
  DEBUG(10, "DF UPDATE, proteinatoms" << topo.disfield_restraints().proteinatoms);
  for(unsigned int a =0; a <= topo.disfield_restraints().proteinatoms; a++){
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

  va_i = topo.disfield_restraints().v1;
  va_j = topo.disfield_restraints().v2;

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

int interaction::Distance_Field_Interaction::init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os,
		     bool quiet) 
{
  DEBUG(10, "Distance field INIT");
  // it could be that this interaction was created but that we only need
  // the perturbed version.
  if(!topo.disfield_restraints().on){
    return 0;
  }
  
  SPLIT_BOUNDARY(_init_grid, conf, topo, sim, va_i, va_j);

  // the first update is called from the interaction function

  if(!quiet) {
    os << "DISTANCEFIELD\n"
       << "\tDistance field restraints\n"
       << "\t\tvirtual atom i (type " << topo.disfield_restraints().v1.type() << "): ";
    for(int i=0; i< va_i.size(); i++)
      os << va_i.atom(i)+1 << " ";
    os << "\n\t\tvirtual atom j (type " << topo.disfield_restraints().v2.type() << "): ";
    for(int i=0; i< va_j.size(); i++)
      os << va_j.atom(i)+1 << " ";
    os << "\n\tForce constant:\t\t" << topo.disfield_restraints().K << "\n"
       << "\tReference distance:\t" << topo.disfield_restraints().r0 << "\n"
       << "\tGrid spacing:\t\t" << sim.param().distancefield.grid << "\n"
       << "\tNumber of grid points:\t" << conf.special().distancefield.distance.size() << "\n"
       << "\tGrid is updated every\t" << sim.param().distancefield.update << " steps\n"
       << "\tGrid points within " << sim.param().distancefield.proteincutoff << " of atoms 1.." 
       << topo.disfield_restraints().proteinatoms+1 << " get a distance offset of " 
       << sim.param().distancefield.proteinoffset << "\n"
       << "END\n";
  }

  return 0;
  
}

