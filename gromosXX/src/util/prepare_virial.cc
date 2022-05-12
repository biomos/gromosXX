/**
 * prepare for virial calculation.
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"

#include "../math/periodicity.h"

#include "prepare_virial.h"

#include "../util/template_split.h"
#include "../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util


/**
 * @TODO should be changed according to phil's plan
 * of following diffusive particles.
 */
template<math::boundary_enum b>
static void _centre_of_mass(topology::Atom_Iterator start, 
			    topology::Atom_Iterator end,
			    topology::Topology const & topo,
			    configuration::Configuration const & conf,
			    math::Vec &com_pos, 
			    math::Matrix &com_e_kin,
			    math::Periodicity<b> const & periodicity,
                               simulation::Simulation const & sim)
{

  com_pos = 0.0;
  double m = 0.0;
  double tot_mass = 0.0;

  math::Vec p;
  math::Vec prev;
  math::Vec v(0.0);

  prev = conf.current().pos(*start);
  //scale the velocity for the adiabatic decoupling
  double scale_vel=1;
  int addc_index = 0;
  for( ; start != end; ++start){

    assert(unsigned(topo.mass().size()) > *start &&
           unsigned(conf.current().pos.size()) > *start);
    addc_index = sim.param().addecouple.check_index_adc(*start);
    scale_vel=1;
    if(addc_index!=-1)
      scale_vel=sqrt(1/
              sim.param().addecouple.adc_index()[addc_index].st);
    
    m = topo.mass()(*start);
    tot_mass += m;
    periodicity.nearest_image(conf.current().pos(*start), prev, p);
    com_pos += m * (p + prev);
    v += m * conf.current().vel(*start)*scale_vel;
    prev += p;
  }

  com_pos /= tot_mass;  

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      com_e_kin(i,j) = 0.5 * v(i) * v(j) / tot_mass;

}


template<math::boundary_enum b>
static void _prepare_virial(topology::Topology const & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation const & sim)
{
  if (sim.param().pcouple.virial == math::molecular_virial){

    DEBUG(10, "lambda = " << topo.lambda());
    
    math::Periodicity<b> periodicity(conf.current().box);

    topology::Pressuregroup_Iterator
      pg_it = topo.pressure_group_begin(),
      pg_to = topo.pressure_group_end();

    math::Vec com_pos;
    math::Matrix com_ekin;

    conf.current().kinetic_energy_tensor = 0.0;

    for( ; pg_it != pg_to; ++pg_it){
      _centre_of_mass(pg_it.begin(),
		      pg_it.end(),
		      topo, conf,
		      com_pos, com_ekin,
		      periodicity, sim);

      conf.current().kinetic_energy_tensor += com_ekin;

      /*
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  conf.current().kinetic_energy_tensor(i,j) += com_ekin(i,j);
      */

    }
  }
  
  else if (sim.param().pcouple.virial == math::atomic_virial){

    conf.current().kinetic_energy_tensor = 0.0;
    //scale the virial for adiabatic decoupling
    double scale_vel = 0.0;
    int addc_index = 0;
    for(unsigned int i=0; i < topo.num_atoms(); ++i){
      addc_index = sim.param().addecouple.check_index_adc(i);
      scale_vel=1;
      if(addc_index!=-1)
       scale_vel=1/sim.param().addecouple.adc_index()[addc_index].st;
      
      for(int a=0; a<3; ++a){
	for(int bb=0; bb<3; ++bb){
	  conf.current().kinetic_energy_tensor(a, bb) +=
	    0.5 * topo.mass()(i) *
	    conf.current().vel(i)(a) * 
	    conf.current().vel(i)(bb) * scale_vel;
	}
      }
    }

    // system().molecular_kinetic_energy() *= 0.5;
  }

}


void util::prepare_virial(topology::Topology const & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation const & sim)
{

  if (conf.boundary_type == math::vacuum) return;

  SPLIT_BOUNDARY(_prepare_virial, topo, conf, sim);
    
}

template<math::boundary_enum boundary>
static void _atomic_to_molecular_virial(topology::Topology const & topo,
					configuration::Configuration & conf,
					simulation::Simulation const & sim)
{
  // this should be done after bonded and nonbonded forces have been calculated
  // but before any special forces are added (because they don't contribute to
  // the virial)

  if (sim.param().pcouple.virial == math::molecular_virial){

    DEBUG(7, "recovering molecular virial from atomic virial");
    DEBUG(10, "lambda = " << topo.lambda());

    
    
    math::Periodicity<boundary> periodicity(conf.current().box);
    math::VArray const &pos = conf.current().pos;
    math::Vec r;

    math::Matrix corrP(0.0);
    // EDS
    const int numstates = conf.special().eds.virial_tensor_endstates.size();
    std::vector<math::Matrix> corrPendstates(numstates,corrP);
    
    topology::Pressuregroup_Iterator
      pg_it = topo.pressure_group_begin(),
      pg_to = topo.pressure_group_end();

    math::Vec com_pos;
    math::Matrix com_ekin;

    for( ; pg_it != pg_to; ++pg_it){
      _centre_of_mass(pg_it.begin(),
		      pg_it.end(),
		      topo, conf,
		      com_pos, com_ekin,
		      periodicity, sim);

      topology::Atom_Iterator a_it = pg_it.begin(),
	a_to = pg_it.end();

      for( ; a_it != a_to; ++a_it){

	// this should be changed to
	// the vector between the free floating atom position
	// and the free floating virial group centre of mass position
	periodicity.nearest_image(pos(*a_it), com_pos, r);

	for(int a=0; a<3; ++a){
	  for(int b=0; b<3; ++b){

	    // conf.current().virial_tensor(b, a) +=
	    corrP(b, a)  +=
	      conf.current().force(*a_it)(a) * r(b);
            
            // EDS
            for(int state = 0; state < numstates; state++){
              corrPendstates[state](b, a) +=
                      conf.special().eds.force_endstates[state](*a_it)(a) * r(b);
            }
             
            
	  }
	}

      } // loop over virial group
    }

    conf.current().virial_tensor -= corrP;
    
    // EDS
    for(int state = 0; state < numstates; state++){
      conf.special().eds.virial_tensor_endstates[state] -= corrPendstates[state];
    }

  } // molecular virial

}


void util::atomic_to_molecular_virial(topology::Topology const & topo,
				      configuration::Configuration & conf,
				      simulation::Simulation const & sim)
{
  
  if (conf.boundary_type == math::vacuum) return;

  SPLIT_BOUNDARY(_atomic_to_molecular_virial, topo, conf, sim);
  
}

template<math::boundary_enum b>
void _centre_of_mass_loop(topology::Topology const & topo,
			  configuration::Configuration & conf,
			  std::vector<math::Vec> & com_pos,
			  std::vector<math::Matrix> & com_ekin, 
                             simulation::Simulation const & sim)
{
  math::Periodicity<b> periodicity(conf.current().box);
  
  const size_t num = topo.pressure_groups().size();
  com_pos.assign(num, math::Vec(0.0));
  com_ekin.assign(num, math::Matrix(0.0));
  
  topology::Pressuregroup_Iterator
    pg_it = topo.pressure_group_begin(),
    pg_to = topo.pressure_group_end();
  
  for(int i=0; pg_it != pg_to; ++pg_it, ++i){
    _centre_of_mass(pg_it.begin(), pg_it.end(),
		    topo, conf,
		    com_pos[i], com_ekin[i],
		    periodicity, sim);
  }
}

void util::centre_of_mass(topology::Topology const & topo,
			  configuration::Configuration & conf,
			  std::vector<math::Vec> & com_pos,
			  std::vector<math::Matrix> & com_ekin, 
                             simulation::Simulation const & sim)
{
  SPLIT_BOUNDARY(_centre_of_mass_loop, topo, conf, com_pos, com_ekin, sim);
}
