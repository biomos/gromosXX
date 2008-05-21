/**
 * prepare for virial calculation.
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include <math/periodicity.h>

#include "prepare_virial.h"

#include <util/template_split.h>
#include <util/debug.h>

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
			    math::Periodicity<b> const & periodicity)
{

  com_pos = 0.0;
  double m;
  double tot_mass = 0.0;

  math::Vec p;
  math::Vec prev;
  math::Vec v(0.0);

  prev = conf.current().pos(*start);

  for( ; start != end; ++start){

    assert(unsigned(topo.mass().size()) > *start &&
           unsigned(conf.current().pos.size()) > *start);

    m = topo.mass()(*start);
    tot_mass += m;
    periodicity.nearest_image(conf.current().pos(*start), prev, p);
    com_pos += m * (p + prev);
    v += m * conf.current().vel(*start);
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

    topology::Molecule_Iterator
      m_it = topo.molecule_begin(),
      m_to = topo.molecule_end();

    math::Vec com_pos;
    math::Matrix com_ekin;

    conf.current().kinetic_energy_tensor = 0.0;

    for( ; m_it != m_to; ++m_it){
      _centre_of_mass(m_it.begin(),
		      m_it.end(),
		      topo, conf,
		      com_pos, com_ekin,
		      periodicity);

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

    for(unsigned int i=0; i < topo.num_atoms(); ++i){
      for(int a=0; a<3; ++a){
	for(int bb=0; bb<3; ++bb){
	  conf.current().kinetic_energy_tensor(a, bb) +=
	    0.5 * topo.mass()(i) *
	    conf.current().vel(i)(a) * 
	    conf.current().vel(i)(bb);
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
    
    topology::Molecule_Iterator
      m_it = topo.molecule_begin(),
      m_to = topo.molecule_end();

    math::Vec com_pos;
    math::Matrix com_ekin;

    for( ; m_it != m_to; ++m_it){
      _centre_of_mass(m_it.begin(),
		      m_it.end(),
		      topo, conf,
		      com_pos, com_ekin,
		      periodicity);

      topology::Atom_Iterator a_it = m_it.begin(),
	a_to = m_it.end();

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
			  std::vector<math::Matrix> & com_ekin)
{
  math::Periodicity<b> periodicity(conf.current().box);
  
  const size_t num = topo.molecules().size();
  com_pos.assign(num, math::Vec(0.0));
  com_ekin.assign(num, math::Matrix(0.0));
  
  topology::Molecule_Iterator
    m_it = topo.molecule_begin(),
    m_to = topo.molecule_end();
  
  for(int i=0; m_it != m_to; ++m_it, ++i){
    _centre_of_mass(m_it.begin(), m_it.end(),
		    topo, conf,
		    com_pos[i], com_ekin[i],
		    periodicity);
  }
}

void util::centre_of_mass(topology::Topology const & topo,
			  configuration::Configuration & conf,
			  std::vector<math::Vec> & com_pos,
			  std::vector<math::Matrix> & com_ekin)
{
  SPLIT_BOUNDARY(_centre_of_mass_loop, topo, conf, com_pos, com_ekin);
}
