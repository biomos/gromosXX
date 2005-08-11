/**
 * @file virtual_grain.cc
 * update positions / forces of virtual atoms (grains)
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <configuration/configuration.h>

#include "virtual_grain.h"

#include "debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util


/**
 * update positions of virtual atoms (grains) from
 * real atom positions
 */
void util::update_virtual_pos(topology::Topology & cg_topo,
			      configuration::Configuration & cg_conf,
			      topology::Topology & topo,
			      configuration::Configuration & conf)
{
  std::cerr << "update virtual pos" << std::endl;

  cg_conf.current().box = conf.current().box;
  
  for(unsigned int i=0; i<cg_topo.virtual_grains().size(); ++i){
    
    assert(i < cg_topo.virtual_grains().size());
    assert(cg_conf.current().pos.size() > unsigned(cg_topo.virtual_grains()[i].i));

    DEBUG(10, "virtual pos " << cg_topo.virtual_grains()[i].i
	  << " = " << math::v2s(cg_topo.virtual_grains()[i].atom.pos(conf)));

    // std::cerr << "virtual pos [" << i << "] : " << cg_topo.virtual_grains()[i].i << std::endl;
    // std::cerr << "\t(" << cg_topo.virtual_grains()[i].atom.size() << ")";
    // std::cerr << std::endl;
    
    cg_conf.current().pos(cg_topo.virtual_grains()[i].i)
      = cg_topo.virtual_grains()[i].atom.pos(conf);
  }
}

/**
 * distribute forces of virtual atoms (grains) on
 * real atoms
 */
void util::update_virtual_force(topology::Topology & cg_topo,
				configuration::Configuration & cg_conf,
				topology::Topology & topo,
				configuration::Configuration & conf)
{
  cg_conf.current().energies.calculate_totals();
  conf.current().energies.external_total += cg_conf.current().energies.potential_total;
  
  conf.current().virial_tensor += cg_conf.current().virial_tensor;
  
  for(unsigned int i=0; i<cg_topo.virtual_grains().size(); ++i){

    DEBUG(10, "virtual force " << cg_topo.virtual_grains()[i].i
	  << " = " << math::v2s(cg_conf.current().force(cg_topo.virtual_grains()[i].i)));
    
    cg_topo.virtual_grains()[i].atom.force(conf, cg_conf.current().force(cg_topo.virtual_grains()[i].i));

  }
}

