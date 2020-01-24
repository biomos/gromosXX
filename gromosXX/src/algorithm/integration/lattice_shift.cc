/**
 * @file lattice_shift.cc
 * implementation of the lattice shift tracking
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../math/periodicity.h"
#include "../../util/template_split.h"

#include "lattice_shift.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

int algorithm::Lattice_Shift_Tracker::
init(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim,
        std::ostream &os,
        bool quiet) {
  
  // initialize the lattice shifts if not read from configuration
  if (!sim.param().start.read_lattice_shifts) {
    conf.special().lattice_shifts = 0.0;
  }
  
  if (!quiet) { // write some stuff to the output file
    os << "LATTICESHIFTS" << std::endl
       << "    keeping track of lattice shifts." << std::endl;
    
    if (sim.param().start.read_lattice_shifts)
      os << "    reading initial shifts from configuration.";
    else
      os << "    setting initial shifts to zero.";
    
    os << std::endl
       << "END" << std::endl;
  }
  
  return 0;
}

int algorithm::Lattice_Shift_Tracker::
apply(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  DEBUG(6, "keeping track of lattice shifts");
  m_timer.start();
  SPLIT_BOUNDARY(_apply, topo, conf, sim);
  m_timer.stop();
  return 0;
}

template<math::boundary_enum b>
void algorithm::Lattice_Shift_Tracker::
_apply(topology::Topology &topo,
        configuration::Configuration &conf,
        simulation::Simulation &sim) {
  math::Periodicity<b> p(conf.current().box);
  // just call the function from the periodicity
  p.put_chargegroups_into_box_saving_shifts(conf, topo);
}



