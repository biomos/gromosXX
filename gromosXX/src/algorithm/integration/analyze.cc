/**
 * @file analyze.cc
 * contains the implementation
 * for trajectory analysis
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../io/blockinput.h"
#include "../../io/instream.h"
#include "../../io/configuration/inframe.h"
#include "../../io/configuration/in_configuration.h"

#include "analyze.h"

#include "../../util/error.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

algorithm::Analyze_Step::Analyze_Step(std::string trajectory)
  : Algorithm("AnalyzeStep")
{
  DEBUG(7, "reading configuration");
  m_trajectory_file.open(trajectory.c_str());

  if (!m_trajectory_file.is_open()){
    std::cout << "\n\ncould not open " << trajectory << "!\n" << std::endl;
    io::messages.add("opening configuration failed", "read_input",
		     io::message::critical);
  }

  io::messages.add("analyzation trajectory read from " + trajectory,
		   "read input",
		   io::message::notice);
  
  // and associate
  m_trajectory.stream(m_trajectory_file);

}


/**
 * analyzation step init
 */
int algorithm::Analyze_Step
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream & os,
       bool quiet)
{
  if (!quiet)
    os << "TRAJECTORY ANALYZATION\n";
  

  if (sim.param().analyze.copy_pos){
    os << "\tcopying current position to old\n"
       << "\tvelocities are set to 0.\n";
  }
  
  if (!quiet)
    os << "END\n";
  
  // save initial configuration as old() state
  conf.exchange_state();

  return 0;
}

/**
 * trajectory analyzation step.
 */
int algorithm::Analyze_Step
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  // just overwrite current conf
  DEBUG(8, "analyze: reading next frame!");
  if (m_trajectory.read_next(topo, conf, sim)){

    if (sim.param().analyze.copy_pos){
      DEBUG(9, "analyze: copy current() pos to old()");
      conf.old().pos = conf.current().pos;
    }

    conf.old().vel = 0;
    conf.current().vel = 0;
    
    return 0;
  }
  
  return E_UNSPECIFIED;
}

algorithm::Analyze_Exchange::Analyze_Exchange()
  : Algorithm("AnalyzeExchange")
{
}

int algorithm::Analyze_Exchange
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream & os,
       bool quiet)
{
  return 0;
}

/**
 * Here the states are exchanged.
 */
int algorithm::Analyze_Exchange
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
        simulation::Simulation &sim) {
  conf.exchange_state();
  return 0;
}

