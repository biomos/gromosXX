/**
 * @file check_forcefield.cc
 */


#include <util/stdheader.h>

#include <topology/core/core.h>

#include <topology/solute.h>
#include <topology/solvent.h>
#include <topology/perturbed_atom.h>
#include <topology/perturbed_solute.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <simulation/simulation.h>
#include <configuration/energy.h>
#include <configuration/energy_average.h>
#include <configuration/configuration.h>
#include <algorithm/algorithm.h>
#include <algorithm/algorithm/algorithm_sequence.h>
#include <interaction/interaction.h>
#include <interaction/forcefield/forcefield.h>

#include <io/argument.h>
#include <util/parse_verbosity.h>
#include <util/error.h>

#include <simulation/parameter.h>
#include <interaction/interaction_types.h>
#include <io/instream.h>
#include <util/parse_tcouple.h>
#include <io/blockinput.h>
#include <io/topology/in_topology.h>

#include <algorithm/integration/leap_frog.h>
#include <algorithm/temperature/temperature_calculation.h>
#include <algorithm/temperature/berendsen_thermostat.h>
#include <algorithm/pressure/pressure_calculation.h>
#include <algorithm/pressure/berendsen_barostat.h>

#include <interaction/forcefield/forcefield.h>
#include <interaction/forcefield/create_forcefield.h>

#include <util/create_simulation.h>
#include <algorithm/create_md_sequence.h>

#include <math/volume.h>

#include <time.h>

#include <config.h>

#include "check.h"
#include "check_state.h"

using namespace std;


void scale_positions(topology::Topology & topo,
		     configuration::Configuration & conf,
		     math::Vec const scale)
{
  for(size_t i = 0; i< topo.num_atoms(); ++i){
    conf.current().pos(i) *= scale;
  }
  conf.current().box(0) *=scale;
  conf.current().box(1) *=scale;
  conf.current().box(2) *=scale;
  
}

int check::check_state(topology::Topology & topo,
		       configuration::Configuration & conf,
		       simulation::Simulation & sim,
		       interaction::Forcefield & ff)
{
  const double epsilon = 0.000001;
  int res=0, total=0;

  CHECKING("Virial (finite diff)",res);

  conf.exchange_state();
  conf.current().pos = conf.old().pos;
  conf.current().box = conf.old().box;
  
  
  math::Matrix finP;
  finP=0;
  
  for(int i=0; i < 3; i++){
    math::Vec s1=1;
    s1(i) += epsilon;
      
    scale_positions(topo, conf, s1);
    ff.apply(topo, conf, sim);
    conf.current().energies.calculate_totals();
    double e1 = conf.current().energies.nonbonded_total;
      
    conf.current().pos = conf.old().pos;
    conf.current().box = conf.old().box;

    math::Vec s2=1;
    s2(i) -= epsilon;

    scale_positions(topo, conf, s2);
    ff.apply(topo, conf, sim);
    conf.current().energies.calculate_totals();
    double e2 = conf.current().energies.nonbonded_total;

    conf.current().pos = conf.old().pos;
    conf.current().box = conf.old().box;

    for(int j=0; j < 3; j++){
      finP(i,j)=0.0;
      for(int k=0; k< 3; k++){
	finP(i,j) += conf.current().box(j)(k) * (e2-e1)/(2*epsilon*conf.current().box(i)(k));
      }
      
    }
  }
    
  // calculate the energy
  ff.apply(topo, conf, sim);

  conf.exchange_state();
  
  algorithm::Pressure_Calculation * pcalc =
    new algorithm::Pressure_Calculation;

  pcalc->apply(topo, conf, sim);
  
  for(int i=0; i < 3; i++){
    for(int j=0; j < 3; j++){
      CHECK_APPROX_EQUAL(conf.old().virial_tensor(i,j), finP(i,j), 0.0000000001, res);
    }
  }

  RESULT(res,total);

  return total;
}
