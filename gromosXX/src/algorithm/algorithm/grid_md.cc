#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include <debug.h>
#include <timing.h>

#include <math/gmath.h>
#include <io/message.h>
#include <simulation/core.h>
#include <math/periodicity.h>
#include <simulation/simulation.h>
#include <simulation/perturbation.h>
#include <interaction/interaction.h>
#include <io/io.h>
#include <algorithm/algorithm.h>

// special includes
// #include <algorithm/integration/runge_kutta.h>

using namespace math;

namespace algorithm
{
  int grid_scaled_md(io::Argument &args, io::InInput &input)
  {
    algorithm::Perturbation_MD<
      algorithm::perturbed_MD_spec<interaction::no_virial>,
      algorithm::Interaction_spec<
      algorithm::perturbed_MD_spec<interaction::no_virial>::simulation_type,
      // perturbation
      false,
      // virial
      interaction::no_virial,
      // atomic cutoff
      false,
      // scaling
      false,
      // bekker
      true
      >
      > 
      the_MD;
    
    return the_MD.do_md(args, input);
  }
}

