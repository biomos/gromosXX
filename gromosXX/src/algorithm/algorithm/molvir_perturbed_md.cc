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
  int molvir_perturbed_md(io::Argument &args, io::InInput &input)
  {
    algorithm::Perturbation_MD<
      algorithm::perturbed_MD_spec<interaction::molecular_virial>,
      algorithm::Interaction_spec<
      algorithm::perturbed_MD_spec<interaction::molecular_virial>::simulation_type,
      // perturbation
      true,
      // virial
      interaction::molecular_virial,
      // atomic cutoff
      false,
      // scaling
      false,
      // bekker
      false
      >
      > 
      the_MD;
    
    return the_MD.do_md(args, input);
  }
}

