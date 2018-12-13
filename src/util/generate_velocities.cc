/**
 * @file generate_velocities.cc
 * generate velocities from Maxwell distribution.
 */

#include <stdheader.h>
#include <math/random.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>

#include "generate_velocities.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

void util::generate_velocities(const simulation::Parameter &param,
                               double const temp, math::SArray const &mass,
			       math::VArray & vel, math::VArray & old_vel,
			       unsigned int const seed,
			       std::ostream & os, bool quiet)
{

  if (!quiet) {
    os << "\n\tgenerating initial velocities\n"
       << "\t\ttemperature         = " << temp << "\n"
       << "\t\trandom number seed  = " << seed << "\n";
  }
  
  std::ostringstream stringseed; stringseed << seed;
  math::RandomGenerator* rng = math::RandomGenerator::create(param, stringseed.str());
  if (!quiet)
    os << "\t\trandom number generator: " << rng->description() << "\n";
  
  for(unsigned int i=0; i<vel.size(); ++i){
    const double sd = sqrt(math::k_Boltzmann * temp / mass(i));
    rng->stddev(sd);
    for(int d=0; d<3; ++d){
      // old_vel(i)(d) = sd * n.random();
      old_vel(i)(d) = rng->get_gauss();
      vel(i)(d) = old_vel(i)(d);
    }
  }

  delete rng;
}
