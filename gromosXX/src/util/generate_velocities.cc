/**
 * @file generate_velocities.cc
 * generate velocities from Maxwell distribution.
 */

#include <stdheader.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "generate_velocities.h"

void util::generate_velocities(double const temp, math::SArray const &mass,
			       math::VArray & vel, math::VArray & old_vel,
			       unsigned int const seed)
{

  std::cout << "\n\tgenerating initial velocities\n"
	    << "\t\ttemperature         = " << temp << "\n"
	    << "\t\trandom number seed  = " << seed
	    << "\n";
  
  // get a random number generator type
  const gsl_rng_type *rng_type;
	
  // enable control via environment variables
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;

  // get the rundom number generator
  gsl_rng *rng = gsl_rng_alloc(rng_type);
	
  // print a comment
  std::cout << "\t\trandom number generator: " << gsl_rng_name(rng)
	    << "\n\t\tdefault seed: " << gsl_rng_default_seed << "\n\n";

  for(unsigned int i=0; i<vel.size(); ++i){
    const double sd = sqrt(math::k_Boltzmann * temp / mass(i));
    for(int d=0; d<3; ++d){
      // old_vel(i)(d) = sd * n.random();
      old_vel(i)(d) = gsl_ran_gaussian(rng, sd);
      vel(i)(d) = old_vel(i)(d);
    }
  }
}
