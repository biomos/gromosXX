/**
 * @file generate_velocities.h
 * generate initial velocities from Maxwell distribution.
 */

#ifndef INCLUDED_GENERATE_VELOCITIES_H
#define INCLUDED_GENERATE_VELOCITIES_H

namespace simulation {
  class Parameter;
}

namespace util
{
  /**
   * generate initial velocities.
   */
  void generate_velocities(const simulation::Parameter &param,
                           double const temp, math::SArray const &mass,
			   math::VArray &vel, math::VArray &old_vel,
			   unsigned int const seed,
			   std::ostream & os = std::cout, 
                           bool quiet = false);
  
}

#endif
