/**
 * @file generate_velocities.cc
 * generate velocities from Maxwell distribution.
 */

#include <util/stdheader.h>
#include <random/normal.h>

#include "generate_velocities.h"

using namespace ranlib;
BZ_USING_NAMESPACE(ranlib)

void util::generate_velocities(double const temp, math::SArray const &mass,
			       math::VArray & vel, math::VArray & old_vel,
			       unsigned int const seed)
{

  NormalUnit<double> n;
  n.seed(seed);

  for(int i=0; i<vel.size(); ++i){
    const double sd = sqrt(math::k_Boltzmann * temp / mass(i));
    for(int d=0; d<3; ++d){
      old_vel(i)(d) = sd * n.random();
      vel(i)(d) = old_vel(i)(d);
    }
  }
}
