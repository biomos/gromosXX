/**
 * @file umbrella_weight.cc
 * implementation of umbrella weighting
 */
#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../util/debug.h"

#include "umbrella_weight.h"

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE leus

std::ostream & operator<<(std::ostream & os, const util::Umbrella_Weight & w) {
  w.write(os);
  return os;
}
std::istream & operator>>(std::istream & is, util::Umbrella_Weight & w) {
  w.read(is);
  return is;
}

void util::Number_Of_Visits_Umbrella_Weight::write(std::ostream & os) const {
  os << std::setw(10) << weight;
}

util::Umbrella_Weight_Factory::~Umbrella_Weight_Factory() {
  // delete all the created instaces
  for(std::vector<util::Umbrella_Weight*>::iterator it = instances.begin(),
          to = instances.end(); it != to; ++it)
    delete *it;
}

