/**
 * @file timing.cc
 */

#include <time.h>

double now()
{
  return double(clock()) / CLOCKS_PER_SEC;
}

void print_timing(std::ostream &os)
{
  os << "TIMING\n"
     << "simulations lasting longer than 30' are not timed correctly!\n\n"
     << std::setw(20) << "PROCEDURE" << std::setw(15) 
     << "TIME" << std::setw(10) << "CALLS\n"
     << std::setw(20) << "total" << std::setw(15) << timing.total << "\n"
    //     << std::setw(20) << "startup" << std::setw(15) << timing.startup << "\n"
     << std::setw(20) << "integration" << std::setw(15) << timing.integration
     << std::setw(10) << timing.count_integration << "\n"
     << std::setw(20) << "constraints" << std::setw(15) << timing.constraints
     << std::setw(10) << timing.count_constraints << "\n"
     << std::setw(20) << "force" << std::setw(15) << timing.force
     << std::setw(10) << timing.count_force << "\n"
     << std::setw(20) << "pairlist" << std::setw(15) << timing.pairlist
     << std::setw(10) << timing.count_pairlist << "\n"
     << std::setw(20) << "shortrange" << std::setw(15) << timing.shortrange
     << std::setw(10) << timing.count_shortrange << "\n"
    // << std::setw(20) << "output" << std::setw(15) << timing.output << "\n"
     << "END\n";

}

timing_struct timing;
