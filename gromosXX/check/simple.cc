/**
 * @file simple.cc
 * simple checks.
 */

#include <config.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

// check macros and functions
#include "check.h"

// include all headers here (once and for all...)
#include "../src/debug.h"

#include <math/gmath.h>
#include <io/message.h>
#include <simulation/core.h>
#include <math/periodicity.h>
#include <simulation/simulation.h>
#include <interaction/interaction.h>
#include <io/io.h>
#include <algorithm/algorithm.h>

// sppecial includes
#include <algorithm/integration/runge_kutta.h>

// gloabl variables for debug
#include "../src/debug.cc"

using namespace math;

#include "math.t.cc"
#include "simulation.t.cc"

int main(int argc, char *argv[])
{
  int result = 0;
  
  { // arguments
    
    char *knowns[] = 
      {
	"verb"
      };
    
    int nknowns = 1;
    
    string usage = argv[0];
    usage += "\t@verb    <[module:][submodule:]level>\n";
    
    io::Argument args(argc, argv, nknowns, knowns, usage);

    // parse the verbosity flag and set debug levels
    parse_verbosity(args);
  }
  
  CHECK_MODULE(math_periodicity, result);
  CHECK_MODULE(simulation_check, result);
  
  std::cout << "\n";
  
  return result;
}
