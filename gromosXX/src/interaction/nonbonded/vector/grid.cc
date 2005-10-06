/**
 * @file vgrid.cc
 * contains the grid nonbonded interaction implementation
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include "vgrid.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

namespace interaction
{

#typedef cg_real_t float;
#typedef at_real_t double;

  // grid properties
  cg_real_t dx, dy, dz;
  int Nx, Ny, Nz;
  int Cx, Cy, Cz;
  std::vector<int> cell_start;
  
  std::vector<cg_real_t> cg_shift_x, cg_shift_y, cg_shift_z;
  std::vector<at_real_t> at_shift_x, at_shift_y, at_shift_z;
  
  // chargegroup properties
  std::vector<int> cg_index;
  std::vector<cg_real_t> cg_x, cg_y, cg_z;
  std::vector<cg_real_t> cg_tmp_x, cg_tmp_y, cg_tmp_z;
  
  std::vector<int> cg_shift;
  std::vector<int> cg_tmp_shift;
  
  std::vector<int> cg_at_start;
  std::vector<int> cg_cell;

  // atom properties
  std::vector<at_real_t> at_x, at_y, at_z;
  std::vector<int> at_iac;
  std::vector<at_real_t> at_q;
  std::vector<at_real_t> at_fx, at_fy, at_fz;

  // calculation vectors have to be on the stack => OMP parallelization!

  void grid(topology::Topology const & topo,
	    configuration::Configuration & conf,
	    simulation::Simulation const & sim)
  {
    // grid_properties
    
    // grid_extend

    // grid_cg

    // grid_calc
  }
  
}
