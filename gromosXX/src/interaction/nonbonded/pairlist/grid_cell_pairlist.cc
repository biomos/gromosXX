/**
 * @file grid_cell_pairlist.cc
 * Implementation of a fast grid-cell algorithm
 */

#include "../../../stdheader.h"
#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../math/gmath.h"
#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"
#include "../../../interaction/nonbonded/pairlist/grid_cell_pairlist.h"
#include "../../../interaction/nonbonded/pairlist/standard_pairlist_algorithm.h"


#include "../../../util/debug.h"
#include "../../../util/template_split.h"

#include <algorithm>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

/**
 * Constructor
 */
interaction::Grid_Cell_Pairlist::Grid_Cell_Pairlist(const topology::Topology & topo,
            const simulation::Simulation &sim) :
        Failing_Pairlist_Algorithm(), is_vacuum(false) {
  DEBUG(10, "Grid_Cell : Constructor");
}

/**
 * Destructor
 */
interaction::Grid_Cell_Pairlist::~Grid_Cell_Pairlist() {
  DEBUG(10, "Grid_Cell : Destructor");
  if (fallback_algorithm != NULL)
    delete fallback_algorithm;
}

/**
 * Initialize the pairlist
 */
int interaction::Grid_Cell_Pairlist::init(topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim,
        std::ostream & os,
        bool quiet) {
  DEBUG(5, "Grid_Cell : init");
  mytopo = &topo;
  myconf = &conf;
  mysim = &sim;

  math::boundary_enum b = conf.boundary_type;
  switch (b) {
    case math::vacuum: {
      is_vacuum = true;
      irregular_shape = false;
      create_vacuum_box();
      break;
    }
    case math::rectangular: {
      irregular_shape = false;
      break;
    }
    case math::triclinic:
    case math::truncoct: {
      irregular_shape = true;
      break;
    }
    default: {
      io::messages.add("No boundary condition given!",
              "Grid Cell", io::message::error);
      return 1;
    }
  } // switch

  // Values for the irregular shape
  if (irregular_shape) {
    DEBUG(10, "Grid cell algorithm with irregular box shapes");
    
    const math::Vec & a = conf.current().box(0);
    const math::Vec & b = conf.current().box(1);
    const math::Vec & c = conf.current().box(2);
    math::Vec d = math::cross(a, b);
    d = math::dot(d.norm(), c) * d.norm();
    math::Vec e = math::cross(d, a);
    e = math::dot(e.norm(), b) * e.norm();
    length_x = math::abs(a);
    length_y = math::abs(e);
    length_z = math::abs(d);

    cases[0] = 0; // place holder for fortran indexing
    cases[1] = 2;
    cases[2] = 1;
    cases[3] = 3;
    cases[4] = 4;
    cases[5] = 2;
    cases[6] = 1;
    cases[7] = 6;
    cases[8] = 5;
    cases[9] = 1;
  }
  else {
    length_x = math::abs(conf.current().box(0));
    length_y = math::abs(conf.current().box(1));
    length_z = math::abs(conf.current().box(2));
  }

  const double grid_size = sim.param().pairlist.grid_size;
  // For having a reasonable amount of boxes
  const unsigned int factor = 1;
  num_x = int (length_x / grid_size) * factor;
  num_y = int (length_y / grid_size) * factor;
  num_z = int (length_z / grid_size) * factor;
  DEBUG(5, "num_x = " << num_x << " num_y = " << num_y << " num_z = " << num_z);
  num_y_minus_1 = num_y - 1;
  num_z_minus_1 = num_z - 1;
  num_x_half = num_x / 2.0;
  num_y_half = num_y / 2.0;
  num_z_half = num_z / 2.0;
  
  cutoff_lr = sim.param().pairlist.cutoff_long;
  cutoff_sr = sim.param().pairlist.cutoff_short;
  cutoff_lr2 = cutoff_lr * cutoff_lr;
  cutoff_sr2 = cutoff_sr * cutoff_sr;
  DEBUG(15, "cutoff_lr = " << cutoff_lr);
  DEBUG(15, "cutoff_lr2 = " << cutoff_lr2);

  // Guess Number of stripes;  eq.(4)
  int num_s =  int(cutoff_lr * cutoff_lr * num_y * num_z / (length_y * length_z));

  // Reserve some space for the mask pointer
  DEBUG(15, "Grid Cell : Reserve Space for the mask pointer")
  mask_pointer.reserve(4 * num_s);

  // Create a mask and the mask pointer, if there is no pressure scaling
  if (sim.param().pcouple.scale == math::pcouple_off){
    unsigned int err = calc_par();
    if(irregular_shape)
      err += make_mask_pointer<Irr_Mask > ();
    else
      err += make_mask_pointer<Reg_Mask > ();
    
    if (err) {
      io::messages.add("Could not create mask!",
              "Grid Cell", io::message::error);
      return 1;
    }
  }

  // Reserve some space for the cell array
  DEBUG(15, "Grid Cell : Reserve Space for the cell")
  cell.reserve(topo.num_chargegroups() + 1);


  // Is cuda on?
  cuda = (sim.param().innerloop.method == simulation::sla_cuda);

  if (is_vacuum)
    restore_vacuum_box();

  if (!quiet) {
    os << "\tgrid cell pairlist algorithm\n"
            << "\t\tcells           : " << num_x << "x"  << num_y << "x"  << num_z << "\n";
    os << "\t\tusing mask routines for " << (irregular_shape ? "triclinic" : "rectangular") << " box shapes.\n";
  }

  fallback_algorithm = new Standard_Pairlist_Algorithm();
  fallback_algorithm->init(topo, conf, sim, os, true);

  return 0;
}

/**
 * calculate center of geometries
 */
int interaction::Grid_Cell_Pairlist::
prepare(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{

  timer().start("pairlist prepare");

  is_vacuum = false;
  if (is_vacuum) {
    create_vacuum_box();
  }

  DEBUG(6, "grid cell pairlist algorithm : prepare");
  // first put the chargegroups into the box. We also do this for atomic cutoffs
  SPLIT_BOUNDARY(_prepare_cog, conf, topo);
  math::VArray const &pos = conf.current().pos;
  if (!sim.param().pairlist.atomic_cutoff){
    // calculate cg cog's
    DEBUG(10, "calculating cg cog (" << topo.num_solute_chargegroups() << ")");
    m_cg_cog.resize(topo.num_chargegroups());
    DEBUG(10, "pos.size() = " << pos.size());

    // calculate solute center of geometries
    topology::Chargegroup_Iterator
      cg1 =   topo.chargegroup_begin();

    // Add solute charge groups
    unsigned int i = 0, num_cg = topo.num_solute_chargegroups();
    first_solvent = num_cg;
    for(i=0; i < num_cg; ++cg1, ++i){
      cg1.cog(pos, m_cg_cog(i));
    }
    // Add solvent charge groups
    num_cg = topo.num_chargegroups();
    for(; i < num_cg; ++cg1, ++i){
      m_cg_cog(i) = pos(*(cg1.begin()));
    }
  } else { // atomic cutoff
    first_solvent = topo.num_solute_atoms();
    num_atoms_per_solvent = topo.solvent(0).num_atoms();
    m_cg_cog = pos;
  }

  failed = false;
  // If pressure is constant, make mask again
  timer().start("pairlist mask");
  unsigned int err = 0;
  if (sim.param().pcouple.scale != math::pcouple_off) {
    DEBUG(7, "Grid Cell : Pressure is constant");
    err += calc_par();

    if (irregular_shape) {
      err += make_mask_pointer<Irr_Mask > ();
    }
    else
      err += make_mask_pointer<Reg_Mask > ();

    if (err) {
      io::messages.add("Could not create mask!",
              "Grid Cell", io::message::error);
      return 1;
    }
  }
  timer().stop("pairlist mask");
  if (err) {
    std::ostringstream msg;
    msg << "At step " << sim.steps() << ": Could not prepare grid. "
            "Falling back to standard algoritm for this step.";
    io::messages.add(msg.str(), "Grid_Cell_Pairlist", io::message::notice);
    failed = true;
    return fallback_algorithm->prepare(topo, conf, sim);
  }

  timer().start("pairlist cell");
  SPLIT_BOUNDARY(make_cell, topo, conf, sim);
  timer().stop("pairlist cell");

  timer().stop("pairlist prepare");
  return 0;
}

/**
 * put the chargegroups into the box
 */
template<math::boundary_enum b>
void interaction::Grid_Cell_Pairlist::_prepare_cog(configuration::Configuration & conf,
        topology::Topology & topo) {
  DEBUG(8, "putting chargegroups into box");
  math::Periodicity<b> periodicity(conf.current().box);
  periodicity.put_chargegroups_into_box(conf, topo);
}

/**
 * update the pairlist
 */
void interaction::Grid_Cell_Pairlist::update(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim,
        interaction::PairlistContainer &pairlist,
        unsigned int begin, unsigned int end,
        unsigned int stride) {

  DEBUG(5, "Grid Cell : Update pairlist");
  if (failed) {
    fallback_algorithm->update(topo, conf, sim, pairlist, begin, end, stride);
    return;
  }

  if (begin == 0) // master
    timer().start("pairlist");

  {
    pairlist.clear();
#ifdef OMP
    #pragma omp barrier
#endif
  }
  // _pairlist(pairlist)
  if (!sim.param().pairlist.atomic_cutoff) {
    SPLIT_BOUNDARY(_pairlist, pairlist, pairlist, begin, stride, cg_cutoff(), no_perturbation());
  } else {
    SPLIT_BOUNDARY(_pairlist, pairlist, pairlist, begin, stride, atomic_cutoff(), no_perturbation());
  }

  if (begin == 0) // master
    timer().stop("pairlist");
#ifdef OMP
  #pragma omp barrier
#endif

  if (begin == 0 && is_vacuum) {
    restore_vacuum_box();
  }
}

/**
 * update the (perturbed) pairlist
 */
void interaction::Grid_Cell_Pairlist::update_perturbed(
        topology::Topology & topo,
        configuration::Configuration & conf,
        simulation::Simulation &sim,
        interaction::PairlistContainer & pairlist,
        interaction::PairlistContainer & perturbed_pairlist,
        unsigned int begin, unsigned int end,
        unsigned int stride) {
  DEBUG(5, "Grid Cell : Update pairlist & pertubred pairlist");
  if (failed) {
    fallback_algorithm->update_perturbed(topo, conf, sim, pairlist, perturbed_pairlist, begin, end, stride);
    return;
  }

  if (begin == 0) // master
    timer().start("perturbed pairlist");

  {
    pairlist.clear();
    perturbed_pairlist.clear();
    #ifdef OMP
    #pragma omp barrier
    #endif
  }
  // _pairlist(pairlist)
  if (!sim.param().pairlist.atomic_cutoff) {
    SPLIT_BOUNDARY(_pairlist, pairlist, perturbed_pairlist, begin, stride, cg_cutoff(), do_perturbation());
  } else {
    SPLIT_BOUNDARY(_pairlist, pairlist, perturbed_pairlist, begin, stride, atomic_cutoff(), do_perturbation());
  }

  if (begin == 0) // master
    timer().stop("perturbed pairlist");

  #ifdef OMP
  #pragma omp barrier
  #endif

  if (begin == 0 && is_vacuum) {
    restore_vacuum_box();
  }
}

/**
 * Calculate parameters like dim_m, length_i etc.
 */
int interaction::Grid_Cell_Pairlist::calc_par(){

  DEBUG(10, "Grid cell : Calculate parameters");
  if (irregular_shape) {
    math::Vec a = myconf->current().box(0);
    const math::Vec & b = myconf->current().box(1);
    const math::Vec & c = myconf->current().box(2);

    double BNX = math::abs(a);
    double BOXX = BNX;

    double FAC1 = a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
    double UX = b(0) - FAC1 / (BOXX * BOXX) * a(0);
    double UY = b(1) - FAC1 / (BOXX * BOXX) * a(1);
    double UZ = b(2) - FAC1 / (BOXX * BOXX) * a(2);
    double BNY = sqrt(UX * UX + UY * UY + UZ * UZ);

    double VX = a(1) * b(2) - a(2) * b(1);
    double VY = a(2) * b(0) - a(0) * b(2);
    double VZ = a(0) * b(1) - a(1) * b(0);
    double FAC2 = sqrt(VX * VX + VY * VY + VZ * VZ);
    double BNZ = (VX * c(0) + VY * c(1) + VZ * c(2)) / FAC2;

    double BSXY = FAC1 / BOXX;

    double WX = c(0) - BNZ * VX / FAC2;
    double WY = c(1) - BNZ * VY / FAC2;
    double WZ = c(2) - BNZ * VZ / FAC2;
    double BSXZ = (WX * a(0) + WY * a(1) + WZ * a(2)) / BOXX;

    double BSYZ = (WX * UX + WY * UY + WZ * UZ) / BNY;

    length_x = BNX;
    length_y = BNY;
    length_z = BNZ;

    double BHX = 0.5 * BNX;
    double BHY = 0.5 * BNY;

    while (BSYZ > BHY) {
      BSYZ = BSYZ - BNY;
      BSXZ = BSXZ - BSXY;
    }

    while (BSYZ < -BHY) {
      BSYZ = BSYZ + BNY;
      BSXZ = BSXZ + BSXY;
    }

    while (BSXZ > BHX) {
      BSXZ = BSXZ - BNX;
    }
    while (BSXZ < -BHX) {
      BSXZ = BSXZ + BNX;
    }

    while (BSXY > BHX) {
      BSXY = BSXY - BNX;
    }

    while (BSXY < -BHX) {
      BSXY = BSXY + BNX;
    }

    delta_l_xy = BSXY;
    delta_l_xz = BSXZ;
    delta_l_yz = BSYZ;

    lambda_xy = delta_l_xy * num_x / length_x;
    lambda_xz = delta_l_xz * num_x / length_x;
    lambda_yz = delta_l_yz * num_y / length_y;

    // B.7
    if (lambda_xy + 0.5 >= 0.0)
      padding = int(lambda_xy + 0.5 - 1.0E-10);
    else
      padding = int(num_x + lambda_xy + 0.5 + 1.0E-10);

    dim_m = (num_x * num_y + padding) * num_z - padding;
    nxnyp = num_x * num_y + padding;

    make_code_table();
  } else {
    // constants needed for A.5
    length_x = math::abs(myconf->current().box(0));
    length_y = math::abs(myconf->current().box(1));
    length_z = math::abs(myconf->current().box(2));

    dim_m = num_x * num_y * num_z;

    padding = 0;
  }

  const double l_x = length_x / num_x;
  const double l_y = length_y / num_y;
  const double l_z = length_z / num_z;

  //Check, if cutoff is too big
  const double double_cutoff = 2 * cutoff_lr;
  if (double_cutoff > length_x - l_x || double_cutoff > length_y - l_y ||
          double_cutoff > length_z - l_z){
    io::messages.add("Cutoff or cell size too big!", "Grid Cell",
            io::message::error);
    return 1;
  }

  l_x2 = l_x * l_x;
  l_y2 = l_y * l_y;
  l_z2 = l_z * l_z;

  return 0;
}

/**
 * Create the mask pointer
 */
template<typename trait_type>
int interaction::Grid_Cell_Pairlist::make_mask_pointer(){

  DEBUG(5, "Grid Cell : Make mask pointer");
  // Create trait type
  trait_type t;

  mask_pointer.clear();

  for (int m = 1; m < dim_m; m++) {
    if (make_mask(m,t)) {
      DEBUG(10, "Start of mask = " << m);
      mask_pointer.push_back(m);

      while (make_mask(++m, t)){
        DEBUG(10, "m = " << m << " true");
        if (m == dim_m)
          break;
      }
      mask_pointer.push_back(m - 1);
      DEBUG(10, "\tEnd of mask = " << m-1);
    }
  }

  stripes = mask_pointer.size() / 2;
  return 0; 
}

/**
 * Creates the mask for rectangular shapes (as described in the appendix A)
 */
bool interaction::Grid_Cell_Pairlist::make_mask(unsigned int delta_m,
        Reg_Mask &t) {

  DEBUG(15, "Grid Cell : Make the mask for delta_m = " << delta_m);

  // A.7 - A.9
  const int delta_m_x = delta_m % num_x;
  const int delta_m_y = (delta_m % (num_x * num_y)) / num_x;
  const int delta_m_z = delta_m / (num_x * num_y);

  // A.14
  const int delta_n_x = abs(minIm(delta_m_x, num_x));
  // A.15
  int delta_n_y = 0;
  if ((delta_m_x == 0) ||
          ((delta_m_y == num_y_minus_1) && (delta_m_z == num_z_minus_1)))
    delta_n_y = abs(minIm(delta_m_y, num_y));
  else
    delta_n_y = std::min(abs(minIm(delta_m_y, num_y)),
          abs(minIm(delta_m_y + 1, num_y)));
  // A.16
  int delta_n_z = 0;
  if ((delta_m_z == num_z_minus_1) || ((delta_m_x == 0) && (delta_m_y == 0)))
    delta_n_z = abs(minIm(delta_m_z, num_z));
  else
    delta_n_z = std::min(abs(minIm(delta_m_z, num_z)),
          abs(minIm(delta_m_z + 1, num_z)));

  // A.5
  double p = (max(delta_n_x) - 1);
  p = p * p * l_x2;
  double s = p;
  p = (max(delta_n_y) - 1);
  p = p * p * l_y2;
  s += p;
  p = (max(delta_n_z) - 1);
  p = p * p * l_z2;
  s += p;

  return (s <= cutoff_lr2);
}

/**
 * Creates the mask for irregular shapes (as described in the appendix B)
 */
bool interaction::Grid_Cell_Pairlist::make_mask(int delta_m,
        Irr_Mask &t) {

  DEBUG(15, "Grid Cell : Make the mask for delta_m = " << delta_m);
  int delta_m_z = int(delta_m) / int(nxnyp);
  int j = int(delta_m) - int(delta_m_z) * int(nxnyp);
  int delta_m_y = j / int(num_x);
  int delta_m_x = j - int(delta_m_y) * int(num_x);
  int mz = delta_m_z;
  int my = delta_m_y;
  int mx = delta_m_x;

  int case_min = 0, case_max = 0;
  if (mz == num_z - 1) {
    if (my != num_y - 1 && mx != 0) {
      case_min = 1;
      case_max = 2;
    } else {
      case_min = 2;
      case_max = 2;
    }
  } else {
    if (my == 0) {
      if (mx != 0) {
        if (mx > padding) {
          case_min = 4;
          case_max = 6;
        } else {
          case_min = 1;
          case_max = 2;
        }
      } else {
        case_min = 2;
        case_max = 2;
      }
    } else if (my == 1) {
      if (mx != 0) {
        if (mx > padding) {
          case_min = 1;
          case_max = 4;
        } else {
          case_min = 1;
          case_max = 3;
        }
      } else {

        if (padding != 0) {
          case_min = 8;
          case_max = 9;
        } else {
          case_min = 2;
          case_max = 3;
        }
      }
    } else if (my == num_y - 1) {
      if (mx != 0) {
        if (mx > padding) {
          case_min = 2;
          case_max = 4;
        } else if (mx == padding) {
          case_min = 2;
          case_max = 3;
        } else {
          case_min = 6;
          case_max = 8;
        }
      } else {
        if (padding != 0) {
          case_min = 6;
          case_max = 8;
        } else {
          case_min = 2;
          case_max = 3;
        }
      }
    } else if (my == num_y) {
      if (mx > padding) {
        case_min = 3;
        case_max = 4;
      } else {
        case_min = 7;
        case_max = 8;
      }
    } else {
      if (mx != 0) {
        if (mx > padding) {
          case_min = 1;
          case_max = 4;
        } else if (mx == padding) {
          case_min = 1;
          case_max = 3;
        } else {
          case_min = 5;
          case_max = 8;
        }
      } else {
        if (padding != 0) {
          case_min = 6;
          case_max = 8;
        } else {
          case_min = 2;
          case_max = 3;
        }
      }
    }
  }

  //loop over cases
  for (int i = case_min; i <= case_max; i++) {
    int k = cases[i];
    //calculate the corresponding vector connecting the cells

    double d_mx = mx + a_code[k];
    double d_my = my + b_code[k];
    double d_mz = mz + c_code[k];

    //convert it to its minimum-image
    if (d_mz > num_z_half) {
      d_mz = d_mz - num_z;
      d_my = d_my - lambda_yz;
      d_mx = d_mx - lambda_xz;
    }

    if (d_my > num_y_half) {
      d_my = d_my - num_y;
      d_mx = d_mx - lambda_xy;
    } else if (d_my <= -num_y_half) {
      d_my = d_my + num_y;
      d_mx = d_mx + lambda_xy;
    }

    if (d_mx > num_x_half) {
      d_mx = d_mx - num_x;
      if (d_mx > num_x_half) {
        d_mx = d_mx - num_x;
      }
    } else if (d_mx <= -num_x_half) {
      d_mx = d_mx + num_x;
      if (d_mx <= -num_x_half) {
        d_mx = d_mx + num_x;
      }
    }

    d_mx = std::abs(d_mx);
    d_my = std::abs(d_my);
    d_mz = std::abs(d_mz);

    //for checking against the minimal possible distance for interaction

    if (d_mx < 1.0) {
      d_mx = 1.0;
    }

    if (d_my < 1.0) {
      d_my = 1.0;
    }

    if (d_mz < 1.0) {
      d_mz = 1.0;
    }

    d_mx = d_mx - 1.0;
    d_my = d_my - 1.0;
    d_mz = d_mz - 1.0;

    // for checking against the maximal possible distance for interaction
    double r2 = d_mx * d_mx * l_x2
            + d_my * d_my * l_y2
            + d_mz * d_mz * l_z2;

    if (r2 <= cutoff_lr2)
      return true;

  } // loop
  return false;
}

/**
 * make the cell array and the cell pointer
 */
template<math::boundary_enum b>
int interaction::Grid_Cell_Pairlist::make_cell(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim) {

  DEBUG(5, "Grid Cell : Make the cell");

  cell.clear();
  cell_pointer.clear();

  // Some constants for making it faster
  const int nxny = num_x * num_y;
  const double nx_p_lx = num_x / length_x;
  const double ny_p_ly = num_y / length_y;
  const double nz_p_lz = num_z / length_z;

  math::Periodicity<b> periodicity(conf.current().box);
  math::Vec v;
  const unsigned int num_atoms = m_cg_cog.size();
  cell_element ce;
  for (ce.i = 0; ce.i < num_atoms; ++ce.i) {
    v = m_cg_cog(ce.i);
    periodicity.put_into_positive_box(v);
    if (irregular_shape) {
      put_into_brickwall(v);

      // B.7 without "+1"
      ce.m = nxnyp * int (nz_p_lz * v(2))
              + num_x * int (ny_p_ly * v(1))
              + int (nx_p_lx * v(0));
    } else {
      // Eq. (5) without "+1"
      ce.m = nxny * int(nz_p_lz * v(2))
              + num_x * int(ny_p_ly * v(1))
              + int(nx_p_lx * v(0));
    }

    cell.push_back(ce);
  } // for all charge groups/atoms

  // sort the vector
  DEBUG(15, "Grid Cell : Sort cell")
  std::sort(cell.begin(), cell.end());

  // Make the cell pointer
  DEBUG(5, "Grid Cell : Make the cell pointer")

  int last_m = 0;
  int current_m = 0;
  cell_pointer.push_back(0);

  for (unsigned int j = 0; j < cell.size();){
    DEBUG(15, "j = " << j << " Size of cell = " << cell.size());
    current_m = cell[j].m;
    DEBUG(15, "Current m = " << current_m << " last m = " << last_m);
    while(last_m < current_m){
      last_m++;
      cell_pointer.push_back(j);
    }

    DEBUG(15, "After first while");

    j++;         
    while(j < cell.size() && cell[j].m == last_m)
      j++;
    DEBUG(15, "After second while");
  } // for all charge groups

  // Add an additional entry to the cell pointer (see paper)
  DEBUG(10, "Grid Cell : Add end element")
  const int diff = dim_m - cell_pointer.size() + 1;
  if (diff <= 0)
    return 1;

  cell_element end;
  end.m = dim_m; 
  DEBUG(13, "num_cg = " << num_atoms);
  end.i = num_atoms;
  cell.push_back(end);
  cell_pointer.insert(cell_pointer.end(), diff, num_atoms);
  DEBUG(10, "Size cell = " << cell.size() << " Size cell pointer = " << cell_pointer.size());

  /*
  DEBUG(1, "Show me the cell")
  std::vector<interaction::Grid_Cell_Pairlist::cell_element>::iterator itc = cell.begin();
  std::vector<interaction::Grid_Cell_Pairlist::cell_element>::iterator toc = cell.end();
  for (unsigned int i = 0; itc != toc; ++itc, ++i){
    DEBUG(1, "i = " << i << " m = " << cell[i].m << " atom : " << cell[i].i)
  }

  DEBUG(1, "Show me the cell pointer")
  std::vector<unsigned int>::iterator it = cell_pointer.begin();
  std::vector<unsigned int>::iterator to = cell_pointer.end();
  for (unsigned int i = 0; it != to; ++it, ++i){
    DEBUG(1, i << "  " << cell[cell_pointer[i]].i)
  } */

  return 0;
}

/**
 * make the pairlist
 */
template<math::boundary_enum b, class cutoff_trait, class perturbation_trait>
int interaction::Grid_Cell_Pairlist::_pairlist(
        interaction::PairlistContainer & pairlist,
        interaction::PairlistContainer & perturbed_pairlist,
        unsigned int offset, unsigned int stride,
        const cutoff_trait & cutoff,
        const perturbation_trait & perturbation){

  DEBUG(12, "Grid Cell : Go through the mask and pair the atoms")
  DEBUG(12, "Number of Charge groupes : " << mytopo->num_chargegroups())
  math::Periodicity<b> periodicity(myconf->current().box);

  // for all cells
  DEBUG(12, "Size of the cell pointer : " << cell_pointer.size())
  DEBUG(12, "Size of the cell : " << cell.size());
  for (int m = offset; m < dim_m; m += stride) {
    DEBUG(15, "m: " << m);

    // for the primary atoms in the cell
    const int n1_start = cell_pointer[m];
    const int n_stop = cell_pointer[m + 1];
    DEBUG(15, "N start: " << n1_start << " end: " << n_stop);
    for (int n1 = n1_start; n1 < n_stop; n1++) {
      for (int n2 = n1; n2 < n_stop; n2++){
        DEBUG(15, "\tFirst pairing")
        DEBUG(15, "f.p n1 = " << n1 << " n2 = " << n2)
        pair<b>(n1, n2, pairlist, perturbed_pairlist, periodicity, cutoff, perturbation); // pair them
      }
      // for all stripes
      DEBUG(15, "Number of stripes = " << stripes);
      for (int s = 0; s < stripes; s++) {
        const int m1 = m + mask_pointer[2 * s];
        if (m1 < dim_m) {
          int min = m + mask_pointer[2 * s + 1];
          int m2 = min < dim_m ? (min + 1) : dim_m;
          const int n2_start = cell_pointer[m1];
          const int n2_end = cell_pointer[m2];
          DEBUG(15, "\tSecond Pairing");
          DEBUG(15, "m1 : " << m1);
          DEBUG(15, "m2 : " << m2);
          DEBUG(15, "m : " << m);
          DEBUG(15, "s : " << s << " of " << stripes);
          DEBUG(15, "mask_pointer[2 * s] : " << mask_pointer[2 * s]);
          DEBUG(15, "mask_pointer[2 * s + 1] : " << mask_pointer[2 * s + 1]);
          DEBUG(15, "dim_m : " << dim_m);
          DEBUG(15, "size of cellpointer : " << cell_pointer.size())
          DEBUG(15, "cell_pointer[m2] : " << cell_pointer[m2]);
          DEBUG(15, "n2 end : " << n2_end);
          for (int n2 = n2_start; n2 < n2_end; n2++) {
            DEBUG(15, "s.p n1 = " << n1 << " n2 = " << n2);
            pair<b>(n1, n2, pairlist, perturbed_pairlist, periodicity, cutoff, perturbation);
          }
        } // }
      } // for all stripes
    } // for the primary atoms in the cell
  } // for all cells

  return 0;
}

/**
 * put the  charge groups inside the cell into the pairlist
 */
template<math::boundary_enum b, class cutoff_trait, class perturbation_trait>
inline int interaction::Grid_Cell_Pairlist::pair(
        unsigned int n1, unsigned int n2,
        interaction::PairlistContainer & pairlist,
        interaction::PairlistContainer & perturbed_pairlist,
        const math::Periodicity<b> &periodicity, 
        const cutoff_trait & cutoff,
        const perturbation_trait & perturbation) {
  unsigned int first = cell[n1].i;
  unsigned int second = cell[n2].i;

  const simulation::qmmm_enum qmmm = mysim->param().qmmm.qmmm;

  DEBUG(15, "First : " << first << " second : " << second);
  if (first > second) {
    std::swap(first, second);
  }

  DEBUG(12, "Grid Cell : Pair " << first << " and " << second);

  // Are they shortrange or not?
  math::Vec r;
  DEBUG(15, "Size of m_cg_cog = " << m_cg_cog.size() << " Num CG : " <<
          mytopo->num_chargegroups());
  
  // Solvent - solvent
  if (first >= first_solvent) {
    DEBUG(15, "Solvent - Solvent");
#ifdef HAVE_LIBCUDART
    if (cuda)
      return 0;
#endif
    periodicity.nearest_image(m_cg_cog(first), m_cg_cog(second), r);
    const double d = math::abs2(r);
    if (d <= cutoff_sr2)
      pair_solvent(first, second, pairlist.solvent_short, cutoff);
    else if (d <= cutoff_lr2)
      pair_solvent(first, second, pairlist.solvent_long, cutoff);
  }// Solvent - solute
  else if (second >= first_solvent) {
    DEBUG(15, "Solute - Solvent");
    if (t_qm_excluded<cutoff_trait>(*mytopo, qmmm, first, second)) {
			DEBUG(9, "Skipping pair: " << first << "-" << second);
      return 0;
    }
    periodicity.nearest_image(m_cg_cog(first), m_cg_cog(second), r);
    const double d = math::abs2(r);
    if (d <= cutoff_sr2)
      pair_solute_solvent(first, second, pairlist.solute_short,
            perturbed_pairlist.solute_short, cutoff, perturbation);
    else if (d <= cutoff_lr2)
      pair_solute_solvent(first, second, pairlist.solute_long,
            perturbed_pairlist.solute_long, cutoff, perturbation);
  }// Solute - solute
  else {
    DEBUG(15, "Solute - Solute");
    if (t_qm_excluded<cutoff_trait>(*mytopo, qmmm, first, second)) {
			DEBUG(9, "Skipping pair: " << first << "-" << second);
      return 0;
    }
    periodicity.nearest_image(m_cg_cog(first), m_cg_cog(second), r);
    // the distance
    const double d = math::abs2(r);
    if (d <= cutoff_sr2) 
      pair_solute(first, second, pairlist.solute_short,
            perturbed_pairlist.solute_short, cutoff, perturbation);
    else if (d <= cutoff_lr2)
      pair_solute(first, second, pairlist.solute_long,
            perturbed_pairlist.solute_long, cutoff, perturbation);
    }

  return 0;
}

/**
 * put the atoms of a chargegroup into the pairlist
 */
template<class perturbation_trait>
inline int interaction::Grid_Cell_Pairlist::pair_solute_solvent(const unsigned int first,
        const unsigned int second,
        interaction::Pairlist & pairlist,
        interaction::Pairlist & perturbed_pairlist,
        const cg_cutoff & cutoff,
        const perturbation_trait & perturbation) {

  assert(first <= second);
  DEBUG(15, "Different charge group");
  for (int a1 = mytopo->chargegroup(first),
          a1_to = mytopo->chargegroup(first + 1);
          a1 < a1_to; ++a1) {
    for (int a2 = mytopo->chargegroup(second),
            a2_to = mytopo->chargegroup(second + 1); a2 < a2_to; ++a2) {


      assert(int(pairlist.size()) > a1);
      insert_pair(pairlist, perturbed_pairlist, a1, a2, perturbation);
    }
  }

  return 0;
}

/**
 * put the atoms into the pairlist
 */
template<class perturbation_trait>
inline int interaction::Grid_Cell_Pairlist::pair_solute_solvent(const unsigned int first,
        const unsigned int second,
        interaction::Pairlist & pairlist,
        interaction::Pairlist & perturbed_pairlist,
        const atomic_cutoff & cutoff,
        const perturbation_trait & perturbation) {
  assert(first <= second);
  assert(pairlist.size() > first);
  insert_pair(pairlist, perturbed_pairlist, first, second, perturbation);
  return 1;
}

/**
 * put the atoms of a chargegroupe into the pairlist
 */
template<class perturbation_trait>
inline int interaction::Grid_Cell_Pairlist::pair_solute(const unsigned int first,
        const unsigned int second,
        interaction::Pairlist & pairlist,
        interaction::Pairlist & perturbed_pairlist,
        const cg_cutoff & cutoff,
        const perturbation_trait & perturbation) {

  assert(first <= second);
  // If they are from the same charge group
  if (first == second) {
    DEBUG(15, "Same charge group");
    for (int a1 = mytopo->chargegroup(first),
            a_to = mytopo->chargegroup(first + 1);
            a1 < a_to; ++a1) {
      for (int a2 = a1 + 1; a2 < a_to; ++a2) {
          // check it is not excluded
          if (excluded_solute_pair(*mytopo, a1, a2))
            continue;
        assert(int(pairlist.size()) > a1);
        insert_pair(pairlist, perturbed_pairlist, a1, a2, perturbation);
      }
    }
  } else { // different charge groups
    DEBUG(15, "Different charge group");
    for (int a1 = mytopo->chargegroup(first),
            a1_to = mytopo->chargegroup(first + 1);
            a1 < a1_to; ++a1) {
      for (int a2 = mytopo->chargegroup(second),
              a2_to = mytopo->chargegroup(second + 1); a2 < a2_to; ++a2) {
          // check if is not excluded
        if (excluded_solute_pair(*mytopo, a1, a2))
          continue;

        assert(int(pairlist.size()) > a1);
        insert_pair(pairlist, perturbed_pairlist, a1, a2, perturbation);
      }
    }
  }

  return 0;
}

/**
 * put the atoms into the pairlist
 */
template<class perturbation_trait>
inline int interaction::Grid_Cell_Pairlist::pair_solute(const unsigned int first,
        const unsigned int second,
        interaction::Pairlist & pairlist,
        interaction::Pairlist & perturbed_pairlist,
        const atomic_cutoff & cutoff,
        const perturbation_trait & perturbation) {
  assert(first <= second);
  // check if is not excluded
  if (first == second || excluded_solute_pair(*mytopo, first, second))
    return 0;

  assert(pairlist.size() > first);
  insert_pair(pairlist, perturbed_pairlist, first, second, perturbation);

  return 1;
}

inline int interaction::Grid_Cell_Pairlist::pair_solvent(const unsigned int first,
        const unsigned int second,
        interaction::Pairlist & pairlist,
        const cg_cutoff & cutoff) {

  assert(first <= second);
  assert(first >= first_solvent);
  // If they are from the same charge group
  if (first == second)
    return 0;
  DEBUG(15, "Different charge group");
  for (int a1 = mytopo->chargegroup(first),
          a1_to = mytopo->chargegroup(first + 1);
          a1 < a1_to; ++a1) {
    for (int a2 = mytopo->chargegroup(second),
            a2_to = mytopo->chargegroup(second + 1); a2 < a2_to; ++a2) {

      assert(int(pairlist.size()) > a1);
      pairlist[a1].push_back(a2);
    }
  }
  return 1;
}

inline int interaction::Grid_Cell_Pairlist::pair_solvent(const unsigned int first,
        const unsigned int second,
        interaction::Pairlist &pairlist, const atomic_cutoff & cutoff) {

  assert(first <= second);
  assert(first >= first_solvent);
  // same solvent mol?
  if ((first - first_solvent) / num_atoms_per_solvent == (second - first_solvent) / num_atoms_per_solvent)
    return 0;

  assert(pairlist.size() > first);
  pairlist[first].push_back(second);

  return 1;
}

inline void interaction::Grid_Cell_Pairlist::insert_pair(
            interaction::Pairlist & pairlist,
            interaction::Pairlist & perturbed_pairlist,
            int first, int second, const no_perturbation & perturbation) {
  pairlist[first].push_back(second);
}

inline void interaction::Grid_Cell_Pairlist::insert_pair(
            interaction::Pairlist & pairlist,
            interaction::Pairlist & perturbed_pairlist,
            int first, int second, const do_perturbation & perturbation) {
  assert(first < second);
  
  if (mytopo->is_perturbed(first) || mytopo->is_eds_perturbed(first))
    perturbed_pairlist[first].push_back(second);
  else if (mytopo->is_perturbed(second) || mytopo->is_eds_perturbed(second))
    perturbed_pairlist[second].push_back(first);
  else
    pairlist[first].push_back(second);
}

/**
 * Check, if two atoms are exclude from the pairlist
 */
inline bool interaction::Grid_Cell_Pairlist::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  return topo.all_exclusion(i).is_excluded(j);
}

/**
 * Makes the combination table of the codes for irregular shapes aswell
 * as the corresponding code tables for the parameters
 * Table 2, p. 1484
 */
inline void interaction::Grid_Cell_Pairlist::make_code_table(){
  // Table 2
  a_code[0] = 0; // place holder for fortran
  a_code[1] = 0;
  a_code[2] = 0;
  a_code[3] = -padding;
  a_code[4] = -padding;
  a_code[5] = num_x - padding;
  a_code[6] = num_x - padding;

  b_code[0] = 0; // place holder for fortran
  b_code[1] = 0;
  b_code[2] = 1;
  b_code[3] = -num_y;
  b_code[4] = 1 - num_y;
  b_code[5] = -num_y;
  b_code[6] = -1 - num_y;

  c_code[0] = 0; // place holder for fortran
  c_code[1] = 0;
  c_code[2] = 0;
  c_code[3] = 1;
  c_code[4] = 1;
  c_code[5] = 1;
  c_code[6] = 1;
}

/**
 * returns the maximum of n or 1
 */
template<typename t>
inline t interaction::Grid_Cell_Pairlist::max(t n){
  return n > 1 ? n : 1;
}

/**
 * returns the minimum image value of its first argument based on the
 * periodicity defined by its second argument.
 */
inline int interaction::Grid_Cell_Pairlist::minIm(int n, int num) {
   const int num_half = num / 2;
   if (n < -num_half)
     return n + num;
   else if (n > num_half)
     return n - num;
   else
     return n;
}

void interaction::Grid_Cell_Pairlist::put_into_brickwall(math::Vec &v) {
  double X1 = v(0);
  double X2 = v(1);
  double X3 = v(2);

  const double & BNZ = length_z;
  const double & BNY = length_y;
  const double & BNX = length_x;
  const double & BSXY = delta_l_xy;
  const double & BSXZ = delta_l_xz;
  const double & BSYZ = delta_l_yz;

  while (X3 >= BNZ) {
    X3 = X3 - BNZ;
    X2 = X2 - BSYZ;
    X1 = X1 - BSXZ;
  }

  while (X3 <= 0.0) {
    X3 = X3 + BNZ;
    X2 = X2 + BSYZ;
    X1 = X1 + BSXZ;
  }

  while (X2 >= BNY) {
    X2 = X2 - BNY;
    X1 = X1 - BSXY;
  }

  while (X2 <= 0.0) {
    X2 = X2 + BNY;
    X1 = X1 + BSXY;
  }

  while (X1 >= BNX) {
    X1 = X1 - BNX;
  }

  while (X1 <= 0.0) {
    X1 = X1 + BNX;
  }

  v(0) = X1;
  v(1) = X2;
  v(2) = X3;
}

void interaction::Grid_Cell_Pairlist::create_vacuum_box() {
  const math::VArray & pos = myconf->current().pos;
  math::Vec mn(pos(0));
  math::Vec mx(pos(0));
  for(unsigned int i = 1; i < mytopo->num_atoms(); ++i) {
    for(unsigned int c = 0; c < 3; ++c) {
      mn(c) = std::min(mn(c), pos(i)(c));
      mx(c) = std::max(mx(c), pos(i)(c));
    }
  }
  math::Box new_box(
          math::Vec(mx(0) - mn(0) + 2.0 * mysim->param().pairlist.cutoff_long, 0.0, 0.0),
          math::Vec(0.0, mx(1) - mn(1) + 2.0 *  mysim->param().pairlist.cutoff_long, 0.0),
          math::Vec(0.0, 0.0, mx(2) - mn(2) + 2.0 * mysim->param().pairlist.cutoff_long));
  box_backup = myconf->current().box;
  myconf->current().box = new_box;
  myconf->boundary_type = math::rectangular;
}

void interaction::Grid_Cell_Pairlist::restore_vacuum_box() {
  myconf->current().box = box_backup;
  myconf->boundary_type = math::vacuum;
}

template<class cutoff_trait>
inline bool interaction::Grid_Cell_Pairlist::t_qm_excluded(const topology::Topology& topo
                                                         , const simulation::qmmm_enum qmmm
                                                         , unsigned first) {
  if (std::is_same<cutoff_trait, Grid_Cell_Pairlist::cg_cutoff>::value) {
    first = topo.chargegroup(first);
  } else if (! std::is_same<cutoff_trait, Grid_Cell_Pairlist::atomic_cutoff>::value)
    io::messages.add("Cutoff trait not implemented",
              "Grid Cell", io::message::critical);
  return Pairlist_Algorithm::qm_excluded(topo, qmmm, first);
}

template<class cutoff_trait>
inline bool interaction::Grid_Cell_Pairlist::t_qm_excluded(const topology::Topology& topo
                                                         , const simulation::qmmm_enum qmmm
                                                         , unsigned first
                                                         , unsigned second){
  if (std::is_same<cutoff_trait, Grid_Cell_Pairlist::cg_cutoff>::value) {
    first = topo.chargegroup(first);
    second = topo.chargegroup(second);
  } else if (! std::is_same<cutoff_trait, Grid_Cell_Pairlist::atomic_cutoff>::value)
    io::messages.add("Cutoff trait not implemented",
              "Grid Cell", io::message::critical);
  return Pairlist_Algorithm::qm_excluded(topo, qmmm, first, second);
}
