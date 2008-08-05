/**
 * @file latticesum.cc
 * implementation of lattice sum methods
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>
#include <interaction/nonbonded/interaction/storage.h>
#include <math/periodicity.h>
#include <math/volume.h>

#include <interaction/nonbonded/interaction/latticesum.h>

#include <util/error.h>
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE nonbonded

void interaction::Lattice_Sum::calculate_charge_density(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const math::VArray & r) {
  
  configuration::Mesh & charge_density = conf.lattice_sum().charge_density;
  charge_density.zero();
  
  const unsigned int Nx = charge_density.x();
  const unsigned int Ny = charge_density.y();
  const unsigned int Nz = charge_density.z();
  
  // the H matrix to transform the grid into coordinates
  const math::Box &box = conf.current().box;
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  // and coordinates into the grid
  math::Matrix H_inv=math::inverse(H);
  math::Matrix H_trans = math::transpose(H);
  // the volume of a grid cell (Vg) MD05.32 eq. 61
  const double cell_volume = math::det(H_trans); 
  DEBUG(11, "grid cell volume: " << cell_volume);
  math::Matrix H_inv_trans = math::transpose(H_inv);
 
  const unsigned int num_atoms = topo.num_atoms();
  
  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  bool assignment_odd = assignment_function_order % 2 == 1;
  
  DEBUG(10,"\t starting to assign charge density to grid ... ");
  
  // do integer divion to calculate how far you have to go away from the first grid point
  const int upper_bound = assignment_function_order / 2;
  const int lower_bound = - (assignment_function_order - 1) / 2;
  DEBUG(10, "Upper bound: " << upper_bound << " lower bound: " << lower_bound);
  
  // loop over all charges
  for (unsigned int i = 0; i < num_atoms; ++i) {
    const double charge = topo.charge(i);
    if (charge == 0.0)
      continue;
    
    DEBUG(12, "\tAssigning charge " << i + 1 << " (" << charge << ") to grid.");

    // we have a charge so let's calculate were it is on the grid.
    math::Vec grid_r = math::product(H_inv, r(i));
    DEBUG(15, "\t\tcoordinates on grid: " << math::v2s(grid_r));
    math::GenericVec<int> nearest_grid_point;
    
    // getting the nearest grid point is dependend on the order of the
    // charge assignment function
    if (assignment_odd) { // assignment order is odd
      // round to the nearest integer to get the nearest grid point
      nearest_grid_point = math::GenericVec<int>(int(rint(grid_r(0))),
               int(rint(grid_r(1))), int(rint(grid_r(2))));
    } else { // assignment order is even
      // round to the nearest integer to get the nearest grid centre
      // we have to subtract math::epsilon (a small number) from the coorinates
      nearest_grid_point = math::GenericVec<int>(int(grid_r(0)),
                                                 int(grid_r(1)), 
                                                 int(grid_r(2)));
    }
    DEBUG(15, "\tnearest grid point: " << math::v2s(nearest_grid_point));
    math::GenericVec<int> point;
    math::Vec point_to_charge;
    double w_x, w_y, w_z;
    for (int dx = lower_bound; dx <= upper_bound; ++dx) {
      point(0) = nearest_grid_point(0) + dx;
      point_to_charge(0) = grid_r(0) - point(0);
      w_x = interaction::Lattice_Sum::charge_assignment_1d(assignment_function_order, point_to_charge(0));
      for (int dy = lower_bound; dy <= upper_bound; ++dy) {
        point(1) = nearest_grid_point(1) + dy;
        point_to_charge(1) = grid_r(1) - point(1);
        w_y = interaction::Lattice_Sum::charge_assignment_1d(assignment_function_order, point_to_charge(1));
        for (int dz = lower_bound; dz <= upper_bound; ++dz) {
          point(2) = nearest_grid_point(2) + dz;
          point_to_charge(2) =  grid_r(2) - point(2);
          w_z = interaction::Lattice_Sum::charge_assignment_1d(assignment_function_order, point_to_charge(2));

          double assignment_function = w_x * w_y * w_z / cell_volume;
          charge_density(point) += assignment_function * charge;
          DEBUG(15, "\t\tgrid point: [" << point(0) << "," << point(1) << "," << point(2) << 
                  "] assignment function: " << assignment_function);
          // output multiplied by 11.78708615 for comparison with promd
          DEBUG(15,"\t\tassignment_function * charge =" << assignment_function * charge* 11.78708615);
          DEBUG(15, "\t\t\tcharge density: " << charge_density(point)* 11.78708615);
        }
      }
    } // loop over grid cells
  } // loop over charges

#ifndef NDEBUG
  for (int z = 0; z < int(Nz); ++z) {
    for (int y = 0; y < int(Ny); ++y) {
      for (int x = 0; x < int(Nx); ++x) {
        DEBUG(20, "charge density (" << x << "," << y << "," << z << ") = " << charge_density(x,y,z).real());
        std::cout.precision(8);
        // output multiplied by 11.78708615 for comparison with promd
        DEBUG(25, x+1 << " " << y+1 << " " << z+1 << std::setw(15) << charge_density(x,y,z).real() * 11.78708615);
      }
    }
  }
#endif 
}
void interaction::Lattice_Sum::calculate_potential_and_energy(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        interaction::Storage & storage) {
  
  
  configuration::Mesh & charge_density = conf.lattice_sum().charge_density;
  configuration::Influence_Function & influence_function = conf.lattice_sum().influence_function;
  configuration::Mesh & potential = conf.lattice_sum().potential;
  
  //assert(charge_density.x() == influence_function.x() && charge_density.x() == potential.x());
  //assert(charge_density.y() == influence_function.y() && charge_density.y() == potential.y());
  //assert(charge_density.z() == influence_function.z() && charge_density.z() == potential.z());
  
  const int Nx = charge_density.x();
  const int Ny = charge_density.y();
  const int Nz = charge_density.z();
  const int grid_volume = Nx * Ny * Nz;
  
  const math::Box &box = conf.current().box;
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  math::Matrix H_trans = math::transpose(H);
  // the volume of a grid cell (Vg) MD05.32 eq. 61
  const double cell_volume = math::det(H_trans); 
  const double sqrt_grid_volume = sqrt(double(grid_volume));
  // now calculate the energy MD05.32 eq (72)
  double energy = 0.0;
  // and the potential MD05.32 eq. 68
  // loop over grid
  for (int x = 0; x < Nx; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        // calculate the potential by MD05.32 eq. 68
        DEBUG(15, "x: " << x + 1 << " y: " << y + 1 << " z: " << z + 1);
        DEBUG(15, "influence_function(hat): " << influence_function(x, y, z));
        DEBUG(15, "charge density(hat): " << charge_density(x, y, z)  / sqrt_grid_volume);
        potential(x, y, z) = math::eps0_i * influence_function(x, y, z) * charge_density(x, y, z) / sqrt_grid_volume;
        // factors in DEBUG statements are for better comparison with promd
        DEBUG(15, "\tpotential (hat)(" << x << "," << y << "," << z << "): " << potential(x,y,z) / math::eps0_i * 11.78708615);
        
        // calculate the energy by MD05.32 eq. 72     
        const std::complex<double> density = charge_density(x, y, z);
        const double abs2_charge_density = density.real() * density.real() + density.imag() * density.imag();
        energy += influence_function(x, y, z) * abs2_charge_density / double(grid_volume);
        DEBUG(15, "\t\t abs_charge_density = " << abs2_charge_density / sqrt_grid_volume
                << ", energy = " << energy);
      }
    }
  } // loop over grid
  
  // See MD05.32 eq. 64, 65 and kpppm.F about this prefactor issue.
  // the 1/(total volume) prefactor is obtained by the FFT automtically
  storage.energies.ls_kspace_total = 0.5 * energy * cell_volume * math::eps0_i;
  DEBUG(10,"ls_kspace_total = " << storage.energies.ls_kspace_total);

}

void interaction::Lattice_Sum::calculate_electric_field(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim) {
  configuration::Mesh & electric_field_x = conf.lattice_sum().electric_field.x;
  configuration::Mesh & electric_field_y = conf.lattice_sum().electric_field.y;
  configuration::Mesh & electric_field_z = conf.lattice_sum().electric_field.z;
    
  configuration::Mesh & potential = conf.lattice_sum().potential;

  const int Nx = potential.x();
  const int Ny = potential.y();
  const int Nz = potential.z();
  const int grid_volume = Nx * Ny * Nz;
  
  const int finite_difference_order = sim.param().nonbonded.p3m_finite_differences_operator;
  std::vector<double> finite_difference_constants = interaction::Lattice_Sum::finite_difference_constants(finite_difference_order);

  const int Nx_2 = Nx / 2;
  const int Ny_2 = Ny / 2;
  const int Nz_2 = Nz / 2;
  const math::Box &box = conf.current().box;
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  const double H_x = math::abs(box(0) / Nx);
  const double H_y = math::abs(box(1) / Ny);
  const double H_z = math::abs(box(2) / Nz);
  DEBUG(15,"H_x = " << H_x << ", H_y = " << H_y << ", H_z = " << H_z);
  math::Matrix H_trans = math::transpose(H);
  // the volume of a grid cell (Vg) MD05.32 eq. 61
  //const double cell_volume = math::det(H_trans);
  const double sqrt_grid_volume = sqrt(double(grid_volume));
  // this is MD.05.32 eq. 35 and 36
  const double volume = math::volume(box, sim.param().boundary.boundary);
  
  const double volume_i_2pi = 2.0 * math::Pi / volume;
  math::Matrix tL_inv_2pi(
          math::cross(box(1), box(2)) * volume_i_2pi,
          math::cross(box(2), box(0)) * volume_i_2pi,
          math::cross(box(0), box(1)) * volume_i_2pi, true /* column wise */);


if (finite_difference_order == 0) {
    DEBUG(10, "Starting ik-differentiation");

    std::complex<double> minus_im(0, -1);
    // loop over grid
    // MD05.32 eq. 76
    math::GenericVec<int> l;
    for (int l1 = 1 - Nx_2; l1 <= Nx_2; ++l1) {
      l(0) = l1;
      for (int l2 = 1 - Ny_2; l2 <= Ny_2; ++l2) {
        l(1) = l2;
        for (int l3 = 1 - Nz_2; l3 <= Nz_2; ++l3) {
          l(2) = l3;
          // MD05.32 eq. 76
          DEBUG(16, "\tl: " << math::v2s(l));                  
          const math::Vec k_l = math::product(tL_inv_2pi, l);
          DEBUG(16, "\tk_l: " << math::v2s(k_l));
          // here we have to divide through sqrt_grid_volume because of numerical
          // FFT algorithms. The electric field (hat) will be scaled but after application
          // of the backward FFT we obtain the correct electrostatic field.
          const std::complex<double> minus_im_pot = minus_im * potential(l) / sqrt_grid_volume;
          electric_field_x(l) = k_l(0) * minus_im_pot;
          electric_field_y(l) = k_l(1) * minus_im_pot;
          electric_field_z(l) = k_l(2) * minus_im_pot;         
        }
      }
    }
#ifndef NDEBUG
    for (int z = 0; z < int(Nz); ++z) {
      for (int y = 0; y < int(Ny); ++y) {
        for (int x = 0; x < int(Nx); ++x) {
          std::cout.precision(5);
          // factors for comparison with promd
          DEBUG(25, " " << x + 1 << " " << y + 1 << " " << z + 1 <<
                  std::setw(12) << electric_field_x(x, y, z).real() * 11.78708615 * sqrt_grid_volume / math::eps0_i <<
                  std::setw(12) << electric_field_x(x, y, z).imag() * 11.78708615 * sqrt_grid_volume / math::eps0_i  <<
                  std::setw(12) << electric_field_y(x, y, z).real() * 11.78708615 * sqrt_grid_volume / math::eps0_i  <<
                  std::setw(12) << electric_field_y(x, y, z).imag() * 11.78708615 * sqrt_grid_volume / math::eps0_i  <<
                  std::setw(12) << electric_field_z(x, y, z).real() * 11.78708615 * sqrt_grid_volume / math::eps0_i  <<
                  std::setw(12) << electric_field_z(x, y, z).imag() * 11.78708615 * sqrt_grid_volume / math::eps0_i );
        }
      }
    }
#endif    
    DEBUG(11, "Starting FFT of electric field...");
    electric_field_x.fft(configuration::Mesh::fft_backward);
    electric_field_y.fft(configuration::Mesh::fft_backward);
    electric_field_z.fft(configuration::Mesh::fft_backward);
    DEBUG(11, "\tdone");
  } else {
    DEBUG(10, "Starting finite-differentiation");
    // backward fft the potential
    potential.fft(configuration::Mesh::fft_backward);
    
    // now loop over the grid
    // running indices MUST BE integers and NOT unsigned integers! Oh I hate it!!
    for (int x = 0; x < Nx; ++x) {
     for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          DEBUG(15, "\tpotential (" << x << "," << y << "," << z << "): " << potential(x,y,z));
          // loop over neighboring grid cells
          math::Vec electric_field(0.0);
          for (int q = 1; q <= finite_difference_order; ++q) {
            // get the q-th neighbor grid points
            const math::GenericVec<int> neighbors[6] = {
              math::GenericVec<int>(x - q, y, z), math::GenericVec<int>(x + q, y, z), // x-q will fail if x is unsigned!!!!!
              math::GenericVec<int>(x, y - q, z), math::GenericVec<int>(x, y + q, z),
              math::GenericVec<int>(x, y, z - q), math::GenericVec<int>(x, y, z + q)
            };
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[0]) << std::setw(15) << potential(neighbors[0]).real());
            DEBUG(18, "\t\t2: " << math::v2s(neighbors[1]) << std::setw(15) << potential(neighbors[1]).real());
            DEBUG(18, "\t\t3: " << math::v2s(neighbors[2]) << std::setw(15) << potential(neighbors[2]).real());
            DEBUG(18, "\t\t4: " << math::v2s(neighbors[3]) << std::setw(15) << potential(neighbors[3]).real());
            DEBUG(18, "\t\t5: " << math::v2s(neighbors[4]) << std::setw(15) << potential(neighbors[4]).real());
            DEBUG(18, "\t\t6: " << math::v2s(neighbors[5]) << std::setw(15) << potential(neighbors[5]).real());

            // calculate Delta Potential for all different access and divide it
            // through the distance between the grid points
            math::Vec delta_pot(
                    (potential(neighbors[0]) - potential(neighbors[1])).real() / (2.0 * H_x * q),
                    (potential(neighbors[2]) - potential(neighbors[3])).real() / (2.0 * H_y * q),
                    (potential(neighbors[4]) - potential(neighbors[5])).real() / (2.0 * H_z * q));      
            electric_field += (finite_difference_constants[q-1] / sqrt_grid_volume) * delta_pot;
            
            DEBUG(14,"delta_pot = " << math::v2s(delta_pot * (finite_difference_constants[q-1] / sqrt_grid_volume)));
            DEBUG(14,"finite_difference_constants[" << q << "] = " << finite_difference_constants[q-1]);

          } // loop over finite difference order
          electric_field_x(x, y, z) = electric_field(0);
          electric_field_y(x, y, z) = electric_field(1);
          electric_field_z(x, y, z) = electric_field(2);
        } 
      }
    } // loop over grid
  }

#ifndef NDEBUG
  for (int z = 0; z < Nz; ++z) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = 0; x < Nx; ++x) {
        std::cout.precision(5);
        DEBUG(20, x+1 << " " << y+1 << " " << z+1 <<
                std::setw(12) << electric_field_x(x, y, z).real() / sqrt(4.0 * math::Pi * math::eps0_i) << 
                std::setw(12) << electric_field_y(x, y, z).real() / sqrt(4.0 * math::Pi * math::eps0_i) << 
                std::setw(12) << electric_field_z(x, y, z).real() / sqrt(4.0 * math::Pi * math::eps0_i));
      }
    }
  }
#endif
}

void interaction::Lattice_Sum::calculate_force(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        Storage & storage,
        const math::VArray & r) {

  configuration::Mesh & electric_field_x = conf.lattice_sum().electric_field.x;
  configuration::Mesh & electric_field_y = conf.lattice_sum().electric_field.y;
  configuration::Mesh & electric_field_z = conf.lattice_sum().electric_field.z;
  
  assert(electric_field_x.x() == electric_field_y.x() && electric_field_x.x() == electric_field_z.x());
  assert(electric_field_x.x() == electric_field_y.y() && electric_field_x.y() == electric_field_z.y());
  assert(electric_field_x.x() == electric_field_y.z() && electric_field_x.z() == electric_field_z.z());
  
  const int Nx = electric_field_x.x();
  const int Ny = electric_field_x.y();
  const int Nz = electric_field_x.z();
  
  const math::Box &box = conf.current().box;
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  // and coordinates into the grid
  math::Matrix H_inv=math::inverse(H);
  
  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  bool assignment_odd = assignment_function_order % 2 == 1;
  const int upper_bound = assignment_function_order / 2;
  const int lower_bound = - (assignment_function_order - 1) / 2;
  DEBUG(10, "Upper bound: " << upper_bound << " lower bound: " << lower_bound);
  
  DEBUG(10, "Calculate the force");
  const unsigned int num_atoms = topo.num_atoms();
  // loop over all charges - again
  for (unsigned int i = 0; i < num_atoms; ++i) {
    const double charge = topo.charge(i);
    if (charge == 0.0)
      continue;

    DEBUG(10, "\tcalculating force for charge " << i + 1 << " (" << charge << ") on grid.");

    // we have a charge so let's calculate were it is on the grid.
    math::Vec grid_r = math::product(H_inv, r(i));
    math::GenericVec<int> nearest_grid_point;

    // getting the nearest grid point is dependend on the order of the
    // charge assignment function
    if (assignment_odd) { // assignment order is odd
      // round to the nearest integer to get the nearest grid point
      nearest_grid_point = math::GenericVec<int>(int(rint(grid_r(0))),
              int(rint(grid_r(1))), int(rint(grid_r(2))));
    } else { // assignment order is even
      // round off to the integer to get the nearest grid centre
      nearest_grid_point = math::GenericVec<int>(int(grid_r(0)),
              int(grid_r(1)),
              int(grid_r(2)));
    }

    math::GenericVec<int> point;
    math::Vec point_to_charge, electric_field_i(0.0);
    double w_x, w_y, w_z;
    for (int dx = lower_bound; dx <= upper_bound; ++dx) {
      point(0) = nearest_grid_point(0) + dx;
      point_to_charge(0) = point(0) - grid_r(0);
      w_x = interaction::Lattice_Sum::charge_assignment_1d(assignment_function_order, point_to_charge(0));
      for (int dy = lower_bound; dy <= upper_bound; ++dy) {
        point(1) = nearest_grid_point(1) + dy;
        point_to_charge(1) = point(1) - grid_r(1);
        w_y = interaction::Lattice_Sum::charge_assignment_1d(assignment_function_order, point_to_charge(1));
        for (int dz = lower_bound; dz <= upper_bound; ++dz) {
          point(2) = nearest_grid_point(2) + dz;
          point_to_charge(2) = point(2) - grid_r(2);
          w_z = interaction::Lattice_Sum::charge_assignment_1d(assignment_function_order, point_to_charge(2));

          double assignment_function = w_x * w_y * w_z;
          // MD05.32 eq. 74
          electric_field_i += assignment_function * math::Vec(
                  electric_field_x(point).real(),
                  electric_field_y(point).real(),
                  electric_field_z(point).real()
                  );
          DEBUG(20, "\t\tgrid point: " << math::v2s(point) << " assignment function: " << assignment_function);
          DEBUG(20, "\t\telectric_field = " << electric_field_x(point).real()
                  << " , " << electric_field_y(point).real()
                  << " , " << electric_field_z(point).real());
        }
      }
    } // loop over grid cells
    DEBUG(12, "electric_field_i = " << math::v2s(electric_field_i));
    // MD05.32 eq. 73
    const math::Vec force = charge * electric_field_i;
    DEBUG(12, "\t\tforce: " << math::v2s(force));
    // safe the force
    storage.force(i) += force;
  } // loop over charges

}

