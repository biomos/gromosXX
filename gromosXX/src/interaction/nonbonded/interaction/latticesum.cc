/**
 * @file latticesum.cc
 * implementation of lattice sum methods
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"
#include "../../../interaction/nonbonded/interaction/storage.h"
#include "../../../math/periodicity.h"
#include "../../../math/volume.h"

#include "../../../interaction/nonbonded/interaction/latticesum.h"

#include "../../../util/error.h"
#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE latticesum

template<class MeshType>
void interaction::Lattice_Sum::calculate_charge_density(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain,
        const math::VArray & r) {
  
  MeshType & charge_density = *reinterpret_cast<MeshType*>(conf.lattice_sum().charge_density);
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
  
  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  bool assignment_odd = assignment_function_order % 2 == 1;
  
  DEBUG(10,"\t starting to assign charge density to grid ... ");
  
  // do integer divion to calculate how far you have to go away from the first grid point
  const int upper_bound = assignment_function_order / 2;
  const int lower_bound = - (assignment_function_order - 1) / 2;
  DEBUG(10, "Upper bound: " << upper_bound << " lower bound: " << lower_bound);
  
  // loop over all charges in the domain
  std::vector<int>::const_iterator it = domain.begin(),
          to = domain.end();
  
  for (;it != to; ++it) {
    const unsigned int i = *it;
    const double charge = topo.charge(i);
    
    DEBUG(10, "\tAssigning charge " << i + 1 << " (" << charge << ") to grid.");

    // we have a charge so let's calculate were it is on the grid.
    math::Vec grid_r = math::product(H_inv, r(i));
    DEBUG(10, "\t\tcoordinates on grid: " << math::v2s(grid_r));
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
    double w_x = 0.0, w_y = 0.0, w_z = 0.0;
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
          DEBUG(15, "\t\tgrid point: [" << point(0) << "," << point(1) << "," << point(2) << 
                  "] assignment function: " << assignment_function);
          charge_density(point) += assignment_function * charge;
          
          // output multiplied by 11.78708615 for comparison with promd
          DEBUG(15,"\t\tassignment_function * charge =" << assignment_function * charge* 11.78708615);
          DEBUG(15, "\t\t\tcharge density: " << charge_density(point)* 11.78708615);
        }
      }
    } // loop over grid cells
  } // loop over charges

  
  charge_density.add_neighbors_caches();
  
#ifndef NDEBUG
  for (int z = 0; z < int(Nz); ++z) {
    for (int y = 0; y < int(Ny); ++y) {
      for (int x = charge_density.left_boundary(); x < int(charge_density.right_boundary()); ++x) {
        DEBUG(25, "charge density (" << x << "," << y << "," << z << ") = " << charge_density(x,y,z).real());
        std::cout.precision(8);
        // output multiplied by 11.78708615 for comparison with promd
        DEBUG(25, x+1 << " " << y+1 << " " << z+1 << std::setw(15) << charge_density(x,y,z).real() * 11.78708615);
      }
    }
  }
#endif
  
}
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_charge_density<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain,
        const math::VArray & r);
template void interaction::Lattice_Sum::calculate_charge_density<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain,
        const math::VArray & r);

template<class MeshType>
void interaction::Lattice_Sum::calculate_squared_charge_grid(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain,
        const math::VArray & r) {
  
  MeshType & squared_charge = *reinterpret_cast<MeshType*>(conf.lattice_sum().squared_charge);
  squared_charge.zero();
  
  const unsigned int Nx = squared_charge.x();
  const unsigned int Ny = squared_charge.y();
  const unsigned int Nz = squared_charge.z();
  
  // the H matrix to transform the grid into coordinates
  const math::Box &box = conf.current().box;
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  // and coordinates into the grid
  math::Matrix H_inv=math::inverse(H);
  math::Matrix H_trans = math::transpose(H);
  // the volume of a grid cell (Vg) MD05.32 eq. 61
  const double cell_volume = math::det(H_trans); 
  const double cell_volume_i = 1.0 / cell_volume;
  DEBUG(11, "grid cell volume: " << cell_volume);
  
  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  bool assignment_odd = assignment_function_order % 2 == 1;
  
  DEBUG(10,"\t starting to assign squared charge to grid ... ");
  
  // how far you have to go away from the first grid point
  const int upper_bound = assignment_function_order - 1;
  const int lower_bound = - upper_bound;
  DEBUG(10, "Upper bound: " << upper_bound << " lower bound: " << lower_bound);
  
  // loop over all charges in the domain
  std::vector<int>::const_iterator it = domain.begin(),
          to = domain.end();
  
  for (;it != to; ++it) {
    const unsigned int i = *it;
    const double q2 = topo.charge(i) * topo.charge(i);
    
    DEBUG(10, "\tAssigning q2 " << i + 1 << " (" << q2 << ") to grid.");

    // we have a charge so let's calculate were it is on the grid.
    math::Vec grid_r = math::product(H_inv, r(i));
    DEBUG(10, "\t\tcoordinates on grid: " << math::v2s(grid_r));
    math::GenericVec<int> nearest_grid_point;
    math::Vec point_to_charge;
    
    // getting the nearest grid point is dependend on the order of the
    // charge assignment function
    if (assignment_odd) { // assignment order is odd
      // round to the nearest integer to get the nearest grid point
      nearest_grid_point = math::GenericVec<int>(int(rint(grid_r(0))),
               int(rint(grid_r(1))), int(rint(grid_r(2))));
      point_to_charge = grid_r - nearest_grid_point;
    } else { // assignment order is even
      nearest_grid_point = math::GenericVec<int>(int(grid_r(0)),
                                                 int(grid_r(1)), 
                                                 int(grid_r(2)));
      point_to_charge = grid_r - nearest_grid_point - math::Vec(0.5, 0.5, 0.5);
    }
    DEBUG(15, "\tnearest grid point: " << math::v2s(nearest_grid_point));
    math::GenericVec<int> point;
    DEBUG(15, "\ts: " << math::v2s(point_to_charge));
    for (int dx = lower_bound; dx <= upper_bound; ++dx) {
      point(0) = dx;
      const double w_x = interaction::Lattice_Sum::p3m_selfterm_fp(assignment_function_order,
              dx, point_to_charge(0));
      for (int dy = lower_bound; dy <= upper_bound; ++dy) {
        point(1) = dy;
        const double w_y = interaction::Lattice_Sum::p3m_selfterm_fp(assignment_function_order,
              dy, point_to_charge(1));
        for (int dz = lower_bound; dz <= upper_bound; ++dz) {
          point(2) = dz;
          const double w_z = interaction::Lattice_Sum::p3m_selfterm_fp(assignment_function_order,
              dz, point_to_charge(2));
          
          double assignment_function = w_x * w_y * w_z * cell_volume_i;
          DEBUG(15, "\t\tgrid point: [" << point(0) << "," << point(1) << "," << point(2) << 
                  "] assignment function: " << assignment_function);
          squared_charge(point) += assignment_function * q2;
          DEBUG(15,"\t\tsquared_charge(p): " << assignment_function * q2);
        }
      }
    } // loop over grid cells
  } // loop over charges

  
  squared_charge.add_neighbors_caches();
}
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_squared_charge_grid<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain,
        const math::VArray & r);
template void interaction::Lattice_Sum::calculate_squared_charge_grid<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain,
        const math::VArray & r);

template<class MeshType>
void interaction::Lattice_Sum::calculate_averaged_squared_charge_grid(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain) {
  
  MeshType & squared_charge = *reinterpret_cast<MeshType*>(conf.lattice_sum().squared_charge);
  squared_charge.zero();
  
  const unsigned int Nx = squared_charge.x();
  const unsigned int Ny = squared_charge.y();
  const unsigned int Nz = squared_charge.z();
  
  // the H matrix to transform the grid into coordinates
  const math::Box &box = conf.current().box;
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  // and coordinates into the grid
  math::Matrix H_trans = math::transpose(H);
  // the volume of a grid cell (Vg) MD05.32 eq. 61
  const double cell_volume = math::det(H_trans); 
  const double cell_volume_i = 1.0 / cell_volume;
  DEBUG(11, "grid cell volume: " << cell_volume);
  
  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  
  DEBUG(10,"\t starting to assign squared charge to grid ... ");
  
  // how far you have to go away from the first grid point
  const int upper_bound = assignment_function_order - 1;
  const int lower_bound = - upper_bound;
  DEBUG(10, "Upper bound: " << upper_bound << " lower bound: " << lower_bound);
  
  // loop over all charges in the domain
  std::vector<int>::const_iterator it = domain.begin(),
          to = domain.end();
  
  for (;it != to; ++it) {
    const unsigned int i = *it;
    const double q2 = topo.charge(i) * topo.charge(i);
    
    DEBUG(10, "\tAssigning q2 " << i + 1 << " (" << q2 << ") to grid.");

    math::GenericVec<int> point;
    for (int dx = lower_bound; dx <= upper_bound; ++dx) {
      point(0) = dx;
      const double w_x = interaction::Lattice_Sum::p3m_selfterm_fp_avg(assignment_function_order, dx);
      for (int dy = lower_bound; dy <= upper_bound; ++dy) {
        point(1) = dy;
        const double w_y = interaction::Lattice_Sum::p3m_selfterm_fp_avg(assignment_function_order, dy);
        for (int dz = lower_bound; dz <= upper_bound; ++dz) {
          point(2) = dz;
          const double w_z = interaction::Lattice_Sum::p3m_selfterm_fp_avg(assignment_function_order, dz);
          
          double assignment_function = w_x * w_y * w_z * cell_volume_i;
          DEBUG(15, "\t\tgrid point: [" << point(0) << "," << point(1) << "," << point(2) << 
                  "] assignment function: " << assignment_function);
          squared_charge(point) += assignment_function * q2;
          DEBUG(15,"\t\tsquared_charge(p): " << assignment_function * q2);
        }
      }
    } // loop over grid cells
  } // loop over charges

  squared_charge.add_neighbors_caches();
}
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_averaged_squared_charge_grid<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain);
template void interaction::Lattice_Sum::calculate_averaged_squared_charge_grid<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        const std::vector<int> & domain);

template<class MeshType>
void interaction::Lattice_Sum::calculate_potential_and_energy(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        interaction::Storage & storage) {
  
  
  MeshType & charge_density = *reinterpret_cast<MeshType*>(conf.lattice_sum().charge_density);
  configuration::Influence_Function & influence_function = conf.lattice_sum().influence_function;
  MeshType & potential = *reinterpret_cast<MeshType*>(conf.lattice_sum().potential);
    
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
  const bool do_virial = sim.param().pcouple.virial != math::no_virial;
  math::SymmetricMatrix virial(0.0);
  // loop over grid
  DEBUG(8, "bounds: " << charge_density.left_boundary() << "-" << charge_density.right_boundary());
  const int to = charge_density.right_boundary();
  for (int x = charge_density.left_boundary(); x < to; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        // calculate the potential by MD05.32 eq. 68
        DEBUG(15, "x: " << x + 1 << " y: " << y + 1 << " z: " << z + 1);
        const double ghat = influence_function(x, y, z);
        DEBUG(15, "influence_function(hat): " << ghat);
        DEBUG(15, "charge density(hat): " << charge_density(x, y, z)  / sqrt_grid_volume);
        potential(x, y, z) = math::eps0_i * ghat * charge_density(x, y, z) / sqrt_grid_volume;
        // factors in DEBUG statements are for better comparison with promd
        DEBUG(15, "\tpotential (hat)(" << x << "," << y << "," << z << "): " << potential(x,y,z) / math::eps0_i * 11.78708615);
        
        // calculate the energy by MD05.32 eq. 72     
        const std::complex<double> density = charge_density(x, y, z);
        const double abs2_charge_density = (density.real() * density.real() + 
                density.imag() * density.imag()) / double(grid_volume);
        energy += ghat * abs2_charge_density;
        DEBUG(15, "\t\t abs_charge_density = " << abs2_charge_density / sqrt_grid_volume
                << ", energy = " << energy);
        
        // calculate the virial if required
        if (do_virial) {
          math::SymmetricMatrix tensor(influence_function.gammahat_0(x,y,z));
          tensor.add_to_diagonal(ghat);
          virial += tensor * abs2_charge_density;
        }
      }
    }
  } // loop over grid
  
  // See MD05.32 eq. 64, 65 and kpppm.F about this prefactor issue.
  // the 1/(total volume) prefactor is obtained by the FFT automtically
  storage.energies.ls_kspace_total = 0.5 * energy * cell_volume * math::eps0_i;
  DEBUG(6,"ls_kspace_total = " << storage.energies.ls_kspace_total);
  if (do_virial) {
    virial *= -0.25 * math::eps0_i * cell_volume;
// !!! since multiplication with -0.5 in pressure calc, multiply with -2 here
  //  storage.virial_tensor += virial;
    storage.virial_tensor +=  ( -2.0 * virial );
    DEBUG(6, "k space virial: \n" << math::m2s(virial));
  }

}
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_potential_and_energy<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        interaction::Storage & storage);
template void interaction::Lattice_Sum::calculate_potential_and_energy<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        interaction::Storage & storage);

template<class MeshType>
void interaction::Lattice_Sum::calculate_p3m_selfterm(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim) {
  
  
  MeshType & squared_charge = *reinterpret_cast<MeshType*>(conf.lattice_sum().squared_charge);
  configuration::Influence_Function & influence_function = conf.lattice_sum().influence_function;
    
  const int Nx = squared_charge.x();
  const int Ny = squared_charge.y();
  const int Nz = squared_charge.z();
  const int grid_volume = Nx * Ny * Nz;
  
  double sum = 0.0;
  const bool do_virial = sim.param().pcouple.virial != math::no_virial;
  math::SymmetricMatrix sum_derivative(0.0);
  
  // loop over grid
  const int to = squared_charge.right_boundary();
  for (int x = squared_charge.left_boundary(); x < to; ++x) {
    for (int y = 0; y < Ny; ++y) {
      for (int z = 0; z < Nz; ++z) {
        // calculate the potential by MD05.32 eq. 68
        DEBUG(15, "x: " << x + 1 << " y: " << y + 1 << " z: " << z + 1);
        DEBUG(15, "influence_function(hat): " << influence_function(x, y, z));
        DEBUG(15, "correction: " << influence_function(x, y, z) - influence_function.ghat_0(x,y,z));
        DEBUG(15, "squared charge (hat): " << squared_charge(x, y, z).real() / grid_volume);
        const double ghat = influence_function(x, y, z);
        sum += ghat * squared_charge(x,y,z).real();
        if(do_virial) {
          math::SymmetricMatrix tensor = influence_function.gammahat_0(x,y,z);
          tensor.add_to_diagonal(ghat);
          sum_derivative += tensor * squared_charge(x,y,z).real();
        } // virial
      }
    }
  } // loop over grid
  
  double s2_tilde = topo.sum_squared_charges();
  
  // 1/vol factor is obtained automatically via numerical FFT
  conf.lattice_sum().a2_tilde = 4.0 * math::Pi * sum / grid_volume / s2_tilde;
  DEBUG(8, "a2_tilde = " << conf.lattice_sum().a2_tilde);
  if (do_virial) {
    conf.lattice_sum().a2_tilde_derivative = (4.0 * math::Pi / grid_volume / s2_tilde) *
            sum_derivative;
    DEBUG(8, "a2_tilde deriv = " << math::m2s(conf.lattice_sum().a2_tilde_derivative));
  }
}
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_p3m_selfterm<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);
template void interaction::Lattice_Sum::calculate_p3m_selfterm<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);

template<class MeshType>
void interaction::Lattice_Sum::calculate_electric_field(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim) {
  MeshType & electric_field_x = *reinterpret_cast<MeshType*>(conf.lattice_sum().electric_field.x);
  MeshType & electric_field_y = *reinterpret_cast<MeshType*>(conf.lattice_sum().electric_field.y);
  MeshType & electric_field_z = *reinterpret_cast<MeshType*>(conf.lattice_sum().electric_field.z);
    
  MeshType & potential = *reinterpret_cast<MeshType*>(conf.lattice_sum().potential);

  const int Nx = potential.x();
  const int Ny = potential.y();
  const int Nz = potential.z();
  const int grid_volume = Nx * Ny * Nz;
  
  const int x_to = potential.right_boundary();
  const int x_from = potential.left_boundary();
  
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
  //math::Matrix H_trans = math::transpose(H);
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

    for(int x = x_from; x < x_to; ++x) {
      if (x > int(Nx_2))
        l(0) = x - Nx;
      else
        l(0) = x;
      
      for (int l2 = 1 - Ny_2; l2 <= Ny_2; ++l2) {
        l(1) = l2;
        const int y = (l2+Ny)%Ny;
        for (int l3 = 1 - Nz_2; l3 <= Nz_2; ++l3) {
          l(2) = l3;
          const int z = (l3+Nz)%Nz;
          // MD05.32 eq. 76
          DEBUG(16, "\tl: " << math::v2s(l));                  
          const math::Vec k_l = math::product(tL_inv_2pi, l);
          DEBUG(16, "\tk_l: " << math::v2s(k_l));
          // here we have to divide through sqrt_grid_volume because of numerical
          // FFT algorithms. The electric field (hat) will be scaled but after application
          // of the backward FFT we obtain the correct electrostatic field.
          const std::complex<double> minus_im_pot = minus_im * potential(x,y,z) / sqrt_grid_volume;
          electric_field_x(x,y,z) = k_l(0) * minus_im_pot;
          electric_field_y(x,y,z) = k_l(1) * minus_im_pot;
          electric_field_z(x,y,z) = k_l(2) * minus_im_pot;         
        }
      }
    }
#ifndef NDEBUG
    for (int z = 0; z < int(Nz); ++z) {
      for (int y = 0; y < int(Ny); ++y) {
        for (int x = x_from; x < x_to; ++x) {
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
    
    // make sure all processes have fourrier coefficients of the potential
    potential.get_neighbors();
    
    // now loop over the grid
    // running indices MUST BE integers and NOT unsigned integers! Oh I hate it!!
    for (int x = x_from; x < x_to; ++x) {
      for (int y = 0; y < Ny; ++y) {
        for (int z = 0; z < Nz; ++z) {
          DEBUG(16, "\tpotential (" << x << "," << y << "," << z << "): " << potential(x,y,z));
          // loop over neighboring grid cells
          math::Vec electric_field(0.0);
          for (int q = 1; q <= finite_difference_order; ++q) {
            // get the q-th neighbor grid points
            const math::GenericVec<int> neighbors[6] = {
              math::GenericVec<int>(x - q, y, z), math::GenericVec<int>(x + q, y, z), // x-q will fail if x is unsigned!!!!!
              math::GenericVec<int>(x, y - q, z), math::GenericVec<int>(x, y + q, z),
              math::GenericVec<int>(x, y, z - q), math::GenericVec<int>(x, y, z + q)
            };
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[0]));
            DEBUG(18, "\t\t\t"  << std::setw(15) << potential(neighbors[0]).real());
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[1]));
            DEBUG(18, "\t\t\t"  << std::setw(15) << potential(neighbors[1]).real());
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[2]));
            DEBUG(18, "\t\t\t"  << std::setw(15) << potential(neighbors[2]).real());
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[3]));
            DEBUG(18, "\t\t\t"  << std::setw(15) << potential(neighbors[3]).real());
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[4]));
            DEBUG(18, "\t\t\t"  << std::setw(15) << potential(neighbors[4]).real());
            DEBUG(18, "\t\t1: " << math::v2s(neighbors[5]));
            DEBUG(18, "\t\t\t"  << std::setw(15) << potential(neighbors[5]).real());

            // calculate Delta Potential for all different axis and divide it
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
  
  // get the caches from the neighbors
  electric_field_x.get_neighbors();
  electric_field_y.get_neighbors();
  electric_field_z.get_neighbors();

#ifndef NDEBUG
  for (int z = 0; z < Nz; ++z) {
    for (int y = 0; y < Ny; ++y) {
      for (int x = x_from; x < x_to; ++x) {
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
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_electric_field<configuration::Mesh>(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);
template void interaction::Lattice_Sum::calculate_electric_field<configuration::ParallelMesh>(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);

template<class MeshType>
void interaction::Lattice_Sum::calculate_force(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        Storage & storage,
        const math::VArray & r) {

  MeshType & electric_field_x = *reinterpret_cast<MeshType*>(conf.lattice_sum().electric_field.x);
  MeshType & electric_field_y = *reinterpret_cast<MeshType*>(conf.lattice_sum().electric_field.y);
  MeshType & electric_field_z = *reinterpret_cast<MeshType*>(conf.lattice_sum().electric_field.z);
  
  assert(electric_field_x.x() == electric_field_y.x() && electric_field_x.x() == electric_field_z.x());
  assert(electric_field_x.y() == electric_field_y.y() && electric_field_x.y() == electric_field_z.y());
  assert(electric_field_x.z() == electric_field_y.z() && electric_field_x.z() == electric_field_z.z());
  
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
  
  // loop over all charges in the domain
  std::vector<int>::const_iterator it = storage.domain.begin(),
          to = storage.domain.end();
  
  for (;it != to; ++it) {
    const unsigned int i = *it;
    const double charge = topo.charge(i);

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
    double w_x = 0.0, w_y = 0.0, w_z = 0.0;
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
// create instances for both mesh types
template void interaction::Lattice_Sum::calculate_force<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        Storage & storage,
        const math::VArray & r);
template void interaction::Lattice_Sum::calculate_force<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        Storage & storage,
        const math::VArray & r);

template<class MeshType>
void interaction::Lattice_Sum::decompose_into_domains(const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim,
        std::vector<int> & domain,
        const math::VArray & r,
        const unsigned int size) {
  DEBUG(8, "Starting domain decomposition for P3M.");
  if (sim.mpi) {
    MeshType & charge_density = *((MeshType*) conf.lattice_sum().charge_density);

    const unsigned int Nx = charge_density.x();

    const unsigned int left_bound = charge_density.left_boundary();
    const unsigned int right_bound = charge_density.right_boundary();
    DEBUG(12, "boundaries: " << left_bound << ", " << right_bound);

    const math::Box &box = conf.current().box;
    math::Matrix L(box(0), box(1), box(2), true);

    // the following is the same as the Hinv multiplication in the
    // assignment part. But we only need the x component ->
    // faster implementation without y and z
    //
    // grid_x = term_x * x + term_y * y + term_z * z

    const math::Vec & b = box(1);
    const math::Vec & c = box(2);

    const double scale = -double(Nx) / math::det(L);
    const double term_x = (b(2) * c(1) - b(1) * c(2)) * scale;
    const double term_y = (b(0) * c(2) - b(2) * c(0)) * scale;
    const double term_z = (b(1) * c(0) - b(0) * c(1)) * scale;

    const unsigned int num_atoms = topo.num_atoms();

    bool assignment_odd = sim.param().nonbonded.p3m_charge_assignment % 2 == 1;

    DEBUG(10, "\t starting to decompose into domain");
    domain.clear();
    domain.reserve(num_atoms / size);

    // loop over all charges
    for (unsigned int i = 0; i < num_atoms; ++i) {
      const double charge = topo.charge(i);
      if (charge == 0.0)
        continue;
      const math::Vec & r_i = r(i);
      // we have a charge so let's calculate were it is on the grid.
      const double grid_x = term_x * r_i(0) + term_y * r_i(1) + term_z * r_i(2);
      DEBUG(13, "\t\treal coordinates: " << math::v2s(r_i));
      DEBUG(13, "\t\tcoordinates on grid (x comp): " << grid_x);
      unsigned int nearest_grid_point = 0;

      // getting the nearest grid point is dependend on the order of the
      // charge assignment function
      if (assignment_odd) { // assignment order is odd
        // round to the nearest integer to get the nearest grid point
        nearest_grid_point = int(rint(grid_x));
      } else { // assignment order is even
        // round to the nearest integer to get the nearest grid centre
        // we have to subtract math::epsilon (a small number) from the coorinates
        nearest_grid_point = int(grid_x);
      }

      nearest_grid_point %= Nx; // last -> 0

      // add it to the domain
      if (nearest_grid_point >= left_bound && nearest_grid_point < right_bound)
        domain.push_back(i);
    }
  } else { // no mpi
    // we simply add all atoms to the domain.
    const unsigned int num_atoms = topo.num_atoms();
    domain.resize(num_atoms);
    
    for(unsigned int i = 0; i < num_atoms; ++i)
      domain[i] = i;
  }
}
// create instances for both mesh types
template void interaction::Lattice_Sum::decompose_into_domains<configuration::Mesh>(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            std::vector<int> & domain,
            const math::VArray & r,
            const unsigned int size);
template void interaction::Lattice_Sum::decompose_into_domains<configuration::ParallelMesh>(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            std::vector<int> & domain,
            const math::VArray & r,
            const unsigned int size);
