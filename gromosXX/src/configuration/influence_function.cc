/**
 * @file influence_function.cc
 * influence function logic
 */
#ifdef XXMPI
#include <mpi.h>
#endif
#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../configuration/mesh.h"
#include "../simulation/multibath.h"
#include "../simulation/parameter.h"
#include "../configuration/influence_function.h"

#include "../interaction/nonbonded/interaction/latticesum.h"

#include "../math/volume.h"

#include "../util/error.h"
#include "../util/debug.h"

#include "influence_function.h"
#ifdef OMP
#include <omp.h>
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

configuration::Influence_Function::Influence_Function() : do_virial(false),
do_scale(false) {
}

void configuration::Influence_Function::setBox(const math::Box & box) {
  DEBUG(12, "setting to box" << math::m2s(math::Matrix(box)));
  DEBUG(12, "initial transpose box: " << math::m2s(tL_0));
  tL = math::transpose(math::Matrix(box));
  tbox_dev = tL - tL_0;
  DEBUG(12, "box deviation: " << math::m2s(tbox_dev));
}

void configuration::Influence_Function::init(const simulation::Parameter & param) {
  unsigned int x = param.nonbonded.p3m_grid_points_x;
  unsigned int y = param.nonbonded.p3m_grid_points_y;
  unsigned int z = param.nonbonded.p3m_grid_points_z;

  if (param.multicell.multicell) {
    x *= param.multicell.x;
    y *= param.multicell.y;
    z *= param.multicell.z;
  }

  ghat.resize(x, y, z);

  // we need the derivative to calculate the pressure!
  do_virial = param.pcouple.virial != math::no_virial;
  if (do_virial) {
    // reserve space for the derivative
    gammahat.resize(x, y, z);
  }
  do_scale = param.pcouple.scale != math::pcouple_off;
}

template<class MeshType>
void configuration::Influence_Function::calculate(const topology::Topology & topo,
configuration::Configuration & conf,
const simulation::Simulation & sim) {
  DEBUG(8, "\tUpdating influence function");
  const unsigned int Nx = ghat.x();
  const unsigned int Ny = ghat.y();
  const unsigned int Nz = ghat.z();

  const double st2 = topo.sum_squared_charges();

  DEBUG(15, "\tgrid dimensions " << Nx << "x" << Ny << "x" << Nz);

  ghat.zero();
  if (do_virial)
    gammahat.zero();

  const int shape = sim.param().nonbonded.ls_charge_shape;
  const double charge_width = sim.param().nonbonded.ls_charge_shape_width;

  const math::Box &box = conf.current().box;
  // the H matrix to transform the grid into coordinates
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  math::Matrix H_inv = math::inverse(H);
  math::Matrix H_trans = math::transpose(H);
  math::Matrix H_inv_trans = math::transpose(H_inv);

  const int mesh_alias = sim.param().nonbonded.p3m_mesh_alias;

  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  const int finite_difference_order = sim.param().nonbonded.p3m_finite_differences_operator;
  std::vector<double> finite_difference_constants = interaction::Lattice_Sum::finite_difference_constants(finite_difference_order);

  const int Nx_2 = Nx / 2;
  const int Ny_2 = Ny / 2;
  const int Nz_2 = Nz / 2;

  const double volume = math::volume(conf.current().box, conf.boundary_type);
  tL_0 = math::transpose(math::Matrix(conf.current().box));
  tLi_0 = math::inverse(tL_0);
  math::Matrix tL_inv_2pi = tLi_0 * (2.0 * math::Pi);

  MeshType & potential = *reinterpret_cast<MeshType*> (conf.lattice_sum().potential);

  const int x_to = potential.right_boundary();
  const int x_from = potential.left_boundary();

  // loop over reciprocal-space grid G
  DEBUG(10, "\tstarting loop over reciprocal-space grid");

  // quality control: quantiy q to calculate the RMS force error
  double my_q = 0.0;
  for (int x = x_from; x < x_to; ++x) {
    double q_thread = 0.0;
    math::GenericVec<int> l;
    if (x > int(Nx_2))
      l(0) = x - Nx;
    else
      l(0) = x;

    for (int l2 = 1 - Ny_2; l2 <= Ny_2; ++l2) {
      l(1) = l2;
      for (int l3 = 1 - Nz_2; l3 <= Nz_2; ++l3) {
        l(2) = l3;

        if (l(0) == 0 && l(1) == 0 && l(2) == 0) {
          ghat(l) = 0.0;
          if (do_virial)
            gammahat(l) = 0.0;
          continue;
        }

        DEBUG(20, "\t l = " << math::v2s(l));
        // this is MD.05.32 eq 63
        math::Vec k_l = math::product(tL_inv_2pi, l);
        DEBUG(20, "\tk_l: " << math::v2s(k_l));

        // loop over m mesh alias -> eq. 90.
        math::GenericVec<int> m;

        // storage for sums over m in eq. 90
        math::Vec sum_ghat_numerator(0.0);
        double sum_ghat_denominator = 0.0;
        double sum_k2ifourier2 = 0.0;

        // storage for tensor sums
        math::SymmetricMatrix sum_gammahat_numerator(0.0);

        // get the finite difference operator MD05.32 eq. 87
        math::Vec D_hat_g(0.0);
        double abs2_D_hat_g = 0.0;
        if (finite_difference_order == 0) {
          // ik differentiation see MD.05.32 eq. 88
          D_hat_g = k_l;
          abs2_D_hat_g = math::abs2(D_hat_g);
        } else {
          // finite differences as in MD05.32 eq. 87
          math::Vec tH_kl = math::product(H_trans, k_l);
          DEBUG(15, "tH_kl" << math::v2s(tH_kl));
          // loop over finit difference order q
          for (int j = 1; j <= finite_difference_order; ++j) {
            math::Vec sin_term(sin(j * tH_kl(0)), sin(j * tH_kl(1)), sin(j * tH_kl(2)));
            DEBUG(17, "sin_term: " << math::v2s(sin_term));
            D_hat_g += sin_term * (finite_difference_constants[j - 1] / j);
          }
          D_hat_g = math::product(H_inv_trans, D_hat_g);
          abs2_D_hat_g = math::abs2(D_hat_g);
        }
        DEBUG(20, "\t D_hat_g = " << math::v2s(D_hat_g));

        // loop over mesh aliases
        for (int m1 = -mesh_alias; m1 <= mesh_alias; ++m1) {
          m(0) = m1;
          for (int m2 = -mesh_alias; m2 <= mesh_alias; ++m2) {
            m(1) = m2;
            for (int m3 = -mesh_alias; m3 <= mesh_alias; ++m3) {
              m(2) = m3;
              DEBUG(25, "\t\t m = " << math::v2s(m));
              // MD.05.32 eq. 91
              math::Vec k_lm = k_l + (2.0 * math::Pi * math::product(H_inv_trans, m));
              const double k_2 = math::abs2(k_lm);
              const double k_2i = 1.0 / k_2;
              const double k = sqrt(k_2);
              const double ak = charge_width * k;

              double fourier_coeff = 0.0, fourier_coeff_deriv = 0.0;
              if (do_virial) // get the derivative as well
                interaction::Lattice_Sum::charge_shape_fourier(shape, ak, fourier_coeff, &fourier_coeff_deriv);
              else
                interaction::Lattice_Sum::charge_shape_fourier(shape, ak, fourier_coeff);

              DEBUG(25, "\t\t k_lm " << math::v2s(k_lm) << ", fourier_coeff = " << fourier_coeff);
              math::Vec tH_k = math::product(H_trans, k_lm);
              const double P_hat =
                      interaction::Lattice_Sum::fourier_charge_assignment_1d(assignment_function_order, tH_k(0)) *
                      interaction::Lattice_Sum::fourier_charge_assignment_1d(assignment_function_order, tH_k(1)) *
                      interaction::Lattice_Sum::fourier_charge_assignment_1d(assignment_function_order, tH_k(2));
              const double P_hat_2 = P_hat * P_hat;
              DEBUG(25, "\t\t P_hat = " << P_hat);
              // add the terms
              const double k2i_fourier = k_2i * fourier_coeff;
              sum_ghat_numerator += k_lm * (k2i_fourier * P_hat_2);
              sum_ghat_denominator += P_hat_2;
              sum_k2ifourier2 += fourier_coeff * k2i_fourier;

              // let's calculate the terms needed for the derivative
              // see GROMOS05 eq. 92
              if (do_virial) {
                const double k_dot_D = math::dot(k_lm, D_hat_g);
                const math::SymmetricMatrix k_x_D = math::symmetric_tensor_product(k_lm, D_hat_g);
                const math::SymmetricMatrix D_x_k = math::symmetric_tensor_product(D_hat_g, k_lm);
                const math::SymmetricMatrix k_x_k = math::symmetric_tensor_product(k_lm, k_lm);
                const math::SymmetricMatrix D_x_D = math::symmetric_tensor_product(D_hat_g, D_hat_g);

                const double D_2i = 1.0 / abs2_D_hat_g;

                sum_gammahat_numerator += (P_hat_2 * k_2i) * ((k_x_D + D_x_k -
                        (2.0 * k_dot_D) * (k_2i * k_x_k + D_2i * D_x_D)) * fourier_coeff +
                        (ak * fourier_coeff_deriv * k_2i * k_dot_D) * k_x_k);
              } // do the virial
            }
          }
        } // loop over mesh alias
        DEBUG(15, "sum ghat denom: " << sum_ghat_denominator);
        DEBUG(15, "sum ghat num: " << math::v2s(sum_ghat_numerator));
        DEBUG(15, "sum gammahat num:\n\t" << math::m2s(sum_gammahat_numerator));
        if (assignment_function_order == 1) {
          sum_ghat_denominator = 1.0;
        }
        // MD05.32 eq. 90

        // if this is zero we have to set the influence function to zero. 
        // Due to numerical problems!!!
        if (abs2_D_hat_g < math::epsilon) {
          ghat(l) = 0.0;
          if (do_virial)
            gammahat(l) = math::SymmetricMatrix(0.0);

          q_thread += sum_k2ifourier2;

          DEBUG(15, "\t influence function" << math::v2s(l) << " ="
                  "0.0 set to zero (numerics).");
          continue;
        }

        const double numerator = math::dot(D_hat_g, sum_ghat_numerator);
        const double denominator = abs2_D_hat_g * sum_ghat_denominator * sum_ghat_denominator;
        DEBUG(15, "numerator: " << numerator << " denominator: " << denominator
                << " sum_k2ifourier2: " << sum_k2ifourier2);

        const double ghat_l = numerator / denominator;
        ghat(l) = ghat_l;

        q_thread += sum_k2ifourier2 - numerator * ghat_l;

        DEBUG(13, "\t influence function (0)" << math::v2s(l) << " ="
                << ghat(l));

        if (do_virial) {
          gammahat(l) = sum_gammahat_numerator / denominator;
          DEBUG(13, "\t influence function derivative:\n\t" << math::m2s(gammahat(l)));
        }

      }
    }
#ifdef OMP
#pragma omp critical(crit1)
#endif
    my_q += q_thread;
  } // loop over reciprocal space grid
  my_q *= 16 * math::Pi * math::Pi / volume;

#ifdef XXMPI
  if (sim.mpi) {
    MPI_Allreduce(&my_q, &q, 1, MPI::DOUBLE, MPI::SUM, sim.mpiControl().comm);
  } else {
#endif
    q = my_q;
#ifdef XXMPI
  }
#endif
  force_error = math::four_pi_eps_i * st2 * sqrt(q / (volume * topo.num_atoms()));
  DEBUG(10, "q = " << q);

  // symmetrize the influence function
  // this can only be done when a rectangular box is used
  /*
  for (int z = 0; z <= Nz_2; ++z) {
    int zz = (Nz - z) % Nz;
    for (int y = 0; y <= Ny_2; ++y) {
      int yy = (Ny - y) % Ny;
      for (int x = 0; x <= Nx_2; ++x) {
        int xx = (Nx - x) % Nx;
        double & tmp = influence_function(x, y, z);
        DEBUG(1, "\tinfluence_function(" << x << "," << y << "," << z << ") = " << tmp);
          
        influence_function(xx, y, z) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << xx << "," << y << "," << z);
        influence_function(x, yy, z) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << x << "," << yy << "," << z);
        influence_function(xx, yy, z) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << xx << "," << yy << "," << z);
        influence_function(x, y, zz) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << x << "," << y << "," << zz);
        influence_function(xx, y, zz) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << xx << "," << y << "," << zz);
        influence_function(x, yy, zz) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << x << "," << yy << "," << zz);
        influence_function(xx, yy, zz) = tmp;
        DEBUG(1,"\t\t - assigning to grid point " << xx << "," << yy << "," << zz);
      }
    }
  }
   */
}
template void configuration::Influence_Function::calculate<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);
template void configuration::Influence_Function::calculate<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);

template<class MeshType>
void configuration::Influence_Function::evaluate_quality(
const topology::Topology & topo,
configuration::Configuration & conf,
const simulation::Simulation & sim) {
  DEBUG(8, "\tEvaluating the quality of the influence function");
  const unsigned int Nx = ghat.x();
  const unsigned int Ny = ghat.y();
  const unsigned int Nz = ghat.z();

  const double st2 = topo.sum_squared_charges();

  DEBUG(15, "\tgrid dimensions " << Nx << "x" << Ny << "x" << Nz);
  const int shape = sim.param().nonbonded.ls_charge_shape;
  const double charge_width = sim.param().nonbonded.ls_charge_shape_width;

  const math::Box &box = conf.current().box;
  // the H matrix to transform the grid into coordinates
  math::Matrix H(box(0) / Nx, box(1) / Ny, box(2) / Nz);
  math::Matrix H_inv = math::inverse(H);
  math::Matrix H_trans = math::transpose(H);
  math::Matrix H_inv_trans = math::transpose(H_inv);

  const int mesh_alias = sim.param().nonbonded.p3m_mesh_alias;

  const int assignment_function_order = sim.param().nonbonded.p3m_charge_assignment;
  const int finite_difference_order = sim.param().nonbonded.p3m_finite_differences_operator;
  std::vector<double> finite_difference_constants = interaction::Lattice_Sum::finite_difference_constants(finite_difference_order);

  const int Nx_2 = Nx / 2;
  const int Ny_2 = Ny / 2;
  const int Nz_2 = Nz / 2;

  // this is MD.05.32 eq. 35 and 36
  const double volume = math::volume(box, sim.param().boundary.boundary);
  const double volume_i_2pi = 2.0 * math::Pi / volume;
  math::Matrix tL_inv_2pi(
          math::cross(box(1), box(2)) * volume_i_2pi,
          math::cross(box(2), box(0)) * volume_i_2pi,
          math::cross(box(0), box(1)) * volume_i_2pi, true /* column wise */);

  MeshType & potential = *reinterpret_cast<MeshType*> (conf.lattice_sum().potential);

  const int x_to = potential.right_boundary();
  const int x_from = potential.left_boundary();

  // loop over reciprocal-space grid G
  DEBUG(10, "\tstarting loop over reciprocal-space grid");
  math::GenericVec<int> l;

  // quality control: quantiy q to calculate the RMS force error
  double my_q = 0.0;
  for (int x = x_from; x < x_to; ++x) {
    double q_thread = 0.0;
    if (x > int(Nx_2))
      l(0) = x - Nx;
    else
      l(0) = x;

    for (int l2 = 1 - Ny_2; l2 <= Ny_2; ++l2) {
      l(1) = l2;
      for (int l3 = 1 - Nz_2; l3 <= Nz_2; ++l3) {
        l(2) = l3;

        if (l(0) == 0 && l2 == 0 && l3 == 0) {
          continue;
        }

        DEBUG(20, "\t l = " << math::v2s(l));
        // this is MD.05.32 eq 63
        math::Vec k_l = math::product(tL_inv_2pi, l);
        DEBUG(20, "\tk_l: " << math::v2s(k_l));

        // loop over m mesh alias -> eq. 90.
        math::GenericVec<int> m;

        // storage for sums over m in eq. 90
        math::Vec sum_ghat_numerator(0.0);
        double sum_ghat_denominator = 0.0;
        double sum_k2ifourier2 = 0.0;

        for (int m1 = -mesh_alias; m1 <= mesh_alias; ++m1) {
          m(0) = m1;
          for (int m2 = -mesh_alias; m2 <= mesh_alias; ++m2) {
            m(1) = m2;
            for (int m3 = -mesh_alias; m3 <= mesh_alias; ++m3) {
              m(2) = m3;
              DEBUG(25, "\t\t m = " << math::v2s(m));
              // MD.05.32 eq. 91
              math::Vec k_lm = k_l + (2.0 * math::Pi * math::product(H_inv_trans, m));
              const double k_2 = math::abs2(k_lm);
              const double k_2i = 1.0 / k_2;
              const double k = sqrt(k_2);
              double fourier_coeff = 0.0;
              interaction::Lattice_Sum::charge_shape_fourier(shape, charge_width * k, fourier_coeff);
              DEBUG(25, "\t\t k_lm " << math::v2s(k_lm) << ", fourier_coeff = " << fourier_coeff);
              math::Vec tH_k = math::product(H_trans, k_lm);
              const double P_hat =
                      interaction::Lattice_Sum::fourier_charge_assignment_1d(assignment_function_order, tH_k(0)) *
                      interaction::Lattice_Sum::fourier_charge_assignment_1d(assignment_function_order, tH_k(1)) *
                      interaction::Lattice_Sum::fourier_charge_assignment_1d(assignment_function_order, tH_k(2));
              const double P_hat_2 = P_hat * P_hat;
              DEBUG(25, "\t\t P_hat = " << P_hat);
              // add the terms
              const double k2i_fourier = k_2i * fourier_coeff;
              sum_ghat_numerator += k_lm * (k2i_fourier * P_hat_2);
              sum_ghat_denominator += P_hat_2;
              sum_k2ifourier2 += fourier_coeff * k2i_fourier;
            }
          }
        } // loop over mesh alias

        if (assignment_function_order == 1) {
          sum_ghat_denominator = 1.0;
        }

        DEBUG(15, "sum ghat denom: " << sum_ghat_denominator);
        DEBUG(15, "sum ghat num: " << math::v2s(sum_ghat_numerator));
        // get the finite difference operator MD05.32 eq. 87
        math::Vec D_hat_g(0.0);
        double abs2_D_hat_g = 0.0;
        if (finite_difference_order == 0) {
          // ik differentiation see MD.05.32 eq. 88
          D_hat_g = k_l;
          abs2_D_hat_g = math::abs2(D_hat_g);
        } else {
          // finite differences as in MD05.32 eq. 87
          math::Vec tH_kl = math::product(H_trans, k_l);
          DEBUG(15, "tH_kl" << math::v2s(tH_kl));
          // loop over finit difference order q
          for (int j = 1; j <= finite_difference_order; ++j) {
            math::Vec sin_term(sin(j * tH_kl(0)), sin(j * tH_kl(1)), sin(j * tH_kl(2)));
            DEBUG(17, "sin_term: " << math::v2s(sin_term));
            D_hat_g += sin_term * (finite_difference_constants[j - 1] / j);
          }
          D_hat_g = math::product(H_inv_trans, D_hat_g);
          abs2_D_hat_g = math::abs2(D_hat_g);

          // if this is zero we have to set the influence function terms to zero. 
          // Due to numerical problems!!!
          if (abs2_D_hat_g < math::epsilon) {
            q_thread += sum_k2ifourier2;
            continue;
          }
        }
        DEBUG(13, "\t D_hat_g = " << math::v2s(D_hat_g));
        // MD99.32 eq. 203
        double ghat = (*this)(l);
        q_thread += ghat * ghat * abs2_D_hat_g * sum_ghat_denominator * sum_ghat_denominator;
        q_thread += sum_k2ifourier2 - 2 * ghat * math::dot(D_hat_g, sum_ghat_numerator);
        DEBUG(20, "running q = " << q);
      }
    }
#ifdef OMP
#pragma omp critical(crit2)
#endif
    my_q += q_thread;
  } // loop over reciprocal space grid
  my_q *= 16 * math::Pi * math::Pi / volume;
#ifdef XXMPI
  if (sim.mpi) {
    MPI_Allreduce(&my_q, &q, 1, MPI::DOUBLE, MPI::SUM, sim.mpiControl().comm);
  } else {
#endif
    q = my_q;
#ifdef XXMPI
  }
#endif

  force_error = math::four_pi_eps_i * st2 * sqrt(q / (volume * topo.num_atoms()));

  DEBUG(8, "q = " << q);
}
template void configuration::Influence_Function::evaluate_quality<configuration::Mesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);
template void configuration::Influence_Function::evaluate_quality<configuration::ParallelMesh>(
        const topology::Topology & topo,
        configuration::Configuration & conf,
        const simulation::Simulation & sim);

