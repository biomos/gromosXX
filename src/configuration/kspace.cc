/**
 * @file kspace.cc
 * dealing with k space
 */

#include "../stdheader.h"

#include "../algorithm/algorithm.h"
#include "../topology/topology.h"
#include "../simulation/simulation.h"
#include "../configuration/configuration.h"
#include "../math/periodicity.h"
#include "../math/volume.h"

#include "../interaction/nonbonded/interaction/latticesum.h"

#include "../util/error.h"
#include "../util/debug.h"

#include "kspace.h"
#include "configuration.h"

#undef MODULE
#undef SUBMODULE
#define MODULE configuration
#define SUBMODULE configuration

void 
configuration::calculate_k_space(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            int rank, int size) {
  
  DEBUG(8, "\tcalculating k-space\n");
  DEBUG(10, "\trank: " << rank << " size: " << size);
  std::vector<configuration::KSpace_Element> & kspace = conf.lattice_sum().kspace;
  kspace.clear();
 
  math::Matrix l_to_k_matrix = configuration::KSpace_Utils::l_to_k_matrix (
          conf.current().box, conf.boundary_type);
  
  math::Vec k(0.0, 0.0, 0.0);
  
  // maximal k components
  const int max_l[] = {
    sim.param().nonbonded.ewald_max_k_x,
    sim.param().nonbonded.ewald_max_k_y,
    sim.param().nonbonded.ewald_max_k_z
  };
  // k cutoff squared for faster comparison
  const double k_cut2 = sim.param().nonbonded.ewald_kspace_cutoff;
  DEBUG(12, "k space cutoff: " <<  k_cut2);
  const double a = sim.param().nonbonded.ls_charge_shape_width;
  
  const double do_virial = sim.param().pcouple.virial != math::no_virial;
  
  // loop over k-space till cutoff is reached
  math::GenericVec<int> l;
  // here we stride
  for(int lx = -max_l[0] + rank; lx <= max_l[0]; lx += size) {
    l(0) = lx;
    for(int ly = -max_l[1]; ly <= max_l[1]; ++ly) {
      l(1) = ly;
      for(int lz = -max_l[2]; lz <= max_l[2]; ++lz) {
        l(2) = lz;
        // G05 formula 37
        k = math::product(l_to_k_matrix, l);
        const double k2 = math::abs2(k);
        if (k2 <= k_cut2 && k2 != 0.0) {
          // calculate some quantities and save them
          DEBUG(10, "\t\tk: " << math::v2s(k));
          configuration::KSpace_Element k_elem;
          k_elem.k = k;
          k_elem.k2 = k2;
          k_elem.abs_k = sqrt(k2);
          k_elem.ak = k_elem.abs_k * a;
          k_elem.k2i = 1.0 / k2;
          DEBUG(10, "\t\tk2i: " << k_elem.k2i);
          
          if (do_virial) {
            interaction::Lattice_Sum::charge_shape_fourier(
                    sim.param().nonbonded.ls_charge_shape,
                    k_elem.ak, k_elem.fourier_coefficient,
                    &k_elem.fourier_coefficient_derivative
                    );
            k_elem.ak_gammahat_prime = k_elem.ak * k_elem.fourier_coefficient_derivative;
          } else {
            interaction::Lattice_Sum::charge_shape_fourier(
                    sim.param().nonbonded.ls_charge_shape,
                    k_elem.ak, k_elem.fourier_coefficient
                    );
          }
          DEBUG(10, "\t\tfourier coefficient: " << k_elem.fourier_coefficient);
          k_elem.k2i_gammahat = k_elem.k2i * k_elem.fourier_coefficient;
          kspace.push_back(k_elem);
        }
      }
    }
  }
  
  if (kspace.empty()) {
    io::messages.add("kspace is empty. Please increase cutoffs.", 
                      "calculate k space", io::message::error);
  }
}

math::Matrix configuration::KSpace_Utils::l_to_k_matrix(const math::Box & box,
        math::boundary_enum b) {
  // this is MD.05.32 eq. 35 and 36
  const double volume_i_2pi = 2.0 * math::Pi / math::volume(box, b);
  return math::Matrix(
          math::cross(box(1), box(2)) * volume_i_2pi,
          math::cross(box(2), box(0)) * volume_i_2pi,
          math::cross(box(0), box(1)) * volume_i_2pi, true /* column wise */);  
}
