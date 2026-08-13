/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file topology_struct.cu
 * Implementation of the light-weight topology struct for GPU
 */


#include "stdheader.h"
#include "topology/topology.h"

#include "topology_struct.h"

gpu::Topology::Topology(const topology::Topology& topo) {
    num_solute_chargegroups     = topo.num_solute_chargegroups();
    num_solute_molecules        = topo.num_solute_molecules();
    num_chargegroups            = topo.num_chargegroups();
    num_atoms                   = topo.num_atoms();
    num_solute_atoms             = topo.num_solute_atoms();

    // allocate each array separately for automatic alignment
    cudaMalloc(&iac,          sizeof(int)   * num_atoms);
    cudaMalloc(&mass,         sizeof(float) * num_atoms);
    cudaMalloc(&inverse_mass, sizeof(float) * num_atoms);
    cudaMalloc(&charge,       sizeof(float) * num_atoms);
    cudaMalloc(&chargegroup,  sizeof(int)   * (num_chargegroups+1));

    update(topo);

    // Build the CSR exclusion list once -- exclusions are static for a
    // normal run (see header comment); update() deliberately does not
    // touch excl_ptr/excl_list.
    std::vector<int> h_excl_ptr(num_solute_atoms + 1);
    std::vector<int> h_excl_list;
    h_excl_list.reserve(num_solute_atoms * 4); // rough guess, just avoids realloc churn
    for (unsigned i = 0; i < num_solute_atoms; ++i) {
        h_excl_ptr[i] = static_cast<int>(h_excl_list.size());
        const topology::Exclusions & excl = topo.all_exclusion(i);
        for (topology::Exclusions::const_iterator it = excl.begin(), to = excl.end(); it != to; ++it) {
            h_excl_list.push_back(*it);
        }
    }
    h_excl_ptr[num_solute_atoms] = static_cast<int>(h_excl_list.size());
    num_exclusion_entries = static_cast<unsigned>(h_excl_list.size());

    cudaMalloc(&excl_ptr,  sizeof(int) * (num_solute_atoms + 1));
    // cudaMalloc(0) is legal but returns a null-ish pointer on some platforms;
    // guard so excl_list is never dereferenced with a zero-sized allocation.
    cudaMalloc(&excl_list, sizeof(int) * (num_exclusion_entries > 0 ? num_exclusion_entries : 1));

    cudaMemcpy(excl_ptr,  h_excl_ptr.data(),  sizeof(int) * (num_solute_atoms + 1), cudaMemcpyHostToDevice);
    if (num_exclusion_entries > 0) {
        cudaMemcpy(excl_list, h_excl_list.data(), sizeof(int) * num_exclusion_entries, cudaMemcpyHostToDevice);
    }
}



gpu::Topology::~Topology() {
    cudaFree(iac);
    cudaFree(mass);
    cudaFree(inverse_mass);
    cudaFree(charge);
    cudaFree(chargegroup);
    cudaFree(excl_ptr);
    cudaFree(excl_list);
}

void gpu::Topology::update(const topology::Topology& topo) {
    num_solute_chargegroups     = topo.num_solute_chargegroups();
    num_solute_molecules        = topo.num_solute_molecules();
    num_chargegroups            = topo.num_chargegroups();
    num_atoms                   = topo.num_atoms();
    // cast doubles to floats
    std::vector<float> h_mass(num_atoms);
    std::vector<float> h_inverse_mass(num_atoms);
    std::vector<float> h_charge(num_atoms);
    
    for (size_t i = 0; i < num_atoms; ++i) {
        h_mass[i]         = static_cast<float>(topo.mass()[i]);
        h_inverse_mass[i] = static_cast<float>(topo.inverse_mass()[i]);
        h_charge[i]       = static_cast<float>(topo.charge()[i]);
    }

    // copy
    cudaMemcpy(iac,          topo.iac().data(),          sizeof(int)   * num_atoms,        cudaMemcpyHostToDevice);
    cudaMemcpy(mass,         h_mass.data(),              sizeof(float) * num_atoms,        cudaMemcpyHostToDevice);
    cudaMemcpy(inverse_mass, h_inverse_mass.data(),      sizeof(float) * num_atoms,        cudaMemcpyHostToDevice);
    cudaMemcpy(charge,       h_charge.data(),            sizeof(float) * num_atoms,        cudaMemcpyHostToDevice);
    cudaMemcpy(chargegroup,  topo.chargegroups().data(), sizeof(int)   * (num_chargegroups+1), cudaMemcpyHostToDevice);
}