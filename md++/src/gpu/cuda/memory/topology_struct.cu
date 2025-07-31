


#include "stdheader.h"
#include "topology/topology.h"

#include "topology_struct.h"

gpu::topology_struct::topology_struct(topology::Topology& topo) {
    num_solute_chargegroups     = topo.num_solute_chargegroups();
    num_solute_molecules        = topo.num_solute_molecules();
    num_chargegroups            = topo.num_chargegroups();
    num_atoms                   = topo.num_atoms();
    size_t total_bytes = 
        sizeof(int) * num_atoms +       // m_iac
        sizeof(float) * num_atoms +     // m_mass - ah, for these we need casting first
        sizeof(float) * num_atoms +     // m_inverse_mass
        sizeof(float) * num_atoms +     // m_charge
        sizeof(int) * num_chargegroups;        // m_chargegroup
    
    // allocate - single strip for better layout
    cudaMalloc(&memory_block, total_bytes);
    char* base                  = reinterpret_cast<char*>(memory_block);
    iac                         = reinterpret_cast<int*>(base);
    mass                        = reinterpret_cast<float*>(base += sizeof(int) * num_atoms);
    inverse_mass                = reinterpret_cast<float*>(base += sizeof(float) * num_atoms);
    charge                      = reinterpret_cast<float*>(base += sizeof(float) * num_atoms);
    chargegroup                 = reinterpret_cast<int*>(base += sizeof(float) * num_atoms);

    update(topo);

}

gpu::topology_struct::~topology_struct() {
    cudaFree(memory_block);
}

void gpu::topology_struct::update(topology::Topology& topo) {
    num_solute_chargegroups     = topo.num_solute_chargegroups();
    num_solute_molecules        = topo.num_solute_molecules();
    num_chargegroups            = topo.num_chargegroups();
    num_atoms                   = topo.num_atoms();
    // cast doubles to floats
    std::vector<float> h_mass, h_inverse_mass, h_charge;
    h_mass.reserve(num_atoms);
    h_inverse_mass.reserve(num_atoms);
    h_charge.reserve(num_atoms);

    for (double d : topo.mass()) {
        h_mass.push_back(static_cast<float>(d));
    }

    for (double d : topo.inverse_mass()) {
        h_inverse_mass.push_back(static_cast<float>(d));
    }

    for (double d : topo.charge()) {
        h_charge.push_back(static_cast<float>(d));
    }

    // copy
    cudaMemcpy(iac, topo.iac().data(), sizeof(int) * num_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(mass, h_mass.data(), sizeof(float) * num_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(inverse_mass, h_inverse_mass.data(), sizeof(float) * num_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(charge, h_charge.data(), sizeof(float) * num_atoms, cudaMemcpyHostToDevice);
    cudaMemcpy(chargegroup, topo.chargegroups().data(), sizeof(int) * num_chargegroups, cudaMemcpyHostToDevice);
}