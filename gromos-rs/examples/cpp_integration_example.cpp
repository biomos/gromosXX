/**
 * Example: Integrating Rust kernels into GROMOS C++ code
 *
 * This file demonstrates how to call the Rust nonbonded inner loop
 * from existing GROMOS C++ code.
 */

#include "../gromos_rs.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

// Mock GROMOS data structures (simplified)
struct Configuration {
    std::vector<float> positions;  // Flat: [x0,y0,z0,x1,y1,z1,...]
    std::vector<float> charges;
    std::vector<float> forces;     // Output
};

struct Topology {
    std::vector<uint32_t> iac;  // Integer atom codes (atom types)
    std::vector<std::vector<double>> lj_params;  // [n_types][n_types][2] (c6, c12)
};

struct Pairlist {
    std::vector<uint32_t> pairs;  // Flat: [i0,j0,i1,j1,...]
};

/**
 * Example 1: Simple 2-atom system
 */
void example_simple() {
    std::cout << "Example 1: Simple 2-atom LJ interaction\n";
    std::cout << "========================================\n";

    // Setup: Two atoms 1 nm apart
    Configuration conf;
    conf.positions = {0.0f, 0.0f, 0.0f,   // Atom 0 at origin
                      1.0f, 0.0f, 0.0f};  // Atom 1 at (1,0,0)
    conf.charges = {0.0f, 0.0f};           // No charges for pure LJ
    conf.forces.resize(6, 0.0f);

    Topology topo;
    topo.iac = {0, 0};  // Both atoms are type 0

    // LJ parameters for type 0-0 interaction
    // Typical values for hydrocarbons: C6 ~ 0.001, C12 ~ 0.0001
    std::vector<double> lj_flat = {0.001, 0.0001};  // c6, c12

    Pairlist plist;
    plist.pairs = {0, 1};  // Single pair: atoms 0 and 1

    // CRF parameters (set to zero for pure LJ)
    double crf_cut = 1.4;
    double crf_2cut3i = 0.0;
    double crf_cut3i = 0.0;

    // Box (not used for vacuum boundaries)
    std::vector<double> box = {10.0, 10.0, 10.0};

    // Output arrays
    std::vector<double> energies(2, 0.0);  // [e_lj, e_crf]
    std::vector<double> virial(9, 0.0);    // 3x3 tensor

    // Call Rust kernel
    rust_lj_crf_innerloop(
        conf.positions.data(),
        conf.charges.data(),
        topo.iac.data(),
        2,  // n_atoms
        plist.pairs.data(),
        1,  // n_pairs
        lj_flat.data(),
        1,  // n_types
        box.data(),
        0,  // boundary_type = vacuum
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
        conf.forces.data(),
        energies.data(),
        virial.data()
    );

    // Print results
    std::cout << "Forces on atom 0: (" << conf.forces[0] << ", "
              << conf.forces[1] << ", " << conf.forces[2] << ")\n";
    std::cout << "Forces on atom 1: (" << conf.forces[3] << ", "
              << conf.forces[4] << ", " << conf.forces[5] << ")\n";
    std::cout << "LJ Energy: " << energies[0] << " kJ/mol\n";
    std::cout << "CRF Energy: " << energies[1] << " kJ/mol\n";

    // Verify Newton's third law: F_i = -F_j
    double fx_diff = std::abs(conf.forces[0] + conf.forces[3]);
    assert(fx_diff < 1e-6 && "Newton's third law violated!");

    std::cout << "✓ Newton's third law verified\n\n";
}

/**
 * Example 2: Water box with periodic boundaries
 */
void example_water_box() {
    std::cout << "Example 2: Water molecules with PBC\n";
    std::cout << "====================================\n";

    // Setup: 3 water molecules (9 atoms) in a periodic box
    Configuration conf;

    // Positions (simplified SPC water, box size 3 nm)
    conf.positions = {
        // Water 1
        1.0f, 1.0f, 1.0f,   // O
        1.1f, 1.0f, 1.0f,   // H1
        0.9f, 1.0f, 1.0f,   // H2
        // Water 2
        2.0f, 2.0f, 2.0f,   // O
        2.1f, 2.0f, 2.0f,   // H1
        1.9f, 2.0f, 2.0f,   // H2
        // Water 3 (near boundary)
        0.1f, 0.1f, 0.1f,   // O
        0.2f, 0.1f, 0.1f,   // H1
        0.0f, 0.1f, 0.1f,   // H2
    };

    // Charges (SPC water: O=-0.82e, H=+0.41e)
    conf.charges = {
        -0.82f, 0.41f, 0.41f,  // Water 1
        -0.82f, 0.41f, 0.41f,  // Water 2
        -0.82f, 0.41f, 0.41f,  // Water 3
    };

    conf.forces.resize(27, 0.0f);

    Topology topo;
    topo.iac = {0, 1, 1,  // O, H, H (water 1)
                0, 1, 1,  // O, H, H (water 2)
                0, 1, 1}; // O, H, H (water 3)

    // LJ parameters: type 0 (O), type 1 (H)
    // Matrix: [0-0, 0-1]
    //         [1-0, 1-1]
    std::vector<double> lj_flat = {
        // O-O
        0.00262, 0.0026,
        // O-H (no interaction)
        0.0, 0.0,
        // H-O (no interaction)
        0.0, 0.0,
        // H-H (no interaction)
        0.0, 0.0
    };

    // Pairlist: all O-O interactions
    Pairlist plist;
    plist.pairs = {
        0, 3,  // Water 1-O to Water 2-O
        0, 6,  // Water 1-O to Water 3-O
        3, 6,  // Water 2-O to Water 3-O
    };

    // Rectangular periodic box: 3x3x3 nm
    std::vector<double> box = {3.0, 3.0, 3.0};

    // CRF parameters for water
    double crf_cut = 1.4;
    double crf_2cut3i = 2.0 / (crf_cut * crf_cut * crf_cut);
    double crf_cut3i = 0.5 * crf_2cut3i;

    // Output arrays
    std::vector<double> energies(2, 0.0);
    std::vector<double> virial(9, 0.0);

    // Call Rust kernel with rectangular PBC
    rust_lj_crf_innerloop(
        conf.positions.data(),
        conf.charges.data(),
        topo.iac.data(),
        9,  // n_atoms
        plist.pairs.data(),
        3,  // n_pairs
        lj_flat.data(),
        2,  // n_types (O and H)
        box.data(),
        1,  // boundary_type = rectangular
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
        conf.forces.data(),
        energies.data(),
        virial.data()
    );

    std::cout << "Total LJ Energy: " << energies[0] << " kJ/mol\n";
    std::cout << "Total CRF Energy: " << energies[1] << " kJ/mol\n";
    std::cout << "Pressure (from virial): " << virial[0] + virial[4] + virial[8] << "\n";

    std::cout << "✓ Water box calculation complete\n\n";
}

/**
 * Example 3: Replacing GROMOS C++ innerloop
 *
 * This shows how to modify the actual GROMOS code to call Rust
 */
void example_gromos_integration() {
    std::cout << "Example 3: GROMOS Integration Pattern\n";
    std::cout << "======================================\n";

    std::cout << R"(
In md++/src/interaction/nonbonded/interaction/nonbonded_innerloop.cc,
replace the innerloop implementation:

// Original C++ code:
template<typename t_nonbonded_spec>
void Nonbonded_Innerloop<t_nonbonded_spec>::lj_crf_innerloop(...) {
    for (pairlist iteration) {
        // C++ implementation
    }
}

// New Rust-accelerated version:
#ifdef USE_RUST_KERNELS
    // Prepare data in C-compatible format
    std::vector<float> pos_flat;
    for (auto& p : conf.current().pos) {
        pos_flat.push_back(p(0));
        pos_flat.push_back(p(1));
        pos_flat.push_back(p(2));
    }

    std::vector<uint32_t> pairlist_flat;
    for (auto& pair : pairlist) {
        pairlist_flat.push_back(pair.i);
        pairlist_flat.push_back(pair.j);
    }

    // Flatten LJ parameters
    std::vector<double> lj_flat;
    for (int i = 0; i < n_types; i++) {
        for (int j = 0; j < n_types; j++) {
            lj_flat.push_back(m_param->lj_parameter(i, j).c6);
            lj_flat.push_back(m_param->lj_parameter(i, j).c12);
        }
    }

    // Prepare outputs
    std::vector<float> forces_flat(3 * n_atoms, 0.0f);
    std::vector<double> energies(2, 0.0);
    std::vector<double> virial(9, 0.0);

    // Call Rust kernel
    rust_lj_crf_innerloop(
        pos_flat.data(),
        charges.data(),
        iac.data(),
        n_atoms,
        pairlist_flat.data(),
        n_pairs,
        lj_flat.data(),
        n_types,
        box_data.data(),
        boundary_type,
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
        forces_flat.data(),
        energies.data(),
        virial.data()
    );

    // Copy results back to GROMOS data structures
    for (size_t i = 0; i < n_atoms; i++) {
        storage.force(i)(0) = forces_flat[i * 3 + 0];
        storage.force(i)(1) = forces_flat[i * 3 + 1];
        storage.force(i)(2) = forces_flat[i * 3 + 2];
    }

    storage.energies.lj_energy = energies[0];
    storage.energies.crf_energy = energies[1];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            storage.virial_tensor(i, j) = virial[i * 3 + j];
        }
    }
#else
    // Fallback to C++ implementation
    cpp_lj_crf_innerloop(...);
#endif

)";

    std::cout << "✓ Integration pattern documented\n\n";
}

int main() {
    std::cout << "GROMOS-RS C++ Integration Examples\n";
    std::cout << "===================================\n\n";

    example_simple();
    example_water_box();
    example_gromos_integration();

    std::cout << "All examples completed successfully!\n";

    return 0;
}

/**
 * Compilation:
 *
 * # Build Rust library first
 * cd gromos-rs
 * cargo build --release
 *
 * # Compile this example
 * g++ -std=c++11 -O3 \
 *     examples/cpp_integration_example.cpp \
 *     -I. \
 *     -L./target/release \
 *     -lgromos_rs -lpthread -ldl -lm \
 *     -o cpp_example
 *
 * # Run
 * ./cpp_example
 */
