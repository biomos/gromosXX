#!/usr/bin/env python3
"""
Hands-On Demo: GROMOS-RS Python Bindings
=========================================

This script provides interactive demonstrations of the key concepts.
Run each demo to see the bindings in action!
"""

import sys
import time
import numpy as np

def try_import():
    """Try to import gromos"""
    try:
        import gromos
        return gromos
    except ImportError:
        print("❌ Error: gromos module not found")
        print("\nPlease build the package first:")
        print("  cd py-gromos")
        print("  pip install maturin")
        print("  maturin develop --release")
        sys.exit(1)

def demo_1_basic_vectors(gromos):
    """Demo 1: Basic vector operations"""
    print("\n" + "="*70)
    print("DEMO 1: Basic Vector Operations (SIMD-accelerated)")
    print("="*70)

    print("\n1. Creating vectors:")
    v1 = gromos.Vec3(1.0, 0.0, 0.0)
    v2 = gromos.Vec3(0.0, 1.0, 0.0)
    v3 = gromos.Vec3(3.0, 4.0, 0.0)

    print(f"  v1 = {v1}")
    print(f"  v2 = {v2}")
    print(f"  v3 = {v3}")

    print("\n2. Vector operations (all SIMD-accelerated):")
    print(f"  v1 + v2 = {v1 + v2}")
    print(f"  v1 - v2 = {v1 - v2}")
    print(f"  v1 * 2.5 = {v1 * 2.5}")

    print("\n3. Vector mathematics:")
    dot = v1.dot(v2)
    cross = v1.cross(v2)
    length = v3.length()
    distance = v1.distance(v2)

    print(f"  v1 · v2 = {dot} (orthogonal!)")
    print(f"  v1 × v2 = {cross}")
    print(f"  |v3| = {length:.4f} (should be 5.0)")
    print(f"  distance(v1, v2) = {distance:.4f}")

    print("\n4. Normalization:")
    v3_norm = v3.normalize()
    print(f"  normalize(v3) = {v3_norm}")
    print(f"  |normalize(v3)| = {v3_norm.length():.4f} (should be 1.0)")

    input("\n✓ Demo complete. Press Enter to continue...")

def demo_2_numpy_integration(gromos):
    """Demo 2: NumPy integration"""
    print("\n" + "="*70)
    print("DEMO 2: Zero-Copy NumPy Integration")
    print("="*70)

    print("\n1. Converting Vec3 to NumPy:")
    v = gromos.Vec3(1.5, 2.5, 3.5)
    arr = v.to_numpy()
    print(f"  Vec3: {v}")
    print(f"  NumPy: {arr}")
    print(f"  Type: {type(arr)}")
    print(f"  Dtype: {arr.dtype}")
    print(f"  Shape: {arr.shape}")

    print("\n2. Creating Vec3 from NumPy:")
    np_vec = np.array([4.0, 5.0, 6.0], dtype=np.float32)
    v2 = gromos.Vec3.from_numpy(np_vec)
    print(f"  NumPy: {np_vec}")
    print(f"  Vec3: {v2}")

    print("\n3. Working with large arrays:")
    n = 10000
    print(f"  Creating {n:,} random positions...")
    positions = np.random.rand(n, 3).astype(np.float32) * 10.0

    print(f"  Creating State...")
    state = gromos.State(num_atoms=n, num_temp_groups=1, num_energy_groups=1)

    print(f"  Setting positions...")
    start = time.time()
    state.set_positions(positions)
    set_time = time.time() - start

    print(f"  Getting positions back...")
    start = time.time()
    retrieved = state.positions()
    get_time = time.time() - start

    print(f"\n  Timings:")
    print(f"    set_positions: {set_time*1000:.2f} ms")
    print(f"    positions():   {get_time*1000:.2f} ms (zero-copy!)")
    print(f"  Array info:")
    print(f"    Shape: {retrieved.shape}")
    print(f"    Memory: {retrieved.nbytes/1024:.2f} KB")
    print(f"    Values match: {np.allclose(positions, retrieved)}")

    input("\n✓ Demo complete. Press Enter to continue...")

def demo_3_performance(gromos):
    """Demo 3: Performance comparison"""
    print("\n" + "="*70)
    print("DEMO 3: Performance Comparison (SIMD vs NumPy)")
    print("="*70)

    n_iter = 50000
    print(f"\nRunning {n_iter:,} iterations...")

    # Create test vectors
    v1 = gromos.Vec3(1.0, 2.0, 3.0)
    v2 = gromos.Vec3(4.0, 5.0, 6.0)
    arr1 = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    arr2 = np.array([4.0, 5.0, 6.0], dtype=np.float32)

    # Benchmark 1: Dot product
    print("\n1. Dot product:")
    print("  Vec3 (SIMD)...", end=" ", flush=True)
    start = time.time()
    for _ in range(n_iter):
        result = v1.dot(v2)
    vec3_dot_time = time.time() - start
    print(f"{vec3_dot_time:.4f}s")

    print("  NumPy.......", end=" ", flush=True)
    start = time.time()
    for _ in range(n_iter):
        result = np.dot(arr1, arr2)
    numpy_dot_time = time.time() - start
    print(f"{numpy_dot_time:.4f}s")
    print(f"  Speedup: {numpy_dot_time/vec3_dot_time:.2f}×")

    # Benchmark 2: Length calculation
    print("\n2. Length calculation:")
    print("  Vec3 (SIMD)...", end=" ", flush=True)
    start = time.time()
    for _ in range(n_iter):
        result = v1.length()
    vec3_len_time = time.time() - start
    print(f"{vec3_len_time:.4f}s")

    print("  NumPy.......", end=" ", flush=True)
    start = time.time()
    for _ in range(n_iter):
        result = np.linalg.norm(arr1)
    numpy_len_time = time.time() - start
    print(f"{numpy_len_time:.4f}s")
    print(f"  Speedup: {numpy_len_time/vec3_len_time:.2f}×")

    # Benchmark 3: Distance calculation
    print("\n3. Distance calculation:")
    print("  Vec3 (SIMD)...", end=" ", flush=True)
    start = time.time()
    for _ in range(n_iter):
        result = v1.distance(v2)
    vec3_dist_time = time.time() - start
    print(f"{vec3_dist_time:.4f}s")

    print("  NumPy.......", end=" ", flush=True)
    start = time.time()
    for _ in range(n_iter):
        result = np.linalg.norm(arr1 - arr2)
    numpy_dist_time = time.time() - start
    print(f"{numpy_dist_time:.4f}s")
    print(f"  Speedup: {numpy_dist_time/vec3_dist_time:.2f}×")

    print("\n4. Summary:")
    avg_speedup = (
        numpy_dot_time/vec3_dot_time +
        numpy_len_time/vec3_len_time +
        numpy_dist_time/vec3_dist_time
    ) / 3
    print(f"  Average speedup: {avg_speedup:.2f}×")
    print("  Why? SIMD vectorization + compiled Rust code!")

    input("\n✓ Demo complete. Press Enter to continue...")

def demo_4_energy_tracking(gromos):
    """Demo 4: Energy tracking"""
    print("\n" + "="*70)
    print("DEMO 4: Energy Tracking")
    print("="*70)

    print("\n1. Creating energy tracker:")
    energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)
    print(f"  {energy}")

    print("\n2. Energy components (kJ/mol):")
    print(f"  Kinetic:   {energy.kinetic:>10.2f}")
    print(f"  Potential: {energy.potential:>10.2f}")
    print(f"  Bond:      {energy.bond:>10.2f}")
    print(f"  Angle:     {energy.angle:>10.2f}")
    print(f"  Dihedral:  {energy.dihedral:>10.2f}")
    print(f"  LJ:        {energy.lj:>10.2f}")
    print(f"  Coulomb:   {energy.coulomb:>10.2f}")
    print(f"  " + "-"*30)
    print(f"  Total:     {energy.total():>10.2f}")

    print("\n3. Getting energy as dictionary:")
    energy_dict = energy.to_dict()
    print("  Keys:", list(energy_dict.keys()))

    print("\n4. Calculating kinetic energy from velocities:")
    n = 100
    temperature = 300.0  # K
    kb = 0.00831446261815324  # kJ/(mol·K)
    mass = 40.0  # g/mol

    vel_std = np.sqrt(kb * temperature / mass)
    velocities = np.random.randn(n, 3).astype(np.float32) * vel_std

    v_squared = (velocities ** 2).sum(axis=1)
    ke = 0.5 * mass * v_squared.sum() / 1000.0

    expected_ke = 1.5 * n * kb * temperature

    print(f"  Calculated KE: {ke:.2f} kJ/mol")
    print(f"  Expected KE:   {expected_ke:.2f} kJ/mol")
    print(f"  Difference:    {abs(ke - expected_ke):.2f} kJ/mol")
    print(f"  Match: {np.isclose(ke, expected_ke, rtol=0.1)}")

    input("\n✓ Demo complete. Press Enter to continue...")

def demo_5_simulation_boxes(gromos):
    """Demo 5: Simulation boxes"""
    print("\n" + "="*70)
    print("DEMO 5: Simulation Boxes & Periodic Boundaries")
    print("="*70)

    print("\n1. Vacuum box (no periodicity):")
    box_vacuum = gromos.Box.vacuum()
    print(f"  {box_vacuum}")

    print("\n2. Rectangular box:")
    box_rect = gromos.Box.rectangular(5.0, 5.0, 5.0)
    print(f"  {box_rect}")
    print(f"  Dimensions: {box_rect.dimensions()}")
    print(f"  Volume: {box_rect.volume():.2f} nm³")

    print("\n3. Different sizes:")
    sizes = [(3, 3, 3), (5, 5, 5), (10, 10, 10), (20, 20, 20)]
    print(f"  {'Size (nm)':<15} {'Volume (nm³)':<15} {'Density (atoms/nm³)':<20}")
    print("  " + "-"*50)
    for lx, ly, lz in sizes:
        box = gromos.Box.rectangular(lx, ly, lz)
        volume = box.volume()
        n_atoms = 1000
        density = n_atoms / volume
        print(f"  {lx}×{ly}×{lz:<10} {volume:<15.1f} {density:<20.3f}")

    print("\n4. Triclinic box (non-orthogonal):")
    v1 = gromos.Vec3(5.0, 0.0, 0.0)
    v2 = gromos.Vec3(1.0, 5.0, 0.0)  # Tilted
    v3 = gromos.Vec3(0.0, 0.0, 5.0)
    mat = gromos.Mat3.from_cols(v1, v2, v3)
    box_triclinic = gromos.Box.triclinic(mat)
    print(f"  {box_triclinic}")
    print(f"  Volume: {box_triclinic.volume():.2f} nm³")

    input("\n✓ Demo complete. Press Enter to continue...")

def demo_6_complete_system(gromos):
    """Demo 6: Complete system setup"""
    print("\n" + "="*70)
    print("DEMO 6: Complete Molecular System")
    print("="*70)

    n_atoms = 256
    box_size = 4.0

    print(f"\n1. System parameters:")
    print(f"  N atoms:      {n_atoms}")
    print(f"  Box size:     {box_size} nm")
    print(f"  Temperature:  300 K")

    print(f"\n2. Creating components...")
    state = gromos.State(num_atoms=n_atoms, num_temp_groups=1, num_energy_groups=1)
    box = gromos.Box.rectangular(box_size, box_size, box_size)
    energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)
    topology = gromos.Topology()
    config = gromos.Configuration(num_atoms=n_atoms, num_temp_groups=1, num_energy_groups=1)

    print("  ✓ State")
    print("  ✓ Box")
    print("  ✓ Energy")
    print("  ✓ Topology")
    print("  ✓ Configuration")

    print(f"\n3. Initializing positions (FCC lattice)...")
    n_side = int(np.ceil(n_atoms ** (1/3)))
    a = box_size / n_side
    positions = []
    for i in range(n_atoms):
        ix = i % n_side
        iy = (i // n_side) % n_side
        iz = i // (n_side * n_side)
        x = (ix + 0.5) * a
        y = (iy + 0.5) * a
        z = (iz + 0.5) * a
        positions.append([x, y, z])
    positions = np.array(positions, dtype=np.float32)
    state.set_positions(positions)
    print(f"  ✓ Positions initialized ({n_atoms} atoms on {n_side}³ lattice)")

    print(f"\n4. Initializing velocities (Maxwell-Boltzmann at 300 K)...")
    T = 300.0
    kb = 0.00831446261815324
    mass = 40.0
    vel_std = np.sqrt(kb * T / mass)
    velocities = np.random.randn(n_atoms, 3).astype(np.float32) * vel_std
    velocities -= velocities.mean(axis=0)  # Remove COM motion
    state.set_velocities(velocities)
    print(f"  ✓ Velocities initialized (σ = {vel_std:.4f} nm/ps)")

    print(f"\n5. System summary:")
    print(f"  State:      {state}")
    print(f"  Box:        {box}")
    print(f"  Energy:     {energy}")
    print(f"  Topology:   {topology}")
    print(f"  Config:     {config}")

    print(f"\n6. Checking system properties:")
    pos = state.positions()
    vel = state.velocities()
    print(f"  Position range: [{pos.min():.2f}, {pos.max():.2f}] nm")
    print(f"  Velocity range: [{vel.min():.2f}, {vel.max():.2f}] nm/ps")
    print(f"  RMS velocity:   {np.sqrt((vel**2).mean()):.4f} nm/ps")
    print(f"  Density:        {n_atoms/box.volume():.3f} atoms/nm³")

    print("\n7. System is ready for simulation!")
    print("   (Use GROMOS-RS MD binary for actual simulation)")

    input("\n✓ Demo complete. Press Enter to continue...")

def main():
    """Main function"""
    print("""
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║               GROMOS-RS Python Bindings                          ║
║               Hands-On Interactive Demo                          ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝

This demo will walk you through the key features of the Python bindings.
Each demo is interactive - press Enter to proceed through them.
""")

    # Try to import
    gromos = try_import()

    print(f"\n✓ Successfully imported gromos v{gromos.__version__}")
    print(f"✓ Available classes: {len(gromos.__all__)} types")

    try:
        input("\nPress Enter to start the demos...")

        demo_1_basic_vectors(gromos)
        demo_2_numpy_integration(gromos)
        demo_3_performance(gromos)
        demo_4_energy_tracking(gromos)
        demo_5_simulation_boxes(gromos)
        demo_6_complete_system(gromos)

        print("\n" + "="*70)
        print("ALL DEMOS COMPLETE!")
        print("="*70)

        print("""
You've seen:
✓ SIMD-accelerated vector operations
✓ Zero-copy NumPy integration
✓ Performance comparisons
✓ Energy tracking
✓ Simulation boxes
✓ Complete system setup

Next steps:
- Try the examples in py-gromos/examples/
- Read the notebooks in py-gromos/notebooks/
- Check the API reference: py-gromos/API_REFERENCE.md
- Run actual MD simulations with the GROMOS-RS binaries

Happy simulating!
""")

    except KeyboardInterrupt:
        print("\n\n❌ Interrupted by user.")
        sys.exit(1)

if __name__ == "__main__":
    main()
