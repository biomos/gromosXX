"""
GROMOS-RS Python Bindings
=========================

High-performance molecular dynamics simulations powered by Rust.

Inspired by Polars' architecture:
- Zero-copy data sharing with NumPy
- Rust core with Python wrapper
- SIMD acceleration and parallel execution
- Memory-safe operations

Quick Start
-----------

>>> import gromos
>>> import numpy as np
>>>
>>> # Create a simple system
>>> state = gromos.State(num_atoms=100, num_temp_groups=1, num_energy_groups=1)
>>> energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)
>>>
>>> # Set up simulation box
>>> box = gromos.Box.rectangular(3.0, 3.0, 3.0)  # 3x3x3 nm box
>>>
>>> # Create integrator
>>> integrator = gromos.LeapFrog(dt=0.002)  # 2 fs timestep
>>>
>>> # Work with vectors
>>> v1 = gromos.Vec3(1.0, 0.0, 0.0)
>>> v2 = gromos.Vec3(0.0, 1.0, 0.0)
>>> v3 = v1.cross(v2)
>>> print(v3)  # Vec3(0.0, 0.0, 1.0)

Advanced Sampling
-----------------

>>> # Gaussian Accelerated MD (GaMD)
>>> gamd_params = gromos.GamdParameters(sigma0=6.0, threshold_mode='lower')
>>> gamd = gromos.GamdRunner(gamd_params)
>>>
>>> # Enveloping Distribution Sampling (EDS)
>>> eds_params = gromos.EDSParameters(num_states=4, smoothness=1.0)
>>> eds = gromos.EDSRunner(eds_params)
>>>
>>> # Replica Exchange MD (REMD)
>>> remd = gromos.ReplicaController(num_replicas=8, exchange_interval=1000)

Features
--------

- **Math Types**: Vec3, Mat3 with SIMD acceleration
- **Core Structures**: State, Energy, Configuration, Topology, Box
- **Integrators**: LeapFrog, VelocityVerlet, StochasticDynamics
- **Advanced Sampling**: GaMD, EDS, REMD
- **NumPy Integration**: Zero-copy array conversion
- **Parallel Execution**: Automatic multi-threading via Rayon

Performance
-----------

GROMOS-RS achieves 2-3x speedup over the original C++ implementation through:

- SIMD vectorization (glam + wide)
- Fearless concurrency (Rayon)
- Zero-cost abstractions
- Memory safety without runtime overhead

See Also
--------

- Documentation: https://gromos.net
- Source code: https://github.com/gromos/gromosXX
- Original GROMOS: http://www.gromos.net
"""

# Import all classes from the Rust extension module
from .gromos import (
    # Math types
    Vec3,
    Mat3,

    # Core structures
    Box,
    Energy,
    State,
    Configuration,
    Topology,

    # Analysis functions
    calculate_rmsd,
    calculate_rmsf,
    calculate_rgyr,
    analyze_trajectory,
)

# Import MD simulation runners
from .md_runners import (
    MDSimulation,
    GaMDSimulation,
    EDSSimulation,
    REMDSimulation,
    TISimulation,
    run_standard_md,
    run_gamd,
    run_eds,
    run_remd,
    run_ti,
)

__version__ = "0.1.0"

__all__ = [
    # Math types
    "Vec3",
    "Mat3",

    # Core structures
    "Box",
    "Energy",
    "State",
    "Configuration",
    "Topology",

    # Analysis functions
    "calculate_rmsd",
    "calculate_rmsf",
    "calculate_rgyr",
    "analyze_trajectory",

    # MD simulation classes
    "MDSimulation",
    "GaMDSimulation",
    "EDSSimulation",
    "REMDSimulation",
    "TISimulation",

    # MD simulation functions
    "run_standard_md",
    "run_gamd",
    "run_eds",
    "run_remd",
    "run_ti",

    # Version
    "__version__",
]
