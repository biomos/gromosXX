"""
Example 3: MD Integrators
==========================

Demonstrates different integration algorithms for molecular dynamics.
"""

import gromos

print("=" * 60)
print("GROMOS-RS: MD Integrators")
print("=" * 60)

# Leap-Frog integrator (velocity Verlet variant)
print("\n--- Leap-Frog Integrator ---")
dt = 0.002  # 2 fs timestep
leapfrog = gromos.LeapFrog(dt=dt)
print(f"{leapfrog}")
print(f"Timestep: {leapfrog.timestep()} ps")
print(f"Timestep: {leapfrog.timestep() * 1000:.2f} fs")

print("\nCharacteristics:")
print("  - Fast and stable")
print("  - Velocities offset by dt/2 from positions")
print("  - Good energy conservation")
print("  - Most commonly used in MD")

# Velocity Verlet integrator
print("\n--- Velocity Verlet Integrator ---")
dt = 0.001  # 1 fs timestep
vv = gromos.VelocityVerlet(dt=dt)
print(f"{vv}")
print(f"Timestep: {vv.timestep()} ps")
print(f"Timestep: {vv.timestep() * 1000:.2f} fs")

print("\nCharacteristics:")
print("  - Higher accuracy than Leap-Frog")
print("  - Positions and velocities at same time")
print("  - Better for analysis")
print("  - Slightly more expensive")

# Stochastic Dynamics (Langevin)
print("\n--- Stochastic Dynamics (Langevin) ---")
dt = 0.002  # 2 fs
gamma = 0.1  # friction coefficient (1/ps)
temperature = 300.0  # K
sd = gromos.StochasticDynamics(dt=dt, gamma=gamma, temperature=temperature)
print(f"{sd}")
print(f"Timestep: {sd.timestep()} ps")
print(f"Friction: {gamma} ps⁻¹")
print(f"Temperature: {temperature} K")

print("\nCharacteristics:")
print("  - Implicit solvent representation")
print("  - Friction + random forces (Brownian motion)")
print("  - Natural temperature control")
print("  - No need for separate thermostat")
print("  - Good for small systems or coarse-grained models")

# Typical simulation parameters
print("\n" + "=" * 60)
print("Typical Simulation Parameters")
print("=" * 60)

configs = [
    ("Small molecule (vacuum)", 0.0005, "LeapFrog"),
    ("Protein in water (explicit)", 0.002, "LeapFrog"),
    ("Coarse-grained system", 0.010, "VelocityVerlet"),
    ("Implicit solvent", 0.002, "StochasticDynamics"),
]

print("\n{:<30} {:<15} {:<20}".format("System", "Timestep", "Integrator"))
print("-" * 65)
for system, dt, integrator in configs:
    print("{:<30} {:<15} {:<20}".format(
        system,
        f"{dt} ps ({dt*1000:.1f} fs)",
        integrator
    ))

# Performance considerations
print("\n" + "=" * 60)
print("Performance & Stability Trade-offs")
print("=" * 60)

print("\nTimestep Selection:")
print("  - Larger dt = faster simulation")
print("  - Smaller dt = more stable, better energy conservation")
print("  - Rule of thumb: dt ≤ 1/10 of fastest motion period")
print("  - Bonds: ~10 fs period → dt ≤ 1 fs (or use constraints)")
print("  - Angles: ~20-50 fs → dt ≤ 2-5 fs")
print("  - With SHAKE/LINCS: dt = 2 fs typical")

print("\nIntegrator Choice:")
print("  LeapFrog:")
print("    + Fast (1 force evaluation per step)")
print("    + Excellent energy conservation")
print("    - Velocities at half-steps (needs correction for analysis)")
print("  ")
print("  VelocityVerlet:")
print("    + Positions/velocities synchronized")
print("    + Better for certain analyses")
print("    = Similar performance to LeapFrog")
print("  ")
print("  StochasticDynamics:")
print("    + Built-in temperature control")
print("    + Good for implicit solvent")
print("    - Random forces → energy not conserved")
print("    - Slower dynamics (friction)")

print("\n" + "=" * 60)
print("Integrator examples completed successfully!")
print("=" * 60)
