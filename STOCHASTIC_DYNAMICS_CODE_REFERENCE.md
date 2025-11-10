# Stochastic Dynamics - Code Reference and Key Algorithms

## Key C++ Snippets for Rust Translation

### 1. Friction Coefficient Calculation (Analytical Path)

**From**: `md++/src/algorithm/integration/stochastic.cc:116-139`

```cpp
// When |gamma * dt| > 0.05, use analytical formulas

const double gg = fabs(topo.stochastic().gamma(i));
const double gdt = gg * sim.time_step_size();
const double gdth = gdt * 0.5;
const double kToverM = sqrt(math::k_Boltzmann * sim.param().stochastic.temp / topo.mass(i));

const double emdth = exp(-gdth);      // exp(-gdt/2)
const double epdth = exp(+gdth);      // exp(+gdt/2)
const double emdt  = emdth * emdth;   // exp(-gdt)
const double epdt  = epdth * epdth;   // exp(+gdt)
const double omdt  = 1.0 - emdt;      // 1 - exp(-gdt)
const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;
const double ddth  = 2.0 - epdth - emdth;
const double bpdth = gdt * (epdt - 1.0) - 4.0 * (epdth - 1.0) * (epdth - 1.0);
const double bmdth = gdt * (1.0 - emdt) - 4.0 * (emdth - 1.0) * (emdth - 1.0);

// Store coefficients
topo.stochastic().c1(i) = emdt;
topo.stochastic().c2(i) = omdt / gdt;
topo.stochastic().c3(i) = sqrt(fabs(omdt));
topo.stochastic().c4(i) = sqrt(fabs(bpdth/cdth));
topo.stochastic().c5(i) = gg * ddth/cdth;
topo.stochastic().c6(i) = (epdth - emdth) / gdt;
topo.stochastic().c7(i) = sqrt(fabs(cdth)) / gg;
topo.stochastic().c8(i) = sqrt(fabs(bmdth/omdt)) / gg;
topo.stochastic().c9(i) = -ddth/(gg * omdt);
```

### 2. Friction Coefficient Calculation (Power Series Path)

**From**: `md++/src/algorithm/integration/stochastic.cc:144-199`

```cpp
// When |gamma * dt| <= 0.05, use power series expansion

const double gdth2 = gdth * gdth;
const double gdth3 = gdth2 * gdth;
const double gdth4 = gdth2 * gdth2;
const double gdth5 = gdth2 * gdth3;
const double gdth6 = gdth3 * gdth3;
const double gdth7 = gdth4 * gdth3;

topo.stochastic().c1(i) = exp(-gdt);

topo.stochastic().c2(i) = 1.0 - gdth + gdth2 * 2.0/3.0 - gdth3/3.0 
    + gdth4 * 2.0/15.0 - gdth5 * 2.0/45.0 + gdth6 * 4.0/315.0;

topo.stochastic().c3(i) = sqrt(fabs(topo.stochastic().c2(i)) * 2.0 * gdth);

topo.stochastic().c4(i) = sqrt(fabs(gdth/2.0 + gdth2 * 7.0/8.0 + gdth3 * 367.0/480.0 
    + gdth4 * 857.0/1920.0 + gdth5 * 52813.0/268800.0 
    + gdth6 * 224881.0/3225600.0 + gdth7 * 1341523.0/64512000.0));

topo.stochastic().c5(i) = -2.0 / sim.time_step_size() *
    (1.5 + gdth * 9.0/8.0 + gdth2 * 71.0/160.0 + gdth3 * 81.0/640.0 
     + gdth4 * 7807.0/268800.0 + gdth5 * 1971.0/358400.0 
     + gdth6 * 56417.0/64512000.0);

topo.stochastic().c6(i) = 1.0 + gdth2/6.0 + gdth4/120.0 + gdth6/5040.0;

topo.stochastic().c7(i) = sim.time_step_size() * 0.5 * 
    sqrt(fabs(gdth * 2.0/3.0 - gdth2/2.0 + gdth3 * 7.0/30.0 
        - gdth4/12.0 + gdth5 * 31.0/1260.0 - gdth6/160.0 
        + gdth7 * 127.0/90720.0));

topo.stochastic().c8(i) = sim.time_step_size() * 0.5 * 
    sqrt(fabs(gdth/6.0 - gdth3/60.0 + gdth5 * 17.0/10080.0 
        - gdth7 * 31.0/181440.0));

topo.stochastic().c9(i) = sim.time_step_size() * 0.5 *
    (0.5 + gdth/2.0 + gdth2 * 5.0/24.0 + gdth3/24.0 
     + gdth4/240.0 + gdth5/720.0 + gdth6 * 5.0/8064.0);
```

### 3. SASA Friction Calculation (Neighbor Counting)

**From**: `md++/src/algorithm/integration/stochastic.cc:65-96`

```cpp
// When NTFR == 3, calculate friction from solvent accessibility

std::vector<unsigned int> neigh(size, 0);
double const cutoff2 = sim.param().stochastic.rcutf * sim.param().stochastic.rcutf;

// Count neighbors within cutoff
for (unsigned int i = 0; i < size; ++i) {
    for (unsigned int j = i+1; j < size; ++j) {
        // Calculate distance with periodic boundary conditions
        periodicity.nearest_image(conf.current().pos(i), conf.current().pos(j), rij);
        const double r2 = math::abs2(rij);
        
        if (r2 < cutoff2) {
            neigh[i] += 1;
            neigh[j] += 1;
        }
    }
}

// Determine gamma values
for (unsigned int i = 0; i < size; ++i) {
    double xh = 0.0;
    switch(sim.param().stochastic.ntfr) {
        case 0: xh = 0.0;                                           break;
        case 1: xh = 1.0;                                           break;
        case 2: xh = topo.stochastic().gamma(i);                    break;
        case 3: xh = std::max(0.0, 1.0 - neigh[i]/double(sim.param().stochastic.nbref)); break;
    }
    topo.stochastic().gamma(i) = sim.param().stochastic.cfric * xh;
}
```

### 4. Velocity Update (Equation 2.11.2.2)

**From**: `md++/src/algorithm/integration/stochastic.cc:330-415`

```cpp
// Core velocity update loop

for (unsigned int i=0; i < topo.num_atoms(); ++i) {
    if(topo.stochastic().gamma(i) != 0.0) {
        
        const double kToverM = sqrt(math::k_Boltzmann * 
                                   sim.param().stochastic.temp / topo.mass(i));
        
        // Scaling factors for random vectors
        double cf = sim.time_step_size() * topo.stochastic().c2(i) / topo.mass(i);
        double sd1 = topo.stochastic().c3(i) * kToverM;
        double sd2 = kToverM * topo.stochastic().c4(i);
        
        // Generate random vectors (Gaussian distributed)
        m_rng->stddev(sd1);
        m_vrand.vrand1(i) = m_rng->get_gaussian_vec();  // For new integral
        m_rng->stddev(sd2);
        m_vrand.vrand2(i) = m_rng->get_gaussian_vec();  // For svh
        
        // Calculate svh term
        math::Vec svh = conf.old().stochastic_integral(i) * 
                        topo.stochastic().c5(i) + m_vrand.vrand2(i);
        
        // Store vrand1 as new stochastic integral
        conf.current().stochastic_integral(i) = m_vrand.vrand1(i);
        
        // Update velocity (equation 2.11.2.2)
        conf.current().vel(i) = (conf.old().vel(i) - svh) * topo.stochastic().c1(i)
            + conf.old().force(i) * cf
            + m_vrand.vrand1(i);
        
        // Generate random vectors for position update (Vel2/Pos2)
        double sd3 = kToverM * topo.stochastic().c7(i);
        double sd4 = kToverM * topo.stochastic().c8(i);
        
        m_rng->stddev(sd3);
        m_vrand.vrand3(i) = m_rng->get_gaussian_vec();
        m_rng->stddev(sd4);
        m_vrand.vrand4(i) = m_rng->get_gaussian_vec();
        
    } else {
        // No friction: standard Verlet step
        conf.current().vel(i) = conf.old().vel(i) + 
                                conf.old().force(i) * sim.time_step_size() / topo.mass(i);
        m_vrand.vrand1(i) = math::Vec(0.0, 0.0, 0.0);
        m_vrand.vrand2(i) = math::Vec(0.0, 0.0, 0.0);
        m_vrand.vrand3(i) = math::Vec(0.0, 0.0, 0.0);
        m_vrand.vrand4(i) = math::Vec(0.0, 0.0, 0.0);
    }
}
```

### 5. Position Update (Vel1 → Pos1)

**From**: `md++/src/algorithm/integration/stochastic.cc:436-442`

```cpp
// Simple position update after velocity calculation

for (unsigned int i=0; i < topo.num_atoms(); ++i) {
    conf.current().pos(i) = conf.old().pos(i) 
                          + conf.current().vel(i) * sim.time_step_size() * 
                            topo.stochastic().c6(i);
}
```

### 6. Velocity Correction for Constraints (Vel2)

**From**: `md++/src/algorithm/integration/stochastic.cc:481-487`

```cpp
// After SHAKE constraint solving, correct velocities

if (sim.param().constraint.solute.algorithm == simulation::constr_shake) {
    for (unsigned int i=0; i < topo.num_atoms(); ++i) {
        double cinv = 1.0 / (topo.stochastic().c6(i) * sim.time_step_size());
        math::Vec r;
        periodicity.nearest_image(conf.current().pos(i), conf.old().pos(i), r);
        conf.current().vel(i) = r * cinv;  // Eq. 2.11.2.24
    }
}
```

### 7. Position Correction for Constraints (Pos2)

**From**: `md++/src/algorithm/integration/stochastic.cc:514-522`

```cpp
// Additional position correction for constrained systems

for (unsigned int i=0; i < topo.num_atoms(); ++i) {
    if(topo.stochastic().gamma(i) != 0.0) {
        // Equation 2.11.2.25
        math::Vec sxh = conf.current().stochastic_integral(i) * topo.stochastic().c9(i)
                      + m_vrand->vrand4(i);
        conf.current().stochastic_integral(i) = m_vrand->vrand3(i);
        
        // Equation 2.11.2.26
        conf.current().pos(i) += m_vrand->vrand3(i) - sxh;
    }
}
```

## Rust Translation Strategy

### Module Structure
```rust
// src/algorithm/stochastic_dynamics.rs

pub mod friction {
    /// Calculate friction coefficients (c1-c9)
    pub fn calc_coefficients_analytical(...) -> StochasticCoefficients;
    pub fn calc_coefficients_power_series(...) -> StochasticCoefficients;
    pub fn calc_friction_sasa(...) -> Vec<f64>;  // gamma values
}

pub mod integration {
    pub struct StochasticDynamics1 { /* Vel1 + Pos1 */ }
    pub struct StochasticDynamics2 { /* Vel2 + Pos2 */ }
    
    impl Integrator for StochasticDynamics1 { ... }
    impl Integrator for StochasticDynamics2 { ... }
}

pub mod rng {
    // Use rand crate or implement G96 algorithm
}

pub struct StochasticState {
    gamma: Vec<f64>,
    c1_c9: Vec<[f64; 9]>,  // Store all coefficients
    stochastic_integral: Vec<Vec3>,
    rng_seed: String,
}
```

### Critical Numerical Considerations

1. **Stability check**: Always verify |gdt| threshold before choosing algorithm path
2. **Sqrt of negative numbers**: Use `abs()` before sqrt (defensive programming)
3. **Division by gdt**: Handle case when gdt approaches 0
4. **kToverM calculation**: Ensure temperature > 0 before sqrt
5. **Power series coefficients**: Use exact rational numbers from C++ code

### Memory Layout Optimization

**Current (Array-of-Structs)**:
```rust
// Not ideal
struct Atom {
    gamma: f64,
    c1: f64, c2: f64, ..., c9: f64,
    stochastic_integral: Vec3,
}
```

**Optimized (Structure-of-Arrays)**:
```rust
struct StochasticState {
    gamma: Vec<f64>,
    coefficients: StochasticCoefficients {
        c1: Vec<f64>,
        c2: Vec<f64>,
        ...,
        c9: Vec<f64>,
    },
    stochastic_integral: Vec<Vec3>,
}
```

### Testing Requirements

1. **Coefficient calculation tests**: Compare analytical vs power series paths
2. **Accuracy tests**: Verify against reference C++ outputs
3. **Energy conservation tests**: Run short NVE simulations
4. **Random number reproducibility**: Verify same seed produces same trajectory
5. **SASA calculation tests**: Verify neighbor counting matches C++

