# GROMOS++ Integration Strategy for gromos-rs

## Executive Summary

**Don't reimplement - integrate!** GROMOS++ has 111 battle-tested analysis tools. gromos-rs should focus on **simulation**, GROMOS++ handles **analysis**.

## Division of Responsibilities

### gromos-rs (Rust MD Engine)
**Role**: Run simulations and **write** trajectory files

✅ **Already implemented**:
- TRC writer (positions, velocities, forces)
- TRF writer (detailed forces)
- TRE writer (energies - simplified)

**Focus**: Core MD algorithms, performance, safety

### GROMOS++ (C++ Analysis Tools)
**Role**: Read and analyze trajectory files

**111 programs available** including:
- `ener_ana` - Energy trajectory analysis
- `frameout` - Extract frames
- `trs_ana` - Trajectory statistics
- `rmsd`, `rmsf` - Structural analysis
- `hbond` - Hydrogen bonds
- `cluster` - Conformational clustering
- `rdf` - Radial distribution functions
- And 100+ more...

## Why This Strategy?

### ✅ Advantages
1. **Avoid duplication** - Don't rewrite 20+ years of development
2. **Proven code** - GROMOS++ tools are extensively tested
3. **Focus resources** - gromos-rs team focuses on MD engine
4. **User familiarity** - Existing GROMOS users know these tools
5. **Maintenance** - Don't have to maintain analysis code

### ❌ What we DON'T do
- Reimplement trajectory parsers in Rust
- Rewrite analysis algorithms
- Maintain duplicate codebases

## Workflow

```
┌─────────────┐     writes      ┌──────────────┐
│ gromos-rs   │ ──────────────> │  trajectory  │
│ (simulate)  │                 │  files       │
└─────────────┘                 │  .trc .tre   │
                                │  .trf        │
                                └──────┬───────┘
                                       │
                                  reads│
                                       v
                        ┌──────────────────────┐
                        │  GROMOS++ tools      │
                        │  (analyze)           │
                        │                      │
                        │  ener_ana, frameout  │
                        │  rmsd, hbond, etc.   │
                        └──────────────────────┘
```

## Building GROMOS++

GROMOS++ source is in `gromosPlusPlus/gromos++/`

### Prerequisites
```bash
# Install GSL (GNU Scientific Library)
sudo apt-get install libgsl-dev

# Install build tools
sudo apt-get install build-essential autoconf automake libtool
```

### Build Steps
```bash
cd gromosPlusPlus/gromos++

# Generate configure script
./Config.sh

# Configure build
./configure

# Compile (uses all cores)
make -j$(nproc)

# Optional: Install to system
# make install
```

### Binaries Location
After building:
- Binaries in: `gromosPlusPlus/gromos++/programs/`
- Or installed to: `$PREFIX/bin/` (if you ran `make install`)

## Usage Pattern

### 1. Run Simulation with gromos-rs
```rust
use gromos_rs::*;

// Run MD simulation
let mut engine = MdEngine::new(topology, config);
let mut traj_writer = TrajectoryWriter::new("output.trc", "My simulation", true, false)?;

for step in 0..num_steps {
    engine.step(dt);
    if step % output_freq == 0 {
        traj_writer.write_frame(step, step as f64 * dt, engine.config())?;
    }
}
```

### 2. Analyze with GROMOS++
```bash
# Extract energy components
ener_ana @traj output.tre @en_files ene.dat

# Calculate RMSD
rmsd @topo system.top @pbc r @traj output.trc @ref reference.cnf > rmsd.dat

# Hydrogen bond analysis
hbond @topo system.top @pbc r @traj output.trc @cut 0.25 0.135 > hbonds.dat

# Extract specific frames
frameout @topo system.top @traj output.trc @frames 0,100,200 @outformat pdb
```

## Optional: Rust Wrappers

If convenient Rust API is desired, create lightweight wrappers:

```rust
// gromos-rs/src/analysis/mod.rs (optional)
use std::process::Command;
use std::path::Path;

pub struct EnergyAnalyzer {
    gromos_bin_path: PathBuf,
}

impl EnergyAnalyzer {
    pub fn analyze(&self, tre_file: &Path) -> Result<EnergyData> {
        let output = Command::new(self.gromos_bin_path.join("ener_ana"))
            .arg("@traj").arg(tre_file)
            .arg("@en_files").arg("/tmp/ene.dat")
            .output()?;

        // Parse output
        parse_energy_file(&output.stdout)
    }
}

pub fn extract_frames(
    trc_file: &Path,
    frames: &[usize],
    gromos_bin: &Path,
) -> Result<()> {
    Command::new(gromos_bin.join("frameout"))
        .arg("@traj").arg(trc_file)
        .arg("@frames").arg(frames.iter().map(|f| f.to_string()).collect::<Vec<_>>().join(","))
        .status()?;
    Ok(())
}
```

But this is **optional** - users can call tools directly.

## GROMOS++ Tool Categories

### Trajectory Manipulation
| Tool | Purpose |
|------|---------|
| `frameout` | Extract specific frames |
| `gathtraj` | Concatenate trajectories |
| `trs_ana` | Trajectory statistics |
| `filter` | Filter frames by criteria |

### Energy Analysis
| Tool | Purpose |
|------|---------|
| `ener_ana` | Complete energy analysis |
| `ener` | Energy calculations |
| `int_ener` | Interaction energies |
| `dg_ener` | Free energy derivatives |
| `ext_ti_ana` | Thermodynamic integration |

### Structural Analysis
| Tool | Purpose |
|------|---------|
| `rmsd` | Root mean square deviation |
| `rmsf` | Root mean square fluctuation |
| `rgyr` | Radius of gyration |
| `dssp` | Secondary structure (Kabsch-Sander) |
| `sasa` | Solvent accessible surface area |
| `cry` | Crystal/lattice analysis |

### Interactions
| Tool | Purpose |
|------|---------|
| `hbond` | Hydrogen bond analysis |
| `close_pair` | Close contacts |
| `ion` | Ion pairing |
| `dipole` | Dipole moments |
| `rdf` | Radial distribution functions |

### Dynamics
| Tool | Purpose |
|------|---------|
| `diffus` | Diffusion coefficients |
| `ditrans` | Dihedral transitions |
| `epath` | Energy path analysis |

### Clustering & Sampling
| Tool | Purpose |
|------|---------|
| `cluster` | Conformational clustering |
| `follow` | Follow specific atoms |

### Free Energy
| Tool | Purpose |
|------|---------|
| `bar` | Bennett acceptance ratio |
| `ext_ti_ana` | Extended TI analysis |
| `dg_ener` | Free energy analysis |
| `eds_update_1/2` | EDS parameter optimization |

## Installation Instructions for Users

Add to `README.md` or documentation:

```markdown
## Analysis Tools

gromos-rs focuses on running simulations. For trajectory analysis, use GROMOS++:

### Installing GROMOS++

1. Install dependencies:
   ```bash
   sudo apt-get install libgsl-dev build-essential
   ```

2. Build GROMOS++:
   ```bash
   cd gromosPlusPlus/gromos++
   ./Config.sh
   ./configure
   make -j$(nproc)
   ```

3. Add to PATH (optional):
   ```bash
   export PATH=$PATH:$(pwd)/programs
   ```

### Example Analysis Workflow

```bash
# 1. Run simulation with gromos-rs
cargo run --release --bin md_sim -- --input sim.imd

# 2. Analyze with GROMOS++
ener_ana @traj output.tre @en_files energies.dat
rmsd @topo system.top @traj output.trc @ref start.cnf
hbond @topo system.top @traj output.trc
```

See [GROMOS++ documentation](https://www.gromos.net) for full tool reference.
```

## Summary

**gromos-rs**: Simulation engine (write trajectories)
**GROMOS++**: Analysis tools (read trajectories)

This division:
- ✅ Avoids code duplication
- ✅ Leverages existing expertise
- ✅ Maintains compatibility
- ✅ Focuses development efforts

**Status**: Strategy documented, GROMOS++ available but needs compilation with GSL.
