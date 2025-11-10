# GROMOS Tutorial Compatibility Plan

**Goal**: Make gromos-rs compatible with the official GROMOS tutorials from biomos/gromos_tutorial_livecoms

**Target**: Tutorial 01 - Basic MD Simulation Workflow

---

## Tutorial Structure

**Repository**: https://github.com/biomos/gromos_tutorial_livecoms
**Target Tutorial**: `review_process/published_2.0/04_tutorial_01/tutorial_01.tex`

The tutorial teaches the complete MD workflow:
1. System preparation (topology, coordinates)
2. Solvation and ionization
3. Energy minimization
4. Equilibration (NVT, NPT)
5. Production MD
6. Analysis

---

## Required Programs (From Tutorial)

### System Preparation
| Program | Function | Source | Status |
|---------|----------|--------|--------|
| **pdb2g96** | Convert PDB to GROMOS coordinates | gromosPlusPlus/programs/pdb2g96.cc | ✅ **IMPLEMENTED** (gromos-rs/src/bin/pdb2g96.rs) |
| **make_top** | Generate topology from molecular specs | gromosPlusPlus/programs/make_top.cc | ⏳ **DEFERRED** (requires .mtb parser) |
| **com_top** | Combine multiple topology files | gromosPlusPlus/programs/com_top.cc | ✅ **IMPLEMENTED** (gromos-rs/src/bin/com_top.rs) |
| **check_top** | Validate topology consistency | gromosPlusPlus/programs/check_top.cc | ✅ **IMPLEMENTED** (gromos-rs/src/bin/check_top.rs) |
| **gch** | Add hydrogens (X-ray structures) | gromosPlusPlus/programs/gch.cc | ❌ Not implemented |

### Solvation & Ions
| Program | Function | Source | Status |
|---------|----------|--------|--------|
| **sim_box** | Solvate protein, create box | gromosPlusPlus/programs/sim_box.cc | ❌ Not implemented |
| **ion** | Place counter ions | gromosPlusPlus/programs/ion.cc | ❌ Not implemented |

### Simulation
| Program | Function | Source | Status |
|---------|----------|--------|--------|
| **mk_script** | Generate MD input files (IMD) | gromosPlusPlus/programs/mk_script.cc | ✅ **IMPLEMENTED** (gromos-rs/src/bin/mk_script.rs) |
| **md** | Run MD simulation | md++/src/ | ⚠️ **Core engine in gromos-rs** (needs main binary) |

### Analysis
| Program | Function | Source | Status |
|---------|----------|--------|--------|
| **ene_ana** | Analyze energy trajectories | gromosPlusPlus/programs/ene_ana.cc | ❌ Not implemented |
| **rmsd** | Calculate RMSD | gromosPlusPlus/programs/rmsd.cc | ❌ Not implemented |
| **rmsf** | Calculate RMS fluctuations | gromosPlusPlus/programs/rmsf.cc | ❌ Not implemented |
| **nhoparam** | Calculate N-H order parameters | gromosPlusPlus/programs/nhoparam.cc | ❌ Not implemented |

---

## Implementation Priority

### Phase 1: System Preparation (2-3 weeks)
**Goal**: Enable users to prepare systems for MD

1. **pdb2g96** (High Priority)
   - Parse PDB files
   - Convert to GROMOS G96 format
   - Handle atom nomenclature differences
   - ~500-800 lines

2. **make_top** (High Priority)
   - Read molecular topology specifications
   - Generate GROMOS topology files
   - Support building blocks (amino acids, nucleotides)
   - ~1,000-1,500 lines

3. **com_top** (Medium Priority)
   - Merge multiple topology files
   - Renumber atoms/molecules
   - ~300-500 lines

4. **check_top** (Medium Priority)
   - Validate topology consistency
   - 34 logical checks (as per tutorial)
   - ~400-600 lines

### Phase 2: Solvation & Ionization (1-2 weeks)
**Goal**: Solvate proteins and add ions

5. **sim_box** (High Priority)
   - Define simulation box
   - Solvate with water models
   - Handle periodic boundaries
   - ~600-900 lines

6. **ion** (Medium Priority)
   - Calculate electrostatic potential
   - Place counter ions
   - Replace water molecules
   - ~400-600 lines

### Phase 3: Main MD Binary (1 week)
**Goal**: Create executable MD engine

7. **md** (Critical)
   - Create main MD driver using gromos-rs library
   - Read IMD input files
   - Run MD loop (integrate, calculate forces, I/O)
   - Handle trajectories (.trc, .tre, .trf)
   - ~800-1,200 lines

### Phase 4: Analysis Tools (2-3 weeks)
**Goal**: Analyze MD results

8. **ene_ana** (High Priority)
   - Read energy trajectory (.tre files)
   - Calculate time averages
   - Generate time series
   - ~600-800 lines

9. **rmsd** (High Priority)
   - Read coordinate trajectories
   - Align structures
   - Calculate RMSD
   - ~500-700 lines

10. **rmsf** (Medium Priority)
    - Calculate atomic fluctuations
    - Per-residue averages
    - ~400-600 lines

---

## Technical Requirements

### I/O Formats Needed

**Already Implemented:**
- ✅ IMD (Input MD parameters) - read/write in gromos-rs/src/io/

**Need to Implement:**
| Format | Extension | Purpose | Status |
|--------|-----------|---------|--------|
| **Topology** | `.top` | System structure | ⚠️ Partial (can read, need write) |
| **Coordinates** | `.cnf, .g96, .pdb` | Positions/velocities | ⚠️ Partial (can read CNF) |
| **Trajectory** | `.trc` | Coordinate time series | ❌ Not implemented |
| **Energy** | `.tre` | Energy time series | ❌ Not implemented |
| **Forces** | `.trf` | Force output | ❌ Not implemented |
| **Building Blocks** | `.mtb, .ifp` | Molecular building blocks | ❌ Not implemented |

### Dependencies

**Rust Crates:**
- `nom` or `pest` - Parsing topology/input files
- `nalgebra` - Matrix operations for RMSD
- `plotters` (optional) - Plotting for analysis tools

**Force Field Data:**
- 54A7 force field parameters (already in gromosPlusPlus)
- Building block definitions (.mtb files)
- Atom type definitions

---

## Validation Strategy

### Test Against Tutorial Data

1. **Download tutorial files**:
   ```bash
   git submodule add https://github.com/biomos/gromos_tutorial_livecoms.git
   ```

2. **Run tutorial_01 workflow**:
   - Convert PDB → G96
   - Generate topology
   - Solvate system
   - Run minimization + equilibration
   - Analyze results

3. **Compare outputs**:
   - Topology files (atom counts, bonds, parameters)
   - Energy values (bonded, nonbonded, total)
   - RMSD values
   - Trajectory statistics

### Success Criteria

- [ ] Successfully complete Tutorial 01 workflow
- [ ] Energy values match GROMOS++ within 0.1%
- [ ] RMSD calculations match within 0.01 nm
- [ ] Topologies validate correctly

---

## Timeline Estimates

| Phase | Duration | Deliverables |
|-------|----------|--------------|
| **Phase 1** | 2-3 weeks | pdb2g96, make_top, com_top, check_top |
| **Phase 2** | 1-2 weeks | sim_box, ion |
| **Phase 3** | 1 week | md binary |
| **Phase 4** | 2-3 weeks | ene_ana, rmsd, rmsf |
| **Testing** | 1 week | Tutorial validation |

**Total**: ~7-10 weeks for full tutorial compatibility

---

## Current Status

**Tier 1 MD Engine**: ✅ 95% complete
**Tier 2 Enhanced Methods**: ✅ 71% complete
**Tutorial Tools**: ✅ 25% complete (4/16 programs)

**Recent Additions** (Session 2):
- ✅ pdb2g96: PDB → GROMOS96 converter (130 lines + I/O modules)
- ✅ com_top: Topology file combiner (200 lines)
- ✅ check_top: Topology validator with 34+ checks (540 lines)
- ✅ Topology writing support (175 lines)
- ✅ PDB and G96 file I/O modules (440 lines)

**Next Steps**:
1. ⏳ Defer make_top (requires .mtb building block parser - complex)
2. Implement md binary (main driver) - **HIGH PRIORITY**
3. Implement sim_box (solvation) - needed for Tutorial 01
4. Implement ene_ana (energy analysis) - for validation
5. Test with Tutorial 01

---

## Notes

- Many analysis tools (rmsd, rmsf, nhoparam) can be implemented in parallel
- System preparation tools (pdb2g96, make_top) are critical path items
- The md binary is relatively simple once the library is feature-complete
- gromos-rs library already has most MD capabilities; tools are mostly I/O wrappers

**Key Insight**: Most "programs" are thin wrappers around the gromos-rs library. The hard work (MD engine, force calculations) is already done. We mainly need:
1. File format parsers (topology, PDB, building blocks)
2. Trajectory I/O (.trc, .tre, .trf)
3. Command-line interfaces for each tool
