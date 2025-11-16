//! Pairlist (neighbor list) generation
//!
//! Direct Rust translation of GROMOS pairlist algorithms from:
//! - md++/src/interaction/nonbonded/pairlist/pairlist.h
//! - md++/src/interaction/nonbonded/pairlist/standard_pairlist_algorithm.cc
//! - md++/src/interaction/nonbonded/pairlist/grid_cell_pairlist.cc

use crate::math::{Vec3, BoundaryCondition};
use crate::topology::Topology;
use crate::configuration::Configuration;
use rayon::prelude::*;

/// Pairlist for nonbonded interactions
///
/// Stores atom pairs (i, j) where i < j that are within cutoff distance.
/// Updated periodically (typically every 5-10 MD steps).
pub type Pairlist = Vec<(usize, usize)>;

/// Complete pairlist container
///
/// Separates short-range and long-range interactions,
/// and solute vs. solvent for optimization.
#[derive(Debug, Clone)]
pub struct PairlistContainer {
    pub solute_short: Pairlist,
    pub solute_long: Pairlist,
    pub solvent_short: Pairlist,
    pub solvent_long: Pairlist,

    // Pairlist parameters
    pub short_range_cutoff: f64,
    pub long_range_cutoff: f64,
    pub skin: f64,  // Extra distance for pairlist stability

    // Tracking for update frequency
    pub steps_since_update: usize,
    pub update_frequency: usize,
}

impl PairlistContainer {
    pub fn new(short_cutoff: f64, long_cutoff: f64, skin: f64) -> Self {
        Self {
            solute_short: Vec::new(),
            solute_long: Vec::new(),
            solvent_short: Vec::new(),
            solvent_long: Vec::new(),
            short_range_cutoff: short_cutoff,
            long_range_cutoff: long_cutoff,
            skin,
            steps_since_update: 0,
            update_frequency: 5,  // Default: update every 5 steps
        }
    }

    /// Check if pairlist needs updating
    pub fn needs_update(&self) -> bool {
        self.steps_since_update >= self.update_frequency
    }

    /// Increment step counter
    pub fn step(&mut self) {
        self.steps_since_update += 1;
    }

    /// Reset step counter after update
    pub fn reset_counter(&mut self) {
        self.steps_since_update = 0;
    }

    /// Total number of pairs
    pub fn total_pairs(&self) -> usize {
        self.solute_short.len() +
        self.solute_long.len() +
        self.solvent_short.len() +
        self.solvent_long.len()
    }
}

/// Standard pairlist algorithm using chargegroups
///
/// Translation of md++/src/interaction/nonbonded/pairlist/standard_pairlist_algorithm.cc
pub struct StandardPairlistAlgorithm {
    pub use_chargegroups: bool,
}

impl StandardPairlistAlgorithm {
    pub fn new(use_chargegroups: bool) -> Self {
        Self { use_chargegroups }
    }

    /// Update pairlist based on current configuration
    ///
    /// Translated from standard_pairlist_algorithm.cc:168-250 (chargegroup version)
    pub fn update<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        pairlist.solute_short.clear();
        pairlist.solute_long.clear();

        if self.use_chargegroups {
            self.update_chargegroup_based(topo, conf, pairlist, periodicity);
        } else {
            self.update_atom_based(topo, conf, pairlist, periodicity);
        }

        pairlist.reset_counter();
    }

    /// Chargegroup-based pairlist generation (GROMOS default)
    ///
    /// Direct translation of lines 168-250 in standard_pairlist_algorithm.cc
    fn update_chargegroup_based<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        let n_chargegroups = topo.chargegroups.len();
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        // Calculate chargegroup centers
        let cg_centers: Vec<Vec3> = topo.chargegroups.iter()
            .map(|cg| cg.center_of_geometry(&conf.current().pos))
            .collect();

        // Double loop over chargegroups
        for i in 0..n_chargegroups {
            for j in (i + 1)..n_chargegroups {
                // Calculate distance between chargegroup centers
                let r = periodicity.nearest_image(cg_centers[i], cg_centers[j]);
                let dist2 = r.length_squared() as f64;

                // Check if within cutoff
                let in_short = dist2 < cutoff2_short;
                let in_long = dist2 < cutoff2_long;

                if !in_short && !in_long {
                    continue;  // Skip if beyond both cutoffs
                }

                // Add all atom pairs from these chargegroups
                for &atom_i in &topo.chargegroups[i].atoms {
                    for &atom_j in &topo.chargegroups[j].atoms {
                        // Check exclusions
                        if topo.is_excluded(atom_i, atom_j) {
                            continue;
                        }

                        // Add to appropriate pairlist
                        if in_short {
                            pairlist.solute_short.push((atom_i, atom_j));
                        } else if in_long {
                            pairlist.solute_long.push((atom_i, atom_j));
                        }
                    }
                }
            }
        }
    }

    /// Atom-based pairlist generation (simpler, no chargegroups)
    fn update_atom_based<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        let n_atoms = topo.num_atoms();
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        for i in 0..n_atoms {
            for j in (i + 1)..n_atoms {
                // Check exclusions
                if topo.is_excluded(i, j) {
                    continue;
                }

                // Calculate distance
                let r = periodicity.nearest_image(
                    conf.current().pos[i],
                    conf.current().pos[j],
                );
                let dist2 = r.length_squared() as f64;

                // Add to appropriate pairlist
                if dist2 < cutoff2_short {
                    pairlist.solute_short.push((i, j));
                } else if dist2 < cutoff2_long {
                    pairlist.solute_long.push((i, j));
                }
            }
        }
    }
}

/// Grid-based pairlist algorithm for large systems
///
/// Translation of md++/src/interaction/nonbonded/pairlist/grid_cell_pairlist.cc
///
/// Uses spatial decomposition into grid cells for O(N) scaling instead of O(N²).
pub struct GridCellPairlistAlgorithm {
    pub cell_size: f32,
    cells: Vec<Vec<usize>>,  // cells[cell_idx] = list of atoms
    grid_dim: [usize; 3],
}

impl GridCellPairlistAlgorithm {
    /// Create grid-based pairlist algorithm
    ///
    /// `cell_size` should be at least cutoff + skin
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cells: Vec::new(),
            grid_dim: [0, 0, 0],
        }
    }

    /// Update pairlist using grid-based algorithm
    pub fn update<BC: BoundaryCondition>(
        &mut self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        pairlist.solute_short.clear();
        pairlist.solute_long.clear();

        // Initialize grid
        self.initialize_grid(conf);

        // Assign atoms to cells
        self.assign_atoms_to_cells(topo, conf);

        // Generate pairs from grid
        self.generate_pairs_from_grid(topo, conf, pairlist, periodicity);

        pairlist.reset_counter();
    }

    /// Initialize grid based on box dimensions
    fn initialize_grid(&mut self, conf: &Configuration) {
        let box_dims = conf.current().box_config.dimensions();

        self.grid_dim = [
            ((box_dims.x / self.cell_size).ceil() as usize).max(1),
            ((box_dims.y / self.cell_size).ceil() as usize).max(1),
            ((box_dims.z / self.cell_size).ceil() as usize).max(1),
        ];

        let total_cells = self.grid_dim[0] * self.grid_dim[1] * self.grid_dim[2];
        self.cells = vec![Vec::new(); total_cells];

        // Clear existing cells
        for cell in &mut self.cells {
            cell.clear();
        }
    }

    /// Assign atoms to grid cells
    fn assign_atoms_to_cells(&mut self, topo: &Topology, conf: &Configuration) {
        let n_atoms = topo.num_atoms();

        for i in 0..n_atoms {
            let pos = conf.current().pos[i];
            let cell_idx = self.position_to_cell(pos);
            self.cells[cell_idx].push(i);
        }
    }

    /// Convert position to cell index
    fn position_to_cell(&self, pos: Vec3) -> usize {
        let ix = ((pos.x / self.cell_size).floor() as usize).min(self.grid_dim[0] - 1);
        let iy = ((pos.y / self.cell_size).floor() as usize).min(self.grid_dim[1] - 1);
        let iz = ((pos.z / self.cell_size).floor() as usize).min(self.grid_dim[2] - 1);

        ix + iy * self.grid_dim[0] + iz * self.grid_dim[0] * self.grid_dim[1]
    }

    /// Generate pairs from grid cells
    fn generate_pairs_from_grid<BC: BoundaryCondition>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) {
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        // Loop over all cells
        for cell_idx in 0..self.cells.len() {
            let (cx, cy, cz) = self.cell_index_to_coords(cell_idx);

            // Check neighboring cells (including self)
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let nx = (cx as isize + dx).rem_euclid(self.grid_dim[0] as isize) as usize;
                        let ny = (cy as isize + dy).rem_euclid(self.grid_dim[1] as isize) as usize;
                        let nz = (cz as isize + dz).rem_euclid(self.grid_dim[2] as isize) as usize;

                        let neighbor_idx = nx + ny * self.grid_dim[0] + nz * self.grid_dim[0] * self.grid_dim[1];

                        // Add pairs between current cell and neighbor
                        self.add_cell_pairs(
                            cell_idx,
                            neighbor_idx,
                            topo,
                            conf,
                            pairlist,
                            periodicity,
                            cutoff2_short,
                            cutoff2_long,
                        );
                    }
                }
            }
        }
    }

    /// Convert 1D cell index to 3D coordinates
    fn cell_index_to_coords(&self, idx: usize) -> (usize, usize, usize) {
        let z = idx / (self.grid_dim[0] * self.grid_dim[1]);
        let remainder = idx % (self.grid_dim[0] * self.grid_dim[1]);
        let y = remainder / self.grid_dim[0];
        let x = remainder % self.grid_dim[0];
        (x, y, z)
    }

    /// Add pairs between two cells
    fn add_cell_pairs<BC: BoundaryCondition>(
        &self,
        cell1: usize,
        cell2: usize,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
        cutoff2_short: f64,
        cutoff2_long: f64,
    ) {
        for &i in &self.cells[cell1] {
            for &j in &self.cells[cell2] {
                // Skip if same cell and j <= i (avoid double counting)
                if cell1 == cell2 && j <= i {
                    continue;
                }

                // Check exclusions
                if topo.is_excluded(i, j) {
                    continue;
                }

                // Calculate distance
                let r = periodicity.nearest_image(
                    conf.current().pos[i],
                    conf.current().pos[j],
                );
                let dist2 = r.length_squared() as f64;

                // Add to appropriate pairlist
                if dist2 < cutoff2_short {
                    pairlist.solute_short.push((i, j));
                } else if dist2 < cutoff2_long {
                    pairlist.solute_long.push((i, j));
                }
            }
        }
    }
}

/// Parallel pairlist generation using Rayon
pub struct ParallelPairlistAlgorithm {
    base_algorithm: StandardPairlistAlgorithm,
}

impl ParallelPairlistAlgorithm {
    pub fn new(use_chargegroups: bool) -> Self {
        Self {
            base_algorithm: StandardPairlistAlgorithm::new(use_chargegroups),
        }
    }

    /// Update pairlist in parallel
    pub fn update<BC: BoundaryCondition + Sync>(
        &self,
        topo: &Topology,
        conf: &Configuration,
        pairlist: &mut PairlistContainer,
        periodicity: &BC,
    ) where
        BC: Clone,
    {
        let n_atoms = topo.num_atoms();
        let cutoff2_short = (pairlist.short_range_cutoff + pairlist.skin).powi(2);
        let cutoff2_long = (pairlist.long_range_cutoff + pairlist.skin).powi(2);

        // Parallel generation of pairs
        let pairs: Vec<(usize, usize, bool)> = (0..n_atoms)
            .into_par_iter()
            .flat_map(|i| {
                let mut local_pairs = Vec::new();
                let periodicity = periodicity.clone();

                for j in (i + 1)..n_atoms {
                    if topo.is_excluded(i, j) {
                        continue;
                    }

                    let r = periodicity.nearest_image(
                        conf.current().pos[i],
                        conf.current().pos[j],
                    );
                    let dist2 = r.length_squared() as f64;

                    if dist2 < cutoff2_short {
                        local_pairs.push((i, j, true));  // true = short range
                    } else if dist2 < cutoff2_long {
                        local_pairs.push((i, j, false)); // false = long range
                    }
                }

                local_pairs
            })
            .collect();

        // Separate into short and long range
        pairlist.solute_short.clear();
        pairlist.solute_long.clear();

        for (i, j, is_short) in pairs {
            if is_short {
                pairlist.solute_short.push((i, j));
            } else {
                pairlist.solute_long.push((i, j));
            }
        }

        pairlist.reset_counter();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Rectangular;
    use crate::topology::ChargeGroup;

    #[test]
    fn test_pairlist_container() {
        let pairlist = PairlistContainer::new(1.0, 1.4, 0.2);

        assert_eq!(pairlist.short_range_cutoff, 1.0);
        assert_eq!(pairlist.long_range_cutoff, 1.4);
        assert!(!pairlist.needs_update());
    }

    #[test]
    fn test_pairlist_update_frequency() {
        let mut pairlist = PairlistContainer::new(1.0, 1.4, 0.2);
        pairlist.update_frequency = 3;

        pairlist.step();
        assert!(!pairlist.needs_update());

        pairlist.step();
        assert!(!pairlist.needs_update());

        pairlist.step();
        assert!(pairlist.needs_update());

        pairlist.reset_counter();
        assert!(!pairlist.needs_update());
    }

    #[test]
    fn test_standard_pairlist() {
        let mut topo = Topology::new();
        topo.iac = vec![0; 4];
        topo.mass = vec![1.0; 4];
        topo.charge = vec![0.0; 4];
        topo.resize_atom_arrays();

        let mut conf = Configuration::new(4, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(0.0, 0.0, 0.0);
        conf.current_mut().pos[1] = Vec3::new(0.5, 0.0, 0.0);  // Close
        conf.current_mut().pos[2] = Vec3::new(2.0, 0.0, 0.0);  // Far
        conf.current_mut().pos[3] = Vec3::new(0.0, 0.5, 0.0);  // Close

        let mut pairlist = PairlistContainer::new(1.0, 2.5, 0.2);
        let algorithm = StandardPairlistAlgorithm::new(false);
        let periodicity = crate::math::Vacuum;

        algorithm.update(&topo, &conf, &mut pairlist, &periodicity);

        // Should have pairs (0,1), (0,3), and (1,3) in short range
        // and (0,2) in long range
        assert!(pairlist.solute_short.len() >= 2);
        assert!(pairlist.solute_long.len() >= 1);
    }

    #[test]
    fn test_grid_cell_algorithm() {
        use crate::topology::Atom;

        let mut topo = Topology::new();

        // Add 10 solute atoms
        for i in 0..10 {
            topo.solute.atoms.push(Atom {
                name: format!("C{}", i),
                residue_nr: 1,
                residue_name: "MOL".to_string(),
                iac: 0,
                mass: 1.0,
                charge: 0.0,
                is_perturbed: false,
                is_polarisable: false,
                is_coarse_grained: false,
            });
        }

        topo.iac = vec![0; 10];
        topo.mass = vec![1.0; 10];
        topo.charge = vec![0.0; 10];
        topo.resize_atom_arrays();

        let mut conf = Configuration::new(10, 1, 1);
        conf.current_mut().box_config = crate::configuration::Box::rectangular(5.0, 5.0, 5.0);

        // Place atoms randomly in box
        for (i, pos) in conf.current_mut().pos.iter_mut().enumerate() {
            *pos = Vec3::new(
                (i as f32 % 5.0),
                ((i / 5) as f32 % 5.0),
                0.0,
            );
        }

        let mut pairlist = PairlistContainer::new(1.5, 2.5, 0.2);
        let mut algorithm = GridCellPairlistAlgorithm::new(2.0);
        let periodicity = Rectangular::new(Vec3::new(5.0, 5.0, 5.0));

        algorithm.update(&topo, &conf, &mut pairlist, &periodicity);

        // Should find some pairs
        assert!(pairlist.total_pairs() > 0);
    }
}
