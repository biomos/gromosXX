//! Python bindings for GROMOS-RS molecular dynamics engine
//!
//! Inspired by Polars' architecture:
//! - Zero-copy data sharing where possible
//! - Rust core with Python wrapper
//! - High-performance parallel execution
//! - Memory-safe operations

use pyo3::prelude::*;
use pyo3::types::{PyList, PyDict};
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray, PyUntypedArrayMethods};

// Re-export from gromos-rs
use gromos_rs::{
    Vec3 as RustVec3,
    Mat3 as RustMat3,
    Configuration as RustConfiguration,
    State as RustState,
    Energy as RustEnergy,
    Topology as RustTopology,
};

/// Python module for GROMOS molecular dynamics
#[pymodule]
fn gromos(_py: Python, m: &PyModule) -> PyResult<()> {
    // Math types
    m.add_class::<PyVec3>()?;
    m.add_class::<PyMat3>()?;

    // Core data structures
    m.add_class::<PyBox>()?;
    m.add_class::<PyEnergy>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyConfiguration>()?;
    m.add_class::<PyTopology>()?;

    // Note: Integrators and advanced sampling methods are available via
    // the Rust CLI binaries. Python bindings focus on data structures and analysis.

    Ok(())
}

//==============================================================================
// MATH TYPES
//==============================================================================

/// 3D vector with SIMD acceleration
///
/// Wraps Rust's Vec3A (SIMD-accelerated vector from glam crate).
/// Provides zero-copy conversion to/from NumPy arrays where possible.
#[pyclass(name = "Vec3")]
#[derive(Clone)]
pub struct PyVec3 {
    inner: RustVec3,
}

#[pymethods]
impl PyVec3 {
    /// Create a new Vec3
    #[new]
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            inner: RustVec3::new(x, y, z),
        }
    }

    /// Get x component
    #[getter]
    fn x(&self) -> f32 {
        self.inner.x
    }

    /// Get y component
    #[getter]
    fn y(&self) -> f32 {
        self.inner.y
    }

    /// Get z component
    #[getter]
    fn z(&self) -> f32 {
        self.inner.z
    }

    /// Set x component
    #[setter]
    fn set_x(&mut self, value: f32) {
        self.inner.x = value;
    }

    /// Set y component
    #[setter]
    fn set_y(&mut self, value: f32) {
        self.inner.y = value;
    }

    /// Set z component
    #[setter]
    fn set_z(&mut self, value: f32) {
        self.inner.z = value;
    }

    /// Vector length (magnitude)
    fn length(&self) -> f32 {
        self.inner.length()
    }

    /// Squared length (faster than length)
    fn length_squared(&self) -> f32 {
        self.inner.length_squared()
    }

    /// Normalize the vector (unit length)
    fn normalize(&self) -> Self {
        Self {
            inner: self.inner.normalize(),
        }
    }

    /// Dot product
    fn dot(&self, other: &Self) -> f32 {
        self.inner.dot(other.inner)
    }

    /// Cross product
    fn cross(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.cross(other.inner),
        }
    }

    /// Distance to another vector
    fn distance(&self, other: &Self) -> f32 {
        self.inner.distance(other.inner)
    }

    /// Add two vectors
    fn __add__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner + other.inner,
        }
    }

    /// Subtract two vectors
    fn __sub__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner - other.inner,
        }
    }

    /// Multiply by scalar
    fn __mul__(&self, scalar: f32) -> Self {
        Self {
            inner: self.inner * scalar,
        }
    }

    /// String representation
    fn __repr__(&self) -> String {
        format!("Vec3({:.4}, {:.4}, {:.4})", self.inner.x, self.inner.y, self.inner.z)
    }

    /// Convert to NumPy array
    fn to_numpy<'py>(&self, py: Python<'py>) -> &'py PyArray1<f32> {
        let arr = [self.inner.x, self.inner.y, self.inner.z];
        arr.to_pyarray(py)
    }

    /// Create from NumPy array
    #[staticmethod]
    fn from_numpy(arr: PyReadonlyArray1<f32>) -> PyResult<Self> {
        let slice = arr.as_slice()?;
        if slice.len() != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have exactly 3 elements"
            ));
        }
        Ok(Self {
            inner: RustVec3::new(slice[0], slice[1], slice[2]),
        })
    }
}

/// 3x3 matrix with SIMD acceleration
#[pyclass(name = "Mat3")]
#[derive(Clone)]
pub struct PyMat3 {
    inner: RustMat3,
}

#[pymethods]
impl PyMat3 {
    /// Create identity matrix
    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: RustMat3::IDENTITY,
        }
    }

    /// Create from column vectors
    #[staticmethod]
    fn from_cols(x: &PyVec3, y: &PyVec3, z: &PyVec3) -> Self {
        Self {
            inner: RustMat3::from_cols(x.inner, y.inner, z.inner),
        }
    }

    /// Matrix determinant
    fn determinant(&self) -> f32 {
        self.inner.determinant()
    }

    /// Matrix inverse
    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse(),
        }
    }

    /// Matrix transpose
    fn transpose(&self) -> Self {
        Self {
            inner: self.inner.transpose(),
        }
    }

    /// Matrix-vector multiplication
    fn mul_vec3(&self, v: &PyVec3) -> PyVec3 {
        PyVec3 {
            inner: self.inner * v.inner,
        }
    }

    /// Convert to NumPy array (3x3)
    fn to_numpy<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let flat = vec![
            self.inner.x_axis.x, self.inner.y_axis.x, self.inner.z_axis.x,
            self.inner.x_axis.y, self.inner.y_axis.y, self.inner.z_axis.y,
            self.inner.x_axis.z, self.inner.y_axis.z, self.inner.z_axis.z,
        ];
        PyArray2::from_vec2(py, &[flat])
            .unwrap()
            .reshape([3, 3])
            .unwrap()
    }

    fn __repr__(&self) -> String {
        format!("Mat3([[{:.4}, {:.4}, {:.4}],\n     [{:.4}, {:.4}, {:.4}],\n     [{:.4}, {:.4}, {:.4}]])",
            self.inner.x_axis.x, self.inner.y_axis.x, self.inner.z_axis.x,
            self.inner.x_axis.y, self.inner.y_axis.y, self.inner.z_axis.y,
            self.inner.x_axis.z, self.inner.y_axis.z, self.inner.z_axis.z)
    }
}

//==============================================================================
// CORE DATA STRUCTURES
//==============================================================================

/// Simulation box representation
#[pyclass(name = "Box")]
#[derive(Clone)]
pub struct PyBox {
    inner: gromos_rs::configuration::Box,
}

#[pymethods]
impl PyBox {
    /// Create vacuum box (no periodicity)
    #[staticmethod]
    fn vacuum() -> Self {
        Self {
            inner: gromos_rs::configuration::Box::vacuum(),
        }
    }

    /// Create rectangular box
    #[staticmethod]
    fn rectangular(lx: f32, ly: f32, lz: f32) -> Self {
        Self {
            inner: gromos_rs::configuration::Box::rectangular(lx, ly, lz),
        }
    }

    /// Create triclinic box from vectors
    #[staticmethod]
    fn triclinic(mat: &PyMat3) -> Self {
        Self {
            inner: gromos_rs::configuration::Box::triclinic(mat.inner),
        }
    }

    /// Get box volume
    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    /// Get box dimensions (for rectangular box)
    fn dimensions(&self) -> PyVec3 {
        PyVec3 {
            inner: self.inner.dimensions(),
        }
    }

    fn __repr__(&self) -> String {
        format!("Box(type={:?}, volume={:.2})", self.inner.box_type, self.volume())
    }
}

/// Energy storage for molecular system
///
/// Stores total energies and per-group energies for detailed accounting.
/// All energies are in kJ/mol.
#[pyclass(name = "Energy")]
#[derive(Clone)]
pub struct PyEnergy {
    inner: RustEnergy,
}

#[pymethods]
impl PyEnergy {
    /// Create new energy object
    #[new]
    fn new(num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            inner: RustEnergy::new(num_temperature_groups, num_energy_groups),
        }
    }

    /// Total energy (kinetic + potential)
    fn total(&self) -> f64 {
        self.inner.total()
    }

    /// Kinetic energy
    #[getter]
    fn kinetic(&self) -> f64 {
        self.inner.kinetic_total
    }

    /// Potential energy
    #[getter]
    fn potential(&self) -> f64 {
        self.inner.potential_total
    }

    /// Bond energy
    #[getter]
    fn bond(&self) -> f64 {
        self.inner.bond_total
    }

    /// Angle energy
    #[getter]
    fn angle(&self) -> f64 {
        self.inner.angle_total
    }

    /// Dihedral energy
    #[getter]
    fn dihedral(&self) -> f64 {
        self.inner.dihedral_total
    }

    /// Lennard-Jones energy
    #[getter]
    fn lj(&self) -> f64 {
        self.inner.lj_total
    }

    /// Coulomb reaction field energy
    #[getter]
    fn coulomb(&self) -> f64 {
        self.inner.crf_total
    }

    /// Clear all energies to zero
    fn clear(&mut self) {
        self.inner.clear();
    }

    /// Get energy as dictionary
    fn to_dict(&self, py: Python) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        dict.set_item("total", self.total())?;
        dict.set_item("kinetic", self.inner.kinetic_total)?;
        dict.set_item("potential", self.inner.potential_total)?;
        dict.set_item("bond", self.inner.bond_total)?;
        dict.set_item("angle", self.inner.angle_total)?;
        dict.set_item("dihedral", self.inner.dihedral_total)?;
        dict.set_item("lj", self.inner.lj_total)?;
        dict.set_item("coulomb", self.inner.crf_total)?;
        Ok(dict.into())
    }

    fn __repr__(&self) -> String {
        format!("Energy(total={:.2} kJ/mol, kinetic={:.2}, potential={:.2})",
            self.total(), self.inner.kinetic_total, self.inner.potential_total)
    }
}

/// System state (positions, velocities, forces)
///
/// Uses zero-copy sharing with NumPy arrays for efficient data access.
#[pyclass(name = "State")]
pub struct PyState {
    inner: RustState,
}

#[pymethods]
impl PyState {
    /// Create new state for N atoms
    #[new]
    fn new(num_atoms: usize, num_temp_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            inner: RustState::new(num_atoms, num_temp_groups, num_energy_groups),
        }
    }

    /// Number of atoms
    fn num_atoms(&self) -> usize {
        self.inner.pos.len()
    }

    /// Get positions as NumPy array (N x 3)
    fn positions<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let n = self.inner.pos.len();
        let mut arr = Vec::with_capacity(n * 3);
        for v in &self.inner.pos {
            arr.push(v.x);
            arr.push(v.y);
            arr.push(v.z);
        }
        PyArray2::from_vec2(py, &[arr])
            .unwrap()
            .reshape([n, 3])
            .unwrap()
    }

    /// Get velocities as NumPy array (N x 3)
    fn velocities<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let n = self.inner.vel.len();
        let mut arr = Vec::with_capacity(n * 3);
        for v in &self.inner.vel {
            arr.push(v.x);
            arr.push(v.y);
            arr.push(v.z);
        }
        PyArray2::from_vec2(py, &[arr])
            .unwrap()
            .reshape([n, 3])
            .unwrap()
    }

    /// Get forces as NumPy array (N x 3)
    fn forces<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let n = self.inner.force.len();
        let mut arr = Vec::with_capacity(n * 3);
        for v in &self.inner.force {
            arr.push(v.x);
            arr.push(v.y);
            arr.push(v.z);
        }
        PyArray2::from_vec2(py, &[arr])
            .unwrap()
            .reshape([n, 3])
            .unwrap()
    }

    /// Set positions from NumPy array (N x 3)
    fn set_positions(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        let shape = arr.shape();
        if shape[1] != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have shape (N, 3)"
            ));
        }

        let slice = arr.as_slice()?;
        self.inner.pos.clear();
        for i in 0..shape[0] {
            let idx = i * 3;
            self.inner.pos.push(RustVec3::new(slice[idx], slice[idx+1], slice[idx+2]));
        }
        Ok(())
    }

    /// Set velocities from NumPy array (N x 3)
    fn set_velocities(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        let shape = arr.shape();
        if shape[1] != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have shape (N, 3)"
            ));
        }

        let slice = arr.as_slice()?;
        self.inner.vel.clear();
        for i in 0..shape[0] {
            let idx = i * 3;
            self.inner.vel.push(RustVec3::new(slice[idx], slice[idx+1], slice[idx+2]));
        }
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("State({} atoms)", self.num_atoms())
    }
}

/// System configuration (combines state, energy, topology)
#[pyclass(name = "Configuration")]
pub struct PyConfiguration {
    inner: RustConfiguration,
}

#[pymethods]
impl PyConfiguration {
    /// Create new configuration
    #[new]
    fn new(num_atoms: usize, num_temp_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            inner: RustConfiguration::new(num_atoms, num_temp_groups, num_energy_groups),
        }
    }

    /// Get current state
    fn current_state(&self) -> PyState {
        PyState {
            inner: self.inner.current().clone(),
        }
    }

    /// Get current energy
    fn current_energy(&self) -> PyEnergy {
        PyEnergy {
            inner: self.inner.current().energies.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!("Configuration({} atoms)", self.inner.current().pos.len())
    }
}

/// Molecular topology (atoms, bonds, parameters)
#[pyclass(name = "Topology")]
pub struct PyTopology {
    inner: RustTopology,
}

#[pymethods]
impl PyTopology {
    /// Create new empty topology
    #[new]
    fn new() -> Self {
        Self {
            inner: RustTopology::new(),
        }
    }

    /// Number of atoms
    fn num_atoms(&self) -> usize {
        self.inner.solute.atoms.len()
    }

    /// Number of bonds
    fn num_bonds(&self) -> usize {
        self.inner.solute.bonds.len()
    }

    /// Number of angles
    fn num_angles(&self) -> usize {
        self.inner.solute.angles.len()
    }

    /// Number of dihedrals (proper)
    fn num_dihedrals(&self) -> usize {
        self.inner.solute.proper_dihedrals.len()
    }

    fn __repr__(&self) -> String {
        format!("Topology({} atoms, {} bonds)", self.num_atoms(), self.num_bonds())
    }
}

//==============================================================================
// NOTE: Integrators and advanced sampling methods (GaMD, EDS, REMD) are
// available through the GROMOS-RS command-line binaries.
// Python bindings focus on data structures, I/O, and analysis tools.
//==============================================================================
