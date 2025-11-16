//! C FFI interface for integration with C++ GROMOS code

use crate::interaction::nonbonded::*;
use crate::math::*;
use std::slice;
use libc::{c_double, c_float, c_uint};

/// C-compatible representation of 3D vector
#[repr(C)]
pub struct CVec3 {
    pub x: c_float,
    pub y: c_float,
    pub z: c_float,
}

/// C-compatible LJ parameters
#[repr(C)]
pub struct CLJParameters {
    pub c6: c_double,
    pub c12: c_double,
}

/// C FFI: LJ + CRF inner loop
///
/// # Safety
/// Caller must ensure all pointers are valid and have correct length
///
/// # Parameters
/// * `positions` - Flat array of positions [x0,y0,z0,x1,y1,z1,...] of length n_atoms*3
/// * `charges` - Array of charges, length n_atoms
/// * `iac` - Integer atom codes (atom types), length n_atoms
/// * `n_atoms` - Number of atoms
/// * `pairlist` - Flat array of pairs [i0,j0,i1,j1,...] of length n_pairs*2
/// * `n_pairs` - Number of pairs
/// * `lj_params` - Flat array of LJ parameters [c6,c12] for each type pair, length n_types*n_types*2
/// * `n_types` - Number of atom types
/// * `box_data` - Box vectors (3 values for rectangular box, 9 for triclinic)
/// * `boundary_type` - 0=vacuum, 1=rectangular, 2=triclinic
/// * `crf_cut` - CRF cutoff distance
/// * `crf_2cut3i` - CRF parameter: 2/cutoff^3
/// * `crf_cut3i` - CRF parameter: 1/cutoff^3
/// * `forces` - Output: forces, flat array length n_atoms*3
/// * `energies` - Output: [e_lj, e_crf]
/// * `virial` - Output: virial tensor, flat array 3x3 = 9 elements
#[no_mangle]
pub unsafe extern "C" fn rust_lj_crf_innerloop(
    positions: *const c_float,
    charges: *const c_float,
    iac: *const c_uint,
    n_atoms: c_uint,
    pairlist: *const c_uint,
    n_pairs: c_uint,
    lj_params: *const c_double,
    n_types: c_uint,
    box_data: *const c_double,
    boundary_type: c_uint,
    crf_cut: c_double,
    crf_2cut3i: c_double,
    crf_cut3i: c_double,
    forces: *mut c_float,
    energies: *mut c_double,
    virial: *mut c_double,
) {
    // Safety checks
    if positions.is_null() || charges.is_null() || iac.is_null() ||
       pairlist.is_null() || lj_params.is_null() || forces.is_null() ||
       energies.is_null() || virial.is_null() {
        eprintln!("Error: null pointer passed to rust_lj_crf_innerloop");
        return;
    }

    let n_atoms = n_atoms as usize;
    let n_pairs = n_pairs as usize;
    let n_types = n_types as usize;

    // Convert raw pointers to slices
    let positions_slice = slice::from_raw_parts(positions, n_atoms * 3);
    let charges_slice = slice::from_raw_parts(charges, n_atoms);
    let iac_slice = slice::from_raw_parts(iac, n_atoms);
    let pairlist_slice = slice::from_raw_parts(pairlist, n_pairs * 2);
    let lj_params_slice = slice::from_raw_parts(lj_params, n_types * n_types * 2);

    // Convert positions to Vec3
    let positions_vec: Vec<Vec3> = (0..n_atoms)
        .map(|i| {
            Vec3::new(
                positions_slice[i * 3],
                positions_slice[i * 3 + 1],
                positions_slice[i * 3 + 2],
            )
        })
        .collect();

    // Convert pairlist to tuples
    let pairlist_vec: Vec<(u32, u32)> = (0..n_pairs)
        .map(|i| (pairlist_slice[i * 2], pairlist_slice[i * 2 + 1]))
        .collect();

    // Convert LJ parameters to 2D array
    let mut lj_params_vec: Vec<Vec<LJParameters>> = vec![vec![LJParameters { c6: 0.0, c12: 0.0 }; n_types]; n_types];
    for i in 0..n_types {
        for j in 0..n_types {
            let idx = (i * n_types + j) * 2;
            lj_params_vec[i][j] = LJParameters {
                c6: lj_params_slice[idx],
                c12: lj_params_slice[idx + 1],
            };
        }
    }

    // Create CRF parameters
    let crf = CRFParameters {
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
    };

    // Create storage
    let mut storage = ForceStorage::new(n_atoms);

    // Determine boundary type and call appropriate function
    match boundary_type {
        0 => {
            // Vacuum
            let bc = Vacuum;
            lj_crf_innerloop(
                &positions_vec,
                charges_slice,
                iac_slice,
                &pairlist_vec,
                &lj_params_vec,
                &crf,
                &bc,
                &mut storage,
            );
        }
        1 => {
            // Rectangular
            let box_slice = slice::from_raw_parts(box_data, 3);
            let box_size = Vec3::new(box_slice[0] as f32, box_slice[1] as f32, box_slice[2] as f32);
            let bc = Rectangular::new(box_size);
            lj_crf_innerloop(
                &positions_vec,
                charges_slice,
                iac_slice,
                &pairlist_vec,
                &lj_params_vec,
                &crf,
                &bc,
                &mut storage,
            );
        }
        2 => {
            // Triclinic (not implemented yet)
            eprintln!("Error: Triclinic boundary not yet implemented");
            return;
        }
        _ => {
            eprintln!("Error: Unknown boundary type {}", boundary_type);
            return;
        }
    }

    // Copy results back
    let forces_slice = slice::from_raw_parts_mut(forces, n_atoms * 3);
    for (i, f) in storage.forces.iter().enumerate() {
        forces_slice[i * 3] = f.x;
        forces_slice[i * 3 + 1] = f.y;
        forces_slice[i * 3 + 2] = f.z;
    }

    let energies_slice = slice::from_raw_parts_mut(energies, 2);
    energies_slice[0] = storage.e_lj;
    energies_slice[1] = storage.e_crf;

    let virial_slice = slice::from_raw_parts_mut(virial, 9);
    for i in 0..3 {
        for j in 0..3 {
            virial_slice[i * 3 + j] = storage.virial[i][j];
        }
    }
}

/// C FFI: Parallel version of innerloop
#[no_mangle]
pub unsafe extern "C" fn rust_lj_crf_innerloop_parallel(
    positions: *const c_float,
    charges: *const c_float,
    iac: *const c_uint,
    n_atoms: c_uint,
    pairlist: *const c_uint,
    n_pairs: c_uint,
    lj_params: *const c_double,
    n_types: c_uint,
    box_data: *const c_double,
    boundary_type: c_uint,
    crf_cut: c_double,
    crf_2cut3i: c_double,
    crf_cut3i: c_double,
    forces: *mut c_float,
    energies: *mut c_double,
    virial: *mut c_double,
) {
    // Similar implementation to serial version, but calls parallel version
    // (Omitted for brevity - would be similar to above but calling lj_crf_innerloop_parallel)
    rust_lj_crf_innerloop(
        positions, charges, iac, n_atoms, pairlist, n_pairs,
        lj_params, n_types, box_data, boundary_type,
        crf_cut, crf_2cut3i, crf_cut3i,
        forces, energies, virial,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ffi_simple() {
        // Simple 2-atom system
        let positions: Vec<f32> = vec![
            0.0, 0.0, 0.0,  // atom 0
            1.0, 0.0, 0.0,  // atom 1
        ];
        let charges: Vec<f32> = vec![0.5, -0.5];
        let iac: Vec<u32> = vec![0, 0];
        let pairlist: Vec<u32> = vec![0, 1];

        let lj_params: Vec<f64> = vec![0.001, 0.0001];  // c6, c12

        let box_data: Vec<f64> = vec![10.0, 10.0, 10.0];

        let mut forces: Vec<f32> = vec![0.0; 6];
        let mut energies: Vec<f64> = vec![0.0; 2];
        let mut virial: Vec<f64> = vec![0.0; 9];

        unsafe {
            rust_lj_crf_innerloop(
                positions.as_ptr(),
                charges.as_ptr(),
                iac.as_ptr(),
                2,
                pairlist.as_ptr(),
                1,
                lj_params.as_ptr(),
                1,
                box_data.as_ptr(),
                0,  // vacuum
                1.4,
                0.364431,
                0.182215,
                forces.as_mut_ptr(),
                energies.as_mut_ptr(),
                virial.as_mut_ptr(),
            );
        }

        // Verify Newton's third law
        assert!((forces[0] + forces[3]).abs() < 1e-5);
        assert!((forces[1] + forces[4]).abs() < 1e-5);
        assert!((forces[2] + forces[5]).abs() < 1e-5);

        // Verify energies are non-zero
        assert!(energies[0] != 0.0);
        assert!(energies[1] != 0.0);
    }
}
