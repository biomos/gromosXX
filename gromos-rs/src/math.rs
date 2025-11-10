//! Mathematical primitives for molecular dynamics

use glam::{Vec3A, Mat3A};

/// 3D vector with SIMD acceleration (16-byte aligned)
pub type Vec3 = Vec3A;

/// 3x3 matrix with SIMD acceleration
pub type Mat3 = Mat3A;

/// Boundary condition types for periodic systems
pub trait BoundaryCondition: Send + Sync {
    /// Calculate nearest image vector between two positions
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3;

    /// Apply periodic boundary conditions to a position
    fn put_into_box(&self, pos: Vec3) -> Vec3;
}

/// Vacuum boundary conditions (no periodicity)
#[derive(Debug, Clone, Copy)]
pub struct Vacuum;

impl BoundaryCondition for Vacuum {
    #[inline(always)]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        rj - ri
    }

    #[inline(always)]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        pos
    }
}

/// Rectangular periodic boundary conditions
#[derive(Debug, Clone, Copy)]
pub struct Rectangular {
    pub box_size: Vec3,
    pub half_box: Vec3,
    pub inv_box: Vec3,
}

impl Rectangular {
    pub fn new(box_size: Vec3) -> Self {
        Self {
            box_size,
            half_box: box_size * 0.5,
            inv_box: box_size.recip(),
        }
    }
}

impl BoundaryCondition for Rectangular {
    #[inline(always)]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        let mut r = rj - ri;

        // Apply minimum image convention
        // Branchless version using SIMD
        r = r - self.box_size * (r * self.inv_box).round();

        r
    }

    #[inline(always)]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        pos - self.box_size * (pos * self.inv_box).floor()
    }
}

/// Triclinic periodic boundary conditions
#[derive(Debug, Clone, Copy)]
pub struct Triclinic {
    pub box_matrix: Mat3,      // Box vectors as columns
    pub inv_box_matrix: Mat3,  // Inverse for wrapping
}

impl Triclinic {
    pub fn new(box_matrix: Mat3) -> Self {
        Self {
            box_matrix,
            inv_box_matrix: box_matrix.inverse(),
        }
    }
}

impl BoundaryCondition for Triclinic {
    #[inline(always)]
    fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        let r = rj - ri;

        // Transform to fractional coordinates
        let frac = self.inv_box_matrix * r;

        // Apply minimum image
        let frac_nearest = frac - frac.round();

        // Transform back to Cartesian
        self.box_matrix * frac_nearest
    }

    #[inline(always)]
    fn put_into_box(&self, pos: Vec3) -> Vec3 {
        let frac = self.inv_box_matrix * pos;
        let frac_wrapped = frac - frac.floor();
        self.box_matrix * frac_wrapped
    }
}

/// Generic periodicity wrapper
#[derive(Debug, Clone)]
pub enum Periodicity {
    Vacuum(Vacuum),
    Rectangular(Rectangular),
    Triclinic(Triclinic),
}

impl Periodicity {
    #[inline]
    pub fn nearest_image(&self, ri: Vec3, rj: Vec3) -> Vec3 {
        match self {
            Periodicity::Vacuum(bc) => bc.nearest_image(ri, rj),
            Periodicity::Rectangular(bc) => bc.nearest_image(ri, rj),
            Periodicity::Triclinic(bc) => bc.nearest_image(ri, rj),
        }
    }

    #[inline]
    pub fn put_into_box(&self, pos: Vec3) -> Vec3 {
        match self {
            Periodicity::Vacuum(bc) => bc.put_into_box(pos),
            Periodicity::Rectangular(bc) => bc.put_into_box(pos),
            Periodicity::Triclinic(bc) => bc.put_into_box(pos),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vacuum_boundary() {
        let bc = Vacuum;
        let ri = Vec3::new(0.0, 0.0, 0.0);
        let rj = Vec3::new(1.0, 2.0, 3.0);

        let r = bc.nearest_image(ri, rj);
        assert_relative_eq!(r.x, 1.0);
        assert_relative_eq!(r.y, 2.0);
        assert_relative_eq!(r.z, 3.0);
    }

    #[test]
    fn test_rectangular_boundary() {
        let bc = Rectangular::new(Vec3::splat(10.0));
        let ri = Vec3::new(9.5, 0.0, 0.0);
        let rj = Vec3::new(0.5, 0.0, 0.0);

        let r = bc.nearest_image(ri, rj);
        // Should wrap: 0.5 - 9.5 = -9.0, nearest = 1.0
        assert_relative_eq!(r.x, 1.0, epsilon = 1e-6);
    }
}
