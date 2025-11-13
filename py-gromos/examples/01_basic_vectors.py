"""
Example 1: Basic Vector Operations
===================================

Demonstrates SIMD-accelerated vector math with zero-copy NumPy integration.
"""

import gromos
import numpy as np

print("=" * 60)
print("GROMOS-RS: Basic Vector Operations")
print("=" * 60)

# Create vectors
v1 = gromos.Vec3(1.0, 0.0, 0.0)
v2 = gromos.Vec3(0.0, 1.0, 0.0)
v3 = gromos.Vec3(0.0, 0.0, 1.0)

print(f"\nv1 = {v1}")
print(f"v2 = {v2}")
print(f"v3 = {v3}")

# Vector operations
print("\n--- Vector Operations ---")
v4 = v1 + v2
print(f"v1 + v2 = {v4}")

v5 = v1 * 2.5
print(f"v1 * 2.5 = {v5}")

# Dot product
dot = v1.dot(v2)
print(f"v1 · v2 = {dot}")

# Cross product
cross = v1.cross(v2)
print(f"v1 × v2 = {cross}")

# Length
length = v1.length()
print(f"|v1| = {length}")

# Distance between vectors
distance = v1.distance(v2)
print(f"distance(v1, v2) = {distance:.4f}")

# Normalization
v6 = gromos.Vec3(3.0, 4.0, 0.0)
print(f"\nv6 = {v6}")
print(f"|v6| = {v6.length()}")
v6_norm = v6.normalize()
print(f"normalize(v6) = {v6_norm}")
print(f"|normalize(v6)| = {v6_norm.length()}")

# NumPy integration (zero-copy where possible)
print("\n--- NumPy Integration ---")
np_array = v1.to_numpy()
print(f"v1 as NumPy array: {np_array}")
print(f"Array dtype: {np_array.dtype}")
print(f"Array shape: {np_array.shape}")

# Create from NumPy
np_vec = np.array([1.5, 2.5, 3.5], dtype=np.float32)
v7 = gromos.Vec3.from_numpy(np_vec)
print(f"\nCreate from NumPy {np_vec}: {v7}")

# Matrix operations
print("\n--- Matrix Operations ---")
identity = gromos.Mat3.identity()
print(f"Identity matrix:\n{identity}")

mat = gromos.Mat3.from_cols(v1, v2, v3)
print(f"\nMatrix from columns:\n{mat}")

mat_det = mat.determinant()
print(f"Determinant: {mat_det}")

mat_inv = mat.inverse()
print(f"Inverse:\n{mat_inv}")

# Matrix-vector multiplication
v_transformed = mat.mul_vec3(v1)
print(f"\nmat * v1 = {v_transformed}")

print("\n" + "=" * 60)
print("All operations completed successfully!")
print("=" * 60)
