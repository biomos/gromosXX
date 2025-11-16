"""
Basic tests for GROMOS-RS Python bindings
"""

import pytest
import numpy as np


def test_import():
    """Test that the module can be imported"""
    import gromos
    assert gromos.__version__ == "0.1.0"


def test_vec3_creation():
    """Test Vec3 creation and basic operations"""
    import gromos

    v = gromos.Vec3(1.0, 2.0, 3.0)
    assert v.x == 1.0
    assert v.y == 2.0
    assert v.z == 3.0


def test_vec3_operations():
    """Test Vec3 mathematical operations"""
    import gromos

    v1 = gromos.Vec3(1.0, 0.0, 0.0)
    v2 = gromos.Vec3(0.0, 1.0, 0.0)

    # Addition
    v3 = v1 + v2
    assert v3.x == 1.0
    assert v3.y == 1.0
    assert v3.z == 0.0

    # Scalar multiplication
    v4 = v1 * 2.0
    assert v4.x == 2.0

    # Dot product
    dot = v1.dot(v2)
    assert dot == 0.0

    # Cross product
    cross = v1.cross(v2)
    assert abs(cross.z - 1.0) < 1e-6


def test_vec3_numpy():
    """Test Vec3 NumPy integration"""
    import gromos

    v = gromos.Vec3(1.0, 2.0, 3.0)
    arr = v.to_numpy()

    assert isinstance(arr, np.ndarray)
    assert arr.shape == (3,)
    assert arr.dtype == np.float32
    assert np.allclose(arr, [1.0, 2.0, 3.0])

    # Create from NumPy
    np_arr = np.array([4.0, 5.0, 6.0], dtype=np.float32)
    v2 = gromos.Vec3.from_numpy(np_arr)
    assert v2.x == 4.0
    assert v2.y == 5.0
    assert v2.z == 6.0


def test_mat3():
    """Test Mat3 creation and operations"""
    import gromos

    mat = gromos.Mat3.identity()
    det = mat.determinant()
    assert abs(det - 1.0) < 1e-6


def test_box():
    """Test Box creation"""
    import gromos

    # Vacuum
    box = gromos.Box.vacuum()
    assert box.volume() >= 0

    # Rectangular
    box = gromos.Box.rectangular(3.0, 3.0, 3.0)
    assert abs(box.volume() - 27.0) < 1e-6

    dims = box.dimensions()
    assert abs(dims.x - 3.0) < 1e-6


def test_energy():
    """Test Energy creation and access"""
    import gromos

    energy = gromos.Energy(num_temperature_groups=1, num_energy_groups=1)

    assert energy.kinetic == 0.0
    assert energy.potential == 0.0
    assert energy.total() == 0.0

    # Get as dictionary
    d = energy.to_dict()
    assert 'total' in d
    assert 'kinetic' in d
    assert 'potential' in d


def test_state():
    """Test State creation"""
    import gromos

    state = gromos.State(num_atoms=10, num_temp_groups=1, num_energy_groups=1)
    assert state.num_atoms() == 10

    # Get arrays
    pos = state.positions()
    assert pos.shape == (10, 3)

    vel = state.velocities()
    assert vel.shape == (10, 3)

    forces = state.forces()
    assert forces.shape == (10, 3)


def test_state_set():
    """Test setting State data from NumPy"""
    import gromos

    state = gromos.State(num_atoms=5, num_temp_groups=1, num_energy_groups=1)

    # Set positions
    pos = np.random.rand(5, 3).astype(np.float32)
    state.set_positions(pos)

    # Set velocities
    vel = np.random.rand(5, 3).astype(np.float32)
    state.set_velocities(vel)


def test_configuration():
    """Test Configuration creation"""
    import gromos

    config = gromos.Configuration(
        num_atoms=100,
        num_temp_groups=1,
        num_energy_groups=1
    )

    state = config.current_state()
    assert state.num_atoms() == 100

    energy = config.current_energy()
    assert energy.total() == 0.0


def test_topology():
    """Test Topology creation"""
    import gromos

    topo = gromos.Topology()
    assert topo.num_atoms() >= 0
    assert topo.num_bonds() >= 0


def test_leapfrog():
    """Test LeapFrog integrator"""
    import gromos

    lf = gromos.LeapFrog(dt=0.002)
    assert lf.timestep() == 0.002


def test_velocity_verlet():
    """Test VelocityVerlet integrator"""
    import gromos

    vv = gromos.VelocityVerlet(dt=0.001)
    assert vv.timestep() == 0.001


def test_stochastic_dynamics():
    """Test StochasticDynamics integrator"""
    import gromos

    sd = gromos.StochasticDynamics(dt=0.002, gamma=0.1, temperature=300.0)
    assert sd.timestep() == 0.002


def test_gamd():
    """Test GaMD parameters and runner"""
    import gromos

    params = gromos.GamdParameters(sigma0=6.0, threshold_mode='lower')
    runner = gromos.GamdRunner(params)

    assert runner is not None


def test_gamd_invalid_mode():
    """Test GaMD with invalid threshold mode"""
    import gromos

    with pytest.raises(Exception):
        gromos.GamdParameters(sigma0=6.0, threshold_mode='invalid')


def test_eds():
    """Test EDS parameters and runner"""
    import gromos

    params = gromos.EDSParameters(num_states=4, smoothness=1.0)
    runner = gromos.EDSRunner(params)

    assert runner is not None


def test_remd():
    """Test REMD controller"""
    import gromos

    remd = gromos.ReplicaController(num_replicas=8, exchange_interval=1000)
    assert remd.num_replicas() == 8


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
