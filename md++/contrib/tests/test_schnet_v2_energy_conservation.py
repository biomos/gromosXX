import sys
import inspect
from pathlib import Path
from types import SimpleNamespace

import ase
import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import schnetpack as spk
from schnet_v2 import Pert_SchNet_V2_Calculator, SchNet_V2_Calculator


class _FakeConverter:
    def __call__(self, system):
        positions = torch.tensor(system.positions, dtype=torch.float64)
        return {spk.properties.R: positions}


class _FakeOutputModule:
    model_outputs = [
        "energy",
        "charges_vac",
        "charges",
        "qeq_embedding_energy",
    ]

    def __init__(self, energy_scale=1.0, charge_scale=1.0):
        self.energy_scale = energy_scale
        self.charge_scale = charge_scale

    def __call__(self, inputs):
        r_qeq_A = inputs[spk.properties.R]
        phi = inputs["phi_static"]

        e_mlp = self.energy_scale * (
            0.7 * torch.sum(r_qeq_A**2) + 0.03 * torch.sum(torch.sin(r_qeq_A))
        )

        raw_q0 = self.charge_scale * 0.2 * r_qeq_A[:, 0]
        q0 = raw_q0 - torch.mean(raw_q0)
        raw_qphi = raw_q0 - self.charge_scale * 0.1 * phi
        qphi = raw_qphi - torch.mean(raw_qphi)
        if not bool(inputs["qeq_polarize"].item()):
            qphi = q0
        embedding = 0.5 * torch.sum((q0 + qphi) * phi)

        inputs["energy"] = e_mlp.reshape(1)
        inputs["charges_vac"] = q0
        inputs["charges"] = qphi
        inputs["qeq_embedding_energy"] = embedding.reshape(1)
        return inputs


class _FakeModel:
    input_modules = []

    def __init__(self, energy_scale=1.0, charge_scale=1.0):
        self.output_modules = [_FakeOutputModule(energy_scale, charge_scale)]

    def representation(self, inputs):
        return inputs


def _fake_calculator(energy_scale=1.0, charge_scale=1.0):
    return SimpleNamespace(
        converter=_FakeConverter(),
        model=_FakeModel(energy_scale, charge_scale),
        energy_key="energy",
        force_key="forces",
    )


def _evaluate(calc, r_qeq_A, r_or_A, q_or, n_link_atoms):
    system = ase.Atoms(numbers=[1] * r_qeq_A.shape[0], positions=r_qeq_A.detach().numpy())
    return calc._b2_eval_model(
        calculator=_fake_calculator(),
        system=system,
        r_qeq_A=r_qeq_A.detach().clone().requires_grad_(True),
        r_or_A=r_or_A.detach().clone().requires_grad_(True),
        q_or=q_or,
        cutoff_nm=1.2,
        n_link_atoms=n_link_atoms,
    )


def test_b2_forces_are_finite_difference_gradient_of_total_energy():
    calc = SchNet_V2_Calculator.__new__(SchNet_V2_Calculator)

    r_qeq_A = torch.tensor(
        [
            [0.00, 0.00, 0.00],
            [10.00, 0.00, 0.00],
            [5.00, 5.00, 5.00],  # trailing link atom
        ],
        dtype=torch.float64,
    )
    r_or_A = torch.tensor(
        [
            [5.02, 5.00, 5.00],
            [8.00, -3.00, 2.50],
        ],
        dtype=torch.float64,
    )
    q_or = torch.tensor([0.8, -0.35], dtype=torch.float64)

    n_link_atoms = 1
    energy, f_qeq, f_or, *_ = _evaluate(calc, r_qeq_A, r_or_A, q_or, n_link_atoms)

    h = 1.0e-4
    for coords, forces, label in ((r_qeq_A, f_qeq, "qeq"), (r_or_A, f_or, "or")):
        for atom in range(coords.shape[0]):
            for dim in range(3):
                plus_qeq = r_qeq_A.clone()
                minus_qeq = r_qeq_A.clone()
                plus_or = r_or_A.clone()
                minus_or = r_or_A.clone()

                if label == "qeq":
                    plus_qeq[atom, dim] += h
                    minus_qeq[atom, dim] -= h
                else:
                    plus_or[atom, dim] += h
                    minus_or[atom, dim] -= h

                e_plus = _evaluate(calc, plus_qeq, plus_or, q_or, n_link_atoms)[0]
                e_minus = _evaluate(calc, minus_qeq, minus_or, q_or, n_link_atoms)[0]
                finite_difference_force = -(e_plus - e_minus) / (2.0 * h)

                assert forces[atom, dim] == pytest.approx(
                    finite_difference_force,
                    rel=2.0e-8,
                    abs=1.0e-3,
                )

    assert np.isfinite(energy)


def test_link_atom_phi_is_zeroed_before_embedding_energy():
    phi = torch.tensor([1.0, 2.0, 999.0])
    zeroed = SchNet_V2_Calculator._zero_link_atom_phi(phi, n_link_atoms=1)

    assert zeroed.tolist() == [1.0, 2.0, 0.0]
    assert phi.tolist() == [1.0, 2.0, 999.0]


def test_runtime_qeq_mode_switches_between_vacqeq_and_polqeq():
    calc = SchNet_V2_Calculator.__new__(SchNet_V2_Calculator)
    r_qeq_A = torch.tensor(
        [[0.0, 0.0, 0.0], [1.3, 0.2, 0.0]], dtype=torch.float64
    )
    r_or_A = torch.tensor([[2.0, 0.4, 0.0]], dtype=torch.float64)
    q_or = torch.tensor([0.2], dtype=torch.float64)

    calc.qeq_mode = "vacuum"
    vacuum = _evaluate(calc, r_qeq_A, r_or_A, q_or, n_link_atoms=0)
    calc.qeq_mode = "polarized"
    polarized = _evaluate(calc, r_qeq_A, r_or_A, q_or, n_link_atoms=0)

    assert not np.allclose(vacuum[3], polarized[3])


def test_vacqeq_is_the_constructor_default():
    signature = inspect.signature(SchNet_V2_Calculator.__init__)
    assert signature.parameters["qeq_mode"].default == "vacuum"
    assert signature.parameters["electrostatic_sigma_A"].default == pytest.approx(0.01)


def test_dynamic_validation_uses_mlp_only_energy_and_forces():
    calc = SchNet_V2_Calculator.__new__(SchNet_V2_Calculator)
    calc.val_calculators = [
        _fake_calculator(energy_scale=1.0, charge_scale=0.8),
        _fake_calculator(energy_scale=1.0, charge_scale=1.2),
    ]

    system = ase.Atoms(
        numbers=[1, 1, 1],
        positions=[
            [0.00, 0.00, 0.00],
            [10.00, 0.00, 0.00],
            [5.00, 5.00, 5.00],
        ],
    )
    r_qeq_A = torch.tensor(system.positions, dtype=torch.float64, requires_grad=True)
    r_or_A = torch.tensor(
        [
            [5.02, 5.00, 5.00],
            [8.00, -3.00, 2.50],
        ],
        dtype=torch.float64,
        requires_grad=True,
    )
    q_or = torch.tensor([0.8, -0.35], dtype=torch.float64)
    n_link_atoms = 1

    prod = _evaluate(calc, r_qeq_A, r_or_A, q_or, n_link_atoms)
    calc.energy = prod[0]
    calc.forces = prod[1]
    calc.energy_mlp = prod[4]
    calc.forces_mlp = prod[5]

    val_energies = calc.validate_prediction(
        system=system,
        dynamic_charges=True,
        r_qeq_A=r_qeq_A,
        r_or_A=r_or_A,
        q_or=q_or,
        cutoff_nm=1.2,
        n_link_atoms=n_link_atoms,
    )

    expected = [
        calc._b2_eval_model(
            calculator=val_calc,
            system=system,
            r_qeq_A=r_qeq_A.detach().clone().requires_grad_(True),
            r_or_A=r_or_A.detach().clone().requires_grad_(True),
            q_or=q_or,
            cutoff_nm=1.2,
            n_link_atoms=n_link_atoms,
        )[4]
        for val_calc in calc.val_calculators
    ]

    assert val_energies == pytest.approx(expected)

    assert len(calc._last_val_embedding_energies) == len(calc.val_calculators)

    mlp_maxf = calc.validate_prediction_maxForceDeviation(
        system=system,
        dynamic_charges=True,
        r_qeq_A=r_qeq_A,
        r_or_A=r_or_A,
        q_or=q_or,
        cutoff_nm=1.2,
        n_link_atoms=n_link_atoms,
    )

    # The fake committee differs only in its charge/embedding head, while its
    # MLP energy surface is identical. MLP-only force validation must therefore
    # remain exactly zero.
    assert mlp_maxf == pytest.approx(0.0)


def test_perturbed_dynamic_qeq_returns_lambda_consistent_or_forces_and_derivative():
    calc = Pert_SchNet_V2_Calculator.__new__(Pert_SchNet_V2_Calculator)
    calc.pred_calculator = _fake_calculator()
    calc.val_calculators = []
    calc.states_idx = {"A": [0], "B": [1]}
    calc.qmzone_size = 2
    calc.num_states = 2
    calc.lam = 0.35
    calc.total_charge = 0
    calc.spin_multiplicity = 1
    calc.torchdevice = torch.device("cpu")
    calc.qeq_mode = "vacuum"

    def vacuum_prediction(system):
        positions = np.asarray(system.positions)
        energy = float(0.4 * np.sum(positions**2))
        forces = -0.8 * positions
        return energy, forces

    calc.predict_energy_and_forces = vacuum_prediction
    atomic_numbers = [1, 1, 1]
    positions = [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.4, 1.0, 0.2]]
    or_positions_nm = [[0.25, -0.04, 0.03]]
    or_charges_e = [0.3]

    calc.calculate_next_step(
        atomic_numbers,
        positions,
        time_step=0,
        dynamic_charges=True,
        or_positions_nm=or_positions_nm,
        or_charges_e=or_charges_e,
        cutoff_nm=1.2,
    )

    state_a = calc.states["A"]
    state_b = calc.states["B"]
    expected_or = (1.0 - calc.lam) * state_a["forces_or"] + calc.lam * state_b["forces_or"]
    expected_derivative = (
        -state_a["energy_burnn"]
        + state_a["energy_vac"]
        + state_b["energy_burnn"]
        - state_b["energy_vac"]
    )
    np.testing.assert_allclose(calc.or_forces, expected_or)
    assert calc.derivative == pytest.approx(expected_derivative)
    expected_charges = np.zeros(len(atomic_numbers))
    for state, weight in (("A", 1.0 - calc.lam), ("B", calc.lam)):
        for local_index, full_index in enumerate(calc.states[state]["burnn_full_indices"]):
            expected_charges[full_index] += (
                weight * calc.states[state]["charges_burnn"][local_index]
            )
    np.testing.assert_allclose(calc.charges, expected_charges)

    # The algebraic lambda derivative must also match a central difference of
    # the fully embedded perturbed energy at fixed coordinates and OR field.
    h = 1.0e-5
    energies = []
    for lam in (calc.lam - h, calc.lam + h):
        calc.lam = lam
        calc.calculate_next_step(
            atomic_numbers,
            positions,
            time_step=1,
            dynamic_charges=True,
            or_positions_nm=or_positions_nm,
            or_charges_e=or_charges_e,
            cutoff_nm=1.2,
        )
        energies.append(calc.energy)
    finite_difference = (energies[1] - energies[0]) / (2.0 * h)
    assert finite_difference == pytest.approx(expected_derivative, rel=1.0e-8, abs=1.0e-5)

    # OR forces must be the coordinate gradient of that same lambda-weighted
    # embedded energy, not a separately mixed point-charge approximation.
    calc.lam = 0.35
    h_nm = 1.0e-5
    displaced_energies = []
    for displacement in (-h_nm, h_nm):
        displaced_or = np.array(or_positions_nm, dtype=float)
        displaced_or[0, 0] += displacement
        calc.calculate_next_step(
            atomic_numbers,
            positions,
            time_step=2,
            dynamic_charges=True,
            or_positions_nm=displaced_or.tolist(),
            or_charges_e=or_charges_e,
            cutoff_nm=1.2,
        )
        displaced_energies.append(calc.energy)
    # Python returns kJ mol^-1 A^-1 while the displacement above is in nm.
    or_force_fd_A = -(
        displaced_energies[1] - displaced_energies[0]
    ) / (2.0 * h_nm * 10.0)
    assert expected_or[0, 0] == pytest.approx(or_force_fd_A, rel=2.0e-3, abs=2.0e-2)
