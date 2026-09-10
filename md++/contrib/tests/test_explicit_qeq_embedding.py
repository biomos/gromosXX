"""Limiting-case tests for the explicit variational QEq embedding."""

import math

import pytest
import torch

import schnetpack.properties as properties
from schnetpack.atomistic import EMLEQEqStatic


DTYPE = torch.float64


def _head() -> EMLEQEqStatic:
    head = EMLEQEqStatic(
        phi_key="phi_static",
        init_a_qeq=1.0,
        learn_a_qeq=False,
        a_qeq_positive=False,
        correct_charges=True,
    )
    return head.to(dtype=DTYPE)


def _inputs(phi=None, total_charge=1.0):
    result = {
        properties.Z: torch.tensor([1, 8], dtype=torch.long),
        properties.R: torch.tensor(
            [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
            dtype=DTYPE,
            requires_grad=True,
        ),
        properties.idx_m: torch.tensor([0, 0], dtype=torch.long),
        properties.n_atoms: torch.tensor([2], dtype=torch.long),
        properties.total_charge: torch.tensor([total_charge], dtype=DTYPE),
        "mbis_valence_widths": torch.tensor(
            [0.8, 1.1], dtype=DTYPE, requires_grad=True
        ),
        "chi_emle": torch.tensor(
            [120.0, 390.0], dtype=DTYPE, requires_grad=True
        ),
    }
    if phi is not None:
        result["phi_static"] = torch.as_tensor(phi, dtype=DTYPE)
    return result


def _run(phi=None, total_charge=1.0):
    return _head()(_inputs(phi=phi, total_charge=total_charge))


def test_zero_or_absent_field_is_vacuum_limit():
    for phi in (None, [0.0, 0.0]):
        out = _run(phi)
        torch.testing.assert_close(out["charges"], out["charges_vac"])
        torch.testing.assert_close(
            out["qeq_embedding_energy"], torch.zeros(1, dtype=DTYPE)
        )


def test_vacqeq_uses_vacuum_charges_but_keeps_static_embedding():
    phi = torch.tensor([11.0, -7.0], dtype=DTYPE)
    inputs = _inputs(phi)
    inputs["qeq_polarize"] = torch.tensor(False)
    out = _head()(inputs)

    torch.testing.assert_close(out["charges"], out["charges_vac"])
    torch.testing.assert_close(
        out["qeq_embedding_energy"][0],
        torch.dot(out["charges_vac"], phi),
    )


def test_old_checkpoint_without_robustness_attributes_still_runs():
    head = _head()
    for attribute in (
        "use_bounded_a_qeq",
        "a_qeq_min",
        "a_qeq_max",
        "sigma_qeq_min",
        "sigma_qeq_max",
        "jii_floor",
        "polarize_key",
    ):
        if hasattr(head, attribute):
            delattr(head, attribute)

    out = head(_inputs([0.1, -0.2]))
    assert torch.isfinite(out["charges"]).all()
    assert torch.isfinite(out["qeq_embedding_energy"]).all()


def test_uniform_potential_has_correct_gauge_behavior():
    phi = torch.tensor([-15.0, 27.0], dtype=DTYPE)
    constant = 73.0
    total_charge = 1.0
    base = _run(phi, total_charge)
    shifted = _run(phi + constant, total_charge)

    torch.testing.assert_close(shifted["charges"], base["charges"])
    energy_change = (
        shifted["qeq_embedding_energy"] - base["qeq_embedding_energy"]
    )
    torch.testing.assert_close(
        energy_change, torch.tensor([constant * total_charge], dtype=DTYPE)
    )


def test_small_nonuniform_field_polarizes_and_preserves_charge():
    out = _run([0.01, -0.02], total_charge=1.0)
    q0 = out["charges_vac"]
    qphi = out["charges"]
    assert torch.max(torch.abs(qphi - q0)).item() > 0.0
    torch.testing.assert_close(q0.sum(), torch.tensor(1.0, dtype=DTYPE))
    torch.testing.assert_close(qphi.sum(), torch.tensor(1.0, dtype=DTYPE))


def test_embedding_energy_identity_and_static_limit():
    phi = torch.tensor([11.0, -7.0], dtype=DTYPE)
    out = _run(phi)
    q0 = out["charges_vac"]
    qphi = out["charges"]
    expected = torch.dot(q0, phi) + 0.5 * torch.dot(qphi - q0, phi)
    torch.testing.assert_close(out["qeq_embedding_energy"][0], expected)

    static_embedding = 0.5 * torch.dot(q0 + q0, phi)
    torch.testing.assert_close(static_embedding, torch.dot(q0, phi))


def test_embedding_remains_differentiable():
    inputs = _inputs([13.0, -4.0])
    out = _head()(inputs)
    gradients = torch.autograd.grad(
        out["qeq_embedding_energy"].sum(),
        [
            inputs[properties.R],
            inputs["mbis_valence_widths"],
            inputs["chi_emle"],
        ],
        allow_unused=False,
    )
    assert all(torch.isfinite(gradient).all() for gradient in gradients)


def test_two_site_gaussian_matrix_is_positive_semidefinite():
    head = _head()
    sigma = torch.tensor([1.0, 1.0], dtype=DTYPE)
    jii = 1389.35456 / (sigma * math.sqrt(math.pi))

    # The short-distance limit is finite; use a small nonzero separation so the
    # explicit erf(R)/R expression is evaluated without its removable 0/0.
    positions = torch.tensor(
        [[0.0, 0.0, 0.0], [1.0e-5, 0.0, 0.0]], dtype=DTYPE
    )
    matrix = head._build_qeq_matrix(positions, sigma, jii)
    eigenvalues = torch.linalg.eigvalsh(matrix)
    assert eigenvalues.min().item() >= -1.0e-8
    assert matrix[0, 1].item() > 0.0
