"""Regression tests for adaptive QM-buffer center and state semantics."""

from pathlib import Path


def _select_buffer(inner, static, adaptive_groups, cutoff):
    """Small reference model of the charge-group selection contract."""
    cutoff2 = cutoff * cutoff
    selected = []
    for group in adaptive_groups:
        center = group[0][1]
        if any(sum((a - b) ** 2 for a, b in zip(qm_pos, center)) < cutoff2
               for _, qm_pos in inner):
            selected.extend(atom for atom, _ in group)
    return sorted([atom for atom, _ in inner] + [atom for atom, _ in static]
                  + selected)


def test_static_buffer_atoms_are_not_adaptive_cutoff_centers():
    inner = [(0, (0.0, 0.0, 0.0))]
    static = [(atom, (1.0 + atom * 0.1, 0.0, 0.0)) for atom in range(1, 4)]
    near_inner = [(4, (0.4, 0.0, 0.0)),
                  (5, (0.41, 0.0, 0.0)),
                  (6, (0.39, 0.0, 0.0))]
    near_static_only = [(7, (1.2, 0.0, 0.0)),
                        (8, (1.21, 0.0, 0.0)),
                        (9, (1.19, 0.0, 0.0))]

    selected = _select_buffer(inner, static,
                              [near_inner, near_static_only], cutoff=0.5)

    assert selected == list(range(7))


def test_adaptive_group_can_leave_and_reenter_while_static_atoms_remain():
    inner = [(0, (0.0, 0.0, 0.0))]
    static = [(1, (1.0, 0.0, 0.0)), (2, (1.1, 0.0, 0.0))]

    def water(x):
        return [(3, (x, 0.0, 0.0)),
                (4, (x + 0.01, 0.0, 0.0)),
                (5, (x - 0.01, 0.0, 0.0))]

    assert _select_buffer(inner, static, [water(0.4)], 0.5) == list(range(6))
    assert _select_buffer(inner, static, [water(0.8)], 0.5) == [0, 1, 2]
    assert _select_buffer(inner, static, [water(0.4)], 0.5) == list(range(6))


def test_cpp_uses_distinct_adaptive_and_or_center_paths():
    source = (Path(__file__).parents[2]
              / "src/interaction/qmmm/qm_zone.cc").read_text()
    buffer_function = source.split(
        "int interaction::QM_Zone::_get_buffer_atoms", 1)[1].split(
        "void interaction::QM_Zone::get_or_cutoff_centers", 1)[0]

    assert "if (topo.is_qm(it->index)) adaptive_centers.push_back" in buffer_function
    assert "gather_chargegroups_from_centers<B>" in buffer_function
    assert "get_or_cutoff_centers" not in buffer_function
