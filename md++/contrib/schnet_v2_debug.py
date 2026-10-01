"""Opt-in, nonperturbed SchNet-v2 force decomposition; see notebooks/README.md.

All forces are kJ/mol/angstrom, energies kJ/mol, charges e, atom IDs one-based.
Region energy partitions are differentiated against *all* coordinates, including
cross-region charge response. This module does not add any force to dynamics.
"""
import os
from pathlib import Path

import numpy as np
import torch
from schnet_v2 import SchNet_V2_Calculator as _ProductionCalculator


class SchNet_V2_Calculator(_ProductionCalculator):
    def __init__(self, *args, **kwargs):
        self._debug_model_path = str(kwargs.get("model_path", args[0] if args else ""))
        super().__init__(*args, **kwargs)

    def configure_debug_regions(self, atom_ids, n_ir, or_atom_ids):
        ids = np.asarray(atom_ids, dtype=np.int64)
        outer = np.asarray(or_atom_ids, dtype=np.int64)
        if ids.ndim != 1 or outer.ndim != 1 or not 0 < n_ir <= len(ids):
            raise ValueError("Specify physical cluster atom IDs, IR count, and OR IDs")
        all_ids = np.concatenate((ids, outer))
        if np.any(all_ids <= 0) or len(np.unique(all_ids)) != len(all_ids):
            raise ValueError("Debug requires unique, positive topology IDs (no polarizable sites)")
        self._debug_regions = (ids, int(n_ir), outer)

    def calculate_next_step(self, *args, **kwargs):
        dynamic = kwargs.get("dynamic_charges", args[3] if len(args) > 3 else False)
        if not dynamic:
            raise ValueError("Force decomposition currently requires dynamic charges")
        return super().calculate_next_step(*args, **kwargs)

    def _observe_force_decomposition(self, E_mlp, E_embedding, q0, qphi, phi,
                                     rq, ro, qo, system, n_link_atoms):
        step = int(getattr(self, "time_step", 0))
        every = int(os.environ.get("SCHNET_DEBUG_EVERY", "1"))
        if every <= 0:
            raise ValueError("SCHNET_DEBUG_EVERY must be positive")
        if step % every:
            return
        if not hasattr(self, "_debug_regions"):
            raise ValueError("Call configure_debug_regions before debug evaluation")
        ids, n_ir, outer_ids = self._debug_regions
        if n_link_atoms:
            raise ValueError("Debug currently supports uncapped IR+BR systems only")
        if len(ids) != len(rq) or len(outer_ids) != len(ro):
            raise ValueError("Debug region metadata does not match coordinate arrays")

        def array(t):
            return t.detach().cpu().numpy().copy()

        def forces(energy):
            if not energy.requires_grad:
                return np.zeros_like(array(rq)), np.zeros_like(array(ro))
            grads = torch.autograd.grad(energy, (rq, ro), retain_graph=True,
                                        allow_unused=True)
            return tuple(np.zeros_like(array(r)) if g is None else -array(g)
                         for g, r in zip(grads, (rq, ro)))

        record = dict(step=np.array(step), atom_ids=ids.copy(), or_atom_ids=outer_ids.copy(),
                      n_ir=np.array(n_ir), positions_A=array(rq), or_positions_A=array(ro),
                      atomic_numbers=np.asarray(system.numbers), or_charges_e=array(qo),
                      charges_vac_e=array(q0), charges_e=array(qphi), phi_kJmol_e=array(phi),
                      qeq_mode=np.array(getattr(self, "qeq_mode", "vacuum")),
                      model_path=np.array(getattr(self, "_debug_model_path", "unknown")),
                      evaluator_file=np.array(str(Path(__file__).resolve())),
                      damping=np.array(getattr(self, "electrostatic_damping", "soft")),
                      sigma_A=np.array(getattr(self, "electrostatic_sigma_A", 0.8)),
                      energy_unit=np.array("kJ/mol"), force_unit=np.array("kJ/mol/angstrom"))

        def add(name, energy):
            record["E_" + name] = array(energy)
            fq, fo = forces(energy)
            record["F_" + name] = fq
            record["F_" + name + "_OR"] = fo
            return fq, fo

        add("BuRNN", E_mlp)
        add("embedding", E_embedding)
        add("total", E_mlp + E_embedding)
        effective_q = 0.5 * (q0.reshape(-1) + qphi.reshape(-1))
        phi = phi.reshape(-1)
        site_energy = effective_q * phi
        partition_sum = site_energy.sum()
        if not torch.allclose(partition_sum, E_embedding, rtol=2e-5, atol=2e-3):
            raise ValueError("Model embedding differs from 0.5*(q0+qphi).phi; cannot partition")
        for label, selection in (("IR_OR", slice(0, n_ir)), ("BR_OR", slice(n_ir, None))):
            full = add(label, site_energy[selection].sum())
            direct = add(label + "_direct", (effective_q.detach()[selection] * phi[selection]).sum())
            # Fixed effective-charge derivative; polarized response includes both
            # geometry dependence and field-induced changes of the QEq solutions.
            for suffix, f, d in zip(("", "_OR"), full, direct):
                record["F_" + label + "_response" + suffix] = f - d

        self.last_force_decomposition = record
        directory = os.environ.get("SCHNET_DEBUG_DIR")
        if directory:
            directory = Path(directory)
            directory.mkdir(parents=True, exist_ok=True)
            # Refuse overwrites so restarts cannot silently destroy prior frames.
            with (directory / f"step_{step:012d}.npz").open("xb") as handle:
                np.savez_compressed(handle, **record)


class Pert_SchNet_V2_Calculator:
    def __init__(self, *args, **kwargs):
        raise ValueError("Force debug is not implemented for perturbed simulations")
