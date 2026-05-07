from pathlib import Path
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import yaml
import ase
import schnetpack as spk
import schnetpack.properties as properties

# Coulomb prefactors
K_E_KJMOL_NM_E2 = 138.935456
K_E_KJMOL_ANG_E2 = K_E_KJMOL_NM_E2 * 10.0
EV_TO_KJMOL = 96.4853321233
K_E_EV_NM_E2 = K_E_KJMOL_NM_E2 / EV_TO_KJMOL
K_E_EV_ANG_E2 = K_E_EV_NM_E2 * 10.0


class YamlParser:
    """Small helper to read SchNetPack training configuration (``config.yaml``)."""

    def __init__(self, yaml_file):
        self.data = self._load_yaml(yaml_file)

    def _load_yaml(self, yaml_file):
        with open(yaml_file, "r") as file:
            return yaml.safe_load(file)

    def get_cutoff(self):
        return self.data["globals"]["cutoff"]

    def get_property_keys(self):
        return list(self.data["data"]["property_units"].keys())

    def get_energy_key(self):
        for key in self.get_property_keys():
            if key.lower().startswith("v") or "ene" in key.lower():
                return key
        raise ValueError("No energy key found in property_units.")

    def get_force_key(self):
        for key in self.get_property_keys():
            if key.lower().startswith("f") or "force" in key.lower():
                return key
        raise ValueError("No force key found in property_units.")

    def get_energy_unit(self):
        return self.data["data"]["property_units"][self.get_energy_key()]

    def get_position_unit(self):
        force_unit = self.data["data"]["property_units"][self.get_force_key()]
        if "/" not in force_unit:
            raise ValueError(f"Cannot infer position unit from force unit '{force_unit}'")
        return force_unit.split("/")[-1]

    def task_output_names(self) -> List[str]:
        outputs = self.data.get("task", {}).get("outputs", [])
        return [out.get("name", "") for out in outputs]

    def has_output(self, name: str) -> bool:
        return name in self.task_output_names()

    def get_charge_prediction(self) -> bool:
        return any("charge" in name.lower() for name in self.task_output_names())

    def get_electronic_embedding(self):
        """
        Identifies if the model was trained with electronic embedding.
        """
        try:
            representation = self.data['model']['representation']
            embeddings = representation.get('electronic_embeddings', [])
            if embeddings:
                return True
            return False
        except Exception:
            return False
class ExtendedConverter(spk.interfaces.AtomsConverter):
    """AtomsConverter that optionally forwards global electronic state inputs."""

    def __call__(self, atoms):
        inputs = super().__call__(atoms)

        if "total_charge" in atoms.info:
            charge = atoms.info.get("total_charge", 0.0)
            inputs["total_charge"] = torch.tensor([charge], dtype=torch.float32, device=self.device)

        if "spin_multiplicity" in atoms.info:
            multiplicity = atoms.info.get("spin_multiplicity", 0.0)
            inputs["spin_multiplicity"] = torch.tensor([multiplicity], dtype=torch.float32, device=self.device)

        return inputs


class EMLEStaticEmbedding:
    r"""
    EMLE-style static embedding against OR/MM point charges.

    ML region:
      - point-like core charges q_core,i
      - valence Slater densities with charges q_val,i and widths s_i

    Runtime QEq field:
      - can be damped by an effective dielectric eps_eff
      - can use Gaussian-smeared MM charges to reduce spill-out
    """

    def __init__(
        self,
        width_min: float = 1.0e-8,
        eps_eff: float = 1.0,
        use_smeared_mm_phi: bool = True,
        sigma_mm_default_A: float = 1.0,
        sigma_mm_table_A: Optional[Dict[int, float]] = None,
    ):
        self.width_min = float(width_min)
        self.eps_eff = float(eps_eff)
        self.use_smeared_mm_phi = bool(use_smeared_mm_phi)
        self.sigma_mm_default_A = float(sigma_mm_default_A)
        self.sigma_mm_table_A = dict(sigma_mm_table_A or {})

    def sigma_mm_from_Z(
        self,
        z_or: torch.Tensor,
        device=None,
        dtype=None,
    ) -> torch.Tensor:
        vals = [
            self.sigma_mm_table_A.get(int(z), self.sigma_mm_default_A)
            for z in z_or.detach().cpu().tolist()
        ]
        return torch.tensor(
            vals,
            device=device or z_or.device,
            dtype=dtype or torch.float32,
        )

    def atomic_potential_on_or(
        self,
        r_mlp_A: torch.Tensor,
        q_core_e: torch.Tensor,
        q_val_e: torch.Tensor,
        s_val_A: torch.Tensor,
        r_or_A: torch.Tensor,
        cutoff_A: Optional[float] = None,
    ) -> torch.Tensor:
        if r_mlp_A.shape[0] == 0 or r_or_A.shape[0] == 0:
            return torch.zeros(
                (r_mlp_A.shape[0], r_or_A.shape[0]),
                device=r_mlp_A.device,
                dtype=r_mlp_A.dtype,
            )

        s_val_A = s_val_A.clamp_min(self.width_min)

        diff = r_mlp_A[:, None, :] - r_or_A[None, :, :]
        dist = torch.norm(diff, dim=-1).clamp_min(1.0e-8)
        inv_r = 1.0 / dist
        exp_term = torch.exp(-dist / s_val_A[:, None])

        v_core = q_core_e[:, None] * inv_r
        v_val = q_val_e[:, None] * (
            inv_r - exp_term * (inv_r + 0.5 / s_val_A[:, None])
        )
        v = K_E_KJMOL_ANG_E2 * (v_core + v_val)

        if cutoff_A is not None:
            v = v * (dist <= float(cutoff_A))

        return v

    def total_energy(
        self,
        r_mlp_A: torch.Tensor,
        q_core_e: torch.Tensor,
        q_val_e: torch.Tensor,
        s_val_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or_e: torch.Tensor,
        cutoff_A: Optional[float] = None,
    ) -> torch.Tensor:
        if r_mlp_A.shape[0] == 0 or r_or_A.shape[0] == 0:
            return torch.zeros((), device=r_mlp_A.device, dtype=r_mlp_A.dtype)

        v_iJ = self.atomic_potential_on_or(
            r_mlp_A=r_mlp_A,
            q_core_e=q_core_e,
            q_val_e=q_val_e,
            s_val_A=s_val_A,
            r_or_A=r_or_A,
            cutoff_A=cutoff_A,
        )
        return torch.sum(v_iJ * q_or_e[None, :])

    def potential_on_mlp_atoms_from_or(
        self,
        r_mlp_A: torch.Tensor,
        s_val_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or_e: torch.Tensor,
        cutoff_A: Optional[float] = None,
    ) -> torch.Tensor:
        """
        External potential on ML atoms due to OR/MM charges.

        Uses width-dependent Slater-like damping.
        DOES NOT depend on ML charges → removes feedback instability.
        """

        if r_or_A.shape[0] == 0:
            return torch.zeros(
                r_mlp_A.shape[0],
                device=r_mlp_A.device,
                dtype=r_mlp_A.dtype,
            )

        s_val_A = s_val_A.clamp_min(self.width_min)

        diff = r_mlp_A[:, None, :] - r_or_A[None, :, :]
        dist = torch.norm(diff, dim=-1).clamp_min(1.0e-8)
        inv_r = 1.0 / dist

        # Slater-like damping kernel
        exp_term = torch.exp(-dist / s_val_A[:, None])
        kernel = inv_r - exp_term * (inv_r + 0.5 / s_val_A[:, None])

        if cutoff_A is not None:
            kernel = kernel * (dist <= float(cutoff_A))

        phi = K_E_KJMOL_ANG_E2 * torch.sum(q_or_e[None, :] * kernel, dim=1)

        return phi


class InducedEmbeddingPlaceholder:
    """Placeholder for the induced EMLE term."""

    def energy(self, device, dtype):
        return torch.zeros((), device=device, dtype=dtype)


class SchNet_V2_Calculator:
    """BuRNN calculator with EMLE-style static embedding and runtime phi injection."""

    @staticmethod
    def model_dir_from_path(model_path: str | Path) -> Path:
        model_path = Path(model_path)
        return model_path.parent if model_path.name == "best_model" else model_path

    @staticmethod
    def model_trained_on_partial_charges(model_path: str | Path) -> bool:
        config_path = SchNet_V2_Calculator.model_dir_from_path(model_path) / "config.yaml"
        return YamlParser(config_path).get_charge_prediction()
    
    @staticmethod
    def model_trained_with_electronic_embedding(model_path: str) -> bool:
        """
        Check if a SchNetPack model was trained with electronic embedding.
        """
        model_dir = Path(model_path).parent if Path(model_path).name == "best_model" else Path(model_path)
        config_path = model_dir / "config.yaml"
        return YamlParser(config_path).get_electronic_embedding()
    def __init__(
        self,
        model_path: str,
        val_model_paths: List[str],
        nn_valid_freq: int,
        write_energy_freq: int,
        spin_multiplicity: int,
        total_charge: int,
        widths_key: str = "mbis_valence_widths",
        charges_key: str = "charges",
        qcore_key: str = "q_core",
        qval_key: str = "q_val",
        phi_key: str = "phi_static",
        phi_scale: float = 1.0,
        phi_clip_kjmol_per_e: Optional[float] = None,
    ) -> None:
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.torchdevice = torch.device(self.device)

        self.spin_multiplicity = spin_multiplicity
        self.total_charge = total_charge
        self.nn_valid_freq = nn_valid_freq
        self.write_energy_freq = write_energy_freq
        self.widths_key = widths_key
        self.charges_key = charges_key
        self.qcore_key = qcore_key
        self.qval_key = qval_key
        self.phi_key = phi_key

        self.phi_scale = float(phi_scale)
        self.phi_clip = phi_clip_kjmol_per_e

        self.pred_calculator = self.get_calculator(Path(model_path))
        self.val_calculators = [self.get_calculator(Path(p)) for p in val_model_paths]

        self.static_embedding = EMLEStaticEmbedding()
        self.induced_embedding = InducedEmbeddingPlaceholder()

        self.energy = None
        self.energy_mlp = None
        self.energy_static = None
        self.energy_induced = 0.0
        self.forces = None
        self.forces_mlp = None
        self.or_forces = np.zeros((0, 3), dtype=float)
        self.charges = None
        self.widths = None
        self.nn_valid_ene = 0.0
        self.nn_valid_maxF = 0.0

    def get_calculator(self, model_path: Path) -> spk.interfaces.SpkCalculator:
        model_args_path = model_path.parent / "config.yaml"
        model_file = os.path.join(model_path.parent, "best_model")

        yaml_parser = YamlParser(model_args_path)
        cutoff = yaml_parser.get_cutoff()
        energy_key = yaml_parser.get_energy_key()
        force_key = yaml_parser.get_force_key()
        energy_unit = yaml_parser.get_energy_unit()
        position_unit = yaml_parser.get_position_unit()
        has_charges_local = yaml_parser.get_charge_prediction()
        charges_key = self.charges_key if has_charges_local else None

        return spk.interfaces.SpkCalculator(
            model_file=model_file,
            dtype=torch.float32,
            converter=ExtendedConverter,
            neighbor_list=spk.transform.ASENeighborList(cutoff=cutoff),
            energy_key=energy_key,
            force_key=force_key,
            charges_key=charges_key,
            energy_unit=energy_unit,
            position_unit=position_unit,
            device=self.torchdevice,
        )

    def _run_output_modules(
        self,
        model,
        inputs: Dict[str, torch.Tensor],
        stop_before_force: bool = True,
    ) -> Dict[str, torch.Tensor]:
        force_key = getattr(self.pred_calculator, "force_key", "forces")
        for mod in model.output_modules:
            outs = getattr(mod, "model_outputs", [])
            if stop_before_force and ((force_key in outs) or ("forces" in outs)):
                continue
            inputs = mod(inputs)
        return inputs

    def _base_inputs(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_mlp_A: torch.Tensor,
    ) -> Dict[str, torch.Tensor]:
        inputs = calculator.converter(system)
        inputs[properties.R] = r_mlp_A

        model = calculator.model
        for mod in getattr(model, "input_modules", []):
            inputs = mod(inputs)
        inputs = model.representation(inputs)
        return inputs

    def _forward_model(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_mlp_A: torch.Tensor,
        extra_inputs: Optional[Dict[str, torch.Tensor]] = None,
    ) -> Dict[str, torch.Tensor]:
        inputs = self._base_inputs(calculator=calculator, system=system, r_mlp_A=r_mlp_A)
        if extra_inputs is not None:
            inputs.update(extra_inputs)
        return self._run_output_modules(calculator.model, inputs, stop_before_force=True)

    def _forward_model_with_runtime_phi(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_mlp_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or_e: torch.Tensor,
        cutoff_nm: Optional[float],
    ) -> Dict[str, torch.Tensor]:
        """
        Two-pass forward for runtime EMLE-QEq with outer-region potential.

        Pass 1:
            obtain widths needed to build the external OR/MM electrostatic potential.
        Pass 2:
            inject phi_static and rerun the output heads so EMLEQEqStatic solves
            the embedded QEq problem.
        """
        cutoff_A = None if cutoff_nm is None else float(cutoff_nm) * 10.0

        # pass 1: in-vacuo
        out0 = self._forward_model(
            calculator=calculator,
            system=system,
            r_mlp_A=r_mlp_A,
            extra_inputs=None,
        )

        if self.widths_key not in out0:
            raise KeyError(
                f"Model output does not contain '{self.widths_key}'. "
                "Runtime EMLE-QEq requires predicted MBIS valence widths."
            )
        if self.qcore_key not in out0:
            raise KeyError(
                f"Model output does not contain '{self.qcore_key}'. "
                "Runtime EMLE-QEq requires q_core from pass 1."
            )
        if self.qval_key not in out0:
            raise KeyError(
                f"Model output does not contain '{self.qval_key}'. "
                "Runtime EMLE-QEq requires q_val from pass 1."
            )

        # Build phi_static from OR/MM point charges
        phi_static = self.static_embedding.potential_on_mlp_atoms_from_or(
            r_mlp_A=r_mlp_A,
            s_val_A=out0[self.widths_key],
            r_or_A=r_or_A,
            q_or_e=q_or_e,
            cutoff_A=cutoff_A,
        )

        # scale field (stability control)
        phi_static = self.phi_scale * phi_static

        # optional clipping (debug safety)
        if self.phi_clip is not None:
            phi_static = torch.clamp(phi_static, -self.phi_clip, self.phi_clip)

        # pass 2: embedded QEq solve
        out1 = self._forward_model(
            calculator=calculator,
            system=system,
            r_mlp_A=r_mlp_A,
            extra_inputs={self.phi_key: phi_static},
        )
        out1[self.phi_key] = phi_static
        return out1

    def _evaluate_with_static_embedding(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_mlp_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or_e: torch.Tensor,
        cutoff_nm: Optional[float],
    ) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, np.ndarray]:
        cutoff_A = None if cutoff_nm is None else float(cutoff_nm) * 10.0

        out = self._forward_model_with_runtime_phi(
            calculator=calculator,
            system=system,
            r_mlp_A=r_mlp_A,
            r_or_A=r_or_A,
            q_or_e=q_or_e,
            cutoff_nm=cutoff_nm,
        )
        
        q = out[self.charges_key].detach()
        phi = out[self.phi_key].detach()

        imax = torch.argmax(torch.abs(q))
        print("\n===== DEBUG QEQ =====")
        print("max |q| atom:", int(imax))
        print("q_max:", float(q[imax]))
        print("phi_max_atom:", float(phi[imax]))
        
        def _stats(name, t):
            t = t.detach()
            return f"{name}: min={float(t.min()):.3e} max={float(t.max()):.3e} mean={float(t.mean()):.3e}"

        print(_stats("phi_static", phi))
        print(_stats("q_total", q))
        print(_stats("mbis_widths", out[self.widths_key]))
        print(_stats("sigma_qeq", out["sigma_qeq"]))
        print(_stats("Jii_qeq", out["Jii_qeq"]))
        print(_stats("q_core", out[self.qcore_key]))
        print(_stats("q_val", out[self.qval_key]))
        print(_stats("a_qeq", out["a_qeq"]))
        
        i = int(imax)

        ri = r_mlp_A[i]                      # (3,)
        r_or = r_or_A.detach()               # (Nor,3)
        q_or_np = q_or_e.detach()

        diff = ri[None, :] - r_or
        dist = torch.norm(diff, dim=1)

        dmin, jmin = torch.min(dist, dim=0)

        print("\n--- Local environment ---")
        print("closest OR atom:", int(jmin))
        print("distance Å:", float(dmin))
        print("q_or of that atom:", float(q_or_np[jmin]))
        
        chi = out["chi_emle"].detach()
        Jii = out["Jii_qeq"].detach()

        print("\n--- QEq terms (problem atom) ---")
        g = chi + phi

        print("chi_emle:", float(chi[i]))
        print("Jii_qeq:", float(Jii[i]))
        print("phi_static:", float(phi[i]))
        print("g = chi + phi:", float(g[i]))
        print("sigma_qeq:", float(out["sigma_qeq"][i].detach()))
        print("a_qeq:", float(out["a_qeq"][i].detach()))

        print("\n--- Charge sensitivity ---")
        print("|g| / Jii:", float(torch.abs(g[i]) / Jii[i]))

        g = chi + phi
        sens = torch.abs(g) / Jii
        imax_sens = torch.argmax(sens)

        print("\n--- Worst sensitivity atom ---")
        print("atom:", int(imax_sens))
        print("chi:", float(chi[imax_sens]))
        print("phi:", float(phi[imax_sens]))
        print("g:", float(g[imax_sens]))
        print("Jii:", float(Jii[imax_sens]))
        print("|g|/Jii:", float(sens[imax_sens]))
        print("q:", float(q[imax_sens]))
        
        if torch.abs(q).max() > 10.0:
            print("⚠️ WARNING: runaway charge detected")

        if self.charges_key not in out:
            raise KeyError(
                f"Model output does not contain '{self.charges_key}'. "
                "Static EMLE embedding requires predicted total atomic charges."
            )
        if self.widths_key not in out:
            raise KeyError(
                f"Model output does not contain '{self.widths_key}'. "
                "Static EMLE embedding requires predicted MBIS valence widths."
            )
        if self.qcore_key not in out:
            raise KeyError(
                f"Model output does not contain '{self.qcore_key}'. "
                "Static EMLE embedding requires q_core from EMLEQEqStatic."
            )
        if self.qval_key not in out:
            raise KeyError(
                f"Model output does not contain '{self.qval_key}'. "
                "Static EMLE embedding requires q_val from EMLEQEqStatic."
            )

        E_mlp = out[calculator.energy_key].sum()      # kJ/mol
        q_total = out[self.charges_key]
        q_core = out[self.qcore_key]
        q_val = out[self.qval_key]
        s_val = out[self.widths_key]

        E_static = self.static_embedding.total_energy(
            r_mlp_A=r_mlp_A,
            q_core_e=q_core,
            q_val_e=q_val,
            s_val_A=s_val,
            r_or_A=r_or_A,
            q_or_e=q_or_e,
            cutoff_A=cutoff_A,
        )

        E_induced = self.induced_embedding.energy(device=r_mlp_A.device, dtype=r_mlp_A.dtype)
        E_total = E_mlp + E_static + E_induced

        dE_dR_mlp, dE_dR_or = torch.autograd.grad(
            E_total,
            [r_mlp_A, r_or_A],
            create_graph=False,
            retain_graph=True,
            allow_unused=True,
        )
        if dE_dR_mlp is None:
            dE_dR_mlp = torch.zeros_like(r_mlp_A)
        if dE_dR_or is None:
            dE_dR_or = torch.zeros_like(r_or_A)

        dEmlp_dR = torch.autograd.grad(
            E_mlp,
            r_mlp_A,
            create_graph=False,
            retain_graph=False,
            allow_unused=True,
        )[0]
        if dEmlp_dR is None:
            dEmlp_dR = torch.zeros_like(r_mlp_A)

        F_total_mlp = (-dE_dR_mlp).detach().cpu().numpy()
        F_total_or = (-dE_dR_or).detach().cpu().numpy()
        F_mlp_only = (-dEmlp_dR).detach().cpu().numpy()
        
        print("E_mlp (kJ/mol)   =", float(E_mlp.detach().cpu()))
        print("E_static (kJ/mol)=", float(E_static.detach().cpu()))
        print("E_total (kJ/mol) =", float(E_total.detach().cpu()))

        print("sum q_total =", float(q_total.sum().detach().cpu()))
        print("sum q_core  =", float(q_core.sum().detach().cpu()))
        print("sum q_val   =", float(q_val.sum().detach().cpu()))

        print("width min/max/mean =",
            float(s_val.min().detach().cpu()),
            float(s_val.max().detach().cpu()),
            float(s_val.mean().detach().cpu()))

        print("q_total min/max =", float(q_total.min().detach().cpu()), float(q_total.max().detach().cpu()))
        print("q_core  min/max =", float(q_core.min().detach().cpu()), float(q_core.max().detach().cpu()))
        print("q_val   min/max =", float(q_val.min().detach().cpu()), float(q_val.max().detach().cpu()))

        if r_or_A.shape[0] > 0:
            d = torch.norm(r_mlp_A[:, None, :] - r_or_A[None, :, :], dim=-1)
            print("min d(IR/BR-OR) Å =", float(d.min().detach().cpu()))
            
        return (
            float(E_total.detach().cpu().item()),
            F_total_mlp,
            F_total_or,
            q_total.detach().cpu().numpy(),
            s_val.detach().cpu().numpy(),
            float(E_mlp.detach().cpu().item()),
            float(E_static.detach().cpu().item()),
            F_mlp_only,
        )

    def _evaluate_plain_mlp(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_mlp_A: torch.Tensor,
    ) -> Tuple[float, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        out = self._forward_model(calculator=calculator, system=system, r_mlp_A=r_mlp_A)
        E_mlp = out[calculator.energy_key].sum()
        dE_dR = torch.autograd.grad(E_mlp, r_mlp_A, create_graph=False, retain_graph=False, allow_unused=True)[0]
        if dE_dR is None:
            dE_dR = torch.zeros_like(r_mlp_A)

        q = out[self.charges_key].detach().cpu().numpy() if self.charges_key in out else None
        s = out[self.widths_key].detach().cpu().numpy() if self.widths_key in out else None
        return float(E_mlp.detach().cpu().item()), (-dE_dR).detach().cpu().numpy(), q, s

    def validate_prediction(
        self,
        system: ase.Atoms,
        use_static_embedding: bool = False,
        r_mlp_A: Optional[torch.Tensor] = None,
        r_or_A: Optional[torch.Tensor] = None,
        q_or_e: Optional[torch.Tensor] = None,
        cutoff_nm: Optional[float] = None,
    ) -> List[float]:
        if len(self.val_calculators) == 0:
            return []

        val_energies = []
        for val_calc in self.val_calculators:
            rq = r_mlp_A.detach().clone().requires_grad_(True) if r_mlp_A is not None else None
            if use_static_embedding:
                ro = r_or_A.detach().clone().requires_grad_(True)
                E_total, *_ = self._evaluate_with_static_embedding(
                    calculator=val_calc,
                    system=system,
                    r_mlp_A=rq,
                    r_or_A=ro,
                    q_or_e=q_or_e,
                    cutoff_nm=cutoff_nm,
                )
                val_energies.append(float(E_total))
            else:
                E_mlp, *_ = self._evaluate_plain_mlp(
                    calculator=val_calc,
                    system=system,
                    r_mlp_A=rq,
                )
                val_energies.append(float(E_mlp))
        return val_energies

    def validate_prediction_maxForceDeviation(
        self,
        system: ase.Atoms,
        use_static_embedding: bool = False,
        r_mlp_A: Optional[torch.Tensor] = None,
        r_or_A: Optional[torch.Tensor] = None,
        q_or_e: Optional[torch.Tensor] = None,
        cutoff_nm: Optional[float] = None,
    ) -> float:
        if len(self.val_calculators) == 0:
            return 0.0

        model_forces = [np.linalg.norm(self.forces, axis=1)]
        for val_calc in self.val_calculators:
            rq = r_mlp_A.detach().clone().requires_grad_(True) if r_mlp_A is not None else None
            if use_static_embedding:
                ro = r_or_A.detach().clone().requires_grad_(True)
                _, F_val, _, *_ = self._evaluate_with_static_embedding(
                    calculator=val_calc,
                    system=system,
                    r_mlp_A=rq,
                    r_or_A=ro,
                    q_or_e=q_or_e,
                    cutoff_nm=cutoff_nm,
                )
                model_forces.append(np.linalg.norm(F_val, axis=1))
            else:
                _, F_val, _, _ = self._evaluate_plain_mlp(
                    calculator=val_calc,
                    system=system,
                    r_mlp_A=rq,
                )
                model_forces.append(np.linalg.norm(F_val, axis=1))

        ensemble = np.mean(model_forces, axis=0)
        sigma = np.sqrt(np.mean([(f - ensemble) ** 2 for f in model_forces], axis=0))
        return float(np.max(sigma))

    def calculate_next_step(
        self,
        atomic_numbers: List[int],
        positions: List[List[float]],
        time_step: int,
        use_static_embedding: bool = False,
        or_positions_nm: Optional[List[List[float]]] = None,
        or_charges_e: Optional[List[float]] = None,
        cutoff_nm: Optional[float] = None,
        or_atomic_numbers: Optional[List[int]] = None,
    ) -> None:
        del or_atomic_numbers  # reserved for future induced embedding

        system = ase.Atoms(numbers=atomic_numbers, positions=positions)
        system.info["total_charge"] = float(self.total_charge)
        system.info["spin_multiplicity"] = float(self.spin_multiplicity)

        r_mlp_A = torch.tensor(positions, dtype=torch.float32, device=self.torchdevice, requires_grad=True)

        if use_static_embedding:
            if or_positions_nm is None or or_charges_e is None or cutoff_nm is None:
                raise ValueError(
                    "use_static_embedding=True requires or_positions_nm, or_charges_e, and cutoff_nm."
                )
            if len(or_positions_nm) > 0:
                r_or_A = torch.tensor(or_positions_nm, dtype=torch.float32, device=self.torchdevice, requires_grad=True) * 10.0
                q_or_e = torch.tensor(or_charges_e, dtype=torch.float32, device=self.torchdevice)
            else:
                r_or_A = torch.zeros((0, 3), dtype=torch.float32, device=self.torchdevice, requires_grad=True)
                q_or_e = torch.zeros((0,), dtype=torch.float32, device=self.torchdevice)

            (
                self.energy,
                self.forces,
                self.or_forces,
                self.charges,
                self.widths,
                self.energy_mlp,
                self.energy_static,
                self.forces_mlp,
            ) = self._evaluate_with_static_embedding(
                calculator=self.pred_calculator,
                system=system,
                r_mlp_A=r_mlp_A,
                r_or_A=r_or_A,
                q_or_e=q_or_e,
                cutoff_nm=cutoff_nm,
            )
            self.energy_induced = 0.0

            if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
                val_energies = self.validate_prediction(
                    system=system,
                    use_static_embedding=True,
                    r_mlp_A=r_mlp_A,
                    r_or_A=r_or_A,
                    q_or_e=q_or_e,
                    cutoff_nm=cutoff_nm,
                )
                self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(
                    system=system,
                    use_static_embedding=True,
                    r_mlp_A=r_mlp_A,
                    r_or_A=r_or_A,
                    q_or_e=q_or_e,
                    cutoff_nm=cutoff_nm,
                )
                if len(val_energies) == 1:
                    self.nn_valid_ene = (self.energy - val_energies[0]) / np.sqrt(2)
                else:
                    all_E = np.array(val_energies + [self.energy], dtype=float)
                    self.nn_valid_ene = all_E.std(ddof=1)
        else:
            self.energy, self.forces, self.charges, self.widths = self._evaluate_plain_mlp(
                calculator=self.pred_calculator,
                system=system,
                r_mlp_A=r_mlp_A,
            )
            self.energy_mlp = self.energy
            self.energy_static = 0.0
            self.energy_induced = 0.0
            self.forces_mlp = self.forces.copy()
            self.or_forces = np.zeros((0, 3), dtype=float)

            if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
                val_energies = self.validate_prediction(
                    system=system,
                    use_static_embedding=False,
                    r_mlp_A=r_mlp_A,
                )
                self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(
                    system=system,
                    use_static_embedding=False,
                    r_mlp_A=r_mlp_A,
                )
                if len(val_energies) == 1:
                    self.nn_valid_ene = (self.energy - val_energies[0]) / np.sqrt(2)
                else:
                    all_E = np.array(val_energies + [self.energy], dtype=float)
                    self.nn_valid_ene = all_E.std(ddof=1)

    def get_energy(self):
        return self.energy

    def get_energy_mlp(self):
        return self.energy_mlp

    def get_energy_static(self):
        return self.energy_static

    def get_energy_induced(self):
        return self.energy_induced

    def get_forces(self):
        return self.forces

    def get_forces_mlp(self):
        return self.forces_mlp

    def get_or_forces(self):
        return self.or_forces

    def get_nn_valid_ene(self):
        return self.nn_valid_ene

    def get_nn_valid_maxF(self):
        return self.nn_valid_maxF

    def get_charges(self):
        return self.charges

    def get_widths(self):
        return self.widths