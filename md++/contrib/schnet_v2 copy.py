import numpy as np
from pathlib import Path
import os
import torch
import yaml
import math
import ase
import schnetpack as spk

# Coulomb prefactor
K_E_KJMOL_NM_E2 = 138.935456  # kJ/mol·nm·e^-2
K_E_kJMOL_ANG_E2 = K_E_KJMOL_NM_E2 * 10

EV_TO_KJMOL = 96.4853321233
K_E_EV_NM_E2 = 138.935456 / EV_TO_KJMOL
K_E_EV_ANG_E2 = K_E_EV_NM_E2 * 10.0
class YamlParser:
    def __init__(self, yaml_file):
        """
        Initializes the parser by loading the YAML file.
        
        :param yaml_file: Path to the YAML file.
        """
        self.data = self._load_yaml(yaml_file)
    
    def _load_yaml(self, yaml_file):
        """
        Loads and parses the YAML file.
        
        :param yaml_file: Path to the YAML file.
        :return: Parsed YAML data as a dictionary.
        """
        with open(yaml_file, 'r') as file:
            return yaml.safe_load(file)
    
    def get_cutoff(self):
        """
        Retrieves the cutoff value from the parsed YAML data.
        
        :return: The cutoff value from the YAML file.
        """
        try:
            return self.data['globals']['cutoff']
        except KeyError as e:
            raise KeyError(f"Missing key in YAML file: {e}")
        
    def get_property_keys(self):
        """
        Returns all property keys defined under data.property_units.
        """
        try:
            return list(self.data["data"]["property_units"].keys())
        except KeyError as e:
            raise KeyError(f"Missing key in YAML file: {e}")

    def get_energy_key(self):
        """
        Attempts to identify the energy property key (usually starts with 'V' or contains 'ene').
        """
        try:
            property_keys = self.get_property_keys()
            # Heuristic: look for a key that starts with 'V' or contains 'ene'
            for key in property_keys:
                if key.lower().startswith("v") or "ene" in key.lower():
                    return key
            raise ValueError("No energy key found in property_units.")
        except Exception as e:
            raise ValueError(f"Error finding energy key: {e}")

    def get_force_key(self):
        """
        Attempts to identify the force property key (usually starts with 'F' or contains 'force').
        """
        try:
            property_keys = self.get_property_keys()
            for key in property_keys:
                if key.lower().startswith("f") or "force" in key.lower():
                    return key
            raise ValueError("No force key found in property_units.")
        except Exception as e:
            raise ValueError(f"Error finding force key: {e}")
        
    def get_energy_unit(self):
        """
        Returns the energy unit defined in data.property_units.
        """
        try:
            energy_key = self.get_energy_key()
            return self.data["data"]["property_units"][energy_key]
        except KeyError as e:
            raise KeyError(f"Energy unit not found in YAML: {e}")
        
    def get_position_unit(self):
        """
        Infers the position unit from the force unit definition.
        Example: 'kJ/mol/Ang' -> 'Ang'
        """
        try:
            force_key = self.get_force_key()
            force_unit = self.data["data"]["property_units"][force_key]

            if "/" in force_unit:
                return force_unit.split("/")[-1]
            else:
                raise ValueError(f"Cannot infer position unit from force unit '{force_unit}'")

        except KeyError as e:
            raise KeyError(f"Position unit not found in YAML: {e}")
    
    def get_charge_prediction(self):
        """
        Identifies if the model was trained to predict charges.
        """
        try:
            outputs = self.data['task']['outputs']
            for out in outputs:
                name = out.get('name', '').lower()
                if 'charge' in name:
                    return True
            return False
        except Exception:
            return False

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
    def __init__(self, *args, **kwargs):
        """
        Initialize the converter.
        Detects which extra keys the model expects.
        """
        super().__init__(*args, **kwargs)
        # Detect which inputs the model actually expects
        self.model_input_keys = getattr(self, "model_input_keys", [])

    def __call__(self, atoms):
        inputs = super().__call__(atoms)

        # Add total_charge only if the model expects it
        if "total_charge" in atoms.info:
            charge = atoms.info.get("total_charge", 0.0)
            inputs["total_charge"] = torch.tensor([charge], dtype=torch.float32, device=self.device)

        # Add spin_multiplicity only if the model expects it
        if "spin_multiplicity" in atoms.info:
            multiplicity = atoms.info.get("spin_multiplicity", 0.0)
            inputs["spin_multiplicity"] = torch.tensor([multiplicity], dtype=torch.float32, device=self.device)
            
        # Add per-atom external potential phi only if present
        if "phi" in atoms.arrays:
            phi = atoms.arrays["phi"]
            phi = torch.tensor(phi, dtype=torch.float32, device=self.device)
            # ensure shape (N, 1) or (N,) depending on what your model expects
            inputs["phi"] = phi

        # Optional ExternalCoulombEmbedding inference inputs.
        # These are injected directly into the model input dict so the
        # output head can add the external electrostatic energy term.
        if "or_positions" in atoms.arrays:
            or_pos = torch.tensor(atoms.arrays["or_positions"], dtype=torch.float32, device=self.device)
            inputs["or_positions"] = or_pos

        if "or_charges" in atoms.arrays:
            or_q = torch.tensor(atoms.arrays["or_charges"], dtype=torch.float32, device=self.device)
            inputs["or_charges"] = or_q

        if "or_mask" in atoms.arrays:
            or_mask = torch.tensor(atoms.arrays["or_mask"], dtype=torch.bool, device=self.device)
            inputs["or_mask"] = or_mask
        return inputs
    

class SchNet_V2_Calculator:
    """
    A calculator class for predicting energy and forces using SchNet models.

    This class provides functionality to load SchNet models, perform energy, 
    force and optionally partial charge predictions for atomic systems
    and validate predictions using additional validation models.

    Attributes:
        device (str): The computation device ('cuda' or 'cpu').
        pred_calculator (spk.interfaces.SpkCalculator): The primary predictive calculator.
        val_calculators (list): A list of validation calculators.
        nn_valid_freq (int): Frequency for performing validation (in terms of time steps).
        write_energy_freq (int): Frequency for writing energy values.
        energy (float): The predicted energy for the current system.
        forces (np.ndarray): The predicted forces for the current system.
        nn_valid_dev (float): Validation deviation calculated during validation.

    Methods:
        model_trained_on_partial_charges(model_path: str) -> bool:
            Checks if a SchNetPack model can predict partial charges
            
        __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int) -> None:
            Initialize the calculator with model paths and configuration.

        get_environment(self, model_args) -> spk.environment.SimpleEnvironmentProvider:
            Retrieve the environment provider for the model.

        get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
            Create a SchnetPack calculator from the given model path.

        predict_energy_and_forces(self, system: ase.Atoms):
            Predict energy and forces for the given atomic system.

        predict_partial_charges(self, system: ase.Atoms):
            Predict partial charges for the given atomic system.

        validate_prediction(self, system: ase.Atoms):
            Validate predictions using the validation calculators.
            
        validate_prediction_maxForceDeviation(self, system: ase.Atoms):
            Calculate Maximum force committee disagreement among all atoms in a structure.

        calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
            Perform energy and force prediction, and validate the predictions if necessary.

        get_energy(self):
            Get the predicted energy for the last atomic system.

        get_forces(self):
            Get the predicted forces for the last atomic system.

        get_nn_valid_ene(self):
            Get the validation deviation calculated during the last validation step.
            
        get_nn_valid_maxF(self):
            Get the maximum force committee disagreement among all atom during the last validation step.
    """
    
    @staticmethod
    def model_trained_on_partial_charges(model_path: str) -> bool:
        """
        Check if a SchNetPack model was trained on partial charges.
        """
        model_dir = Path(model_path).parent if Path(model_path).name == "best_model" else Path(model_path)
        config_path = model_dir / "config.yaml"
        
        return YamlParser(config_path).get_charge_prediction()
        
    @staticmethod
    def model_trained_with_electronic_embedding(model_path: str) -> bool:
        """
        Check if a SchNetPack model was trained with electronic embedding.
        """
        model_dir = Path(model_path).parent if Path(model_path).name == "best_model" else Path(model_path)
        config_path = model_dir / "config.yaml"
        return YamlParser(config_path).get_electronic_embedding()
        

    def __init__(self, model_path: str, val_model_paths: list,
                 nn_valid_freq: int, write_energy_freq: int,
                 spin_multiplicity:int,total_charge:int) -> None:
        """
        Initialize the calculator.
        """
        # Device
        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.torchdevice = torch.device(self.device)

        # Detect if the primary model was trained on charges
        self.has_partial_charges = self.model_trained_on_partial_charges(model_path)

        # Initialize primary calculator
        self.pred_calculator = self.get_calculator(model_path=Path(model_path))

        # Validation calculators
        self.val_calculators = []
        for val_path in val_model_paths:
            self.val_calculators.append(self.get_calculator(model_path=Path(val_path)))

        self.nn_valid_freq = nn_valid_freq
        self.write_energy_freq = write_energy_freq
        self.spin_multiplicity = spin_multiplicity
        self.total_charge = total_charge
    
    def get_calculator(self, model_path) -> spk.interfaces.SpkCalculator:
        """
        Create a SchnetPack calculator from the given model path.

        Args:
            model_path (Path): Path to the SchNet model.

        Returns:
            spk.interfaces.SpkCalculator: A configured SchnetPack calculator.
        """
        model_args_path = model_path.parent / 'config.yaml'
        model_path = os.path.join(model_path.parent,'best_model')

        # get cutoff and properties from configuration file
        yaml_parser = YamlParser(model_args_path)
        cutoff = yaml_parser.get_cutoff()
        energy_key = yaml_parser.get_energy_key()
        force_key = yaml_parser.get_force_key()
        energy_unit = yaml_parser.get_energy_unit()
        position_unit = yaml_parser.get_position_unit()
        
        # Only provide charges_key if model was trained on partial charges
        charges_key = "charges" if self.has_partial_charges else None
        
        calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            converter=ExtendedConverter,
            neighbor_list=spk.transform.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key=energy_key, # name of energy property in model
            force_key=force_key, # name of force property in model
            charges_key=charges_key,
            energy_unit=energy_unit,
            position_unit=position_unit,
            device= self.torchdevice, # device for computation
        )
        
        return calculator
        
    def _energy_forces(
        self,
        calculator,
        system,
        r_mlp_A,
        r_or_A: torch.Tensor | None = None,
        q_or: torch.Tensor | None = None,
        compute_or_forces: bool = False,
    ):
        """
        Evaluate model energy and forces.

        Modes:
        - vacuum / legacy: only r_mlp_A provided
        - external embedding via model head: r_or_A and q_or are injected into inputs
        - runtime fallback: if OR inputs are provided but the model energy does not depend on
          them, add E_ext = sum_i q_i(R) * phi_i(R, OR) outside the checkpoint using the
          predicted vacuum charges. This keeps inference usable even without retraining.

        Returns
        -------
        energy : float
            Total energy in the calculator energy units.
        F_qm : np.ndarray
            Forces on QM atoms, shape (N_qm, 3).
        F_or : np.ndarray
            Forces on OR atoms, shape (N_or, 3). Zero-size in vacuum mode.
        q_pred : np.ndarray | None
            Predicted QM charges if available, else None.
        """
        inputs = calculator.converter(system)
        inputs[spk.properties.R] = r_mlp_A

        has_or = (r_or_A is not None) and (q_or is not None)
        if has_or:
            inputs["or_positions"] = r_or_A
            inputs["or_charges"] = q_or
            inputs["or_mask"] = torch.ones(
                q_or.shape[0], device=q_or.device, dtype=torch.bool
            )

        model = calculator.model

        input_modules = getattr(model, "input_modules", None)
        if input_modules is not None:
            for mod in input_modules:
                inputs = mod(inputs)

        inputs = model.representation(inputs)

        force_key = getattr(calculator, "force_key", "forces")
        for mod in model.output_modules:
            outs = getattr(mod, "model_outputs", [])
            if (force_key in outs) or ("forces" in outs):
                continue
            inputs = mod(inputs)

        out = inputs
        E_mlp = out[calculator.energy_key].sum()
        q_pred = out["charges"] if "charges" in out else None

        # Runtime fallback: if OR inputs are present but the checkpoint does not actually
        # use them in its energy graph, autograd wrt OR coordinates will be None. In that
        # case, add the external Coulomb term explicitly from the predicted vacuum charges.
        if has_or and q_pred is not None:
            phi_or = self._coulomb_phi(
                r_qeq_A=r_mlp_A,
                r_or_A=r_or_A,
                q_or_e=q_or,
                cutoff_A=None,
            )
            
            E_elec_IRBR2OR = torch.sum(q_pred * phi_or)
            E_total = E_mlp + E_elec_IRBR2OR
        else:
            E_total = E_mlp

        grad_targets = [r_mlp_A]
        if has_or and compute_or_forces:
            grad_targets.append(r_or_A)

        grads = torch.autograd.grad(
            E_total,
            grad_targets,
            create_graph=False,
            retain_graph=False,
            allow_unused=True,
        )

        dE_dRqm = grads[0]
        if dE_dRqm is None:
            dE_dRqm = torch.zeros_like(r_mlp_A)
        F_qm = (-dE_dRqm).detach().cpu().numpy()

        if has_or and compute_or_forces:
            dE_dRor = grads[1]
            if dE_dRor is None:
                dE_dRor = torch.zeros_like(r_or_A)
            F_or = (-dE_dRor).detach().cpu().numpy()
        else:
            F_or = np.zeros((0, 3), dtype=float)

        q_np = None if q_pred is None else q_pred.detach().cpu().numpy()
        return float(E_total.detach().cpu().item()), F_qm, F_or, q_np, float(E_mlp.detach().cpu().item())

    def predict_energy_and_forces(
        self,
        system: ase.Atoms,
        r_mlp_A: torch.Tensor,
        r_or_A: torch.Tensor | None = None,
        q_or: torch.Tensor | None = None,
        compute_or_forces: bool = False,
    ):
        """
        Predict energy and forces for the given atomic system.

        Args:
            system (ase.Atoms): Atomic system for which predictions are required.

        Returns:
            tuple: Predicted energy (float) and forces (np.ndarray).
        """
        # Set the predictive calculator for the atomic system.
        system.calc = self.pred_calculator

        # Perform predictions.
        return self._energy_forces(
            system=system,
            calculator=self.pred_calculator,
            r_mlp_A=r_mlp_A,
            r_or_A=r_or_A,
            q_or=q_or,
            compute_or_forces=compute_or_forces,
        )
    
    def predict_partial_charges(self, system: ase.Atoms):
        """
        Predict partial charges for the given atomic system.

        Args:
            system (ase.Atoms): Atomic system for which predictions are required.

        Returns:
            np.ndarray: Partial charges.
        """
        # Set the predictive calculator for the atomic system.
        system.calc = self.pred_calculator

        # Perform predictions.
        partial_charges = system.get_charges()
        return partial_charges
    

    def validate_prediction(
        self,
        system: ase.Atoms,
        dynamic_charges: bool = False,
        r_qeq_A: torch.Tensor | None = None,
        r_or_A: torch.Tensor | None = None,
        q_or: torch.Tensor | None = None,
        cutoff_nm: float | None = None,
    ) -> list[float]:
        """
        Committee energy list.

        - legacy: model ASE energy (kJ/mol) (this is your E_mlp)
        - dynamic_charges (B2): returns E_mlp (NOT E_total)
        """
        if len(self.val_calculators) == 0:
            return []

        if not dynamic_charges:
            val_energies: list[float] = []
            for val_calculator in self.val_calculators:
                r_val_A = torch.tensor(system.get_positions(), dtype=torch.float32, device=self.torchdevice, requires_grad=True)
                e_kj, _, _, _, e_mlp = self._energy_forces(system=system, calculator=val_calculator, r_mlp_A=r_val_A)
                val_energies.append(float(e_kj))
            return val_energies

        assert r_qeq_A is not None and r_or_A is not None and q_or is not None

        val_energies: list[float] = []
        for val_calc in self.val_calculators:
            rq = r_qeq_A.detach().clone().requires_grad_(True)
            ro = r_or_A.detach().clone().requires_grad_(True)
            e_tot, _, _, _, e_mlp = self._energy_forces(
                system=system,
                calculator=val_calc,
                r_mlp_A=rq,
                r_or_A=ro,
                q_or=q_or,
                compute_or_forces=False,
            )
            val_energies.append(float(e_mlp))

        return val_energies
    
    def validate_prediction_maxForceDeviation(
        self,
        system: ase.Atoms,
        dynamic_charges: bool = False,
        r_qeq_A: torch.Tensor | None = None,
        r_or_A: torch.Tensor | None = None,
        q_or: torch.Tensor | None = None,
        cutoff_nm: float | None = None,
    ) -> float:
        """
        Committee max-force disagreement (max over atoms of sigmaF_alpha).

        - legacy: uses ASE forces (MLP forces)
        - dynamic_charges (B2): uses forces from E_mlp only (NOT E_total)
        """
        if len(self.val_calculators) == 0:
            return 0.0

        if not dynamic_charges:
            model_forces = [np.linalg.norm(self.forces, axis=1)]
            for val_calculator in self.val_calculators:
                r_val_A = torch.tensor(system.get_positions(), dtype=torch.float32, device=self.torchdevice, requires_grad=True)
                _, f_kj_A, _, _,_ = self._energy_forces(system=system, calculator=val_calculator, r_mlp_A=r_val_A)
                model_forces.append(np.linalg.norm(f_kj_A, axis=1))
        else:
            assert r_qeq_A is not None and r_or_A is not None and q_or is not None
            assert self.forces is not None

            # production model total forces (including external embedding)
            model_forces = [np.linalg.norm(self.forces, axis=1)]

            for val_calc in self.val_calculators:
                rq = r_qeq_A.detach().clone().requires_grad_(True)
                ro = r_or_A.detach().clone().requires_grad_(True)

                _, F_qm, _, _,_ = self._energy_forces(
                    system=system,
                    calculator=val_calc,
                    r_mlp_A=rq,
                    r_or_A=ro,
                    q_or=q_or,
                    compute_or_forces=False,
                )
                model_forces.append(np.linalg.norm(F_qm, axis=1))

        ensemble = np.mean(model_forces, axis=0)
        sigma = np.sqrt(np.mean([(f - ensemble) ** 2 for f in model_forces], axis=0))
        return float(np.max(sigma))


    def calculate_next_step(
        self,
        atomic_numbers: list,
        positions_A: list,
        time_step: int,
        dynamic_charges: bool = True,
        or_positions_nm: list | None = None,
        or_charges_e: list | None = None,
        cutoff_nm: float | None = None,
        or_atomic_numbers: list | None = None,
        n_caps: int = 0) -> None:
        """
        Unified next-step:
        - dynamic_charges=False: classic SchNetPack ASE energy+forces
        - dynamic_charges=True : B2 energy+forces (E_mlp + q·phi) + OR forces
        Committee validation uses the SAME definition in each regime.
        """
        system = ase.Atoms(numbers=atomic_numbers, positions=positions_A)
        system.info["total_charge"] = float(self.total_charge)
        system.info["spin_multiplicity"] = float(self.spin_multiplicity)

        self.time_step = time_step
        self.or_forces = np.zeros((0, 3), dtype=float)

        # -------------------------
        # Legacy mode
        # -------------------------
        if not dynamic_charges:
            r_mlp_A = torch.tensor(positions_A, dtype=torch.float32, device=self.torchdevice, requires_grad=True)
            self.energy, self.forces, _, outputs = self.predict_energy_and_forces(system=system, r_mlp_A=r_mlp_A)

            if self.has_partial_charges:
                if "charges" in outputs:
                    self.charges = outputs["charges"].detach().cpu().numpy()
                else:
                    self.charges = self.predict_partial_charges(system=system)

            if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
                val_energies = self.validate_prediction(system=system, dynamic_charges=False)

                self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(
                    system=system,
                    dynamic_charges=False,
                )

                if len(val_energies) == 1:
                    self.nn_valid_ene = (self.energy - val_energies[0]) / np.sqrt(2)
                else:
                    all_E = np.array(val_energies + [self.energy], dtype=float)
                    self.nn_valid_ene = all_E.std(ddof=1)

            return None

        # -------------------------
        # Dynamic charges / B2 mode
        # -------------------------
        if or_positions_nm is None or or_charges_e is None or cutoff_nm is None:
            raise ValueError("dynamic_charges=True requires or_positions_nm, or_charges_e, cutoff_nm")

        r_qeq_A = torch.tensor(positions_A, dtype=torch.float32, device=self.torchdevice, requires_grad=True)

        if len(or_positions_nm) > 0:
            r_or_A = torch.tensor(or_positions_nm, dtype=torch.float32, device=self.torchdevice, requires_grad=True) * 10.0
            q_or = torch.tensor(or_charges_e, dtype=torch.float32, device=self.torchdevice)
            z_or = torch.tensor(or_atomic_numbers, dtype=torch.long, device=self.torchdevice)
        else:
            r_or_A = torch.zeros((0, 3), dtype=torch.float32, device=self.torchdevice, requires_grad=True)
            q_or = torch.zeros((0,), dtype=torch.float32, device=self.torchdevice)

        # production model with ExternalCoulombEmbedding at inference
        E_total, F_total_qeq, F_total_or, outputs, E_MLP = self.predict_energy_and_forces(
            system=system,
            r_mlp_A=r_qeq_A,
            r_or_A=r_or_A,
            q_or=q_or,
            compute_or_forces=True,
        )

        self.energy = float(E_total)
        self.forces = F_total_qeq
        self.or_forces = F_total_or

        if "charges" in outputs:
            self.charges = outputs["charges"].detach().cpu().numpy()
        else:
            self.charges = np.zeros((len(positions_A),), dtype=float)

        # committee validation (ExternalCoulombEmbedding-consistent)
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(
                system=system,
                dynamic_charges=True,
                r_qeq_A=r_qeq_A,
                r_or_A=r_or_A,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
            )

            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(
                system=system,
                dynamic_charges=True,
                r_qeq_A=r_qeq_A,
                r_or_A=r_or_A,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
            )
            
            #print("validation energies: ", val_energies)
            #print("E_MLP: ", E_MLP)

            if len(val_energies) == 1:
                self.nn_valid_ene = (E_MLP - val_energies[0]) / np.sqrt(2)
            else:
                all_E = np.array(val_energies + [E_MLP], dtype=float)
                self.nn_valid_ene = all_E.std(ddof=1)

        return None

    def _coulomb_phi(
        self,
        r_qeq_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or_e: torch.Tensor,
        cutoff_A: float | None = None,
        make_zero: bool = False,):
        """
        Screened electrostatic potential at QEq atoms due to OR static charges.

        If sigma_qeq_A is None:
            falls back to bare Coulomb
                phi_i = k_e * sum_j q_j / r_ij

        If sigma_qeq_A is provided:
            uses Gaussian-screened Coulomb
                phi_i = k_e * sum_j q_j * erf(r_ij / sigma_ij) / r_ij
            with
                sigma_ij = sqrt(sigma_i^2 + sigma_or_j^2)

        Units:
            r in Ang
            q in e
            returns phi in eV
        """
        if r_or_A.shape[0] == 0:
            return torch.zeros(
                r_qeq_A.shape[0],
                device=r_qeq_A.device,
                dtype=r_qeq_A.dtype,
            )

        diff = r_qeq_A[:, None, :] - r_or_A[None, :, :]
        dist = torch.norm(diff, dim=-1).clamp_min(1e-8)   # (Nq, Nor)

        # optional cutoff
        mask = None
        if cutoff_A is not None:
            mask = (dist <= float(cutoff_A))

        # bare Coulomb fallback
        kernel = 1.0 / dist

        if mask is not None:
            kernel = kernel * mask

        phi = K_E_kJMOL_ANG_E2 * torch.sum(kernel * q_or_e[None, :], dim=1)

        if make_zero:
            phi = torch.zeros_like(phi)

        return phi

    def get_or_forces(self):
        return self.or_forces

    def get_energy(self):
        """
        Get the predicted energy for the last atomic system.

        Returns:
            float: Predicted energy.
        """
        return self.energy

    def get_forces(self):
        """
        Get the predicted forces for the last atomic system.

        Returns:
            np.ndarray: Predicted forces.
        """
        return self.forces

    def get_nn_valid_ene(self):
        """
        Get the validation energy deviation calculated during the last validation step.

        Returns:
            float: Validation energy deviation.
        """
        return self.nn_valid_ene

    def get_nn_valid_maxF(self):
        """
        Get the maximum force committee disagreement among all atom during the last validation step.

        Returns:
            float: Validation force deviation.
        """
        return self.nn_valid_maxF
     
    def get_charges(self):
        """
        Get the predicted partial charges for the last atomic system.

        Returns:
            np.ndarray: Predicted partial charges.
        """
        return self.charges

class Pert_SchNet_V2_Calculator(SchNet_V2_Calculator):
    """Thermodynamic integration / alchemical perturbation wrapper.

    This class evaluates two endstates (A and B) and interpolates them with λ.

    TI (linear endstate interpolation)
    ---------------------------------
    - Energy:
        ``E(λ) = (1-λ) E_A + λ E_B``
    - Forces:
        ``F(λ) = (1-λ) F_A + λ F_B``
    - Derivative:
        ``dE/dλ = E_B - E_A``

    Dual topology / state labels
    ----------------------------
    The per-atom state label array ``perturbed_qm_states`` is defined for the **IR block**
    (i.e., for the QM zone atoms only; BR and caps are handled separately):

    - ``0``: present/active in both endstates
    - ``1``: active only in A
    - ``2``: active only in B

    For each endstate we build two auxiliary systems:
    - **vac**: active IR subset only (no BR, no caps)
    - **burnn**: active IR subset + all BR + caps

    We then compute a corrected energy:

        ``E_corr = E_burnn - E_vac``

    If ``reference_energies`` are provided (default in this file), we use those instead of
    evaluating ``E_vac`` on the fly.

    Ordering contract
    -----------------
    Input ordering from C++ is:

        ``[ IR ][ BR ][ caps ]``

    burnn ordering (constructed here) is:

        ``[ active IR ][ all BR ][ caps ]``

    To return forces/charges back to the driver in the original ordering, we *scatter*
    burnn outputs back into the full-size arrays, leaving inactive IR atoms as zeros.
    """

    def __init__(
        self,
        model_path: str,
        val_model_paths: list,
        nn_valid_freq: int,
        write_energy_freq: int,
        spin_multiplicity: int,
        total_charge: int,
        lam: float,
        perturbed_qm_states: list,
        perturbed_total_charge: int,
        perturbed_spin_multiplicity: int,
        ref_vacA=None,
        ref_vacB=None,
    ) -> None:
        super().__init__(
            model_path,
            val_model_paths,
            nn_valid_freq,
            write_energy_freq,
            spin_multiplicity,
            total_charge,
        )

        self.lam = float(lam)

        self.ref_vacA = ref_vacA
        self.ref_vacB = ref_vacB
        self.reference_energies = (
            {"A": float(ref_vacA), "B": float(ref_vacB)}
            if ref_vacA is not None and ref_vacB is not None
            else None
        )

        self.perturbed_qm_states = np.asarray(perturbed_qm_states, dtype=int)
        self.qmzone_size = len(self.perturbed_qm_states)

        self.charge_A = int(total_charge)
        self.spin_A = int(spin_multiplicity)
        self.charge_B = int(perturbed_total_charge)
        self.spin_B = int(perturbed_spin_multiplicity)

        self.idx_A = [i for i, s in enumerate(self.perturbed_qm_states) if s in (0, 1)]
        self.idx_B = [i for i, s in enumerate(self.perturbed_qm_states) if s in (0, 2)]

        if len(self.idx_A) == 0 or len(self.idx_B) == 0:
            raise ValueError("Endstate A or B has zero active QM atoms.")

        self.pred_calculators = {
            "A": self.pred_calculator,
            "B": self.get_calculator(Path(model_path)),
        }

        self.val_calculators_A = self.val_calculators
        self.val_calculators_B = [self.get_calculator(Path(p)) for p in val_model_paths]

        self.energy = None
        self.forces = None
        self.or_forces = np.zeros((0, 3), dtype=float)

        self.derivative = None
        self.derivative_mlp = None
        self.derivative_elec = None

        self.q_A = None
        self.q_B = None
        self.charges = None

        self.nn_valid_ene = 0.0
        self.nn_valid_maxF = 0.0

    def _active_ir_indices(self, endstate: str):
        if endstate == "A":
            return self.idx_A, self.charge_A, self.spin_A
        if endstate == "B":
            return self.idx_B, self.charge_B, self.spin_B
        raise ValueError(f"Unknown endstate: {endstate}")

    def _build_endstate_systems(self, atomic_numbers, positions, endstate: str, n_caps: int):
        n_total = len(atomic_numbers)
        n_caps = int(n_caps)

        if n_caps < 0 or n_caps > n_total:
            raise RuntimeError(f"Invalid n_caps={n_caps} for n_total={n_total}")

        n_ir = self.qmzone_size
        cap_start = n_total - n_caps

        if cap_start < n_ir:
            raise RuntimeError(
                f"Ordering mismatch: cap_start={cap_start} < n_ir={n_ir}."
            )

        ir_idx, charge, spin = self._active_ir_indices(endstate)
        br_idx = list(range(n_ir, cap_start))
        qmbr_idx = ir_idx + br_idx

        vac = ase.Atoms(
            numbers=[atomic_numbers[i] for i in ir_idx],
            positions=[positions[i] for i in ir_idx],
        )
        vac.info["total_charge"] = float(charge)
        vac.info["spin_multiplicity"] = float(spin)

        burnn = ase.Atoms(
            numbers=[atomic_numbers[i] for i in qmbr_idx] + list(atomic_numbers[cap_start:]),
            positions=[positions[i] for i in qmbr_idx] + list(positions[cap_start:]),
        )
        burnn.info["total_charge"] = float(charge)
        burnn.info["spin_multiplicity"] = float(spin)

        return vac, burnn, ir_idx, br_idx, cap_start

    def _scatter_burnn_vector_to_full(self, vec_burnn, ir_idx, br_idx, cap_start, n_total):
        vec_burnn = np.asarray(vec_burnn)

        n_ir_active = len(ir_idx)
        n_br = len(br_idx)
        n_caps = n_total - int(cap_start)
        expected = n_ir_active + n_br + n_caps

        if vec_burnn.shape[0] != expected:
            raise RuntimeError(
                f"burnn vector length mismatch: got {vec_burnn.shape[0]}, "
                f"expected {expected}."
            )

        out = np.zeros((n_total,) + vec_burnn.shape[1:], dtype=vec_burnn.dtype)

        out_ir = vec_burnn[:n_ir_active]
        out_br = vec_burnn[n_ir_active:n_ir_active + n_br]
        out_caps = vec_burnn[n_ir_active + n_br:]

        for local_i, global_i in enumerate(ir_idx):
            out[global_i] = out_ir[local_i]

        for local_i, global_i in enumerate(br_idx):
            out[global_i] = out_br[local_i]

        out[cap_start:] = out_caps

        return out

    def _eval_with_calculator(
        self,
        system,
        calculator,
        r_or_A=None,
        q_or=None,
        compute_or_forces=False,
    ):
        r_mlp_A = torch.tensor(
            system.get_positions(),
            dtype=torch.float32,
            device=self.torchdevice,
            requires_grad=True,
        )

        return self._energy_forces(
            calculator=calculator,
            system=system,
            r_mlp_A=r_mlp_A,
            r_or_A=r_or_A,
            q_or=q_or,
            compute_or_forces=compute_or_forces,
        )

    def _evaluate_corrected_endstate(
        self,
        atomic_numbers,
        positions,
        endstate: str,
        n_caps: int,
        dynamic_charges: bool = False,
        r_or_A=None,
        q_or=None,
    ):
        n_total = len(atomic_numbers)

        vac, burnn, ir_idx, br_idx, cap_start = self._build_endstate_systems(
            atomic_numbers, positions, endstate, n_caps
        )

        calc = self.pred_calculators[endstate]

        E_burnn_total, F_burnn, F_or, q_burnn, E_burnn_mlp = self._eval_with_calculator(
            burnn,
            calc,
            r_or_A=r_or_A if dynamic_charges else None,
            q_or=q_or if dynamic_charges else None,
            compute_or_forces=dynamic_charges,
        )

        if self.reference_energies is None:
            E_vac_total, F_vac, _, _, E_vac_mlp = self._eval_with_calculator(
                vac,
                calc,
                r_or_A=None,
                q_or=None,
                compute_or_forces=False,
            )
            E_vac = E_vac_total
        else:
            E_vac = float(self.reference_energies[endstate])
            E_vac_mlp = E_vac
            F_vac = np.zeros((len(vac), 3), dtype=float)

        E_elec = float(E_burnn_total - E_burnn_mlp)
        E_corr_total = float(E_burnn_total - E_vac)
        E_corr_mlp = float(E_burnn_mlp - E_vac_mlp)

        n_ir_active = len(ir_idx)
        n_br = len(br_idx)

        F_ir_corr = F_burnn[:n_ir_active] - F_vac
        F_br = F_burnn[n_ir_active:n_ir_active + n_br]
        F_caps = F_burnn[n_ir_active + n_br:]

        F_compact = np.concatenate([F_ir_corr, F_br, F_caps], axis=0)

        F_full = np.stack(
            [
                self._scatter_burnn_vector_to_full(
                    F_compact[:, k],
                    ir_idx,
                    br_idx,
                    cap_start,
                    n_total,
                )
                for k in range(3)
            ],
            axis=1,
        )

        q_full = None
        if q_burnn is not None:
            q_full = self._scatter_burnn_vector_to_full(
                q_burnn,
                ir_idx,
                br_idx,
                cap_start,
                n_total,
            ).reshape((n_total,))

        return {
            "E_total": E_corr_total,
            "E_mlp": E_corr_mlp,
            "E_elec": E_elec,
            "F_full": np.asarray(F_full, dtype=float),
            "F_or": np.asarray(F_or, dtype=float),
            "q_full": None if q_full is None else np.asarray(q_full, dtype=float),
        }

    def calculate_next_step(
        self,
        atomic_numbers,
        positions,
        time_step: int,
        n_caps: int = 0,
        dynamic_charges: bool = False,
        or_positions_nm=None,
        or_charges_e=None,
        cutoff_nm=None,
        or_atomic_numbers=None,
    ) -> None:
        lam = float(self.lam)

        r_or_A = None
        q_or = None

        if dynamic_charges:
            if or_positions_nm is None or or_charges_e is None:
                raise ValueError(
                    "dynamic_charges=True requires or_positions_nm and or_charges_e."
                )

            if len(or_positions_nm) > 0:
                r_or_A = (
                    torch.tensor(
                        or_positions_nm,
                        dtype=torch.float32,
                        device=self.torchdevice,
                        requires_grad=True,
                    )
                    * 10.0
                )
                q_or = torch.tensor(
                    or_charges_e,
                    dtype=torch.float32,
                    device=self.torchdevice,
                )
            else:
                r_or_A = torch.zeros(
                    (0, 3),
                    dtype=torch.float32,
                    device=self.torchdevice,
                    requires_grad=True,
                )
                q_or = torch.zeros(
                    (0,),
                    dtype=torch.float32,
                    device=self.torchdevice,
                )

        A = self._evaluate_corrected_endstate(
            atomic_numbers,
            positions,
            "A",
            n_caps,
            dynamic_charges=dynamic_charges,
            r_or_A=r_or_A,
            q_or=q_or,
        )

        B = self._evaluate_corrected_endstate(
            atomic_numbers,
            positions,
            "B",
            n_caps,
            dynamic_charges=dynamic_charges,
            r_or_A=r_or_A,
            q_or=q_or,
        )

        self.energy = (1.0 - lam) * A["E_total"] + lam * B["E_total"]
        self.forces = (1.0 - lam) * A["F_full"] + lam * B["F_full"]

        self.derivative_mlp = B["E_mlp"] - A["E_mlp"]
        self.derivative_elec = B["E_elec"] - A["E_elec"]
        self.derivative = self.derivative_mlp + self.derivative_elec
        if self.debug:
            print("\n[DEBUG] ===== Perturbation Step =====")
            print(f"[DEBUG] lambda = {lam: .4f}")

            if self.debug_level >= 1:
                print(f"[DEBUG] E_A_total = {A['E_total']: .6f}")
                print(f"[DEBUG] E_B_total = {B['E_total']: .6f}")
                print(f"[DEBUG] E_A_mlp   = {A['E_mlp']: .6f}")
                print(f"[DEBUG] E_B_mlp   = {B['E_mlp']: .6f}")
                print(f"[DEBUG] E_A_elec  = {A['E_elec']: .6f}")
                print(f"[DEBUG] E_B_elec  = {B['E_elec']: .6f}")

                print(f"[DEBUG] dE/dλ (total) = {self.derivative: .6f}")
                print(f"[DEBUG] dE/dλ (MLP)   = {self.derivative_mlp: .6f}")
                print(f"[DEBUG] dE/dλ (elec)  = {self.derivative_elec: .6f}")

            if self.debug_level >= 2:
                if A["q_full"] is not None:
                    print("[DEBUG] Charges A:")
                    print(A["q_full"])
                    print(f"[DEBUG] sum(q_A) = {np.sum(A['q_full']): .6f}")

                if B["q_full"] is not None:
                    print("[DEBUG] Charges B:")
                    print(B["q_full"])
                    print(f"[DEBUG] sum(q_B) = {np.sum(B['q_full']): .6f}")

            if self.debug_level >= 3:
                print("[DEBUG] Interpolated charges:")
                print(self.charges)
                print(f"[DEBUG] sum(q_interp) = {np.sum(self.charges): .6f}")
                
        if dynamic_charges:
            self.or_forces = (1.0 - lam) * A["F_or"] + lam * B["F_or"]
        else:
            self.or_forces = np.zeros((0, 3), dtype=float)

        self.q_A = A["q_full"]
        self.q_B = B["q_full"]

        if self.q_A is not None and self.q_B is not None:
            self.charges = (1.0 - lam) * self.q_A + lam * self.q_B
        elif self.q_A is not None:
            self.charges = self.q_A
        elif self.q_B is not None:
            self.charges = self.q_B
        else:
            self.charges = np.zeros((len(atomic_numbers),), dtype=float)

        if len(self.val_calculators_A) > 0 and time_step % self.nn_valid_freq == 0:
            E_vals = []

            for val_A, val_B in zip(self.val_calculators_A, self.val_calculators_B):
                A_val = self._evaluate_corrected_endstate_with_calculator(
                    atomic_numbers,
                    positions,
                    "A",
                    val_A,
                    n_caps,
                    dynamic_charges=dynamic_charges,
                    r_or_A=r_or_A,
                    q_or=q_or,
                )

                B_val = self._evaluate_corrected_endstate_with_calculator(
                    atomic_numbers,
                    positions,
                    "B",
                    val_B,
                    n_caps,
                    dynamic_charges=dynamic_charges,
                    r_or_A=r_or_A,
                    q_or=q_or,
                )

                E_vals.append((1.0 - lam) * A_val + lam * B_val)

            if len(E_vals) == 1:
                self.nn_valid_ene = abs(self.energy - E_vals[0]) / np.sqrt(2.0)
            elif len(E_vals) > 1:
                self.nn_valid_ene = np.array(E_vals + [self.energy], dtype=float).std(ddof=1)

        return None

    def _evaluate_corrected_endstate_with_calculator(
        self,
        atomic_numbers,
        positions,
        endstate: str,
        val_calculator,
        n_caps: int,
        dynamic_charges: bool = False,
        r_or_A=None,
        q_or=None,
    ) -> float:
        vac, burnn, _, _, _ = self._build_endstate_systems(
            atomic_numbers,
            positions,
            endstate,
            n_caps,
        )

        E_burnn_total, _, _, _, _ = self._eval_with_calculator(
            burnn,
            val_calculator,
            r_or_A=r_or_A if dynamic_charges else None,
            q_or=q_or if dynamic_charges else None,
            compute_or_forces=False,
        )

        if self.reference_energies is None:
            E_vac_total, _, _, _, _ = self._eval_with_calculator(
                vac,
                val_calculator,
                r_or_A=None,
                q_or=None,
                compute_or_forces=False,
            )
            E_vac = E_vac_total
        else:
            E_vac = float(self.reference_energies[endstate])

        return float(E_burnn_total - E_vac)

    def get_derivative(self) -> float:
        return float(self.derivative)

    def get_mlp_derivative(self) -> float:
        return float(self.derivative_mlp)

    def get_electrostatic_derivative(self) -> float:
        return float(self.derivative_elec)

    def get_charges_A(self):
        if self.q_A is None:
            raise AttributeError("q_A not available.")
        return np.asarray(self.q_A)

    def get_charges_B(self):
        if self.q_B is None:
            raise AttributeError("q_B not available.")
        return np.asarray(self.q_B)

    def get_or_forces(self):
        return self.or_forces