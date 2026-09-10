import numpy as np
from pathlib import Path
import fcntl
import os
import torch
import yaml
import ase
import schnetpack as spk

EV_TO_KJMOL = 96.4853321233  
# Coulomb prefactor in the same distance unit used by the MLP.
K_E_KJMOL_A_E2 = 1389.35456  # kJ/mol·Å·e^-2

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
        if "phi_static" in atoms.arrays:
            phi = atoms.arrays["phi_static"]
            phi = torch.tensor(phi, dtype=torch.float32, device=self.device)
            # ensure shape (N, 1) or (N,) depending on what your model expects
            inputs["phi_static"] = phi
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
                 spin_multiplicity:int, total_charge:int,
                 electrostatic_damping: str = "soft",
                 electrostatic_sigma_A: float = 0.01,
                 qeq_mode: str = "vacuum",
                 ntwse: int = 0,
                 val_thresh: float = 0.0) -> None:
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

        sigma_from_environment = os.environ.get(
            "SCHNET_ELECTROSTATIC_SIGMA_A"
        )
        if sigma_from_environment is not None:
            try:
                electrostatic_sigma_A = float(sigma_from_environment)
            except ValueError as error:
                raise ValueError(
                    "SCHNET_ELECTROSTATIC_SIGMA_A must be a number in angstrom"
                ) from error

        if electrostatic_damping not in {"soft", "erf", "none"}:
            raise ValueError("electrostatic_damping must be 'soft', 'erf', or 'none'")
        if not np.isfinite(electrostatic_sigma_A) or electrostatic_sigma_A <= 0.0:
            raise ValueError(
                "electrostatic_sigma_A must be a finite positive number in angstrom"
            )
        self.electrostatic_damping = electrostatic_damping
        self.electrostatic_sigma_A = float(electrostatic_sigma_A)
        qeq_mode = os.environ.get("SCHNET_QEQ_MODE", qeq_mode).strip().lower()
        qeq_mode_aliases = {
            "vac": "vacuum",
            "vacqeq": "vacuum",
            "vacuum": "vacuum",
            "pol": "polarized",
            "polqeq": "polarized",
            "polarized": "polarized",
        }
        if qeq_mode not in qeq_mode_aliases:
            raise ValueError(
                "qeq_mode/SCHNET_QEQ_MODE must be vacuum (VacQEq) or "
                "polarized (PolQEq)"
            )
        self.qeq_mode = qeq_mode_aliases[qeq_mode]
        print(
            "B2ELEC "
            f"damping={self.electrostatic_damping} "
            f"sigma_A={self.electrostatic_sigma_A:.8g} "
            f"qeq_mode={self.qeq_mode}",
            flush=True,
        )
        self.ntwse = int(ntwse)
        self.val_thresh = float(val_thresh)

        phi_output_file = os.environ.get("SCHNET_PHI_STATIC_FILE")
        self.phi_static_output_file = Path(phi_output_file) if phi_output_file else None
        if self.phi_static_output_file is not None:
            self.phi_static_output_file.parent.mkdir(parents=True, exist_ok=True)

    
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
            energy_unit="kJ/mol",
            position_unit="Angstrom",
            device= self.torchdevice, # device for computation
        )
        
        return calculator
        
    def _energy_forces_kjmol(self, system: ase.Atoms, calculator) -> tuple[float, np.ndarray]:
        system.calc = calculator
        energy_ev = system.get_potential_energy()      # ASE-native
        forces_ev_A = system.get_forces()              # ASE-native
        energy_kj = energy_ev * EV_TO_KJMOL
        forces_kj_A = forces_ev_A * EV_TO_KJMOL
        return energy_kj, forces_kj_A

    def predict_energy_and_forces(self, system: ase.Atoms):
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
        return self._energy_forces_kjmol(system=system, calculator=self.pred_calculator)
    
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
        n_link_atoms: int = 0,
    ) -> list[float]:
        """
        Committee energy list.

        - legacy: model ASE energy (kJ/mol) (this is your E_mlp)
        - dynamic_charges (B2): returns E_mlp only. Embedding-energy committee
          values are retained separately as an opt-in diagnostic.
        """
        self._last_val_embedding_energies = []
        if len(self.val_calculators) == 0:
            return []

        if not dynamic_charges:
            val_energies: list[float] = []
            for val_calculator in self.val_calculators:
                e_kj, _ = self._energy_forces_kjmol(system=system, calculator=val_calculator)
                val_energies.append(float(e_kj))
            return val_energies

        assert r_qeq_A is not None and r_or_A is not None and q_or is not None and cutoff_nm is not None

        val_energies: list[float] = []
        for val_calc in self.val_calculators:
            rq = r_qeq_A.detach().clone().requires_grad_(True)
            ro = r_or_A.detach().clone().requires_grad_(True)

            E_total, _, _, _, E_mlp, _ = self._b2_eval_model(
                calculator=val_calc,
                system=system,
                r_qeq_A=rq,
                r_or_A=ro,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
            )
            val_energies.append(float(E_mlp))
            self._last_val_embedding_energies.append(float(E_total - E_mlp))

        return val_energies
    
    def validate_prediction_maxForceDeviation(
        self,
        system: ase.Atoms,
        dynamic_charges: bool = False,
        r_qeq_A: torch.Tensor | None = None,
        r_or_A: torch.Tensor | None = None,
        q_or: torch.Tensor | None = None,
        cutoff_nm: float | None = None,
        n_link_atoms: int = 0,
    ) -> float:
        """
        Committee max-force disagreement (max over atoms of sigmaF_alpha).

        - legacy: uses ASE forces (MLP forces)
        - dynamic_charges (B2): uses MLP-only forces on QEq sites
        """
        if len(self.val_calculators) == 0:
            return 0.0

        if not dynamic_charges:
            model_forces = [np.linalg.norm(self.forces, axis=1)]
            for val_calculator in self.val_calculators:
                _, f_kj_A = self._energy_forces_kjmol(system=system, calculator=val_calculator)
                model_forces.append(np.linalg.norm(f_kj_A, axis=1))
        else:
            assert r_qeq_A is not None and r_or_A is not None and q_or is not None and cutoff_nm is not None
            model_forces = [np.linalg.norm(self.forces_mlp, axis=1)]

            for val_calc in self.val_calculators:
                rq = r_qeq_A.detach().clone().requires_grad_(True)
                ro = r_or_A.detach().clone().requires_grad_(True)

                _, _, _, _, _, F_mlp_qeq = self._b2_eval_model(
                    calculator=val_calc,
                    system=system,
                    r_qeq_A=rq,
                    r_or_A=ro,
                    q_or=q_or,
                    cutoff_nm=cutoff_nm,
                    n_link_atoms=n_link_atoms,
                )
                model_forces.append(np.linalg.norm(F_mlp_qeq, axis=1))

        ensemble = np.mean(model_forces, axis=0)
        sigma = np.sqrt(np.mean([(f - ensemble) ** 2 for f in model_forces], axis=0))
        return float(np.max(sigma))
    
    def _b2_eval_model(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_qeq_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or: torch.Tensor,
        cutoff_nm: float,
        n_link_atoms: int = 0,
        store_diagnostics: bool = False,
        energy_components: dict | None = None,
    ):
        """
        Evaluate one model in B2.
        Returns both:
        - total: E_total = E_mlp + 0.5*(q0 + q_phi) dot phi and its forces
        - mlp-only: E_mlp and forces from E_mlp

        Returns:
        E_total_kJmol: float
        F_total_qeq_kJmol_A: np.ndarray (Nq,3)
        F_total_or_kJmol_A : np.ndarray (Nor,3)
        q_e               : np.ndarray (Nq,)
        E_mlp_kJmol       : float
        F_mlp_qeq_kJmol_A : np.ndarray (Nq,3)
        """
        inputs = calculator.converter(system)
        inputs[spk.properties.R] = r_qeq_A

        # GROMOS already supplies the pairlisted OR atoms. No second distance
        # cutoff is applied here.
        phi = self._coulomb_phi_kjmol(r_qeq_A, r_or_A, q_or)
        phi = self._zero_link_atom_phi(phi, n_link_atoms)
        inputs["phi_static"] = phi
        inputs["qeq_polarize"] = torch.tensor(
            getattr(self, "qeq_mode", "vacuum") == "polarized",
            device=phi.device,
        )

        # Forward without force modules
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

        # This manual B2 path bypasses AtomisticModel.forward(), so reproduce
        # its saved postprocessing step (including model-specific offsets).
        if getattr(model, "do_postprocessing", False):
            for postprocessor in getattr(model, "postprocessors", []):
                inputs = postprocessor(inputs)

        out = inputs

        # --- energies ---
        E_mlp = out[calculator.energy_key].sum()   # kJ/mol (per your SpkCalculator config)
        q0 = out["charges_vac"]                    # vacuum QEq charges, e
        q_phi = out["charges"]                     # field-polarized charges, e
        E_embedding = out["qeq_embedding_energy"].sum()  # kJ/mol
        E_total = E_mlp + E_embedding

        # Optional scalar diagnostics; do not retain or modify the force graph.
        if energy_components is not None:
            with torch.no_grad():
                energy_components.update(
                    burnn=float(E_mlp.detach().cpu().item()),
                    elecstatic=float(torch.sum(q0 * phi).cpu().item()),
                    elecinduced=float((0.5 * torch.sum((q_phi - q0) * phi)).cpu().item()),
                )

        # Fixed-vacuum-charge reference used only for force decomposition.
        # The physical production force below differentiates the full
        # variational embedding energy without detaching either charge solution.
        if store_diagnostics:
            E_static_direct = torch.sum(q0.detach() * phi)
            dEdirect_dRqeq, dEdirect_dRor = torch.autograd.grad(
                E_mlp + E_static_direct,
                [r_qeq_A, r_or_A],
                create_graph=False,
                retain_graph=True,
                allow_unused=True,
            )
            if dEdirect_dRqeq is None:
                dEdirect_dRqeq = torch.zeros_like(r_qeq_A)
            if dEdirect_dRor is None:
                dEdirect_dRor = torch.zeros_like(r_or_A)

        # --- forces from E_total (what you integrate) ---
        dEt_dRqeq, dEt_dRor = torch.autograd.grad(
            E_total,
            [r_qeq_A, r_or_A],
            create_graph=False,
            retain_graph=True,   # keep graph for E_mlp grad below
            allow_unused=True,
        )
        if dEt_dRqeq is None:
            dEt_dRqeq = torch.zeros_like(r_qeq_A)
        if dEt_dRor is None:
            dEt_dRor = torch.zeros_like(r_or_A)

        F_total_qeq = (-dEt_dRqeq).detach().cpu().numpy()
        F_total_or  = (-dEt_dRor).detach().cpu().numpy()

        # --- forces from E_mlp only (committee metric) ---
        dEm_dRqeq = torch.autograd.grad(
            E_mlp,
            r_qeq_A,
            create_graph=False,
            retain_graph=False,
            allow_unused=True,
        )[0]
        if dEm_dRqeq is None:
            dEm_dRqeq = torch.zeros_like(r_qeq_A)

        F_mlp_qeq = (-dEm_dRqeq).detach().cpu().numpy()

        if store_diagnostics:
            self._last_phi_static = phi.detach().cpu().numpy().astype(np.float32, copy=False)
            F_direct_qeq = (-dEdirect_dRqeq).detach().cpu().numpy()
            F_direct_or = (-dEdirect_dRor).detach().cpu().numpy()
            F_response_qeq = F_total_qeq - F_direct_qeq
            F_response_or = F_total_or - F_direct_or
            self._last_b2_diag = self._collect_b2_diagnostics(
                E_total=E_total,
                E_mlp=E_mlp,
                E_embedding=E_embedding,
                qeq_mode=getattr(self, "qeq_mode", "vacuum"),
                q0=q0,
                q_phi=q_phi,
                phi=phi,
                r_qeq_A=r_qeq_A,
                r_or_A=r_or_A,
                F_total_qeq=F_total_qeq,
                F_total_or=F_total_or,
                F_direct_qeq=F_direct_qeq,
                F_direct_or=F_direct_or,
                F_response_qeq=F_response_qeq,
                F_response_or=F_response_or,
                F_mlp_qeq=F_mlp_qeq,
                n_link_atoms=n_link_atoms,
            )

        return (
            float(E_total.detach().cpu().item()),
            F_total_qeq,
            F_total_or,
            q_phi.detach().cpu().numpy(),
            float(E_mlp.detach().cpu().item()),
            F_mlp_qeq,
        )

    def calculate_next_step(
        self,
        atomic_numbers: list,
        positions_A: list,
        time_step: int,
        dynamic_charges: bool = False,
        or_positions_nm: list | None = None,
        or_charges_e: list | None = None,
        cutoff_nm: float | None = None,
        n_link_atoms: int = 0,
    ) -> None:
        """
        Unified next-step:
        - dynamic_charges=False: classic SchNetPack ASE energy+forces
        - dynamic_charges=True : variational explicit-QEq energy+forces
          (E_mlp + 0.5*(q0 + q_phi) dot phi) + OR forces
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
            self.energy, self.forces = self.predict_energy_and_forces(system=system)

            if self.has_partial_charges:
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
                self.nn_valid_mlp = self.nn_valid_ene
                self.nn_valid_mlp_maxF = self.nn_valid_maxF

            return None

        # -------------------------
        # Dynamic charges / B2 mode
        # -------------------------
        if or_positions_nm is None or or_charges_e is None or cutoff_nm is None:
            raise ValueError("dynamic_charges=True requires or_positions_nm, or_charges_e, cutoff_nm")

        r_qeq_A = torch.tensor(positions_A, dtype=torch.float32, device=self.torchdevice, requires_grad=True)

        if len(or_positions_nm) > 0:
            r_or_A = torch.tensor(np.asarray(or_positions_nm) * 10.0, dtype=torch.float32, device=self.torchdevice, requires_grad=True)
            q_or = torch.tensor(or_charges_e, dtype=torch.float32, device=self.torchdevice)
        else:
            r_or_A = torch.zeros((0, 3), dtype=torch.float32, device=self.torchdevice, requires_grad=True)
            q_or = torch.zeros((0,), dtype=torch.float32, device=self.torchdevice)

        # production model B2
        E_total, F_total_qeq, F_total_or, q_qeq, E_mlp, F_mlp_qeq = self._b2_eval_model(
            calculator=self.pred_calculator,
            system=system,
            r_qeq_A=r_qeq_A,
            r_or_A=r_or_A,
            q_or=q_or,
            cutoff_nm=cutoff_nm,
            n_link_atoms=n_link_atoms,
            store_diagnostics=True,
        )

        # what MD uses
        self.energy = float(E_total)
        self.forces = F_total_qeq
        self.or_forces = F_total_or
        self.charges = q_qeq

        # what validation uses
        self.energy_mlp = float(E_mlp)
        self.forces_mlp = F_mlp_qeq

        if time_step % self.write_energy_freq == 0:
            self._print_b2_diagnostics(time_step)

        fd_freq = int(os.environ.get("SCHNET_B2_FD_CHECK_FREQ", "0"))
        if fd_freq > 0 and time_step % fd_freq == 0:
            self._print_b2_finite_difference_check(
                calculator=self.pred_calculator,
                system=system,
                r_qeq_A=r_qeq_A,
                r_or_A=r_or_A,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
                F_total_qeq=F_total_qeq,
                F_total_or=F_total_or,
            )

        # committee validation (B2-consistent)
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            val_energies = self.validate_prediction(
                system=system,
                dynamic_charges=True,
                r_qeq_A=r_qeq_A,
                r_or_A=r_or_A,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
            )

            self.nn_valid_maxF = self.validate_prediction_maxForceDeviation(
                system=system,
                dynamic_charges=True,
                r_qeq_A=r_qeq_A,
                r_or_A=r_or_A,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
            )
            
            if len(val_energies) == 1:
                self.nn_valid_ene = (self.energy_mlp - val_energies[0]) / np.sqrt(2)
            else:
                all_E = np.array(val_energies + [self.energy_mlp], dtype=float)
                self.nn_valid_ene = all_E.std(ddof=1)

            # The legacy fields/getters now deliberately carry MLP-only
            # disagreement in explicit-QEq mode.
            self.nn_valid_mlp = self.nn_valid_ene
            self.nn_valid_mlp_maxF = self.nn_valid_maxF

            embedding_energies = list(
                getattr(self, "_last_val_embedding_energies", [])
            )
            production_embedding = self.energy - self.energy_mlp
            if len(embedding_energies) == 1:
                self.nn_valid_embedding = (
                    production_embedding - embedding_energies[0]
                ) / np.sqrt(2)
            elif len(embedding_energies) > 1:
                all_embedding = np.array(
                    embedding_energies + [production_embedding], dtype=float
                )
                self.nn_valid_embedding = all_embedding.std(ddof=1)
            else:
                self.nn_valid_embedding = 0.0

            if os.environ.get(
                "SCHNET_NN_VALID_EMBEDDING", ""
            ).strip().lower() in {"1", "true", "yes", "on"}:
                print(
                    "NN_VALIDATION "
                    f"nn_valid_mlp={self.nn_valid_mlp:.10g} "
                    f"nn_valid_mlp_maxF={self.nn_valid_mlp_maxF:.10g} "
                    f"nn_valid_embedding={self.nn_valid_embedding:.10g}",
                    flush=True,
                )

        if self._should_write_phi_static(time_step):
            self._write_phi_static(time_step=time_step, n_mlp_atoms=len(atomic_numbers))

        return None

    def _should_write_phi_static(self, time_step: int) -> bool:
        """Apply the NTWSE/adaptive-sampling policy to potential output."""
        if getattr(self, "phi_static_output_file", None) is None:
            return False

        ntwse = getattr(self, "ntwse", 0)
        if ntwse < 0:
            validation_is_current = (
                len(self.val_calculators) > 0
                and self.nn_valid_freq > 0
                and time_step % self.nn_valid_freq == 0
            )
            if not validation_is_current:
                return False

            return abs(self.nn_valid_ene) > self.val_thresh

        # NTWSE == 0 and NTWSE > 0 both use the forwarded write_val_step.
        return self.nn_valid_freq > 0 and time_step % self.nn_valid_freq == 0

    def _write_phi_static(self, time_step: int, n_mlp_atoms: int) -> None:
        """Append a timestep and its float32 potential to one NumPy stream."""
        output_file = getattr(self, "phi_static_output_file", None)
        if output_file is None:
            return

        phi_static = getattr(self, "_last_phi_static", None)
        if phi_static is None:
            raise RuntimeError("phi_static export requested before the potential was evaluated")

        phi_static = np.asarray(phi_static, dtype=np.float32)
        if phi_static.shape != (n_mlp_atoms,):
            raise ValueError(
                f"phi_static must have shape ({n_mlp_atoms},), got {phi_static.shape}"
            )

        record = np.empty(
            1,
            dtype=[("step", np.int64), ("phi_static", np.float32, (n_mlp_atoms,))],
        )
        record["step"][0] = time_step
        record["phi_static"][0] = phi_static

        # Each np.save block is self-describing. Locking prevents concurrent
        # writers from interleaving blocks in a shared output file.
        with output_file.open("ab") as stream:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
            try:
                np.save(stream, record, allow_pickle=False)
                stream.flush()
                os.fsync(stream.fileno())
            finally:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)

    def _collect_b2_diagnostics(
        self,
        E_total: torch.Tensor,
        E_mlp: torch.Tensor,
        E_embedding: torch.Tensor,
        qeq_mode: str,
        q0: torch.Tensor,
        q_phi: torch.Tensor,
        phi: torch.Tensor,
        r_qeq_A: torch.Tensor,
        r_or_A: torch.Tensor,
        F_total_qeq: np.ndarray,
        F_total_or: np.ndarray,
        F_direct_qeq: np.ndarray,
        F_direct_or: np.ndarray,
        F_response_qeq: np.ndarray,
        F_response_or: np.ndarray,
        F_mlp_qeq: np.ndarray,
        n_link_atoms: int,
    ) -> dict:
        n_qeq = int(q_phi.shape[0])
        n_phys = max(n_qeq - int(n_link_atoms), 0)

        q0_np = q0.detach().cpu().numpy()
        qphi_np = q_phi.detach().cpu().numpy()
        phi_np = phi.detach().cpu().numpy()
        E_static = torch.sum(q0 * phi)
        E_induced = 0.5 * torch.sum((q_phi - q0) * phi)
        identity_error = torch.abs(E_embedding - E_static - E_induced)
        f_qeq_sum = np.sum(F_total_qeq, axis=0) if F_total_qeq.size else np.zeros(3)
        f_or_sum = np.sum(F_total_or, axis=0) if F_total_or.size else np.zeros(3)
        f_mlp_sum = np.sum(F_mlp_qeq, axis=0) if F_mlp_qeq.size else np.zeros(3)

        if r_or_A.shape[0] > 0 and n_phys > 0:
            dist = torch.cdist(r_qeq_A[:n_phys], r_or_A)
            min_qeq_or_A = float(torch.min(dist).detach().cpu().item())
        else:
            min_qeq_or_A = float("nan")

        if n_link_atoms > 0 and r_or_A.shape[0] > 0:
            dist_link = torch.cdist(r_qeq_A[n_phys:], r_or_A)
            min_link_or_A = float(torch.min(dist_link).detach().cpu().item())
        else:
            min_link_or_A = float("nan")

        return {
            "E_total": float(E_total.detach().cpu().item()),
            "E_mlp": float(E_mlp.detach().cpu().item()),
            "E_static": float(E_static.detach().cpu().item()),
            "E_induced": float(E_induced.detach().cpu().item()),
            "E_embedding": float(E_embedding.detach().cpu().item()),
            "qeq_mode": qeq_mode,
            "embedding_identity_error": float(identity_error.detach().cpu().item()),
            "q0_sum": float(np.sum(q0_np)),
            "qphi_sum": float(np.sum(qphi_np)),
            "q_phys_sum": float(np.sum(qphi_np[:n_phys])),
            "q_link_sum": float(np.sum(qphi_np[n_phys:])),
            "q_min": float(np.min(qphi_np)) if qphi_np.size else 0.0,
            "q_max": float(np.max(qphi_np)) if qphi_np.size else 0.0,
            "max_abs_dq": float(np.max(np.abs(qphi_np - q0_np))) if qphi_np.size else 0.0,
            "phi_phys_min": float(np.min(phi_np[:n_phys])) if n_phys > 0 else 0.0,
            "phi_phys_max": float(np.max(phi_np[:n_phys])) if n_phys > 0 else 0.0,
            "phi_link_max_abs": float(np.max(np.abs(phi_np[n_phys:]))) if n_link_atoms > 0 else 0.0,
            "F_qeq_max": float(np.max(np.linalg.norm(F_total_qeq, axis=1))) if F_total_qeq.size else 0.0,
            "F_or_max": float(np.max(np.linalg.norm(F_total_or, axis=1))) if F_total_or.size else 0.0,
            "F_direct_qeq_max": float(np.max(np.linalg.norm(F_direct_qeq, axis=1))) if F_direct_qeq.size else 0.0,
            "F_direct_or_max": float(np.max(np.linalg.norm(F_direct_or, axis=1))) if F_direct_or.size else 0.0,
            "F_response_qeq_max": float(np.max(np.linalg.norm(F_response_qeq, axis=1))) if F_response_qeq.size else 0.0,
            "F_response_or_max": float(np.max(np.linalg.norm(F_response_or, axis=1))) if F_response_or.size else 0.0,
            "F_mlp_qeq_max": float(np.max(np.linalg.norm(F_mlp_qeq, axis=1))) if F_mlp_qeq.size else 0.0,
            "F_total_sum": f_qeq_sum + f_or_sum,
            "F_qeq_sum": f_qeq_sum,
            "F_or_sum": f_or_sum,
            "F_mlp_sum": f_mlp_sum,
            "min_qeq_or_A": min_qeq_or_A,
            "min_link_or_A": min_link_or_A,
            "n_qeq": n_qeq,
            "n_link_atoms": int(n_link_atoms),
            "n_or": int(r_or_A.shape[0]),
        }

    def _print_b2_diagnostics(self, time_step: int) -> None:
        diag = getattr(self, "_last_b2_diag", None)
        if not diag:
            return

        def v2s(vec):
            return "[" + ", ".join(f"{x:.6e}" for x in vec) + "]"

        print(
            "B2DBG "
            f"step={time_step} "
            f"qeq_mode={diag['qeq_mode']} "
            f"n_qeq={diag['n_qeq']} n_link={diag['n_link_atoms']} n_or={diag['n_or']} "
            f"E_total={diag['E_total']:.10g} E_mlp={diag['E_mlp']:.10g} "
            f"E_static={diag['E_static']:.10g} E_induced={diag['E_induced']:.10g} "
            f"E_embedding={diag['E_embedding']:.10g} "
            f"E_identity_err={diag['embedding_identity_error']:.3e} "
            f"q0_sum={diag['q0_sum']:.8g} qphi_sum={diag['qphi_sum']:.8g} "
            f"max_abs_dq={diag['max_abs_dq']:.6g} "
            f"q_phys={diag['q_phys_sum']:.8g} q_link={diag['q_link_sum']:.8g} "
            f"q_range=[{diag['q_min']:.6g},{diag['q_max']:.6g}] "
            f"phi_phys=[{diag['phi_phys_min']:.6g},{diag['phi_phys_max']:.6g}] "
            f"phi_link_max_abs={diag['phi_link_max_abs']:.3e} "
            f"Fmax_qeq={diag['F_qeq_max']:.6g} Fmax_or={diag['F_or_max']:.6g} "
            f"Fmax_mlp={diag['F_mlp_qeq_max']:.6g} "
            f"Fsum_qeq={v2s(diag['F_qeq_sum'])} Fsum_or={v2s(diag['F_or_sum'])} "
            f"Fsum_total={v2s(diag['F_total_sum'])} Fsum_mlp={v2s(diag['F_mlp_sum'])} "
            f"min_qeq_or_A={diag['min_qeq_or_A']:.6g} min_link_or_A={diag['min_link_or_A']:.6g}",
            flush=True,
        )
        print(
            "B2RESP "
            f"step={time_step} "
            f"Fmax_qeq_full={diag['F_qeq_max']:.6g} "
            f"Fmax_qeq_static_direct={diag['F_direct_qeq_max']:.6g} "
            f"Fmax_qeq_polarization={diag['F_response_qeq_max']:.6g} "
            f"Fmax_or_full={diag['F_or_max']:.6g} "
            f"Fmax_or_static_direct={diag['F_direct_or_max']:.6g} "
            f"Fmax_or_polarization={diag['F_response_or_max']:.6g}",
            flush=True,
        )

    def _print_b2_finite_difference_check(
        self,
        calculator: spk.interfaces.SpkCalculator,
        system: ase.Atoms,
        r_qeq_A: torch.Tensor,
        r_or_A: torch.Tensor,
        q_or: torch.Tensor,
        cutoff_nm: float,
        n_link_atoms: int,
        F_total_qeq: np.ndarray,
        F_total_or: np.ndarray,
    ) -> None:
        h_A = float(os.environ.get("SCHNET_B2_FD_STEP_A", "0.001"))

        probes: list[tuple[str, int, int, float]] = []
        if F_total_qeq.size:
            flat = int(np.argmax(np.abs(F_total_qeq)))
            atom, dim = np.unravel_index(flat, F_total_qeq.shape)
            probes.append(("qeq", int(atom), int(dim), float(F_total_qeq[atom, dim])))
        if F_total_or.size:
            flat = int(np.argmax(np.abs(F_total_or)))
            atom, dim = np.unravel_index(flat, F_total_or.shape)
            probes.append(("or", int(atom), int(dim), float(F_total_or[atom, dim])))

        for kind, atom, dim, analytic_force in probes:
            rq_plus = r_qeq_A.detach().clone()
            rq_minus = r_qeq_A.detach().clone()
            ro_plus = r_or_A.detach().clone()
            ro_minus = r_or_A.detach().clone()

            if kind == "qeq":
                rq_plus[atom, dim] += h_A
                rq_minus[atom, dim] -= h_A
            else:
                ro_plus[atom, dim] += h_A
                ro_minus[atom, dim] -= h_A

            e_plus = self._b2_eval_model(
                calculator=calculator,
                system=system,
                r_qeq_A=rq_plus.requires_grad_(True),
                r_or_A=ro_plus.requires_grad_(True),
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
            )[0]
            e_minus = self._b2_eval_model(
                calculator=calculator,
                system=system,
                r_qeq_A=rq_minus.requires_grad_(True),
                r_or_A=ro_minus.requires_grad_(True),
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
            )[0]

            fd_force = -(e_plus - e_minus) / (2.0 * h_A)
            abs_err = abs(analytic_force - fd_force)
            rel_err = abs_err / max(abs(analytic_force), abs(fd_force), 1.0)
            print(
                "B2FD "
                f"kind={kind} atom={atom} dim={dim} h_A={h_A:g} "
                f"F_autograd={analytic_force:.10g} F_fd={fd_force:.10g} "
                f"abs_err={abs_err:.4e} rel_err={rel_err:.4e}",
                flush=True,
            )

    @staticmethod
    def _zero_link_atom_phi(phi: torch.Tensor, n_link_atoms: int) -> torch.Tensor:
        """
        Link atoms are geometric capping sites. They may be present in the model
        input and receive model/QEq quantities, but they should not couple
        directly to the OR electrostatic potential.
        """
        if n_link_atoms <= 0:
            return phi
        if n_link_atoms > phi.shape[0]:
            raise ValueError("n_link_atoms cannot exceed the number of QEq sites")

        phi = phi.clone()
        phi[-n_link_atoms:] = 0.0
        return phi

    def _coulomb_phi_kjmol(
        self,
        r_qeq_A,
        r_or_A,
        q_or_e,
        damping=None,
        sigma_A=None,
    ):
        """
        Compute the damped external electrostatic potential on QEq atoms.

        Units:
            coordinates and sigma in Å
            q in e
            returns phi in kJ/mol/e

        The OR atoms are assumed to have already been pairlisted by GROMOS.
        """
        if r_or_A.shape[0] == 0:
            return torch.zeros(r_qeq_A.shape[0], dtype=r_qeq_A.dtype, device=r_qeq_A.device)

        damping = damping or getattr(self, "electrostatic_damping", "soft")
        sigma_A = float(sigma_A if sigma_A is not None else getattr(self, "electrostatic_sigma_A", 0.8))
        if damping not in {"soft", "erf", "none"}:
            raise ValueError("damping must be 'soft', 'erf', or 'none'")
        if sigma_A <= 0.0:
            raise ValueError("sigma must be positive")

        diff = r_qeq_A[:, None, :] - r_or_A[None, :, :]
        r = torch.linalg.norm(diff, dim=-1).clamp_min(1.0e-8)

        if damping == "soft":
            inv_r = torch.rsqrt(r * r + sigma_A * sigma_A)
        elif damping == "erf":
            inv_r = torch.erf(r / sigma_A) / r
        else:
            inv_r = 1.0 / r

        return K_E_KJMOL_A_E2 * torch.sum(inv_r * q_or_e[None, :], dim=1)

    def _model_forward_no_forces(self, inputs: dict) -> dict:
        """
        Forward through SchNetPack model while skipping any output module that
        computes forces via autograd. Still runs input_modules so PaiNN gets
        neighbor geometry (_Rij, idx_i/j, ...).
        """
        model = self.pred_calculator.model

        # 1) run input modules (neighbor list etc.)
        input_modules = getattr(model, "input_modules", None)
        if input_modules is not None:
            for mod in input_modules:
                inputs = mod(inputs)

        # 2) representation (PaiNN)
        inputs = model.representation(inputs)

        # 3) output modules (skip forces-producing ones)
        force_key = getattr(self.pred_calculator, "force_key", "forces")

        for mod in model.output_modules:
            outs = getattr(mod, "model_outputs", [])
            # Skip anything that outputs forces
            if (force_key in outs) or ("forces" in outs):
                continue
            inputs = mod(inputs)

        return inputs

    def calculate_next_step_b2(
        self,
        atomic_numbers: list,
        positions_A: list,
        time_step: int,
        or_positions_nm: list,
        or_charges_e: list,
        cutoff_nm: float,
        n_link_atoms: int = 0,
    ) -> None:
        """
        Backward-compatible entry point for dynamic-charge/B2 calculations.

        The public GROMOS interface remains in nm; ``calculate_next_step``
        converts OR coordinates and the cutoff once and performs the model and
        electrostatic calculation in Å.
        """
        return self.calculate_next_step(
            atomic_numbers=atomic_numbers,
            positions_A=positions_A,
            time_step=time_step,
            dynamic_charges=True,
            or_positions_nm=or_positions_nm,
            or_charges_e=or_charges_e,
            cutoff_nm=cutoff_nm,
            n_link_atoms=n_link_atoms,
        )

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
        Get the MLP-only validation energy deviation.

        The legacy name is retained for the C++ interface.

        Returns:
            float: Validation energy deviation.
        """
        return self.nn_valid_ene

    def get_nn_valid_maxF(self):
        """
        Get the MLP-only maximum force committee disagreement.

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
    """
    Pert_SchNet_V1_Calculator is a specialized calculator for computing perturbed energies and forces 
    in molecular dynamics simulations using the SchNet model. It extends the SchNet_V1_Calculator class 
    to handle perturbations between different quantum mechanical (QM) end states.

    This class supports the following functionalities:
    - Initialization with model paths, validation frequencies, and perturbation settings.
    - Retrieval of atomic numbers and positions for specified QM states.
    - Calculation of perturbed energies and their derivatives.
    - Validation of perturbed energy calculations using multiple validation models.
    - Calculation of perturbed forces for the system.
    - Prediction of the next simulation step by updating energies, forces, and derivatives.

    Attributes:
        lam (float): Parameter for interpolating between different states.
        perturbed_qm_states (np.ndarray): Array specifying the QM states of each atom.
        qmzone_size (int): Length of the QM zone.
        num_states (int): Number of perturbed QM states.
        states_idx (dict): Dictionary storing indices of atoms in different states.
        states (dict): Dictionary storing atomic states and their properties.
        energy (float): Total perturbed energy.
        derivative (float): Derivative of the perturbed energy.
        forces (np.ndarray): Perturbed forces for the system.
        nn_valid_dev (float): Deviation of the perturbed energy during validation.

    Methods:
        __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int, lam: float, perturbed_qm_states: list) -> None:
            Initialize the Pert_SchNet_V1_Calculator with necessary model paths, validation frequency, and perturbation settings.

        get_z_and_positions_vac(self, atomic_numbers, positions, state_id):
            Retrieve the atomic numbers and positions within the quantum mechanical (QM) zone for a specified state.

        get_z_and_positions_vac_and_burnn(self, atomic_numbers, positions, state_id):
            Retrieve the atomic numbers and positions for both the vacuum and burnn configurations based on a specified state ID.

        get_state(self, atomic_numbers: list, positions: list, state_id: str) -> ase.Atoms:
            Construct two ASE `Atoms` objects representing the system end state in the vacuum and BuRNN environments for the specified perturbed quantum mechanical (QM) state.

        calculate_perturbed_energy_and_derivative(self) -> float:
            Calculate the perturbed energy and its derivative for two states (A and B) in both vacuum and Burnn environments.

        validate_state(self, state: str, calculator: spk.interfaces.SpkCalculator) -> float:
            Validate the state by calculating the potential energy using the specified validation calculator.

        validate_perturbed_energy(self):
            Validate the perturbed energy calculations using the provided end states in both vacuum and BuRNN environments.

        calculate_perturbed_forces(self):
            Calculate perturbed forces for the system.

        calculate_next_step(self, atomic_numbers: list, positions: list, time_step: int) -> None:
            Calculate the next step in the MLP/MM simulation by updating the states, predicting energies and forces, and calculating the perturbed energy and forces.

        get_derivative(self):
            Get the derivative of the perturbed energy.
    """
    def __init__(self, model_path: str, val_model_paths: list, nn_valid_freq: int, write_energy_freq: int,
                 spin_multiplicity:int, total_charge:int, lam: float, perturbed_qm_states: list,
                 endpoint_atomic_numbers_A=None, endpoint_atomic_numbers_B=None,
                 ntwse: int = 0, val_thresh: float = 0.0,
                 endpoint_total_charge_B=None, endpoint_spin_multiplicity_B=None) -> None:
        """
        Initialize the Pert_SchNet_V1_Calculator with necessary model paths, validation frequency, and perturbation settings.

        Args:
            model_path (str): Path to the primary SchNet model for energy and force calculations.
            val_model_paths (list): List of paths to validation models used to validate predictions.
            nn_valid_freq (int): Frequency (in steps) at which to validate the neural network predictions.
            write_energy_freq (int): Frequency (in steps) at which to write the energy predictions to output.
            lam (float): A parameter used for interpolating between different states.
            perturbed_qm_states (list): List specifying the QM states of each atom, where 0 indicates both states, 1 indicates state A, and 2 indicates state B.
            endpoint_total_charge_B: Net QM + buffer charge for endpoint B; defaults to A.
            endpoint_spin_multiplicity_B: Combined QM + buffer multiplicity for B; defaults to A.

        Returns:
            None
        """
        # Call the constructor of the parent class, SchNet_V1_Calculator, with the provided model paths and frequencies.
        super().__init__(
            model_path,
            val_model_paths,
            nn_valid_freq,
            write_energy_freq,
            spin_multiplicity,
            total_charge,
            ntwse=ntwse,
            val_thresh=val_thresh,
        )
        
        # Store the perturbation parameter lambda (lam) for interpolating between states.
        self.lam = lam
        # Keep endpoint electronic states fixed as lambda changes. Omitted B
        # metadata preserves the historical same-charge/spin Python API.
        self.endpoint_total_charge = {
            "A": total_charge,
            "B": total_charge if endpoint_total_charge_B is None else endpoint_total_charge_B,
        }
        self.endpoint_spin_multiplicity = {
            "A": spin_multiplicity,
            "B": spin_multiplicity if endpoint_spin_multiplicity_B is None else endpoint_spin_multiplicity_B,
        }
        
        # Convert the list of perturbed QM states into a numpy array for easier manipulation and indexing.
        self.perturbed_qm_states = np.array(perturbed_qm_states)
        self.endpoint_atomic_numbers = {
            "A": np.asarray(endpoint_atomic_numbers_A or [], dtype=int),
            "B": np.asarray(endpoint_atomic_numbers_B or [], dtype=int),
        }
        # A shared physical atom whose endpoint nuclear charges differ is an
        # alchemical transmutation (for example Mn, Z=25 -> Zn, Z=30).  Such a
        # perturbation interpolates the two complete BuRNN/QEq endpoint
        # Hamiltonians directly; it must not add the other endpoint's vacuum
        # energy, as required by the older separate-region insertion/deletion
        # construction.
        za = self.endpoint_atomic_numbers["A"]
        zb = self.endpoint_atomic_numbers["B"]
        n_endpoint_z = min(len(self.perturbed_qm_states), len(za), len(zb))
        self.direct_endpoint_interpolation = any(
            self.perturbed_qm_states[i] == 0
            and za[i] > 0
            and zb[i] > 0
            and za[i] != zb[i]
            for i in range(n_endpoint_z)
        )
        
        # Length of the QM zone.
        self.qmzone_size = len(self.perturbed_qm_states)

        # Initialize the `num_states` attribute with the number of perturbed QM states.
        self.num_states = len(set([s for s in self.perturbed_qm_states if s != 0]))

        # Create a dictionary to store indices of atoms in different states.
        self.states_idx = {}

        # Iterate over the perturbed QM states and populate the `states_idx` dictionary.
        for i, s in enumerate(self.perturbed_qm_states):
            if s == 0:
                # If state is both (0), add index to all relevant lists.
                try: 
                    self.states_idx['A'].append(i)
                    self.states_idx['B'].append(i)
                    self.states_idx['both'].append(i)
                except KeyError:
                    self.states_idx['A'] = [i]
                    self.states_idx['B'] = [i]
                    self.states_idx['both'] = [i]
            elif s == 1:
                # If state is 1, add index to the 'A' list.
                try: 
                    self.states_idx['A'].append(i)
                except KeyError: 
                    self.states_idx['A'] = [i]
            elif s == 2:
                # If state is 2, add index to the 'B' list.
                try: 
                    self.states_idx['B'].append(i)
                except KeyError: 
                    self.states_idx['B'] = [i]
            else:
                print("Error: Unknown state value:", s)
        
        # Check if both states are present and populate their respective lists.
        if 'both' in self.states_idx.keys():
            # Add indices of atoms only in A to the 'only_A' list.
            self.states_idx['only_A'] = [i for i in self.states_idx['A'] if i not in self.states_idx['both']]
            # Add indices of atoms only in B to the 'only_B' list.
            self.states_idx['only_B'] = [i for i in self.states_idx['B'] if i not in self.states_idx['both']]
        
        return


    def get_z_and_positions_vac(self, atomic_numbers, positions, state_id):
        """
        Retrieve the atomic numbers and positions within the quantum mechanical (QM) zone for a specified state.

        This method filters the atomic numbers and positions based on the provided state ID, returning
        only those that belong to the specified state within the QM zone.

        Args:
            atomic_numbers (list): A list of atomic numbers for the entire system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (str): The state identifier ('A' for state A, 'B' for state B) to filter the QM zone.

        Returns:
            tuple: Two lists containing:
                - atomic_numbers_vac (list): Atomic numbers within the QM zone for the specified state.
                - positions_vac (list): Positions corresponding to the atomic numbers within the QM zone.
        """
        atomic_numbers_vac = []  # Initialize an empty list to store atomic numbers for the specified state.
        positions_vac = []       # Initialize an empty list to store positions for the specified state.

        # Iterate over the indices of atoms in the specified state.
        for i in self.states_idx[state_id]:
            # Add atomic number and position to the respective lists.
            endpoint_z = self.endpoint_atomic_numbers[state_id]
            z = int(endpoint_z[i]) if i < len(endpoint_z) else 0
            atomic_numbers_vac.append(z if z > 0 else atomic_numbers[i])
            positions_vac.append(positions[i])

        return atomic_numbers_vac, positions_vac


    def get_z_and_positions_vac_and_burnn(self, atomic_numbers, positions, state_id):
        """
        Retrieve the atomic numbers and positions for both the vacuum and burnn configurations based on a specified state ID.

        This method separates atomic numbers and positions into two groups:
        - Vacuum (vac): Atoms within the quantum mechanical (QM) zone corresponding to the specified state.
        - BuRNN (burnn): All atoms including the buffer zone.

        Args:
            atomic_numbers (list): A list of atomic numbers for the system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (str): The state identifier ('A' for state A, 'B' for state B) to filter the QM zone.

        Returns:
            tuple: Four lists containing:
                - atomic_numbers_vac (list): Atomic numbers within the QM zone for the specified state.
                - atomic_numbers_burnn (list): Atomic numbers for the entire Inner + Buffer zone, starting with those in the Inner zone.
                - positions_vac (list): Positions corresponding to the atomic numbers within the QM zone.
                - positions_burnn (list): Positions for the entire Inner + Buffer zone, starting with those in the Inner zone.
        """
        # Get atomic numbers and positions for the specified state within the QM zone.
        atomic_numbers_vac, positions_vac = self.get_z_and_positions_vac(
            atomic_numbers[:self.qmzone_size],  # Slice atomic numbers up to the QM zone size.
            positions[:self.qmzone_size],       # Slice positions up to the QM zone size.
            state_id                            # State identifier to filter the QM zone.
        )

        # Create a copy of the vacuum atomic numbers to start the burnn configuration.
        atomic_numbers_burnn = atomic_numbers_vac.copy()

        # Create a copy of the vacuum positions to start the burnn configuration.
        positions_burnn = positions_vac.copy()

        # Iterate over the atomic numbers and positions beyond the QM zone.
        for i, a in enumerate(atomic_numbers[self.qmzone_size:]):
            # Append the atomic number to the burnn configuration list.
            atomic_numbers_burnn.append(a)

            # Append the position to the burnn configuration list.
            positions_burnn.append(positions[i + self.qmzone_size])

        return atomic_numbers_vac, atomic_numbers_burnn, positions_vac, positions_burnn


    def get_state(self, atomic_numbers: list, positions: list, state_id: str) -> ase.Atoms:
        """
        Construct two ASE `Atoms` objects representing the system end state in the vacuum and BuRNN environments 
        for the specified perturbed quantum mechanical (QM) state.

        This method divides the atomic numbers and positions into two states: vacuum and Burnn, based on 
        the provided `state_id`, and then creates corresponding ASE `Atoms` objects for each state.

        Args:
            atomic_numbers (list): A list of atomic numbers for the entire system.
            positions (list): A list of atomic positions corresponding to the atomic numbers.
            state_id (str): The state identifier ('A' for state A, 'B' for state B) used to filter and construct end states.

        Returns:
            tuple: A pair of ASE `Atoms` objects:
                - state_vac (ase.Atoms): The vacuum state representing the Z and positions for the given end state.
                - state_burnn (ase.Atoms): The BuRNN state representing the Z and positions for the given end state, including the buffer region.
        """
        # Retrieve atomic numbers and positions for both vacuum and Burnn states for the specified state ID.
        atomic_numbers_vac, atomic_numbers_burnn, positions_vac, positions_burnn = self.get_z_and_positions_vac_and_burnn(atomic_numbers, positions, state_id)

        # Create ASE Atoms object for the vacuum state.
        state_vac = ase.Atoms(numbers=atomic_numbers_vac, positions=positions_vac)

        # Create ASE Atoms object for the Burnn state.
        state_burnn = ase.Atoms(numbers=atomic_numbers_burnn, positions=positions_burnn)

        for system in (state_vac, state_burnn):
            system.info["total_charge"] = float(self.endpoint_total_charge[state_id])
            system.info["spin_multiplicity"] = float(self.endpoint_spin_multiplicity[state_id])

        return state_vac, state_burnn


    def calculate_perturbed_energy_and_derivative(self) -> float:
        """
        Calculate the perturbed energy and its derivative for two states (A and B) in both vacuum and Burnn environments.

        This method computes a linear combination of the energies for state A and state B in both the vacuum 
        and Burnn conditions, using the lambda (`lam`) parameter. It also calculates the derivative of the 
        perturbed energy with respect to lambda.

        Returns:
            tuple: A tuple containing:
                - perturbed_energy (float): The total perturbed energy.
                - perturbed_energy_derivative (float): The derivative of the perturbed energy with respect to lambda.
        """
        if (
            getattr(self, "direct_endpoint_interpolation", False)
            and "A" in self.states
            and "B" in self.states
        ):
            energy_a = self.states["A"]["energy_burnn"]
            energy_b = self.states["B"]["energy_burnn"]
            self.states["A"]["perturbed_energy"] = (1.0 - self.lam) * energy_a
            self.states["B"]["perturbed_energy"] = self.lam * energy_b
            self.states["A"]["derivative"] = -energy_a
            self.states["B"]["derivative"] = energy_b
            return (
                (1.0 - self.lam) * energy_a + self.lam * energy_b,
                energy_b - energy_a,
            )

        # Older insertion/deletion behavior for distinct state-specific atoms.
        # Iterate over each state to calculate the perturbed energy and its derivative.
        for state in self.states.keys():
            if state == 'A':
                # Calculate perturbed energy and derivative for state A.
                self.states[state]['perturbed_energy'] = (1 - self.lam) * self.states[state]['energy_burnn'] + self.lam * self.states[state]['energy_vac']
                self.states[state]['derivative'] = -self.states[state]['energy_burnn'] + self.states[state]['energy_vac']
            else:
                # Calculate perturbed energy and derivative for state B.
                self.states[state]['perturbed_energy'] = self.lam * self.states[state]['energy_burnn'] + (1 - self.lam) * self.states[state]['energy_vac']
                self.states[state]['derivative'] = self.states[state]['energy_burnn'] - self.states[state]['energy_vac']

        # Return the perturbed energies and derivatives for state A
        if len(self.states.keys()) == 1:
            return self.states['A']['perturbed_energy'], self.states['A']['derivative']
        # Sum the perturbed energies and derivatives of states A and B to get the total perturbed energy and derivative.
        else:
            return self.states['A']['perturbed_energy'] + self.states['B']['perturbed_energy'], self.states['A']['derivative'] + self.states['B']['derivative']

    def validate_state(self, state: str, calculator: spk.interfaces.SpkCalculator) -> float:
        """
        Validate the state by calculating the potential energy using the specified validation calculator.

        Args:
            state (str): The state to validate. Expected values are 'A' or other states.
            calculator (spk.interfaces.SpkCalculator): The calculator to set for the states.

        Returns:
            float: The calculated potential energy based on the state and lambda value.
        """
        self.states[state]['vac'].set_calculator(calculator)
        self.states[state]['burnn'].set_calculator(calculator)

        energy_burnn = self.states[state]['burnn'].get_potential_energy()
        energy_vac = self.states[state]['vac'].get_potential_energy()

        if getattr(self, "direct_endpoint_interpolation", False):
            return (1.0 - self.lam) * energy_burnn if state == 'A' else self.lam * energy_burnn

        if state == 'A':
            return (1 - self.lam) * energy_burnn + self.lam * energy_vac
        else:
            return self.lam * energy_burnn + (1 - self.lam) * energy_vac
    
    def validate_perturbed_energy(self):
        """
        Validate the perturbed energy calculations using the provided end states in both vacuum and BuRNN environments.

        This method applies validation calculators to the provided atomic states to compute the perturbed energies 
        for both states A and B in vacuum and Burnn environments. It then calculates the deviation between the 
        previously computed perturbed energies and the validation energies, using either a single or multiple 
        validation calculators.

        Returns:
            float: The maximum absolute deviation between the computed perturbed energies and the validation energies.
        """
        # Case when only one state (A) is present.
        if len(self.states.keys()) == 1:
            # Case when only one validation calculator is provided.
            if len(self.val_calculators) == 1:
                # Compute validation perturbed energy for state A.
                perturbed_energy_a_val = self.validate_state('A', self.val_calculators[0])

                # Calculate deviation between computed and validation perturbed energy.
                dev = abs(self.states['A']['perturbed_energy'] - perturbed_energy_a_val)
        
            else:
                # Initialize lists to store validation perturbed energies for multiple calculators.
                perturbed_energies_a_val = [self.states['A']['perturbed_energy']]
                for val_calculator in self.val_calculators:
                    # Compute validation perturbed energy for state A.
                    perturbed_energy_a_val = self.validate_state('A', val_calculator)
                    # Append the calculated validation energy to the list.
                    perturbed_energies_a_val.append(perturbed_energy_a_val)
                
                # Convert the list of validation energies to a numpy array.
                perturbed_energies_a_val = np.array(perturbed_energies_a_val)
                # Calculate variance of the validation energies.
                dev = perturbed_energies_a_val.std()
            return dev

        # Case when two states (A and B) are present.
        else:
            # Case when only one validation calculator is provided.
            if len(self.val_calculators) == 1:
                # Compute validation perturbed energies for state A and state B.
                perturbed_energy_a_val = self.validate_state('A', self.val_calculators[0])
                perturbed_energy_b_val = self.validate_state('B', self.val_calculators[0])

                # Calculate deviations between computed and validation perturbed energies.
                dev_a = abs(self.states['A']['perturbed_energy'] - perturbed_energy_a_val)
                dev_b = abs(self.states['B']['perturbed_energy'] - perturbed_energy_b_val)

            else:
                # Initialize lists to store validation perturbed energies for multiple calculators.
                perturbed_energies_a_val = [self.states['A']['perturbed_energy']]
                perturbed_energies_b_val = [self.states['B']['perturbed_energy']]

                # Iterate over each validation calculator and compute perturbed energies.
                for val_calculator in self.val_calculators:
                    # Compute validation perturbed energies for state A and state B.
                    perturbed_energy_a_val = self.validate_state('A', val_calculator)
                    perturbed_energy_b_val = self.validate_state('B', val_calculator)

                    # Append the calculated validation energies to the lists.
                    perturbed_energies_a_val.append(perturbed_energy_a_val)
                    perturbed_energies_b_val.append(perturbed_energy_b_val)

                # Convert the lists of validation energies to numpy arrays.
                perturbed_energies_a_val = np.array(perturbed_energies_a_val)
                perturbed_energies_b_val = np.array(perturbed_energies_b_val)

                # Calculate variance of the validation energies for state A and state B.
                dev_a = perturbed_energies_a_val.var() / (len(perturbed_energies_a_val) - 1)
                dev_b = perturbed_energies_b_val.var() / (len(perturbed_energies_b_val) - 1)
            # Return the maximum absolute deviation between computed and validation energies.
            return np.max([dev_a, dev_b])


    def calculate_perturbed_forces(self):
        """
        Calculate the perturbed forces for the QM zone and buffer zone based on the lambda factor and the end states.
        This method calculates the perturbed forces for the QM zone and buffer zone by interpolating between the forces
        from different states ('A' and 'B') using the lambda factor. It handles cases where only one state ('A') is present
        or both states ('A' and 'B') are present.
        Returns:
            np.ndarray: A numpy array of shape (qm zone size + buffer size, 3) containing the calculated perturbed forces.
        """
        if (
            getattr(self, "direct_endpoint_interpolation", False)
            and "A" in self.states
            and "B" in self.states
            and "burnn_full_indices" in self.states["A"]
        ):
            n_total = max(
                max(self.states["A"]["burnn_full_indices"], default=-1),
                max(self.states["B"]["burnn_full_indices"], default=-1),
            ) + 1
            forces = np.zeros((n_total, 3), dtype=float)
            for state, weight in (("A", 1.0 - self.lam), ("B", self.lam)):
                for local_index, full_index in enumerate(
                    self.states[state]["burnn_full_indices"]
                ):
                    forces[full_index] += (
                        weight * self.states[state]["forces_burnn"][local_index]
                    )
            return forces

        # Initialize the perturbed forces array, shape (qm zone size + buffer size, 3).
        forces = np.zeros((self.qmzone_size + (len(self.states['A']['forces_burnn']) - len(self.states['A']['forces_vac'])), 3), dtype=float)

        def calculate_forces(state, indices, lam_factor)->None:
            """
            Calculate the forces for a given state and set of indices.

            Args:
                state (int): The state for which to calculate the forces.
                indices (list of int): The indices of the particles for which to calculate the forces.
                lam_factor (float): The lambda factor used to interpolate between two sets of forces.

            Returns:
                None
            """
            for i in indices:
                i_state = self.states_idx[state].index(i)
                forces[i] = lam_factor * self.states[state]['forces_burnn'][i_state] + (1 - lam_factor) * self.states[state]['forces_vac'][i_state]
        
        def calculate_forces_both(indices, lam_factor)->None:
            """
            Calculate the combined forces for given indices using a lambda factor.

            This function computes the forces for each index in the provided list of indices
            by combining forces from different states ('A' and 'B') and conditions ('burnn' and 'vac')
            using the specified lambda factor.

            Args:
                indices (list of int): List of indices for which the forces need to be calculated.
                lam_factor (float): Lambda factor used to weight the contributions from different states and conditions.

            Returns:
                None: The function updates the `forces` array in place.
            """
            for i in indices:
                i_a = self.states_idx['A'].index(i)
                i_b = self.states_idx['B'].index(i)
                forces[i] = lam_factor * self.states['A']['forces_burnn'][i_a] + (1 - lam_factor) * self.states['A']['forces_vac'][i_a] + (1 - lam_factor) * self.states['B']['forces_burnn'][i_b] + lam_factor * self.states['B']['forces_vac'][i_b]

        # Case when only one state (A) is present.
        if len(self.states.keys()) == 1:
            # Calculate forces for the QM zone.
            forces[:self.qmzone_size] = (1 - self.lam) * self.states['A']['forces_burnn'][:self.qmzone_size] + self.lam * self.states['A']['forces_vac']
            # Calculate forces for the buffer zone.
            forces[self.qmzone_size:] = (1 - self.lam) * self.states['A']['forces_burnn'][self.qmzone_size:]
            return forces

        # Case when two states (A and B) are present.
        if len(self.states.keys()) == 2:
            if 'both' in self.states_idx.keys():
                # Calculate forces for atoms in both states.
                calculate_forces_both(self.states_idx['both'], 1 - self.lam)
                # Calculate forces for atoms only in state A.
                calculate_forces('A', self.states_idx['only_A'], 1 - self.lam)
                # Calculate forces for atoms only in state B.
                calculate_forces('B', self.states_idx['only_B'], self.lam)
            else:
                # Calculate forces for atoms in state A.
                calculate_forces('A', self.states_idx['A'], 1 - self.lam)
                # Calculate forces for atoms in state B.
                calculate_forces('B', self.states_idx['B'], self.lam)
            # Calculate forces for the buffer zone.
            forces[self.qmzone_size:] = (1 - self.lam) * self.states['A']['forces_burnn'][len(self.states_idx['A']):] + self.lam * self.states['B']['forces_burnn'][len(self.states_idx['B']):]
            return forces

    def _burnn_full_indices(self, state: str, n_total: int) -> list[int]:
        """Map a perturbed-state BuRNN array back to the full GROMOS order."""
        return list(self.states_idx[state]) + list(range(self.qmzone_size, n_total))

    def _print_dynamic_dhdl(self, time_step: int, endpoint_components: dict) -> None:
        """Report endpoint-weight derivatives in kJ/mol (lambda is dimensionless).

        Static = q0 dot phi; induced = 0.5*(q_phi-q0) dot phi.
        Endpoint charge solutions are lambda-independent in this Hamiltonian.
        The residual exposes any mismatch with the model's embedding energy.
        """
        components = dict(burnn=0.0, elecstatic=0.0, elecinduced=0.0, vacuum=0.0)
        direct = (
            getattr(self, "direct_endpoint_interpolation", False)
            and "A" in self.states and "B" in self.states
        )
        for state, energies in endpoint_components.items():
            sign = -1.0 if state == "A" else 1.0
            for name, energy in energies.items():
                components[name] += sign * energy
            if not direct:
                components["vacuum"] -= sign * self.states[state]["energy_vac"]
        total = sum(components.values())
        print(
            f"DHDL_DEBUG step={time_step} lambda={self.lam:.10g} units=kJ/mol "
            + " ".join(f"d{name}/dlambda={value:.10g}" for name, value in components.items())
            + f" sum={total:.10g} dH/dlambda={self.derivative:.10g} "
            + f"residual={self.derivative - total:.10g}",
            flush=True,
        )

    def _calculate_next_step_dynamic(
        self,
        atomic_numbers: list,
        positions: list,
        time_step: int,
        or_positions_nm: list,
        or_charges_e: list,
        cutoff_nm: float,
        n_link_atoms: int,
    ) -> None:
        """Evaluate the lambda Hamiltonian with explicit QEq embedding.

        Only the BuRNN endpoint contains the OR embedding, matching the
        existing delta-model perturbation Hamiltonian.  Consequently the OR
        force and effective dynamic charge use the same endpoint weights as
        the BuRNN energies: (1-lambda) for A and lambda for B.
        """
        # Collect and print the derivative breakdown on NN validation steps.
        debug_dhdl = len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0
        endpoint_components = {}
        self.states = {}
        n_total = len(atomic_numbers)
        r_or_A = torch.tensor(
            np.asarray(or_positions_nm, dtype=float).reshape((-1, 3)) * 10.0,
            dtype=torch.float32,
            device=self.torchdevice,
            requires_grad=True,
        )
        q_or = torch.tensor(or_charges_e, dtype=torch.float32, device=self.torchdevice)

        for state in ("A", "B"):
            if state not in self.states_idx:
                continue
            state_vac, state_burnn = self.get_state(atomic_numbers, positions, state)

            energy_vac, forces_vac = self.predict_energy_and_forces(state_vac)
            r_burnn_A = torch.tensor(
                state_burnn.positions,
                dtype=torch.float32,
                device=self.torchdevice,
                requires_grad=True,
            )
            # Each endpoint needs its own OR coordinate leaf because autograd
            # consumes the graph during the endpoint force evaluation.
            endpoint_r_or_A = r_or_A.detach().clone().requires_grad_(True)
            components = {} if debug_dhdl else None
            (
                energy_burnn,
                forces_burnn,
                forces_or,
                charges_burnn,
                energy_mlp_burnn,
                forces_mlp_burnn,
            ) = self._b2_eval_model(
                calculator=self.pred_calculator,
                system=state_burnn,
                r_qeq_A=r_burnn_A,
                r_or_A=endpoint_r_or_A,
                q_or=q_or,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
                **({"energy_components": components} if debug_dhdl else {}),
            )
            if debug_dhdl:
                endpoint_components[state] = components
            self.states[state] = {
                "vac": state_vac,
                "burnn": state_burnn,
                "energy_vac": energy_vac,
                "forces_vac": forces_vac,
                "energy_burnn": energy_burnn,
                "forces_burnn": forces_burnn,
                "energy_mlp_burnn": energy_mlp_burnn,
                "forces_mlp_burnn": forces_mlp_burnn,
                "forces_or": forces_or,
                "charges_burnn": charges_burnn,
                "burnn_full_indices": self._burnn_full_indices(state, n_total),
            }

        self.energy, self.derivative = self.calculate_perturbed_energy_and_derivative()

        if debug_dhdl:
            self._print_dynamic_dhdl(time_step, endpoint_components)
        self.forces = self.calculate_perturbed_forces()

        # dE/dR_OR and dE/dphi use exactly the BuRNN endpoint coefficients.
        self.or_forces = np.zeros((len(or_charges_e), 3), dtype=float)
        self.charges = np.zeros(n_total, dtype=float)
        endpoint_weights = {"A": 1.0 - self.lam, "B": self.lam}
        for state, data in self.states.items():
            weight = endpoint_weights[state]
            self.or_forces += weight * data["forces_or"]
            for local_index, full_index in enumerate(data["burnn_full_indices"]):
                self.charges[full_index] += weight * data["charges_burnn"][local_index]

        # Keep the established validation output usable. It remains an
        # MLP-only committee diagnostic, consistent with the non-perturbed B2
        # path, and is deliberately separate from the integrated embedding.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            if getattr(self, "direct_endpoint_interpolation", False):
                weights = {"A": 1.0 - self.lam, "B": self.lam}

                def combine_mlp(endpoint_data):
                    energy = 0.0
                    forces = np.zeros((n_total, 3), dtype=float)
                    for endpoint, data in endpoint_data.items():
                        weight = weights[endpoint]
                        energy += weight * data["energy_mlp_burnn"]
                        for local_index, full_index in enumerate(
                            data["burnn_full_indices"]
                        ):
                            forces[full_index] += (
                                weight * data["forces_mlp_burnn"][local_index]
                            )
                    return float(energy), forces

                production_mlp_energy, production_mlp_forces = combine_mlp(
                    self.states
                )
                committee_energies = [production_mlp_energy]
                committee_forces = [
                    np.linalg.norm(production_mlp_forces, axis=1)
                ]

                for val_calculator in self.val_calculators:
                    validation_states = {}
                    for endpoint, data in self.states.items():
                        rq = torch.tensor(
                            data["burnn"].positions,
                            dtype=torch.float32,
                            device=self.torchdevice,
                            requires_grad=True,
                        )
                        ro = r_or_A.detach().clone().requires_grad_(True)
                        _, _, _, _, energy_mlp, forces_mlp = self._b2_eval_model(
                            calculator=val_calculator,
                            system=data["burnn"],
                            r_qeq_A=rq,
                            r_or_A=ro,
                            q_or=q_or,
                            cutoff_nm=cutoff_nm,
                            n_link_atoms=n_link_atoms,
                        )
                        validation_states[endpoint] = {
                            "energy_mlp_burnn": energy_mlp,
                            "forces_mlp_burnn": forces_mlp,
                            "burnn_full_indices": data["burnn_full_indices"],
                        }
                    val_energy, val_forces = combine_mlp(validation_states)
                    committee_energies.append(val_energy)
                    committee_forces.append(np.linalg.norm(val_forces, axis=1))

                if len(committee_energies) == 2:
                    self.nn_valid_ene = (
                        committee_energies[0] - committee_energies[1]
                    ) / np.sqrt(2.0)
                else:
                    self.nn_valid_ene = float(
                        np.asarray(committee_energies).std(ddof=1)
                    )

                mean_force = np.mean(committee_forces, axis=0)
                sigma_force = np.sqrt(
                    np.mean(
                        [(force - mean_force) ** 2 for force in committee_forces],
                        axis=0,
                    )
                )
                self.nn_valid_maxF = float(np.max(sigma_force))
                self.nn_valid_mlp = self.nn_valid_ene
                self.nn_valid_mlp_maxF = self.nn_valid_maxF
            else:
                self.nn_valid_ene = self.validate_perturbed_energy()
                self.nn_valid_maxF = 0.0

    def calculate_next_step(
        self,
        atomic_numbers: list,
        positions: list,
        time_step: int,
        dynamic_charges: bool = False,
        or_positions_nm: list | None = None,
        or_charges_e: list | None = None,
        cutoff_nm: float | None = None,
        n_link_atoms: int = 0,
    ) -> None:
        """
        Calculate the next step in the MLP/MM simulation by updating the states, predicting energies and forces,
        and calculating the perturbed energy and forces.

        Args:
            atomic_numbers (list): List of atomic numbers for the atoms in the system.
            positions (list): List of positions for the atoms in the system.
            time_step (int): The current time step of the simulation.

        Returns:
            None
        """
        if dynamic_charges:
            if or_positions_nm is None or or_charges_e is None or cutoff_nm is None:
                raise ValueError(
                    "dynamic_charges=True requires or_positions_nm, "
                    "or_charges_e, and cutoff_nm"
                )
            return self._calculate_next_step_dynamic(
                atomic_numbers=atomic_numbers,
                positions=positions,
                time_step=time_step,
                or_positions_nm=or_positions_nm,
                or_charges_e=or_charges_e,
                cutoff_nm=cutoff_nm,
                n_link_atoms=n_link_atoms,
            )

        # Initialize the states dictionary to store vacuum and burnn states for each end state.
        self.states = {}

        # Iterate over the number of states (A and B) to populate the states dictionary.
        for i in range(self.num_states):
            if i == 0:
                state_id = 'A'
                self.states[state_id] = {'vac': [], 'burnn': []}
                # Retrieve the vacuum and burnn states for state A.
                self.states[state_id]['vac'], self.states[state_id]['burnn'] = self.get_state(atomic_numbers, positions, state_id=state_id)
            else:
                state_id = 'B'
                self.states[state_id] = {'vac': [], 'burnn': []}
                # Retrieve the vacuum and burnn states for state B.
                self.states[state_id]['vac'], self.states[state_id]['burnn'] = self.get_state(atomic_numbers, positions, state_id=state_id)

        # Predict energies and forces for each state in both vacuum and burnn environments.
        for state in self.states.keys():
            self.states[state]['energy_vac'], self.states[state]['forces_vac'] = self.predict_energy_and_forces(system=self.states[state]['vac'])
            self.states[state]['energy_burnn'], self.states[state]['forces_burnn'] = self.predict_energy_and_forces(self.states[state]['burnn'])
            print(f"{state}: forces Vac:\n{self.states[state]['forces_vac']}")
            print(f"{state}: forces BuRNN:\n{self.states[state]['forces_burnn']}")
        # Calculate the perturbed energy and its derivative using the predicted energies.
        self.energy, self.derivative = self.calculate_perturbed_energy_and_derivative()
        self.forces = self.calculate_perturbed_forces()
        print(self.forces)
        # Validate the perturbed energy if the current time step matches the validation frequency.
        if len(self.val_calculators) > 0 and time_step % self.nn_valid_freq == 0:
            self.nn_valid_ene = self.validate_perturbed_energy()
        return None


    def get_derivative(self):
        """
        Get the derivative of the perturbed energy.

        Returns:
            float: Derivative of the perturbed energy.
        """
        return self.derivative
