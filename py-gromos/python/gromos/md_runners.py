"""
High-level MD simulation runners for Python

This module provides Python interfaces to run molecular dynamics simulations
using different methods: standard MD, GaMD, EDS, and REMD.
"""

import subprocess
import os
import tempfile
import json
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import numpy as np

try:
    import gromos
except ImportError:
    raise ImportError(
        "gromos module not found. Please build with: maturin develop --release"
    )


class MDSimulation:
    """
    High-level interface for running MD simulations.

    This class wraps the GROMOS-RS command-line tools and provides
    a Pythonic interface for setting up and running simulations.
    """

    def __init__(self, topology_file: str, coordinate_file: str,
                 input_file: str, output_prefix: str = "md"):
        """
        Initialize MD simulation.

        Parameters
        ----------
        topology_file : str
            Path to topology file (.top)
        coordinate_file : str
            Path to coordinate file (.cnf)
        input_file : str
            Path to input parameter file (.imd)
        output_prefix : str
            Prefix for output files
        """
        self.topology_file = Path(topology_file)
        self.coordinate_file = Path(coordinate_file)
        self.input_file = Path(input_file)
        self.output_prefix = output_prefix

        # Check files exist
        for f in [self.topology_file, self.coordinate_file, self.input_file]:
            if not f.exists():
                raise FileNotFoundError(f"File not found: {f}")

        # Find GROMOS-RS binary
        self.md_binary = self._find_binary("md")

    def _find_binary(self, name: str) -> Path:
        """Find GROMOS-RS binary"""
        # Check in gromos-rs/target/release
        release_path = Path(__file__).parent.parent.parent / "gromos-rs" / "target" / "release" / name
        if release_path.exists():
            return release_path

        # Check in gromos-rs/target/debug
        debug_path = Path(__file__).parent.parent.parent / "gromos-rs" / "target" / "debug" / name
        if debug_path.exists():
            return debug_path

        # Check in PATH
        import shutil
        binary = shutil.which(name)
        if binary:
            return Path(binary)

        raise FileNotFoundError(
            f"Could not find '{name}' binary. "
            f"Please build gromos-rs first: cd gromos-rs && cargo build --release"
        )

    def run(self, steps: Optional[int] = None, verbose: bool = True) -> Dict[str, Path]:
        """
        Run standard MD simulation.

        Parameters
        ----------
        steps : int, optional
            Number of MD steps (overrides input file)
        verbose : bool
            Print output to console

        Returns
        -------
        dict
            Dictionary with paths to output files
        """
        cmd = [
            str(self.md_binary),
            str(self.topology_file),
            str(self.coordinate_file),
            str(self.input_file),
        ]

        if steps is not None:
            cmd.extend(["--steps", str(steps)])

        print(f"Running MD simulation...")
        print(f"Command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=not verbose,
            text=True,
            check=True
        )

        if not verbose and result.stdout:
            print(result.stdout)

        # Return output file paths
        return {
            'trajectory': Path(f"{self.output_prefix}.trc"),
            'energy': Path(f"{self.output_prefix}.tre"),
            'final_config': Path(f"{self.output_prefix}_final.cnf"),
        }


class GaMDSimulation(MDSimulation):
    """
    Gaussian Accelerated Molecular Dynamics (GaMD) simulation.

    GaMD adds a harmonic boost potential to smooth the energy landscape,
    accelerating sampling of rare events.
    """

    def __init__(self, topology_file: str, coordinate_file: str,
                 input_file: str, output_prefix: str = "gamd",
                 sigma0: float = 6.0, threshold_mode: str = "lower"):
        """
        Initialize GaMD simulation.

        Parameters
        ----------
        sigma0 : float
            Standard deviation of boost potential (kJ/mol)
        threshold_mode : str
            'lower' or 'upper' - when to apply boost
        """
        super().__init__(topology_file, coordinate_file, input_file, output_prefix)
        self.gamd_binary = self._find_binary("gamd")
        self.sigma0 = sigma0
        self.threshold_mode = threshold_mode

    def run(self, equilibration_steps: int = 10000, production_steps: int = 100000,
            verbose: bool = True) -> Dict[str, Path]:
        """
        Run GaMD simulation.

        Parameters
        ----------
        equilibration_steps : int
            Steps for equilibration (collecting statistics)
        production_steps : int
            Steps for production with boost potential
        verbose : bool
            Print output to console

        Returns
        -------
        dict
            Dictionary with paths to output files
        """
        cmd = [
            str(self.gamd_binary),
            str(self.topology_file),
            str(self.coordinate_file),
            str(self.input_file),
            "--sigma0", str(self.sigma0),
            "--threshold", self.threshold_mode,
            "--equilibration", str(equilibration_steps),
            "--production", str(production_steps),
        ]

        print(f"Running GaMD simulation...")
        print(f"  Equilibration: {equilibration_steps} steps")
        print(f"  Production: {production_steps} steps")
        print(f"  Sigma0: {self.sigma0} kJ/mol")
        print(f"Command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=not verbose,
            text=True,
            check=True
        )

        if not verbose and result.stdout:
            print(result.stdout)

        return {
            'trajectory': Path(f"{self.output_prefix}.trc"),
            'energy': Path(f"{self.output_prefix}.tre"),
            'boost': Path(f"{self.output_prefix}_boost.dat"),
            'statistics': Path(f"{self.output_prefix}_stats.dat"),
        }


class EDSSimulation(MDSimulation):
    """
    Enveloping Distribution Sampling (EDS) simulation.

    EDS uses a smoothed energy landscape to enhance sampling of
    multiple end-states simultaneously.
    """

    def __init__(self, topology_file: str, coordinate_file: str,
                 input_file: str, output_prefix: str = "eds",
                 num_states: int = 2, smoothness: float = 1.0):
        """
        Initialize EDS simulation.

        Parameters
        ----------
        num_states : int
            Number of end-states
        smoothness : float
            Smoothness parameter (kJ/mol)
        """
        super().__init__(topology_file, coordinate_file, input_file, output_prefix)
        self.eds_binary = self._find_binary("eds")
        self.num_states = num_states
        self.smoothness = smoothness

    def run(self, steps: int = 100000, verbose: bool = True) -> Dict[str, Path]:
        """
        Run EDS simulation.

        Parameters
        ----------
        steps : int
            Number of simulation steps
        verbose : bool
            Print output to console

        Returns
        -------
        dict
            Dictionary with paths to output files
        """
        cmd = [
            str(self.eds_binary),
            str(self.topology_file),
            str(self.coordinate_file),
            str(self.input_file),
            "--num-states", str(self.num_states),
            "--smoothness", str(self.smoothness),
            "--steps", str(steps),
        ]

        print(f"Running EDS simulation...")
        print(f"  States: {self.num_states}")
        print(f"  Smoothness: {self.smoothness} kJ/mol")
        print(f"  Steps: {steps}")
        print(f"Command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=not verbose,
            text=True,
            check=True
        )

        if not verbose and result.stdout:
            print(result.stdout)

        return {
            'trajectory': Path(f"{self.output_prefix}.trc"),
            'energy': Path(f"{self.output_prefix}.tre"),
            'free_energy': Path(f"{self.output_prefix}.dlg"),
            'state_probabilities': Path(f"{self.output_prefix}_prob.dat"),
        }


class REMDSimulation:
    """
    Replica Exchange Molecular Dynamics (REMD) simulation.

    REMD runs multiple replicas at different temperatures and
    periodically attempts to exchange configurations.
    """

    def __init__(self, topology_file: str, coordinate_file: str,
                 input_file: str, temperatures: List[float],
                 output_prefix: str = "remd"):
        """
        Initialize REMD simulation.

        Parameters
        ----------
        topology_file : str
            Path to topology file
        coordinate_file : str
            Path to coordinate file
        input_file : str
            Path to input parameter file
        temperatures : list of float
            Temperatures for each replica (K)
        output_prefix : str
            Prefix for output files
        """
        self.topology_file = Path(topology_file)
        self.coordinate_file = Path(coordinate_file)
        self.input_file = Path(input_file)
        self.temperatures = temperatures
        self.output_prefix = output_prefix
        self.num_replicas = len(temperatures)

        # Check files exist
        for f in [self.topology_file, self.coordinate_file, self.input_file]:
            if not f.exists():
                raise FileNotFoundError(f"File not found: {f}")

        # Find binary
        self.remd_binary = self._find_binary("remd")

    def _find_binary(self, name: str) -> Path:
        """Find GROMOS-RS binary"""
        release_path = Path(__file__).parent.parent.parent / "gromos-rs" / "target" / "release" / name
        if release_path.exists():
            return release_path

        debug_path = Path(__file__).parent.parent.parent / "gromos-rs" / "target" / "debug" / name
        if debug_path.exists():
            return debug_path

        import shutil
        binary = shutil.which(name)
        if binary:
            return Path(binary)

        raise FileNotFoundError(f"Could not find '{name}' binary")

    def run(self, steps: int = 100000, exchange_interval: int = 1000,
            verbose: bool = True) -> Dict[str, List[Path]]:
        """
        Run REMD simulation.

        Parameters
        ----------
        steps : int
            Number of steps per replica
        exchange_interval : int
            Steps between exchange attempts
        verbose : bool
            Print output to console

        Returns
        -------
        dict
            Dictionary with lists of output files for each replica
        """
        # Create temperature file
        temp_file = Path(f"{self.output_prefix}_temperatures.dat")
        with open(temp_file, 'w') as f:
            for i, T in enumerate(self.temperatures):
                f.write(f"{i} {T}\n")

        cmd = [
            str(self.remd_binary),
            str(self.topology_file),
            str(self.coordinate_file),
            str(self.input_file),
            "--temperatures", str(temp_file),
            "--steps", str(steps),
            "--exchange-interval", str(exchange_interval),
        ]

        print(f"Running REMD simulation...")
        print(f"  Replicas: {self.num_replicas}")
        print(f"  Temperatures: {min(self.temperatures):.1f} - {max(self.temperatures):.1f} K")
        print(f"  Steps: {steps}")
        print(f"  Exchange interval: {exchange_interval}")
        print(f"Command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=not verbose,
            text=True,
            check=True
        )

        if not verbose and result.stdout:
            print(result.stdout)

        # Collect output files
        trajectories = [Path(f"{self.output_prefix}_replica{i}.trc")
                       for i in range(self.num_replicas)]
        energies = [Path(f"{self.output_prefix}_replica{i}.tre")
                   for i in range(self.num_replicas)]

        return {
            'trajectories': trajectories,
            'energies': energies,
            'exchange_log': Path(f"{self.output_prefix}_exchanges.log"),
            'statistics': Path(f"{self.output_prefix}_stats.dat"),
        }


def run_standard_md(topology: str, coordinates: str, input_file: str,
                    steps: int = 10000, output_prefix: str = "md") -> Dict[str, Path]:
    """
    Convenience function to run standard MD simulation.

    Parameters
    ----------
    topology : str
        Path to topology file (.top)
    coordinates : str
        Path to coordinate file (.cnf)
    input_file : str
        Path to input parameter file (.imd)
    steps : int
        Number of MD steps
    output_prefix : str
        Prefix for output files

    Returns
    -------
    dict
        Dictionary with paths to output files

    Examples
    --------
    >>> from gromos.md_runners import run_standard_md
    >>> outputs = run_standard_md(
    ...     topology="system.top",
    ...     coordinates="start.cnf",
    ...     input_file="md.imd",
    ...     steps=10000
    ... )
    >>> print(f"Trajectory: {outputs['trajectory']}")
    """
    sim = MDSimulation(topology, coordinates, input_file, output_prefix)
    return sim.run(steps=steps)


def run_gamd(topology: str, coordinates: str, input_file: str,
             equilibration_steps: int = 10000, production_steps: int = 100000,
             sigma0: float = 6.0, threshold_mode: str = "lower",
             output_prefix: str = "gamd") -> Dict[str, Path]:
    """
    Convenience function to run GaMD simulation.

    Parameters
    ----------
    topology : str
        Path to topology file
    coordinates : str
        Path to coordinate file
    input_file : str
        Path to input parameter file
    equilibration_steps : int
        Steps for equilibration phase
    production_steps : int
        Steps for production phase
    sigma0 : float
        Standard deviation of boost (kJ/mol)
    threshold_mode : str
        'lower' or 'upper'
    output_prefix : str
        Prefix for output files

    Returns
    -------
    dict
        Dictionary with paths to output files

    Examples
    --------
    >>> from gromos.md_runners import run_gamd
    >>> outputs = run_gamd(
    ...     topology="system.top",
    ...     coordinates="start.cnf",
    ...     input_file="gamd.imd",
    ...     equilibration_steps=10000,
    ...     production_steps=100000,
    ...     sigma0=6.0
    ... )
    """
    sim = GaMDSimulation(topology, coordinates, input_file, output_prefix,
                        sigma0=sigma0, threshold_mode=threshold_mode)
    return sim.run(equilibration_steps=equilibration_steps,
                  production_steps=production_steps)


def run_eds(topology: str, coordinates: str, input_file: str,
            num_states: int = 2, smoothness: float = 1.0,
            steps: int = 100000, output_prefix: str = "eds") -> Dict[str, Path]:
    """
    Convenience function to run EDS simulation.

    Parameters
    ----------
    topology : str
        Path to topology file
    coordinates : str
        Path to coordinate file
    input_file : str
        Path to input parameter file
    num_states : int
        Number of end-states
    smoothness : float
        Smoothness parameter (kJ/mol)
    steps : int
        Number of simulation steps
    output_prefix : str
        Prefix for output files

    Returns
    -------
    dict
        Dictionary with paths to output files

    Examples
    --------
    >>> from gromos.md_runners import run_eds
    >>> outputs = run_eds(
    ...     topology="system.top",
    ...     coordinates="start.cnf",
    ...     input_file="eds.imd",
    ...     num_states=4,
    ...     smoothness=1.0,
    ...     steps=100000
    ... )
    """
    sim = EDSSimulation(topology, coordinates, input_file, output_prefix,
                       num_states=num_states, smoothness=smoothness)
    return sim.run(steps=steps)


def run_remd(topology: str, coordinates: str, input_file: str,
             temperatures: List[float], steps: int = 100000,
             exchange_interval: int = 1000,
             output_prefix: str = "remd") -> Dict[str, List[Path]]:
    """
    Convenience function to run REMD simulation.

    Parameters
    ----------
    topology : str
        Path to topology file
    coordinates : str
        Path to coordinate file
    input_file : str
        Path to input parameter file
    temperatures : list of float
        Temperatures for each replica (K)
    steps : int
        Number of steps per replica
    exchange_interval : int
        Steps between exchange attempts
    output_prefix : str
        Prefix for output files

    Returns
    -------
    dict
        Dictionary with lists of output files

    Examples
    --------
    >>> from gromos.md_runners import run_remd
    >>> temperatures = [300, 310, 320, 330, 345, 360, 375, 390]
    >>> outputs = run_remd(
    ...     topology="system.top",
    ...     coordinates="start.cnf",
    ...     input_file="remd.imd",
    ...     temperatures=temperatures,
    ...     steps=100000,
    ...     exchange_interval=1000
    ... )
    """
    sim = REMDSimulation(topology, coordinates, input_file, temperatures, output_prefix)
    return sim.run(steps=steps, exchange_interval=exchange_interval)


class TISimulation:
    """
    Thermodynamic Integration (TI) simulation.

    TI calculates free energy differences by integrating the derivative
    of the Hamiltonian with respect to the coupling parameter λ.
    """

    def __init__(self, topology_file: str, coordinate_file: str,
                 input_file: str, lambda_values: Optional[List[float]] = None,
                 output_prefix: str = "ti", lambda_type: str = "total",
                 soft_core: bool = False, soft_core_alpha: float = 0.5):
        """
        Initialize TI simulation.

        Parameters
        ----------
        topology_file : str
            Path to topology file (should be hybrid topology for TI)
        coordinate_file : str
            Path to coordinate file
        input_file : str
            Path to input parameter file
        lambda_values : list of float, optional
            Lambda values for windows. If None, uses default linear spacing.
        output_prefix : str
            Prefix for output files
        lambda_type : str
            Type of coupling: "total", "coulomb", "lj"
        soft_core : bool
            Use soft-core potentials (recommended for LJ)
        soft_core_alpha : float
            Soft-core parameter alpha
        """
        self.topology_file = Path(topology_file)
        self.coordinate_file = Path(coordinate_file)
        self.input_file = Path(input_file)
        self.output_prefix = output_prefix
        self.lambda_type = lambda_type
        self.soft_core = soft_core
        self.soft_core_alpha = soft_core_alpha

        # Default lambda values if not provided
        if lambda_values is None:
            # Use 11 evenly spaced points
            self.lambda_values = np.linspace(0.0, 1.0, 11)
        else:
            self.lambda_values = np.array(lambda_values)

        self.num_windows = len(self.lambda_values)

        # Check files exist
        for f in [self.topology_file, self.coordinate_file, self.input_file]:
            if not f.exists():
                raise FileNotFoundError(f"File not found: {f}")

        # Find binary
        self.ti_binary = self._find_binary("ti")

    def _find_binary(self, name: str) -> Path:
        """Find GROMOS-RS binary"""
        release_path = Path(__file__).parent.parent.parent / "gromos-rs" / "target" / "release" / name
        if release_path.exists():
            return release_path

        debug_path = Path(__file__).parent.parent.parent / "gromos-rs" / "target" / "debug" / name
        if debug_path.exists():
            return debug_path

        import shutil
        binary = shutil.which(name)
        if binary:
            return Path(binary)

        raise FileNotFoundError(f"Could not find '{name}' binary")

    def run(self, steps_per_window: int = 100000, equilibration_steps: int = 10000,
            verbose: bool = True) -> Dict[str, any]:
        """
        Run TI simulation for all λ windows.

        Parameters
        ----------
        steps_per_window : int
            Number of production steps per λ window
        equilibration_steps : int
            Number of equilibration steps per window
        verbose : bool
            Print output to console

        Returns
        -------
        dict
            Dictionary with lambda values, dH/dλ values, and output files
        """
        print(f"Running TI simulation...")
        print(f"  Windows: {self.num_windows}")
        print(f"  Lambda type: {self.lambda_type}")
        print(f"  Soft-core: {self.soft_core}")
        print(f"  Steps per window: {steps_per_window}")

        dhdl_values = []
        window_outputs = []

        for i, lam in enumerate(self.lambda_values):
            print(f"\n  Window {i+1}/{self.num_windows}: λ = {lam:.3f}")

            # Create lambda-specific input file
            lambda_input = self._create_lambda_input(lam)

            # Run MD for this lambda
            cmd = [
                str(self.ti_binary),
                str(self.topology_file),
                str(self.coordinate_file),
                str(lambda_input),
                "--lambda", str(lam),
                "--lambda-type", self.lambda_type,
                "--steps", str(steps_per_window),
                "--equilibration", str(equilibration_steps),
                "--output", f"{self.output_prefix}_lambda{lam:.3f}",
            ]

            if self.soft_core:
                cmd.extend(["--soft-core", "--soft-core-alpha", str(self.soft_core_alpha)])

            if verbose:
                print(f"    Command: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=not verbose,
                text=True,
                check=True
            )

            if not verbose and result.stdout:
                print(result.stdout)

            # Parse dH/dλ from output
            dhdl_file = Path(f"{self.output_prefix}_lambda{lam:.3f}_dhdl.dat")
            dhdl = self._parse_dhdl(dhdl_file)
            dhdl_values.append(dhdl)

            window_outputs.append({
                'lambda': lam,
                'trajectory': Path(f"{self.output_prefix}_lambda{lam:.3f}.trc"),
                'energy': Path(f"{self.output_prefix}_lambda{lam:.3f}.tre"),
                'dhdl_file': dhdl_file,
                'dhdl_mean': dhdl,
            })

            print(f"    ⟨∂H/∂λ⟩ = {dhdl:.2f} kJ/mol")

        return {
            'lambda_values': self.lambda_values,
            'dhdl_values': np.array(dhdl_values),
            'window_outputs': window_outputs,
        }

    def _create_lambda_input(self, lam: float) -> Path:
        """
        Create input file with specific lambda value.

        Parameters
        ----------
        lam : float
            Lambda value

        Returns
        -------
        Path
            Path to created input file
        """
        # Read base input file
        with open(self.input_file, 'r') as f:
            content = f.read()

        # Modify lambda value (simplified - actual implementation
        # would need to parse GROMOS input format properly)
        # For now, just copy the file
        lambda_input = Path(f"{self.output_prefix}_lambda{lam:.3f}.imd")
        with open(lambda_input, 'w') as f:
            f.write(content)
            f.write(f"\n# Lambda = {lam}\n")

        return lambda_input

    def _parse_dhdl(self, dhdl_file: Path) -> float:
        """
        Parse ∂H/∂λ from output file.

        Parameters
        ----------
        dhdl_file : Path
            Path to dH/dλ file

        Returns
        -------
        float
            Mean ∂H/∂λ value
        """
        if not dhdl_file.exists():
            raise FileNotFoundError(f"dH/dλ file not found: {dhdl_file}")

        # Read data (skip comments and equilibration)
        data = []
        with open(dhdl_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                try:
                    value = float(line.strip().split()[1])  # Assuming 2nd column
                    data.append(value)
                except (ValueError, IndexError):
                    continue

        if not data:
            raise ValueError(f"No valid data in {dhdl_file}")

        # Discard first 20% as equilibration
        n_equil = len(data) // 5
        data = data[n_equil:]

        return np.mean(data)

    def integrate(self, method: str = "trapezoid") -> Tuple[float, float]:
        """
        Integrate ∂H/∂λ to calculate free energy.

        Parameters
        ----------
        method : str
            Integration method: "trapezoid", "simpson", or "spline"

        Returns
        -------
        tuple
            (free_energy, error) in kJ/mol
        """
        # This would be called after run() to integrate the results
        # For now, raise NotImplementedError as it needs the results
        raise NotImplementedError(
            "Call run() first to generate dH/dλ data, then use "
            "numpy.trapz or scipy.integrate to integrate the results"
        )


def run_ti(topology: str, coordinates: str, input_file: str,
           num_windows: int = 11, lambda_spacing: str = "linear",
           steps_per_window: int = 100000,
           lambda_type: str = "total",
           soft_core: bool = False,
           output_prefix: str = "ti") -> Dict[str, any]:
    """
    Convenience function to run TI simulation.

    Parameters
    ----------
    topology : str
        Path to topology file (hybrid topology)
    coordinates : str
        Path to coordinate file
    input_file : str
        Path to input parameter file
    num_windows : int
        Number of λ windows
    lambda_spacing : str
        Spacing type: "linear", "nonlinear", or "optimized"
    steps_per_window : int
        Number of steps per λ window
    lambda_type : str
        Type of coupling: "total", "coulomb", "lj"
    soft_core : bool
        Use soft-core potentials
    output_prefix : str
        Prefix for output files

    Returns
    -------
    dict
        Dictionary with lambda values, dH/dλ values, and output files

    Examples
    --------
    >>> from gromos.md_runners import run_ti
    >>> outputs = run_ti(
    ...     topology="hybrid.top",
    ...     coordinates="start.cnf",
    ...     input_file="ti.imd",
    ...     num_windows=11,
    ...     steps_per_window=100000,
    ...     lambda_type="total"
    ... )
    >>> # Integrate using numpy
    >>> import numpy as np
    >>> dG = np.trapz(outputs['dhdl_values'], outputs['lambda_values'])
    >>> print(f"ΔG = {dG:.2f} kJ/mol")
    """
    # Generate lambda values based on spacing
    if lambda_spacing == "linear":
        lambda_values = np.linspace(0.0, 1.0, num_windows)
    elif lambda_spacing == "nonlinear":
        # More points near endpoints
        half = num_windows // 2
        lambda_values = np.concatenate([
            np.linspace(0.0, 0.5, half, endpoint=False)**2 * 0.5,
            0.5 + (np.linspace(0.0, 0.5, num_windows - half)**2 * 0.5)
        ])
    else:
        # Default to linear
        lambda_values = np.linspace(0.0, 1.0, num_windows)

    sim = TISimulation(
        topology, coordinates, input_file,
        lambda_values=lambda_values.tolist(),
        output_prefix=output_prefix,
        lambda_type=lambda_type,
        soft_core=soft_core
    )

    outputs = sim.run(steps_per_window=steps_per_window)

    # Calculate free energy using trapezoid rule
    dG = np.trapz(outputs['dhdl_values'], outputs['lambda_values'])

    # Simple error estimation (would need proper block averaging in practice)
    # This is a placeholder
    dG_error = np.std(outputs['dhdl_values']) / np.sqrt(len(outputs['dhdl_values']))

    outputs['free_energy'] = dG
    outputs['error'] = dG_error

    return outputs
