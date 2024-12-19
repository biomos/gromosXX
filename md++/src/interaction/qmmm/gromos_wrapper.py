import subprocess as sbp
import pandas as pd
import numpy as np
import os
from glob import glob
from importlib import reload
from shutil import rmtree

import ase
import schnetpack as spk
import schnetpack.transform as trn
import pytorch_lightning as pl
import torch
import torchmetrics
import yaml

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
    

class calculator():
    def __init__(self, 
                 model_directory: str=None,
                 device: str='cuda',
                 ) -> None:

        self.model_dir=model_directory
        self.device=device
        
        # set device
        torchdevice = torch.device(self.device)

        # create ase Atoms object
        #self.atoms = Atoms(numbers=self.iacs, positions=self.positions)
        
        # set up calculator
        cutoff = YamlParser(os.path.join(self.model_dir,'config.yaml')).get_cutoff()
        model_path = os.path.join(self.model_dir,'best_model')
        
        self.calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            neighbor_list=trn.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key='BuRNN_energy', # name of energy property in model
            force_key='BuRNN_forces', # name of force property in model
            energy_unit="eV", # units of energy property predicted by ase
            device=torchdevice, # device for computation
        )
        
        return None
    
    def spk_calculator(self):
        return self.calculator
        
        

class gromos2MLP():
    def __init__(self, 
                 model_directory: str=None,
                 device: str='cuda',
                 positions: list=None,
                 iacs: list=None,
                 
                 ) -> None:

        self.model_dir=model_directory
        self.device=device
        self.positions=positions
        self.iacs=iacs
        self.calculator=None
        
        # set device
        self.torchdevice = torch.device(self.device)

    def ase_atoms(self):
    	pos_arr = np.asarray(self.positions).reshape(int(len(self.positions)/3),3)
        iac_arr = np.asarray(self.iacs)
    	
    	atoms = ase.Atoms(numbers=iac_arr, positions=pos_arr)
    	
    	return atoms
 
    
    def calculate_energies(self):
    
        atoms = self.ase_atoms()
        atoms.set_calculator(self.calculator)
        
        return atoms.get_total_energy()
    
    def calculate_forces(self):
    
        atoms = self.ase_atoms()
        atoms.set_calculator(self.calculator)

        return atoms.get_forces()

        
        
class gromos2schnet_v2(gromos2MLP):
    def __init__(self, 
                 model_directory: str=None,
                 device: str='cuda'
                 ) -> None:

        self.model_dir=model_directory
        self.device=device
        # set device
        self.torchdevice = torch.device(self.device)
        
    def spk_calculator(self):
        # set up calculator
        cutoff = YamlParser(os.path.join(self.model_dir,'config.yaml')).get_cutoff()
        model_path = os.path.join(self.model_dir,'best_model')
        
        self.calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            neighbor_list=trn.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key='BuRNN_energy', # name of energy property in model
            force_key='BuRNN_forces', # name of force property in model
            energy_unit="eV", # units of energy property predicted by ase
            device= self.torchdevice, # device for computation
        )
        
        
    
class gromos2schnet_v1(gromos2MLP):
    def __init__(self, 
                 model_directory: str=None,
                 device: str='cuda',
                 ) -> None:

        self.model_dir=model_directory
        self.device=device
        
        # set device
        self.torchdevice = torch.device(self.device)
        
        
    def spk_calculator(self):
        # set up calculator
        cutoff = YamlParser(os.path.join(self.model_dir,'config.yaml')).get_cutoff()
        model_path = os.path.join(self.model_dir,'best_model')
        
        self.calculator = spk.interfaces.SpkCalculator(
            model_file = model_path, # path to model
            dtype=torch.float32, # 
            neighbor_list=trn.ASENeighborList(cutoff=cutoff), # neighbor list
            energy_key='BuRNN_energy', # name of energy property in model
            force_key='BuRNN_forces', # name of force property in model
            energy_unit="eV", # units of energy property predicted by ase
            device= self.torchdevice, # device for computation
        )   
