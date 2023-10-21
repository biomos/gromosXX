"""
This file is part of GROMOS.

Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
See <https://www.gromos.net> for details.

GROMOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

import os
import stat
import shutil
import yaml
import numpy as np

# default paths for GITLAB
_MD_BIN = '${CI_PROJECT_DIR}/${BIN_PATH}'
_CI_PROJECT_DIR = os.getenv('CI_PROJECT_DIR')
_ENE_ANA_LIB = f'{_CI_PROJECT_DIR}/md++/data/ene_ana.md++.lib'

_run_file_header = """#!/bin/sh

# first we set some variables
NAME=`whoami`
PROGRAM={md_bin}
#output log file
OUNIT={omd}

"""

_md_header = """MDOK=1

${{PROGRAM}} \\
@topo        {TOPO} \\
@conf        {INPUTCRD} \\
@input       {IUNIT} \\
"""

_run_file_tail = """>            ${OUNIT}
grep "finished successfully" ${OUNIT} > /dev/null || MDOK=0

uname -a >> ${OUNIT}"""

class ConfigTest():
    def __init__(self, config_path: str) -> None:
        """
        Creates executable run file from default values.
        -----
        input:
                config_path: str --> path to yaml file
        output:
                md.run: file --> gromos run file
        -----
        """
        self.expected_values = dict()
        self.config_path = config_path
        self.config_path_dir = os.path.dirname(self.config_path)
        self.load_config_yaml()
        self.get_run_string()
        #self.load_expected_values()
        #self.dump_full_config_yaml() #enable for local testing
        #self.run_string()
        #self.write_run_file()

    def __get_input_file_path(self, file_path):
        if os.path.exists(file_path):
            return file_path, os.path.abspath(file_path)
        else:
            file_path = os.path.join(self.input_dir, file_path)
            return file_path, os.path.abspath(file_path)

    def load_config_yaml(self) -> dict:
        """
        Loads content of yaml file into a python dict.
        -----
        input:
                self.config_path: str --> path to yaml file
        output:
                settings_loaded: dict --> params need for run file
                                          in a nested dict.
        -----
        """
        with open(self.config_path, "r") as stream:
            try:
                # load the yaml file
                settings_loaded = yaml.safe_load(stream)
                self.config = settings_loaded
                # input dir
                if self.config['input']['input_dir']:
                    self.input_dir_abs = os.path.abspath(self.config['input']['input_dir'])
                    self.input_dir_rel = self.config['input']['input_dir']
                else:
                    self.input_dir_abs = os.path.abspath(self.config_path_dir)
                    self.input_dir_rel = self.config_path_dir
                self.input_dir = self.input_dir_rel
                # simulation dir
                if self.config['input']['simulation_dir']:
                    self.input_dir = self.input_dir_abs
                    self.sim_dir = self.config['input']['simulation_dir']
                    if os.path.abspath(self.sim_dir) != self.sim_dir:
                        self.sim_dir = os.path.join(self.input_dir, self.sim_dir)
                else:
                    self.sim_dir = self.config_path_dir
                # ene ana lib
                if self.config['input']['ene_ana_lib']:
                    ene_ana_lib = self.__get_input_file_path(self.config['input']['ene_ana_lib'])
                    self.ene_ana_lib = ene_ana_lib[1]
                else:
                    self.ene_ana_lib = _ENE_ANA_LIB
                # get the binary
                if self.config['input']['md_binary']:
                    self.md_bin = self.config['input']['md_binary']
                else:
                    self.md_bin = _MD_BIN
                _pre_comm = os.getenv('PRE_MD_BIN')
                if _pre_comm: # pre_command (e.g. mpirun)
                    self.md_bin = '"${PRE_MD_BIN} ' + self.md_bin + '"'
            except yaml.YAMLError as exc:
                print(exc)
            return settings_loaded

    def dump_full_config_yaml(self):
        """
        Dumps content of dict into yaml file.
        """
        with open(f"{self.output}/used_inp.yaml", 'w') as file:
            yaml.dump(self.config, file)

    def get_run_string(self):
        """
        Inserts params from settings dict into run file string.
        -----
        input:
                self.md_bin: str --> path to md binary. Default path to md in docker image.
                self.sim_dir: str --> path of simulation dir. Default current dir.
                self.input_dir: str --> path to input dir containing input files.
                self.topo: str --> name of topology f.ile
                self.iunit: str --> name of imd file.
                self.icrd: str --> name of input cnf file.
                self.omd: str --> name for gromos output file.
                self.ocrd: str --> name of output cnf file.
                self.trc: str --> name of coordinate trajectory.
                self.tre: str --> name of energy trajectory.
                self.trf: str --> name of force trajectory.
                self.trv: str --> name of velocity trajectory.
        output:
                self.run_str: str --> formated string
        """
        _md_line_format = '@{:<12}{:} \\\n'
        out_files_dict = dict(OUTPUTCRD='fin', OUTPUTTRX='trc', OUTPUTTRE='tre', 
                              OUTPUTTRG='trg', OUTPUTTRF='trf', OUTPUTTRV='trv')
        in_files_dict = dict() # pos res for instance, not implemented yet...
        # get the run file
        # header
        self.omd = self.config['output']['OUNIT'] # out log file (omd)
        self.run_str = _run_file_header.format(md_bin=self.md_bin, sim_dir=self.sim_dir, omd=self.omd)
        # input files
        self.topo = os.path.join(self.input_dir, self.config['input']['TOPO'])
        self.iunit = os.path.join(self.input_dir, self.config['input']['IUNIT'])
        self.icrd = os.path.join(self.input_dir, self.config['input']['INPUTCRD'])
        self.run_str += _md_header.format(TOPO=self.topo, IUNIT=self.iunit, INPUTCRD=self.icrd)
        for var_name, flag in in_files_dict.items():
            if var_name in self.config['input'] and self.config['input'][var_name]:
                temp_file_name = os.path.join(self.input_dir, self.config['input'][var_name])
                self.run_str += _md_line_format.format(flag, temp_file_name)
        # output files
        self.out_files = {}
        for var_name, flag in out_files_dict.items():
            if var_name in self.config['output'] and self.config['output'][var_name]:
                temp_file = self.config['output'][var_name]
                self.run_str += _md_line_format.format(flag, temp_file)
                self.out_files[flag] = temp_file
        self.run_str += _run_file_tail

    def write_run_file(self, run_file="md.run", run_str=None, sim_dir=None):
        """
        Writes out run_string into a executable run file.
        -----
        input:
                run_file: str --> path to a file
                run_str: str --> formated string containing run file structure (taken from self.run)
                sim_dir: str --> simulation directory (also where run_file is made; taken from self.sim_dir)
        output:
                run_file: str  --> path to the gromos run file.
        -----
        """
        if sim_dir is None:
            sim_dir = self.sim_dir
        if sim_dir:
            if not os.path.isdir(sim_dir):
                os.mkdir(sim_dir)
            run_file = os.path.join(sim_dir, run_file)
        if run_str is None:
            run_str = self.run_str
        with open(run_file, 'w+') as f:
            f.write(run_str)
        st = os.stat(run_file)
        os.chmod(run_file, st.st_mode | stat.S_IEXEC)
        return os.path.abspath(run_file)

    def __store_expected_values(self, trj_type, temp_file):
        temp_file = self.__get_input_file_path(temp_file)
        with open(temp_file[1], "rb") as f:
            self.expected_values[trj_type] = np.load(f)
        return

    def load_expected_values(self, *args):
        """
        loads hardcoded expected values (based on the self.conif file paths)
        -----
        input:
                *args: list of file types (e.g. 'tre', 'trg')
        output:
                None (stores the values in self.expected_values)
        -----
        """
        if args:
            for trj_type in args:
                temp_file = self.config['expected_values'][trj_type]
                self.__store_expected_values(trj_type, temp_file)
        else:
            for trj_type, temp_file in self.config['expected_values'].items():
                self.__store_expected_values(trj_type, temp_file)

    def clean_up(self, flag_rm_input_dir=False):
        if os.path.abspath(self.sim_dir) != self.input_dir_abs or flag_rm_input_dir:
            if os.path.isdir(self.sim_dir):
                shutil.rmtree(self.sim_dir)
