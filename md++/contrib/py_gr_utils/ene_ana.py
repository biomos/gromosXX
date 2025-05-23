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
from pathlib import Path
import numpy as np
try:
    import pandas as pd
except:
    pd = None
import re
import gzip
from dataclasses import dataclass

@dataclass
class EnergyTrajectory():
    """
    A class that reads and stores all informations form a (free) energy trajectory.
    """
    lib_path:str=None
    trj_files:list=None
    num_type:type=np.double

    def __post_init__(self):
        # initialize some vars
        self._data_type = None
        self._data_type_dict = dict()
        self._trj = list()
        if self.lib_path:
            self.read_ene_ana_lib(self.lib_path)
        if self.trj_files:
            assert self.lib_path
            if isinstance(self.trj_files, str):
                self.trj_files = [self.trj_files]
            self.read_trj(*self.trj_files, flag_update_trjFs=False)

    def read_ene_ana_lib(self, lib_path):
        self.lib_path = lib_path
        self.variables = dict()
        read_flag = False
        self.param_dict = {}
        self.size_dict = {} 
        with open(self.lib_path, 'r') as inp:
            for line in inp:
                if line.startswith('#'):
                    continue
                elif line.startswith('ENERTRJ'):
                    read_flag = 1
                elif line.startswith('FRENERTRJ'):
                    read_flag = 1
                elif line.startswith('VARIABLES'):
                    read_flag = 2
                if line.startswith('END'):
                    read_flag = False
                elif read_flag==1:
                    parts = line.split()
                    if parts[0] == 'block':
                        self.param_dict[parts[1]] = {}
                        self.size_dict[parts[1]] = {}
                        block = parts[1]
                    if parts[0] == 'subblock':    
                        self.param_dict[block][parts[1]] = []
                        if parts[2].isnumeric():
                            self.param_dict[block][parts[1]].append(int(parts[2]))
                        else:
                            self.param_dict[block][parts[1]].append(str(parts[2]))
                        if parts[3].isnumeric():
                            self.param_dict[block][parts[1]].append(int(parts[3]))
                        else:
                            self.param_dict[block][parts[1]].append(str(parts[3]))
                    elif parts[0] == 'size':
                        self.size_dict[block][parts[1]] = None
                        self.param_dict[block][parts[1]] = None
                elif read_flag==2 and not line.startswith('VARIABLES'):
                    temp_match = re.match('(.*?)=(.*?)\[', line)
                    if temp_match:
                        v_name = temp_match.group(1).strip()
                        subbl_name = temp_match.group(2).strip()
                        idx = tuple(int(i)-1 for i in re.findall('\[(\d*)', line))
                        self.variables[v_name] = (subbl_name, idx)
        self.block_names = list(self.param_dict)

    def __get_values_block(self, block_name, value_list):
        data_all = list()
        shapes_all = list()
        pos = 0
        for subblock, temp_shape in self.param_dict[block_name].items():
            if temp_shape:
                for i, temp_value in enumerate(temp_shape):
                    if isinstance(temp_value, str):
                        temp_shape[i] = self.size_dict[block_name][temp_value]
                if temp_shape[-1] == 1:
                    temp_shape = temp_shape[:-1]
                temp_shape = tuple(temp_shape)
                num_values = 1
                for temp_value in temp_shape:
                    num_values *= temp_value
                new_pos = pos + num_values
                temp_data = np.array(value_list[pos:new_pos]).reshape(temp_shape)
                data_all.append(temp_data)
                if self._data_type is None:
                    shapes_all.append((subblock, self.num_type, temp_shape))
                pos = new_pos
            else:
                self.size_dict[block_name][subblock] = int(value_list[pos])
                if subblock=="NUM_ENERGY_GROUPS":
                    temp_value = self.size_dict[block_name][subblock]
                    self.size_dict[block_name]["matrix_NUM_ENERGY_GROUPS"] = temp_value * (temp_value+1) // 2
                pos += 1
        if self._data_type is None:
            self._data_type_dict[block_name] = np.dtype(shapes_all)
        return data_all

    def __create_dtype(self):
        subblock_data_types = list()
        for block_name, temp_data_type in self._data_type_dict.items():
            subblock_data_types.append((block_name, temp_data_type))
        self._data_type = np.dtype(subblock_data_types)

    def __store_step_data(self, data_step):
        if self._data_type is None:
            self.__create_dtype()
        formated_data = list()
        for i, temp_data in enumerate(data_step):
            temp_data = np.array(tuple(temp_data), dtype=self._data_type[i])
            formated_data.append(temp_data)
        temp_step_data = np.array(tuple(formated_data), dtype=self._data_type)
        self._trj.append(temp_step_data)

    def _trj_line_generator(self, trj_path):
        try:
            with gzip.open(trj_path) as inp:
                for l in inp:
                    line = l.decode('ascii')
                    yield line
        except:
            with open(trj_path) as inp:
                for line in inp:
                    yield line

    def read_trj(self, *trj_files, flag_update_trjFs=True):
        for trj_path in trj_files:
            if flag_update_trjFs:
                self.trj_files.append(trj_path)
            data_step = None
            flag = False
            value_list = list()
            for line in self._trj_line_generator(trj_path):
                if line.startswith("END"):
                    if flag:
                        temp_data = self.__get_values_block(flag, value_list)
                        data_step.append(temp_data)
                    flag=False
                line = line.strip()
                pos = line.find("#")
                if pos!=-1:
                    line = line[:pos]
                if line and flag:
                    value_list += line.split()
                elif line in self.block_names:
                    flag = line
                    value_list = list()
                    if flag=="TIMESTEP":
                        if data_step:
                            self.__store_step_data(data_step)
                        data_step = list()
                    continue
            if data_step:
                self.__store_step_data(data_step)
        self.trj = np.array(self._trj, dtype=self._data_type)
        self.get_subblock_map()

    @staticmethod
    def __fix_ext(fout, ext):
        if not fout.endswith(ext):
            fout += ext
        return fout

    def write_trj_np(self, fout):
        fout = self.__fix_ext(fout, ext='.npz')
        with open(fout, "wb") as f:
            np.savez_compressed(f, trj=self.trj, var_dict=self.variables,
                               size_dict=self.size_dict, param_dict=self.param_dict)
        return

    def load_trj_np(self, fin, flag_load_var=False, allow_pickle=False):
        if flag_load_var:
            allow_pickle=True
        with open(fin, "rb") as f:
            temp = np.load(f, allow_pickle=allow_pickle)
            self.trj = temp['trj']
            self._data_type = self.trj.dtype
            if allow_pickle:
                if 'var_dict' in temp:
                    self.variables = temp['var_dict'].tolist()
                if 'size_dict' in temp:
                    self.size_dict = temp['size_dict'].tolist()
                if 'param_dict' in temp:
                    self.param_dict = temp['param_dict'].tolist()
        self.get_subblock_map()

    def get_subblock_map(self):
        self.subbl_map = {}
        for bl in self._data_type.fields:
            for subbl in self._data_type[bl].fields:
                self.subbl_map[subbl] = bl
        return self.subbl_map

    def __get_values_from_var(self, var):
        subbl, idx = self.variables[var]
        return self.extract_values(subbl, idx)

    def extract_values(self, subbl, idx=()):
        if not isinstance(idx, tuple):
            idx = (idx,)
        bl = self.subbl_map[subbl]
        data = self.trj[bl][subbl][(slice(None, None, None),) + idx]
        return data

    def get_values(self, *var, flag_df=False):
        data = []
        for v in var:
            data.append(self.__get_values_from_var(v))
        data = np.array(data).T
        if flag_df:
            if pd:
                return pd.DataFrame(data, columns = var)
            else:
                print('pandas not available')
        return data

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="read GROMOS format (free) energy trajectory")
    parser.add_argument('-l', '--library', type=str, help="ene_ana library")
    parser.add_argument('--def_lib', default=False, action='store_true', help='use library file from the same git repo as the script')
    parser.add_argument('-e', '--en_files', type=str, nargs='+', help="energy trajectory files (GROMOS format)")
    parser.add_argument('-f', '--fr_files', type=str, nargs='+', help="free energy trajectory files (GROMOS format)")
    parser.add_argument('-d', '--data_t', type=str, default='d', choices=['s', 'd', 'h'],
                        help='data type to use: s-single (default), d-double, h-half')
    parser.add_argument('--np_tre', type=str, help="energy trajectory files (stored as np array)")
    parser.add_argument('--np_trg', type=str, help="free energy trajectory files (stored as np array)")
    parser.add_argument('--np_load_pickle', default=False, action='store_true',
                        help="sets allow_pickle=True when loading to be able to retrieve the variables")
    parser.add_argument('-o', '--out_file', type=str, help="out file to save loaded trajectory in np format")
    args=parser.parse_args()

    num_type_map = dict(d=np.double, s=np.single, h=np.half)
    num_type = num_type_map[args.data_t]

    if args.library is None and args.def_lib:
        fpath = os.path.abspath(__file__)
        md_fd = Path(fpath).parents[2]
        args.library = os.path.join(md_fd, 'data','ene_ana.md++.lib')
    # load GROMOS trajectory files
    if args.en_files:
        assert args.library
        tre = EnergyTrajectory(lib_path=args.library, trj_files=args.en_files, num_type=num_type)
        if args.out_file:
            tre.write_trj_np(args.out_file)
    if args.fr_files:
        assert args.library
        trg = EnergyTrajectory(lib_path=args.library, trj_files=args.fr_files, num_type=num_type)
        if args.out_file:
            trg.write_trj_np(args.out_file)

    # load np trajectory files
    if args.np_tre:
        tre = EnergyTrajectory()
        tre.load_trj_np(args.np_tre, allow_pickle=args.np_load_pickle)
    if args.np_trg:
        trg = EnergyTrajectory()
        trg.load_trj_np(args.np_trg, allow_pickle=args.np_load_pickle)
