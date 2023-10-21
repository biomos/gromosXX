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

import numpy as np

class EnegyTrajectory():
    def __init__(self, lib_path, ene_trj_path=None, fene_trj_path=None):
        """
        A class that reads and stores all informations form the tre.
        """
        
        # initialize some vars
        self._data_type = None
        self._data_type_dict = dict()
        self._tre = list()

        self.lib_path = lib_path
        self.read_ene_ana_lib()
        self.block_names = list(self.param_dict)
        if ene_trj_path:
            self.trj_path = ene_trj_path
            self.read_trj()
        elif fene_trj_path:
            self.trj_path = fene_trj_path
            self.trf = self.read_trj()

    def read_ene_ana_lib(self):  
        read_flag = 0
        self.param_dict = {}
        self.size_dict = {} 
        with open(self.lib_path, 'r') as inp:
            for line in inp:
                if line.startswith('#'):
                    continue
                elif line.startswith('ENERTRJ'):
                    read_flag = 1
                elif line.startswith('FRENERTRJ'):
                    read_flag = 2
                if line.startswith('END'):
                    read_flag = 0
                elif read_flag == 1:
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

    def __get_values_block(self, block_name, value_list):
        data_all = list()
        shapes_all = list()
        pos = 0
        for subblock, temp_shape in self.param_dict[block_name].items():
            if temp_shape:
                for i, temp_value in enumerate(temp_shape):
                    if isinstance(temp_value, str):
                        temp_shape[i] = self.size_dict[block_name][temp_value]
                temp_shape = tuple(temp_shape)
                num_values = 1
                for temp_value in temp_shape:
                    num_values *= temp_value
                new_pos = pos + num_values
                temp_data = np.array(value_list[pos:new_pos]).reshape(temp_shape)
                data_all.append(temp_data)
                if self._data_type is None:
                    shapes_all.append((subblock, float, temp_shape))
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
            subblock_data_types.append((block_name, temp_data_type, (1,)))
        self._data_type = np.dtype(subblock_data_types)

    def __store_step_data(self, data_step):
        if self._data_type is None:
            self.__create_dtype()
        formated_data = list()
        for i, temp_data in enumerate(data_step):
            temp_data = np.array(tuple(temp_data), dtype=self._data_type[i])
            formated_data.append(temp_data)
        temp_step_data = np.array(tuple(formated_data), dtype=self._data_type)
        self._tre.append(temp_step_data)

    def read_trj(self):
        data_step = None
        flag = False
        value_list = list()
        with open(self.trj_path, 'r') as inp:
            for line in inp:
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
            #print(data_step)
            self.__store_step_data(data_step)
        self.tre = np.array(self._tre, dtype=self._data_type)

if __name__ == "__main__":
    tr = Ene_Ana('ene_ana.md++.lib',ene_trj_path='md_peptide_1.tre')
    with open("md_peptide_1.tre.npy", "wb") as f:
        np.save(f, tr._tre)
    with open("md_peptide_1.tre.npy", "rb") as f:
        tre = np.load(f)
    
