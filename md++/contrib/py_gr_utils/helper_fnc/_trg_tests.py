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

import sys
import os
import numpy as np

_current = os.path.dirname(os.path.realpath(__file__))
_current = os.path.dirname(_current)
_current = os.path.dirname(_current)
sys.path.append(_current)

from ene_ana import EnergyTrajectory

# energy trajectory tests
class Base_Test_trg:
    "base class for all test related to energy trajectory"

    @classmethod
    def load_trg(cls):
        # simulated trajectory
        ene_trj_file = cls.test_conf.out_files.get('trg')
        assert os.path.isfile(ene_trj_file)
        assert os.path.isfile(cls.test_conf.ene_ana_lib)
        cls.trg = EnergyTrajectory(cls.test_conf.ene_ana_lib, trj_files=(ene_trj_file,), num_type=np.single)
        # hardcoded expected values
        cls.test_conf.load_expected_values('trg')
        cls.HC_trg = cls.test_conf.expected_values['trg']

    def test_load_trg(self):
        self.load_trg()


class Basic_Tests_trg(Base_Test_trg):
    "class for basic tests related to free energy trajectory"
    precision = 8

    @classmethod
    def __call_assert_almost_equal(cls, trg_data, HC_trg_data):
        np.testing.assert_almost_equal(trg_data, HC_trg_data, decimal=cls.precision) # make this into a variable

    def _test_almost_eq_var(self, *var):
        trg_data = self.trg.get_values(*var)
        HC_trg_data = self.HC_trg.get_values(*var)
        self.__call_assert_almost_equal(trg_data, HC_trg_data)

    def _test_almost_eq_subbl(self, subblock, idx=()):
        trg_data = self.trg.extract_values(subblock, idx=idx)
        HC_trg_data = self.HC_trg.extract_values(subblock, idx=idx)
        self.__call_assert_almost_equal(trg_data, HC_trg_data)

    # test functions
    def test_TIME_block(self):
        """
        Checks if the TIMESTEP in the reference file and the new md are consistent
        """
        self.trg.variables['step'] = ('TIME', (0,))
        self.HC_trg.variables['step'] = ('TIME', (0,))
        self._test_almost_eq_var('step', 'time')

    def test_FREEENER_block(self):
        """
        Checks for energy consistency between ref data and new md.
        """
        self._test_almost_eq_var('totfren', 'totfrpot')

    def test_FREEKINENER_block(self):
        """
        Checks for kinetic energy consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('FREEKINENER')

    def test_FREEBONDED_block(self):
        """
        Checks for bonded terms consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('FREEBONDED')

    def test_FREENONBONDED_block(self):
        """
        Checks for nonbonded terms consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('FREENONBONDED')

class extTI_Tests_trg(Basic_Tests_trg):
    "class for extTI tests related to free energy trajectory"

    def test_extTI_freeenergy(self):
        """
        Checks for extendedTI energy consistency between ref data and new md.
        """
        self._test_almost_eq_subbl('FREEPRECALCLAM')
