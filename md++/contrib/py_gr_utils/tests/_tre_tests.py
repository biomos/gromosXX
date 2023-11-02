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
class Base_Test_tre:
    "base class for all test related to energy trajectory"

    @classmethod
    def load_tre(cls):
        # simulated trajectory
        ene_trj_file = cls.test_conf.out_files.get('tre')
        assert os.path.isfile(ene_trj_file)
        assert os.path.isfile(cls.test_conf.ene_ana_lib)
        cls.tre = EnergyTrajectory(cls.test_conf.ene_ana_lib, trj_files=(ene_trj_file,), num_type=np.single)
        # hardcoded expected values
        cls.test_conf.load_expected_values('tre')
        cls.HC_tre = cls.test_conf.expected_values['tre']

    def test_load_tre(self):
        self.load_tre()


class Basic_Tests_tre(Base_Test_tre):
    "class for basic tests related to energy trajectory"

    def _test_almost_eq_var(self, *var):
        tre_data = self.tre.get_values(*var)
        HC_tre_data = self.HC_tre.get_values(*var)
        np.testing.assert_almost_equal(tre_data, HC_tre_data)

    def _test_almost_eq_subbl(self, subblock, idx=()):
        tre_data = self.tre.extract_values(subblock, idx=idx)
        HC_tre_data = self.HC_tre.extract_values(subblock, idx=idx)
        np.testing.assert_almost_equal(tre_data, HC_tre_data)

    # test functions
    def test_TIME_block(self):
        """
        Checks if the TIMESTEP in the reference file and the new md are consistent
        """
        self.tre.variables['step'] = ('TIME', (0,))
        self.HC_tre.variables['step'] = ('TIME', (0,))
        self._test_almost_eq_var('step', 'time')

    def test_ENER_block(self):
        """
        Checks for energy consistency between ref data and new md.
        """
        self._test_almost_eq_var('totene', 'totkin', 'totpot', 'totcov')
        self._test_almost_eq_var('totbond', 'totangle', 'totimproper', 'totdihedral')
        self._test_almost_eq_var('totnonbonded', 'totlj','totcrf')

    def test_KINENER_block(self):
        """
        Checks for kinetic energy consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('KINENER')

    def test_BONDED_block(self):
        """
        Checks for bonded terms consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('BONDED')

    def test_NONBONDED_block(self):
        """
        Checks for nonbonded terms consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('NONBONDED')

    def test_TEMP_block(self):
        """
        Checks for temperature consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('TEMPERATURE')

    def test_VOLUME_block(self):
        self._test_almost_eq_subbl('VOLUME')

    def test_PRESSURE_block(self):
        """
        Checks for pressure consistency in different timesteps between ref data and new md.
        """
        self._test_almost_eq_subbl('PRESSURE', idx=slice(10)) # first 10 numbers only


class AEDS_Tests_tre(Basic_Tests_tre):
    "class for EDS tests related to energy trajectory"

    def test_AEDS_energy(self):
        """
        Checks for EDS energy consistency between ref data and new md.
        """
        self._test_almost_eq_var('eds_vmix', 'eds_vr')
        self._test_almost_eq_var('eds_emax', 'eds_emin', 'eds_globmin', 'eds_globminfluc')

    def test_AEDS_state_energy(self):
        """
        Checks for EDS energy consistency between ref data and new md.
        """
        self._test_almost_eq_var('e1', 'e2', 'e3', 'e4')
        self._test_almost_eq_var('e1r', 'e2r', 'e3r', 'e4r')
