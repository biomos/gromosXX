import sys
import os
import numpy as np

_current = os.path.dirname(os.path.realpath(__file__))
_current = os.path.dirname(_current)
_current = os.path.dirname(_current)
sys.path.append(_current)

from ene_ana import EnegyTrajectory

# energy trajectory tests
class Test_tre():
    
    @classmethod
    def load_tre(cls):
        # simulated trajectory
        ene_trj_file = cls.test_conf.out_files.get('tre')
        assert os.path.isfile(ene_trj_file)
        assert os.path.isfile(cls.test_conf.ene_ana_lib)
        cls.tr = EnegyTrajectory(cls.test_conf.ene_ana_lib, ene_trj_path=ene_trj_file)
        # hardcoded expected values
        cls.test_conf.load_expected_values('tre')
        cls.tre = cls.test_conf.expected_values['tre']
    
    @staticmethod
    def test_load_tre():
        Test_tre.load_tre()

    # test functions
    def test_TIME_block(self):
        """
        Checks if the 0th, 1st and the last timestep in the reference file and the new md are consistent
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["TIMESTEP"]["TIME"][0,0,:2,0] ,tr.tre["TIMESTEP"]["TIME"][0,0,:2,0], decimal=8)
        np.testing.assert_almost_equal(tre["TIMESTEP"]["TIME"][1,0,:2,0] ,tr.tre["TIMESTEP"]["TIME"][1,0,:2,0], decimal=8)
        np.testing.assert_almost_equal(tre["TIMESTEP"]["TIME"][2,0,:2,0] ,tr.tre["TIMESTEP"]["TIME"][2,0,:2,0], decimal=8)
        np.testing.assert_almost_equal(tre["TIMESTEP"]["TIME"][-1,0,:2,0] ,tr.tre["TIMESTEP"]["TIME"][-1,0,:2,0], decimal=8)

    def test_ENER_block(self):
        """
        Checks for energy consistency in different timesteps between ref data and new md.
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["ENERGY03"]["ENER"][0,0,:5,0] ,tr.tre["ENERGY03"]["ENER"][0,0,:5,0], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["ENER"][1,0,:5,0] ,tr.tre["ENERGY03"]["ENER"][1,0,:5,0], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["ENER"][2,0,:5,0] ,tr.tre["ENERGY03"]["ENER"][2,0,:5,0], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["ENER"][-1,0,:5,0] ,tr.tre["ENERGY03"]["ENER"][-1,0,:5,0], decimal=8)

    def test_KINENER_block(self):
        """
        Checks for kinetic energy consistency in different timesteps between ref data and new md.
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["ENERGY03"]["KINENER"][0,0,:2,:3] ,tr.tre["ENERGY03"]["KINENER"][0,0,:2,:3], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["KINENER"][1,0,:2,:3] ,tr.tre["ENERGY03"]["KINENER"][1,0,:2,:3], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["KINENER"][2,0,:2,:3] ,tr.tre["ENERGY03"]["KINENER"][2,0,:2,:3], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["KINENER"][-1,0,:2,:3] ,tr.tre["ENERGY03"]["KINENER"][-1,0,:2,:3], decimal=8)

    def test_BONDED_block(self):
        """
        Checks for bonded terms consistency in different timesteps between ref data and new md.
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["ENERGY03"]["BONDED"][0,0,:4,:5] ,tr.tre["ENERGY03"]["BONDED"][0,0,:4,:5], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["BONDED"][1,0,:4,:5] ,tr.tre["ENERGY03"]["BONDED"][1,0,:4,:5], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["BONDED"][2,0,:4,:5] ,tr.tre["ENERGY03"]["BONDED"][2,0,:4,:5], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["BONDED"][-1,0,:4,:5] ,tr.tre["ENERGY03"]["BONDED"][-1,0,:4,:5], decimal=8)

    def test_NONBONDED_block(self):
        """
        Checks for nonbonded terms consistency in different timesteps between ref data and new md.
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["ENERGY03"]["NONBONDED"][0,0,:8,:4] ,tr.tre["ENERGY03"]["NONBONDED"][0,0,:8,:4], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["NONBONDED"][1,0,:8,:4] ,tr.tre["ENERGY03"]["NONBONDED"][1,0,:8,:4], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["NONBONDED"][2,0,:8,:4] ,tr.tre["ENERGY03"]["NONBONDED"][2,0,:8,:4], decimal=8)
        np.testing.assert_almost_equal(tre["ENERGY03"]["NONBONDED"][-1,0,:8,:4] ,tr.tre["ENERGY03"]["NONBONDED"][-1,0,:8,:4], decimal=8)

    def test_TEMP_block(self):
        """
        Checks for temperature consistency in different timesteps between ref data and new md.
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["TEMPERATURE"][0,0,:2,:4] ,tr.tre["VOLUMEPRESSURE03"]["TEMPERATURE"][0,0,:2,:4], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["TEMPERATURE"][1,0,:2,:4] ,tr.tre["VOLUMEPRESSURE03"]["TEMPERATURE"][1,0,:2,:4], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["TEMPERATURE"][2,0,:2,:4] ,tr.tre["VOLUMEPRESSURE03"]["TEMPERATURE"][2,0,:2,:4], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["TEMPERATURE"][-1,0,:2,:4] ,tr.tre["VOLUMEPRESSURE03"]["TEMPERATURE"][-1,0,:2,:4], decimal=8)

    def test_VOLUME_block(self):
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["VOLUME"][0,0,:5,0] ,tr.tre["VOLUMEPRESSURE03"]["VOLUME"][0,0,:5,0], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["VOLUME"][1,0,:5,0] ,tr.tre["VOLUMEPRESSURE03"]["VOLUME"][1,0,:5,0], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["VOLUME"][2,0,:5,0] ,tr.tre["VOLUMEPRESSURE03"]["VOLUME"][2,0,:5,0], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["VOLUME"][-1,0,:5,0] ,tr.tre["VOLUMEPRESSURE03"]["VOLUME"][-1,0,:5,0], decimal=8)

    def test_PRESSURE_block(self):
        """
        Checks for pressure consistency in different timesteps between ref data and new md.
        """
        tre, tr = self.tre, self.tr
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["PRESSURE"][0,0,:10,0] ,tr.tre["VOLUMEPRESSURE03"]["PRESSURE"][0,0,:10,0], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["PRESSURE"][1,0,:10,0] ,tr.tre["VOLUMEPRESSURE03"]["PRESSURE"][1,0,:10,0], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["PRESSURE"][2,0,:10,0] ,tr.tre["VOLUMEPRESSURE03"]["PRESSURE"][2,0,:10,0], decimal=8)
        np.testing.assert_almost_equal(tre["VOLUMEPRESSURE03"]["PRESSURE"][-1,0,:10,0] ,tr.tre["VOLUMEPRESSURE03"]["PRESSURE"][-1,0,:10,0], decimal=8)
