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

######## default code ########
_current = os.path.dirname(os.path.realpath(__file__))
_parent = os.path.dirname(_current)
sys.path.append(_parent)
from tests._config_test_params import ConfigTest

 # for local testing
_config_file = os.getenv('ConfigFile', 'default_param.yaml')
if os.path.isfile(_config_file):
    _test_conf = ConfigTest(_config_file)
else:
    # default settings based on the path to the Test_Repo (gromos_test_files)
    _test_repo_path = os.getenv('TEST_REPO', 'gromos_test_files')
    _test_dir = 'tutorial_pept'
    _config_file = 'default_param.yaml'
    _config_file = os.path.join(_test_repo_path, _test_dir, _config_file)
    _test_conf = ConfigTest(_config_file)

# make sure to clear_up the sim_dir before running the tests
if os.path.isdir(_test_conf.sim_dir):_test_conf.clean_up()
######## default code ########

######## specific code for the test ########
# test MD
from tests._md_tests import Test_MD
Test_MD.test_conf = _test_conf

# test ene
from tests._tre_tests import Test_tre
Test_tre.test_conf = _test_conf
######## specific code for the test ########

######## clean-up code, relevant for the local testing ########
flag_clean_up=False # set True for local testing
# or set it up using export FLAG_TEST_CLEAN_UP=1
temp_flag_clean_up = os.getenv('FLAG_TEST_CLEAN_UP')
if temp_flag_clean_up:
    if str(temp_flag_clean_up).lower() not in ('0', 'false'):
        flag_clean_up = True
# or set it up using flag_clean_up: True in the yaml file
if _test_conf.config.get('flag_clean_up'):
    flag_clean_up = True

if flag_clean_up:
    # this has to be declared in this way such that it can be run automatically by pytest
    def test_clean_up():
        _test_conf.clean_up()
