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

######## specific code for the test ########
# set test directory in the repository (gromos_test_files)
_repo_test_dir='AEDS_test'
_config_file='GA_param.yaml'
######## specific code for the test ########

######## default code ########
_current = os.path.dirname(os.path.realpath(__file__))
_parent = os.path.dirname(_current)
sys.path.append(_parent)
from tests._config_test_params import ConfigTest

_config_file = ConfigTest.find_config_file(_repo_test_dir, conf_file=_config_file)
_test_conf = ConfigTest(_config_file)

# make sure to clear_up the sim_dir before running the tests
if os.path.isdir(_test_conf.sim_dir):_test_conf.clean_up()
######## default code ########

######## specific code for the test ########
# test MD
from tests._md_tests import Base_MD_Test
class Test_AEDS_GA_MD(Base_MD_Test):
    test_conf = _test_conf

# test ene
from tests._tre_tests import AEDS_Tests_tre
class Test_AEDS_GA_tre(AEDS_Tests_tre):
    test_conf = _test_conf
######## specific code for the test ########


######## clean-up code, relevant for the local testing ########
if _test_conf.flag_cleanup():
    def test_clean_up():
        _test_conf.clean_up()
