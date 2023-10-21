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
import subprocess

class Test_MD():
    def test_md_run(self):
        """
        Executes a short simulation of a test system checks for output files.
        """
        run_file = self.test_conf.write_run_file()
        if os.path.abspath(self.test_conf.sim_dir) != os.getcwd():
            os.chdir(self.test_conf.sim_dir)

        exe = subprocess.run(['bash', run_file],
                        check=True,
                        capture_output=True,
                        text=True
                    )
        exe.check_returncode()
        assert os.path.isfile(self.test_conf.omd)
        assert os.path.isfile(self.test_conf.out_files['fin'])
        out_files = self.test_conf.out_files
        if 'trc' in out_files: assert os.path.isfile(out_files['trc'])
        if 'tre' in out_files: assert os.path.isfile(out_files['tre'])
        if 'trg' in out_files: assert os.path.isfile(out_files['trg'])
        if 'trv' in out_files: assert os.path.isfile(out_files['trv'])
        if 'trf' in out_files: assert os.path.isfile(out_files['trf'])
        # this is all same as checking all files like that, but it makes a nicer output
        for temp_file in self.test_conf.out_files.values():
            assert os.path.isfile(temp_file)
