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
