import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from test_srp_ads_reconstruct import load_module

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'../../SRP/boxwing_modes_DRO-2_init'


def run_native(*args):
    env=os.environ.copy()
    env['LD_LIBRARY_PATH']=':'.join(map(str,[ROOT/'lib',Path.home()/'local_dace/lib',Path.home()/'local_openblas/lib']))+':'+env.get('LD_LIBRARY_PATH','')
    return subprocess.run(['fpm','run','--target','run_srp_modes_ads','--',*map(str,args)],cwd=ROOT,env=env,text=True,
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)


class NativeModesRunner(unittest.TestCase):
    def test_help_and_required_cli(self):
        result=run_native('--help')
        self.assertEqual(result.returncode,0,result.stdout[-4000:])
        self.assertIn('8D',result.stdout)
        self.assertIn('--hours',result.stdout)
        self.assertIn('--rel-tol',result.stdout)
        self.assertIn('--abs-tol',result.stdout)
        result=run_native('--hours','NaN')
        self.assertNotEqual(result.returncode,0,result.stdout[-2000:])

    def test_rejects_nonpositive_and_nonfinite_integrator_tolerances(self):
        with tempfile.TemporaryDirectory() as folder:
            base=['-opm',ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                  '-cfg',ROOT/'config/config_srp_modes_DRO-2.txt',
                  '-o',Path(folder)/'invalid']
            for option in ('--rel-tol','--abs-tol'):
                for value in ('0','-1e-13','NaN','Infinity'):
                    with self.subTest(option=option,value=value):
                        result=run_native(*base,option,value)
                        self.assertNotEqual(result.returncode,0,result.stdout[-2000:])

    def test_explicit_integrator_tolerances_reach_raw_report(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix=Path(folder)/'tight'
            result=run_native('-opm',ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                '-cfg',ROOT/'config/config_srp_modes_DRO-2.txt','-o',prefix,
                '--hours','0.01','--save-hours','0.01','--rel-tol','1e-13','--abs-tol','1e-13')
            self.assertEqual(result.returncode,0,result.stdout[-5000:])
            report=json.loads(Path(f'{prefix}_report.json').read_text())
            self.assertEqual(report['rel_tol'],1e-13)
            self.assertEqual(report['abs_tol'],1e-13)

    def test_short_raw_run_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix=Path(folder)/"c2 'quoted space"
            args=['-opm',ROOT/'OPM/DRO-2/DRO-2_init.opm.json','-cfg',ROOT/'config/config_srp_modes_DRO-2.txt',
                  '-o',prefix,'--hours','0.01','--save-hours','0.01',
                  ]
            result=run_native(*args)
            self.assertEqual(result.returncode,0,result.stdout[-5000:])
            report=json.loads(Path(f'{prefix}_report.json').read_text())
            self.assertEqual(report['n_variables'],8)
            self.assertFalse(report['formal_result'])
            self.assertEqual(report['rel_tol'],1e-12)
            self.assertEqual(report['abs_tol'],1e-12)
            self.assertEqual(report['parameterization'],'SRP_MODES')
            expected_hash=hashlib.sha256((DATA/'boxwing_modes_meta.json').read_bytes()).hexdigest()
            self.assertEqual(report['mode_metadata_sha256'],expected_hash)
            for snapshot in report['snapshot_prefixes']:
                self.assertEqual(load_module().load_snapshot(snapshot).mode_hashes['mode_metadata_sha256'],expected_hash)
            self.assertLessEqual(report['ads_numerical_error']['max_position_km'],0.1)
            self.assertLessEqual(report['ads_numerical_error']['max_velocity_kms'],1e-6)
            self.assertLessEqual(report['zero_mode_error']['max_position_km'],1e-6)
            self.assertLessEqual(report['zero_mode_error']['max_velocity_kms'],1e-9)
            self.assertEqual(report['da_handles_final'],report['da_handles_baseline'])
            original=Path(f'{prefix}_report.json').read_bytes()
            second=run_native(*args)
            self.assertNotEqual(second.returncode,0)
            self.assertEqual(Path(f'{prefix}_report.json').read_bytes(),original)
            rejected=run_native(*args[:args.index('-o')],'-o',Path(folder)/'bad',
                                '--hours','0.01','--save-hours','0.01','--mode-meta-sha256','0'*64)
            self.assertNotEqual(rejected.returncode,0)
            self.assertIn('metadata hash differs',rejected.stdout)


if __name__=='__main__': unittest.main()
