"""Real-only endpoint propagation and audited publication, independent of DA."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCENARIOS={'Deta','proxy_earth','proxy_moon','boxwing_earth','boxwing_moon'}


class ModesTruth(unittest.TestCase):
    def wrapper(self):
        path=ROOT/'tools/run_srp_modes_truth.py'
        if not path.is_file(): self.fail('Real truth wrapper is not implemented')
        spec=importlib.util.spec_from_file_location('run_srp_modes_truth',path)
        module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        return module

    def native(self,*args):
        return subprocess.run(['fpm','run','--target','run_srp_modes_truth','--',*map(str,args)],
                              cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)

    def test_native_help_and_real_short_profiles(self):
        help_result=self.native('--help')
        self.assertEqual(help_result.returncode,0,help_result.stdout[-3000:])
        self.assertIn('Real-only',help_result.stdout)
        with tempfile.TemporaryDirectory() as folder:
            prefix=Path(folder)/'truth'
            result=self.native('-opm',ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                '-cfg',ROOT/'config/config_srp_modes_DRO-2.txt','-o',prefix,
                '--hours','0.01','--save-hours','0.01')
            self.assertEqual(result.returncode,0,result.stdout[-3000:])
            with Path(f'{prefix}_states.csv').open() as stream: rows=list(csv.DictReader(stream))
            self.assertEqual(len(rows),10)
            self.assertEqual({row['scenario'] for row in rows},SCENARIOS)
            report=json.loads(Path(f'{prefix}_report.json').read_text())
            self.assertEqual(report['da_handles_final'],0)
            self.assertIn('geometry_fallback_count',report)
            self.assertEqual(report['geometry_fallback_count'],0)
            self.assertFalse(report['formal_result'])
            final={row['scenario']:tuple(float(row[key]) for key in ('x_km','y_km','z_km'))
                   for row in rows if float(row['hour'])>0}
            self.assertNotEqual(final['proxy_earth'],final['proxy_moon'])
            self.assertNotEqual(final['boxwing_earth'],final['boxwing_moon'])

    def test_formal_truth_identity_no_overwrite_and_tamper(self):
        module=self.wrapper()
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder)/'formal'
            report=module.run_truth(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                ROOT/'config/config_srp_modes_DRO-2.txt',output,hours=0.01,save_hours=0.01)
            self.assertTrue(report['formal_result'])
            self.assertEqual(report['state_source'],'independent_real_nominal')
            self.assertEqual(report['rel_tol'],1e-13)
            self.assertEqual(report['row_count'],10)
            self.assertEqual(report['da_handles_baseline'],report['da_handles_final'])
            original=(output/'truth_report.json').read_bytes()
            with self.assertRaisesRegex(module.RunError,'exists'):
                module.run_truth(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                    ROOT/'config/config_srp_modes_DRO-2.txt',output,hours=0.01,save_hours=0.01)
            self.assertEqual((output/'truth_report.json').read_bytes(),original)
            path=output/'truth_states.csv'
            with path.open() as stream:
                reader=csv.DictReader(stream); fields=reader.fieldnames; rows=list(reader)
            rows[-1]['x_km']='NaN'
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
            with self.assertRaisesRegex(module.RunError,'nonfinite'):
                module.validate_truth_output(output,report,0.01,0.01,
                    report['epoch0_tdb_s'],[float(rows[0][k]) for k in module.STATE_KEYS])

    def test_physical_geometry_and_bad_tolerances_rejected(self):
        module=self.wrapper()
        with tempfile.TemporaryDirectory() as folder:
            cfg=Path(folder)/'bad.txt'
            cfg.write_text((ROOT/'config/config_srp_modes_DRO-2.txt').read_text().replace(
                'srp_box_dimensions_m = 2.0 3.0 4.0','srp_box_dimensions_m = 9.0 3.0 4.0'))
            with self.assertRaisesRegex(module.RunError,'physical'):
                module.preflight(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',cfg)
            for value in (0,float('nan'),float('inf')):
                with self.subTest(value=value),self.assertRaisesRegex(module.RunError,'invalid'):
                    module.run_truth(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                        ROOT/'config/config_srp_modes_DRO-2.txt',Path(folder)/'bad',rel_tol=value)

    def test_incomplete_integration_preserves_diagnostics(self):
        module=self.wrapper()
        with tempfile.TemporaryDirectory() as folder:
            cfg=Path(folder)/'limited.txt'
            cfg.write_text((ROOT/'config/config_srp_modes_DRO-2.txt').read_text().replace(
                'max_propagation_steps = 10000','max_propagation_steps = 2'))
            output=Path(folder)/'incomplete'
            with self.assertRaisesRegex(module.RunError,'diagnostics preserved'):
                module.run_truth(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',cfg,output,hours=0.01,save_hours=0.01)
            self.assertFalse(output.exists())
            staging=list(Path(folder).glob('.incomplete.staging-*'))
            self.assertEqual(len(staging),1)
            self.assertIn('did not reach', (staging[0]/'run.log').read_text())


if __name__=='__main__': unittest.main()
