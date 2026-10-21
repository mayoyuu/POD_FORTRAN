import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


def load_wrapper():
    spec=importlib.util.spec_from_file_location('run_srp_modes_ads',ROOT/'tools/run_srp_modes_ads.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


class FormalModesRunner(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.wrapper=load_wrapper()

    def test_rejects_malformed_integrator_tolerances(self):
        with tempfile.TemporaryDirectory() as folder:
            for option in ('rel_tol','abs_tol'):
                for value in (0.0,-1e-13,float('nan'),float('inf')):
                    with self.subTest(option=option,value=value):
                        kwargs={option:value}
                        with self.assertRaisesRegex(self.wrapper.RunError,'invalid run limits'):
                            self.wrapper.run_modes(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                                ROOT/'config/config_srp_modes_DRO-2.txt',Path(folder)/'invalid',**kwargs)

    def test_run_settings_reject_tiny_tolerance_identity_mismatch(self):
        report={'rel_tol':1.5e-13,'abs_tol':1.5e-13}
        expected={'rel_tol':1e-13,'abs_tol':1e-13}
        with self.assertRaisesRegex(self.wrapper.RunError,'run setting rel_tol mismatch'):
            self.wrapper.bind_run_settings(report,expected)
        report['rel_tol']=1e-13
        with self.assertRaisesRegex(self.wrapper.RunError,'run setting abs_tol mismatch'):
            self.wrapper.bind_run_settings(report,expected)

    def test_explicit_integrator_tolerances_reach_formal_report(self):
        with tempfile.TemporaryDirectory() as folder:
            report=self.wrapper.run_modes(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                ROOT/'config/config_srp_modes_DRO-2.txt',Path(folder)/'tight',
                hours=0.01,save_hours=0.01,rel_tol=1e-13,abs_tol=1e-13)
            self.assertTrue(report['formal_result'])
            self.assertEqual(report['rel_tol'],1e-13)
            self.assertEqual(report['abs_tol'],1e-13)

    def test_preflight_binds_runtime_opm_and_parameters(self):
        opm=ROOT/'OPM/DRO-2/DRO-2_init.opm.json'
        cfg=ROOT/'config/config_srp_modes_DRO-2.txt'
        self.wrapper.preflight(opm,cfg)
        with tempfile.TemporaryDirectory() as folder:
            bad=Path(folder)/'bad.txt'
            bad.write_text(cfg.read_text().replace('srp_mass_kg = 1200.0','srp_mass_kg = 1201.0'))
            with self.assertRaisesRegex(self.wrapper.RunError,'physical'):
                self.wrapper.preflight(opm,bad)
            bad_opm=Path(folder)/'bad.opm.json'; bad_opm.write_bytes(opm.read_bytes()+b'\n')
            with self.assertRaisesRegex(self.wrapper.RunError,'OPM'):
                self.wrapper.preflight(bad_opm,cfg)

    def test_publication_reconstruction_and_tamper_rejection(self):
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder)/'formal'
            report=self.wrapper.run_modes(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                ROOT/'config/config_srp_modes_DRO-2.txt',output,hours=0.01,save_hours=0.01)
            self.assertTrue(report['formal_result'])
            self.assertTrue(output.is_dir())
            self.assertEqual(report['rel_tol'],1e-12)
            self.assertEqual(report['abs_tol'],1e-12)
            self.assertTrue(report['provenance']['inputs']['mode_metadata']['sha256'])
            self.assertIn('fpm.toml',report['provenance']['source_files_sha256'])
            self.assertEqual(len(report['snapshot_prefixes']),2)
            self.assertTrue(all(item['verification']['status']=='passed' for item in report['snapshot_audits']))
            original=(output/'c2_report.json').read_bytes()
            with self.assertRaisesRegex(self.wrapper.RunError,'exists'):
                self.wrapper.run_modes(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',ROOT/'config/config_srp_modes_DRO-2.txt',output)
            self.assertEqual((output/'c2_report.json').read_bytes(),original)
            path=output/'c2_probes.csv'
            with path.open(newline='') as stream:
                reader=csv.DictReader(stream); fields=reader.fieldnames; rows=list(reader)
            rows[-1]['ads_x']=str(float(rows[-1]['ads_x'])+1)
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(rows)
            raw=dict(report); raw['formal_result']=False
            with self.assertRaisesRegex(self.wrapper.RunError,'probe'):
                self.wrapper.audit_outputs(output,raw,report['mode_csv_sha256'],report['mode_metadata_sha256'])

    def test_failure_preserves_staging_without_publication(self):
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder)/'failed'
            with self.assertRaisesRegex(self.wrapper.RunError,'diagnostics preserved'):
                self.wrapper.run_modes(ROOT/'OPM/DRO-2/DRO-2_init.opm.json',
                    ROOT/'config/config_srp_modes_DRO-2.txt',output,hours=0.01,save_hours=0.01,executable='/usr/bin/false')
            self.assertFalse(output.exists())
            staging=list(Path(folder).glob('.failed.staging-*'))
            self.assertEqual(len(staging),1)
            self.assertTrue((staging[0]/'runner').is_file())
            self.assertTrue((staging[0]/'run.log').is_file())
            self.assertFalse((Path(folder)/'.failed.publish.lock').exists())


if __name__=='__main__': unittest.main()
