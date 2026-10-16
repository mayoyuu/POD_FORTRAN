import csv
import inspect
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
from test_run_srp_modes_ads import load_wrapper
from test_srp_modes_probability_samples import load_generator

ROOT = Path(__file__).resolve().parents[1]
OPM = ROOT/'OPM/DRO-2/DRO-2_init.opm.json'
CFG = ROOT/'config/config_srp_modes_DRO-2.txt'


class ProbabilityRunner(unittest.TestCase):
    def require_feature(self):
        wrapper = load_wrapper()
        self.assertIn('probability_samples', inspect.signature(wrapper.run_modes).parameters,
                      'C2 probability opt-in is not implemented')
        return wrapper

    def test_requires_paired_inputs_and_exact_counts_before_running(self):
        wrapper = self.require_feature()
        with tempfile.TemporaryDirectory() as folder:
            csv, manifest = Path(folder)/'samples.csv', Path(folder)/'samples.json'
            load_generator().write_sample_bundle(csv, manifest, count=32)
            for kwargs in ({'probability_samples':csv}, {'probability_manifest':manifest},
                           {'probability_samples':csv, 'probability_manifest':manifest, 'probability_count':99},
                           {'probability_samples':csv, 'probability_manifest':manifest, 'probability_count':32,
                            'probability_visualization_count':0}):
                with self.subTest(kwargs=kwargs), self.assertRaises(wrapper.RunError):
                    wrapper.run_modes(OPM, CFG, Path(folder)/'rejected', **kwargs)
                self.assertFalse((Path(folder)/'rejected').exists())
                self.assertFalse(list(Path(folder).glob('.rejected.staging-*')))

    def test_short_probability_run_preserves_default_and_statistics(self):
        wrapper = self.require_feature()
        with tempfile.TemporaryDirectory() as folder:
            csv, manifest = Path(folder)/'samples.csv', Path(folder)/'samples.json'
            load_generator().write_sample_bundle(csv, manifest, count=32)
            base = dict(hours=.02, save_hours=.01)
            normal = Path(folder)/'normal'
            wrapper.run_modes(OPM, CFG, normal, **base)
            saved_moments = None
            for v in (5, 11):
                output = Path(folder)/f'probability{v}'
                report = wrapper.run_modes(OPM, CFG, output, **base, probability_samples=csv,
                    probability_manifest=manifest, probability_count=32,
                    probability_visualization_count=v, probability_validation_count=8)
                self.assertTrue(report['formal_result'])
                self.assertEqual(report['shape_count'], 512)
                self.assertEqual(report['validation_count'], 337)
                self.assertEqual(report['probability_outputs']['statistics_sample_count'], 32)
                self.assertEqual(report['probability_outputs']['validation_count'], 8)
                self.assertEqual(len(report['snapshot_prefixes']), 3)
                self.assertFalse(report['probability_outputs']['probability_output_accepted'])
                self.assertEqual((output/'c2_history.csv').read_bytes(), (normal/'c2_history.csv').read_bytes())
                self.assertEqual((output/'c2_shape.csv').read_bytes(), (normal/'c2_shape.csv').read_bytes())
                self.assertFalse((output/'c2_probability_visualization.csv').exists())
                self.assertFalse((output/'c2_probability_keyframes.csv').exists())
                self.assertEqual(json.loads((output/'probability_samples.json').read_text())['sample_count'], 32)
                with h5py.File(output/'c2_probability_points.h5') as data:
                    self.assertEqual(data['visualization/state'].shape, (3, v, 6))
                    self.assertEqual(data['keyframes/state'].shape, (2, 32, 6))
                    current = data['moments/covariance'][:]
                    if saved_moments is not None: np.testing.assert_array_equal(current, saved_moments)
                    saved_moments = current
                self.assertEqual(report['probability_outputs']['snapshot_probe_checks'], 3*8)


if __name__ == '__main__': unittest.main()
