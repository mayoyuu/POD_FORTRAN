import csv
import gzip
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / 'tools' / 'generate_srp_mc_samples.py'


def load_module():
    spec = importlib.util.spec_from_file_location('generate_srp_mc_samples', MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FixedSampleGeneratorTests(unittest.TestCase):
    def test_deterministic_truncated_samples_and_mapping(self):
        mod = load_module()
        first, rejected_first = mod.generate_samples(50_000, 20261005, 3.0)
        second, rejected_second = mod.generate_samples(50_000, 20261005, 3.0)
        self.assertEqual(first.shape, (50_000, 7))
        self.assertTrue(np.array_equal(first, second))
        self.assertEqual(rejected_first, rejected_second)
        self.assertGreater(rejected_first, 0)
        self.assertTrue(np.isfinite(first).all())
        self.assertLessEqual(float(np.max(np.abs(first))), 3.0)
        self.assertFalse(np.array_equal(first, mod.generate_samples(50_000, 20261006, 3.0)[0]))

    def test_bundle_has_stable_ids_coordinates_eta_and_hash(self):
        mod = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp) / 'samples.csv'
            manifest_path = Path(tmp) / 'samples.json'
            npz_path = Path(tmp) / 'samples.npz'
            gzip_path = Path(tmp) / 'samples.csv.gz'
            manifest = mod.write_sample_bundle(
                csv_path, manifest_path, count=256, seed=20261005,
                sigma_limit=3.0, eta_sigma=0.3, npz_path=npz_path,
                gzip_path=gzip_path,
            )
            with csv_path.open(newline='', encoding='utf-8') as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 256)
            self.assertEqual([int(row['sample_id']) for row in rows], list(range(1, 257)))
            for row in rows:
                for i in range(1, 8):
                    self.assertAlmostEqual(float(row[f'u{i}']), float(row[f'z{i}']) / 3.0, places=15)
                self.assertAlmostEqual(float(row['eta_treatment']), 0.3 * float(row['z7']), places=15)
            digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
            loaded = json.loads(manifest_path.read_text(encoding='utf-8'))
            self.assertEqual(manifest, loaded)
            self.assertEqual(loaded['csv_sha256'], digest)
            self.assertEqual(loaded['sample_count'], 256)
            self.assertEqual(loaded['prefix_counts'], [10_000, 20_000, 50_000])
            self.assertEqual(loaded['generator'], 'numpy.random.PCG64DXSM')
            self.assertEqual(loaded['coordinate_rule'], 'u=z/3')
            self.assertEqual(loaded['control_eta'], 0.0)
            self.assertEqual(loaded['npz_sha256'],
                             hashlib.sha256(npz_path.read_bytes()).hexdigest())
            self.assertEqual(loaded['gzip_csv_sha256'],
                             hashlib.sha256(gzip_path.read_bytes()).hexdigest())
            with np.load(npz_path) as bundle:
                self.assertEqual(bundle['z'].shape,(256,7))
                np.testing.assert_allclose(bundle['u'],bundle['z']/3.0)
                np.testing.assert_allclose(bundle['eta_treatment'],0.3*bundle['z'][:,6])
            with gzip.open(gzip_path,'rt',newline='',encoding='utf-8') as handle:
                compressed_rows=list(csv.DictReader(handle))
            self.assertEqual(compressed_rows,rows)

    def test_rejects_invalid_arguments(self):
        mod = load_module()
        for args in ((0, 1, 3.0), (10, 1, 0.0), (10, 1, float('nan'))):
            with self.assertRaises(ValueError):
                mod.generate_samples(*args)


if __name__ == '__main__':
    unittest.main()
