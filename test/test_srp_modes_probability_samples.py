import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load_generator():
    path = ROOT / 'tools/generate_srp_modes_probability_samples.py'
    assert path.is_file(), '8D probability sample generator is not implemented'
    spec = importlib.util.spec_from_file_location('probability_samples', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ProbabilitySamples(unittest.TestCase):
    def test_deterministic_prefixes_and_prior_moments(self):
        generator = load_generator()
        small, _ = generator.generate_samples(73)
        large, metadata = generator.generate_samples(50000)
        np.testing.assert_array_equal(small, large[:73])
        self.assertEqual(large.shape, (50000, 19))
        z, q = large[:, 1:7], large[:, 7:15]
        self.assertTrue(np.all(np.abs(z) < 3))
        np.testing.assert_allclose(q[:, :6], z / 3, rtol=0, atol=0)
        np.testing.assert_allclose(q[:, 6:], 2 * large[:, 15:17] - 1)
        np.testing.assert_allclose(large[:, 17], large[:, 15])
        np.testing.assert_allclose(large[:, 18], (1-large[:, 15])*large[:, 16])
        kappa = 1 - 6*math.exp(-4.5)/math.sqrt(2*math.pi)/math.erf(3/math.sqrt(2))
        np.testing.assert_allclose(z.mean(0), 0, atol=.02)
        np.testing.assert_allclose(z.var(0), kappa, atol=.025)
        np.testing.assert_allclose(large[:, 15:17].mean(0), .5, atol=.007)
        np.testing.assert_allclose(large[:, 15:17].var(0), 1/12, atol=.002)
        self.assertGreater(metadata['rejected_state_rows'], 0)
        self.assertEqual(metadata['state_spawn_key'], [0])
        self.assertEqual(metadata['modes_spawn_key'], [1])

    def test_fixed_visualization_ids(self):
        generator = load_generator()
        ids = generator.choose_ids(50000, 2000)
        self.assertEqual(len(set(ids)), 2000)
        self.assertEqual(ids[[0, -1]].tolist(), [1, 50000])
        self.assertEqual(generator.choose_ids(3, 2000).tolist(), [1, 2, 3])
        self.assertEqual(generator.choose_ids(10, 1).tolist(), [1])
        for n, v in ((0, 2), (10, 0), (2, -1)):
            with self.assertRaises(ValueError): generator.choose_ids(n, v)

    def test_bundle_integrity_and_no_overwrite(self):
        generator = load_generator()
        with tempfile.TemporaryDirectory() as folder:
            csv, manifest = Path(folder)/'samples.csv', Path(folder)/'samples.json'
            generator.write_sample_bundle(csv, manifest, count=37, visualization_count=9)
            samples, metadata = generator.load_sample_bundle(csv, manifest, expected_count=37)
            self.assertEqual(samples.shape, (37, 19))
            self.assertEqual(metadata['visualization_ids'], generator.choose_ids(37, 9).tolist())
            original = csv.read_bytes()
            with self.assertRaises(FileExistsError): generator.write_sample_bundle(csv, manifest, count=37)
            self.assertEqual(csv.read_bytes(), original)
            csv.write_bytes(original+b'\n')
            with self.assertRaisesRegex(ValueError, 'hash'): generator.load_sample_bundle(csv, manifest)

    def test_rejects_bad_schema_ids_nonfinite_and_mapping(self):
        generator = load_generator()
        with tempfile.TemporaryDirectory() as folder:
            csv, manifest = Path(folder)/'samples.csv', Path(folder)/'samples.json'
            generator.write_sample_bundle(csv, manifest, count=8)
            original, meta = csv.read_text(), json.loads(manifest.read_text())
            cases = [original.replace('samples_v1', 'samples_v0'),
                     original.replace('\n2,', '\n1,'),
                     original.replace('\n1,', '\n1,nan,', 1)]
            fields = original.splitlines()
            row = fields[2].split(','); row[17] = '0.999'
            cases.append('\n'.join(fields[:2]+[','.join(row)]+fields[3:])+'\n')
            for changed in cases:
                csv.write_text(changed)
                meta['csv_sha256'] = hashlib.sha256(csv.read_bytes()).hexdigest()
                manifest.write_text(json.dumps(meta))
                with self.subTest(changed=changed[:60]), self.assertRaises(ValueError):
                    generator.load_sample_bundle(csv, manifest)


if __name__ == '__main__': unittest.main()
