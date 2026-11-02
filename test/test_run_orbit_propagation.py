"""Native integration tests; set POD_ORBIT_EXE to an already-built executable."""
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
KEYS = ['X', 'Y', 'Z', 'X_DOT', 'Y_DOT', 'Z_DOT']
COV_KEYS = ['C' + a + '_' + b for i, a in enumerate(KEYS) for b in KEYS[:i+1]]


@unittest.skipUnless(os.environ.get('POD_ORBIT_EXE'), 'set POD_ORBIT_EXE to run native tests')
class OrbitPropagationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.exe = Path(os.environ['POD_ORBIT_EXE']).resolve()
        assert cls.exe.is_file(), 'run_orbit_propagation executable has not been built'

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.initial = json.loads((ROOT/'OPM/DRO-1/DRO-1_init.opm.json').read_text())
        self.initial.update(EPOCH='2027-01-01T00:01:00.000000', MASS=4321.0,
                            REF_FRAME='J2000', CENTER_NAME='EARTH', TIME_SYSTEM='UTC')
        # Diagonal covariance gives independently reproducible UT points.
        for key in COV_KEYS:
            self.initial[key] = 0.0
        self.variances = [1.0, 2.0, 3.0, 1e-10, 2e-10, 3e-10]
        for key, variance in zip(KEYS, self.variances):
            self.initial['C'+key+'_'+key] = variance

    def run_case(self, data, args, name='case', success=True):
        opm = self.directory/(name+'.json')
        opm.write_text(json.dumps(data))  # compact JSON must also work
        prefix = self.directory/name
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        result = subprocess.run([str(self.exe), '--opm', str(opm), '--config',
                                 'config/config.txt', '--output', str(prefix), *args],
                                cwd=ROOT, env=env, text=True, capture_output=True, timeout=90)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            with Path(str(prefix)+'_ephemeris.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            return rows, json.loads(Path(str(prefix)+'_final.opm.json').read_text())
        self.assertNotEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertFalse(Path(str(prefix)+'_final.opm.json').exists())
        return result.stdout+result.stderr

    def without_cov(self, data=None):
        return {k: v for k, v in (data or self.initial).items() if k not in COV_KEYS}

    def test_hourly_grid_and_explicit_endpoint(self):
        rows, final = self.run_case(self.initial, ['--hours', '2.5'])
        self.assertEqual([float(r['elapsed_seconds']) for r in rows], [0, 3600, 7200, 9000])
        other, end = self.run_case(self.initial, ['--end-utc', final['EPOCH']], 'end')
        self.assertEqual(len(other), 4)
        self.assertEqual(final['PROPAGATION_METHOD'], 'UT')
        self.assertEqual(final['MASS'], 4321)
        self.assertNotIn('GMM_COMPONENTS', final)
        self.assertNotIn('RMS', final)
        for key, column in zip(KEYS, ['x_km','y_km','z_km','vx_kms','vy_kms','vz_kms']):
            self.assertAlmostEqual(final[key], float(rows[-1][column]), places=10)
            self.assertAlmostEqual(final[key], end[key], delta=2e-6)
        self.assertTrue(all(math.isfinite(final[k]) for k in COV_KEYS))

    def test_no_covariance_and_zero_duration_roundtrip(self):
        rows, final = self.run_case(self.without_cov(), ['--hours', '0'])
        self.assertEqual(len(rows), 1)
        self.assertEqual(final['PROPAGATION_METHOD'], 'DETERMINISTIC')
        self.assertFalse(any(k in final for k in COV_KEYS))
        for key in KEYS:
            self.assertEqual(final[key], self.initial[key])
        _, roundtrip = self.run_case(final, ['--hours', '0.1'], 'roundtrip')
        self.assertNotEqual(roundtrip['X'], final['X'])

    def test_ut_matches_independent_sigma_trajectories(self):
        _, final = self.run_case(self.initial, ['--hours', '0.1'], 'ut')
        sigma_states = []
        weights = [0.5/6.5]+[1/13]*12
        for index in range(13):
            point = self.without_cov()
            if index:
                axis = (index-1) % 6
                sign = 1 if index <= 6 else -1
                point[KEYS[axis]] += sign*math.sqrt(6.5*self.variances[axis])
            _, propagated = self.run_case(point, ['--hours', '0.1'], 'point'+str(index))
            sigma_states.append([propagated[key] for key in KEYS])
        mean = [sum(w*x[j] for w, x in zip(weights, sigma_states)) for j in range(6)]
        for i, key in enumerate(KEYS):
            self.assertAlmostEqual(final[key], mean[i], delta=1e-9)
            for j in range(i+1):
                cov = sum(w*(x[i]-mean[i])*(x[j]-mean[j]) for w,x in zip(weights, sigma_states))
                self.assertAlmostEqual(final['C'+key+'_'+KEYS[j]], cov, delta=max(1e-18,abs(cov)*1e-7))

    def test_rejects_partial_invalid_covariance_and_time_conflicts(self):
        partial = dict(self.initial)
        del partial['CY_X']
        self.assertIn('covariance', self.run_case(partial, ['--hours','1'], 'partial', False))
        bad = dict(self.initial, CX_X=-1)
        self.assertIn('covariance', self.run_case(bad, ['--hours','1'], 'bad', False))
        self.run_case(self.initial, ['--hours','1','--end-utc','2027-01-02'], 'conflict', False)
        self.run_case(self.initial, ['--hours','-1'], 'negative', False)
        self.run_case(self.initial, ['--end-utc','2026-01-01'], 'past', False)
        self.run_case(dict(self.initial, CENTER_NAME='MOON'), ['--hours','1'], 'frame', False)


if __name__ == '__main__':
    unittest.main()
