import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class BatchShellTests(unittest.TestCase):
    def test_parallel_wrapper_can_preview_without_parallel_installed(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'not-created'
            proc = subprocess.run(['bash', str(ROOT/'.sh/halo_15d_parallel.sh'), str(out)],
                                  env=dict(os.environ, DRY_RUN='1'), text=True, capture_output=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(sum('--cases ' in s for s in proc.stdout.splitlines()), 12)
            self.assertFalse(out.exists())

    def test_dry_run_lists_twelve_360_hour_single_jobs_without_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'not-created'
            env = dict(os.environ, DRY_RUN='1', BACKEND='parallel', JOBS='3')
            proc = subprocess.run(['bash', str(ROOT/'.sh/halo_15d.sh'), str(out)],
                                  env=env, text=True, capture_output=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            commands = [s for s in proc.stdout.splitlines() if '--cases ' in s]
            self.assertEqual(len(commands), 12)
            self.assertTrue(all('--hours 360' in s and '--jobs 1' in s for s in commands))
            self.assertTrue(all('--count 50000' in s and '--save-hours 1' in s for s in commands))
            self.assertFalse(out.exists())

    def test_existing_results_are_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(['bash', str(ROOT/'.sh/halo_15d.sh'), tmp],
                                  env=dict(os.environ, DRY_RUN='1'), text=True, capture_output=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn('exists', proc.stderr)


if __name__ == '__main__':
    unittest.main()
