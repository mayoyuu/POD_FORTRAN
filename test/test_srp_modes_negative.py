"""Run native error-stop contracts in isolated processes; never accept a missing target."""
import os
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ModeForceFailures(unittest.TestCase):
    def test_explicit_failure_contracts(self):
        binaries = list((ROOT / "build").glob("gfortran_*/test/test_srp_modes_force_consistency"))
        self.assertEqual(len(binaries), 1, "build native test first with fpm test test_srp_modes_force_consistency")
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = ":".join(map(str, [ROOT / "lib", Path.home() / "local_dace/lib",
                                                       Path.home() / "local_openblas/lib"])) + ":" + env.get("LD_LIBRARY_PATH", "")
        cases = {
            "invalid-real": "coefficients must be in [0,1]",
            "nan-real": "coefficients must be finite",
            "invalid-index": "invalid SRP mode DA indices",
            "out-of-range-index": "invalid SRP mode DA indices",
            "invalid-span": "DA domain must stay in [0,1]",
            "nan-da": "DA coefficients must be finite",
            "unset-real": "SRP mode coefficients are unset",
            "unset-da": "SRP mode DA coefficients are unset",
            "no-table": "SRP mode table not loaded",
            "mixed-real": "cannot mix eta/ballistic overrides",
            "mixed-da": "cannot mix eta/ballistic overrides",
            "override-real": "cannot mix eta/ballistic overrides",
            "override-da": "cannot mix eta/ballistic overrides",
            "extrapolation": "outside table coverage",
            "init-clears-real": "SRP mode coefficients are unset",
            "init-clears-da": "SRP mode DA coefficients are unset",
        }
        for arg, expected in cases.items():
            with self.subTest(arg=arg):
                result = subprocess.run([str(binaries[0]), arg], cwd=ROOT, env=env, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
                self.assertNotEqual(result.returncode, 0, result.stdout[-4000:])
                self.assertIn(expected, result.stdout, result.stdout[-4000:])


if __name__ == "__main__":
    unittest.main()
