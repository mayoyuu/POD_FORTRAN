import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


sys.dont_write_bytecode = True

SCRIPT = Path(__file__).parents[1] / "tools" / "run_srp_ads_batch.py"
SPEC = importlib.util.spec_from_file_location("run_srp_ads_batch", SCRIPT)
batch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch)


def write_opm(path: Path, pos_rss_km: float, vel_rss_mps: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pos_var = pos_rss_km**2 / 3.0
    vel_var = (vel_rss_mps / 1000.0) ** 2 / 3.0
    data = {
        "CX_X": pos_var,
        "CY_Y": pos_var,
        "CZ_Z": pos_var,
        "CX_DOT_X_DOT": vel_var,
        "CY_DOT_Y_DOT": vel_var,
        "CZ_DOT_Z_DOT": vel_var,
    }
    path.write_text(json.dumps(data), encoding="utf-8")


class BatchDiscoveryTest(unittest.TestCase):
    def test_discovers_mirrors_and_filters_expected_initial_error(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            opm = root / "OPM"
            output = root / "SRP" / "SRP_UCERTAINTY_PROP_WITHINIT"
            good = opm / "DRO-1" / "DRO-1_init.opm.json"
            bad = opm / "nested" / "DRO-1_corrected_init.opm.json"
            ignored = opm / "DRO-2" / "DRO-2_final.opm.json"
            write_opm(good, 10.0, 0.03)
            write_opm(bad, 1.0, 0.01)
            write_opm(ignored, 10.0, 0.03)

            found = batch.find_initial_opms(opm)
            self.assertEqual(found, [good, bad])
            pos_rss, vel_rss = batch.covariance_rss(good)
            self.assertTrue(math.isclose(pos_rss, 10.0, rel_tol=1e-12))
            self.assertTrue(math.isclose(vel_rss, 0.03, rel_tol=1e-12))
            self.assertTrue(batch.matches_expected_error(good, 10.0, 0.03))
            self.assertFalse(batch.matches_expected_error(bad, 10.0, 0.03))
            self.assertEqual(
                batch.output_prefix_for(good, opm, output),
                output / "DRO-1" / "DRO-1_init",
            )


if __name__ == "__main__":
    unittest.main()
