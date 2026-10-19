#!/usr/bin/env python3
import csv
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import run_no_srp_initial_mc_campaign as campaign


FAKE_PROPAGATOR = r"""#!/usr/bin/env python3
import csv,json,sys
from pathlib import Path
args=sys.argv[1:]
def value(name):
    return args[args.index(name)+1]
prefix=Path(value("-o"))
count=int(value("--mc-count"))
visual=int(value("--visualization-count"))
hourly=Path(str(prefix)+"_mc_hourly_stats.csv")
fields=["hour","epoch_et","group","sample_count","actual_sample_count",
        "active_dimensions","patch_count","depth_limited_patches",
        *[f"split_dim{i}" for i in range(1,8)],
        "validation_position_error_km","validation_velocity_error_kms"]
with hourly.open("w",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader()
    for label in (10000,20000,50000):
        writer.writerow({"hour":0,"epoch_et":0,"group":"no_srp",
            "sample_count":label,"actual_sample_count":count,
            "active_dimensions":6,"patch_count":1,"depth_limited_patches":0,
            **{f"split_dim{i}":0 for i in range(1,8)},
            "validation_position_error_km":0,
            "validation_velocity_error_kms":0})
header=["hour","keyframe_kind","sample_id","x_km","y_km","z_km",
        "vx_kms","vy_kms","vz_kms","dx_km","dy_km","dz_km",
        "dvx_kms","dvy_kms","dvz_kms"]
def points(path,n,kind):
    with path.open("w",newline="") as stream:
        stream.write("schema=srp_mc_initial_points_v1\n")
        writer=csv.writer(stream); writer.writerow(header)
        for sample_id in range(1,n+1):
            state=[sample_id+j for j in range(6)]
            error=[0.01*(sample_id+j) for j in range(6)]
            writer.writerow([0,kind,sample_id,*state,*error])
points(Path(str(prefix)+"_mc_visualization_points.csv"),visual,
       "hourly_visualization")
points(Path(str(prefix)+"_mc_keyframe_points.csv"),count,"fixed_start")
Path(str(prefix)+"_split_events.csv").write_text(
    "interval_start_hour,interval_end_hour,event_id\n")
Path(str(prefix)+"_mc_report.json").write_text(json.dumps({
    "output_schema_version":"srp_ads_mc_initial_only_v1",
    "nominal_srp_enabled":False,
    "srp_parameter_uncertainty_enabled":False,
    "sample_count":count,"active_dimensions":6,"dace_storage_dimensions":6,
    "srp_eta":0.0,"srp_model":"cannonball","srp_attitude_mode":"no_srp",
    "geometry_fallback_diagnostics":{"real_count":0,"da_count":0,
        "first_real_message":"","first_da_message":""},
    "saved_keyframe_hours":[0.0],
    "events":{"first_split":{"observed":False,
        "interval_start_hour":None,"hour":None,"full_criterion":None},
        "first_unconfirmed":{"observed":False,"hour":None},
        "first_validation_failure":{"observed":False,"hour":None,
            "previous_passing_hour":None}}}))
Path(str(prefix)+"_patch_s0000000_coeff.csv").write_text("snapshot\n")
"""


class NoSrpCampaignTest(unittest.TestCase):
    def make_args(self):
        return SimpleNamespace(
            mc_count=2,
            visualization_count=2,
            days=0.0,
            save_hours=1.0,
            da_order=4,
            max_depth=8,
            pos_tol_km=0.1,
            vel_tol_kms=1e-6,
            force=False,
            refresh_incompatible=False,
        )

    def make_fake(self, root):
        path = Path(root) / "fake_propagator.py"
        path.write_text(FAKE_PROPAGATOR)
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def test_work_items_cover_each_orbit_once(self):
        opms = {
            "DRO-1": Path("/tmp/DRO-1_init.opm.json"),
            "L1Halo-1": Path("/tmp/L1Halo-1_init.opm.json"),
        }
        config = Path("/tmp/config_no_srp.txt")
        items = campaign.make_work_items(opms, config)
        self.assertEqual(
            [
                (group, orbit, opm, item_config)
                for group, orbit, opm, item_config in items
            ],
            [
                ("no_srp", "DRO-1", opms["DRO-1"], config),
                ("no_srp", "L1Halo-1", opms["L1Halo-1"], config),
            ],
        )

    def test_command_uses_no_srp_config_and_omits_eta(self):
        command = campaign.build_command(
            Path("/bin/fake"),
            Path("/tmp/DRO-1_init.opm.json"),
            Path("/tmp/out"),
            Path("/tmp/samples.csv"),
            Path("/tmp/config_no_srp.txt"),
            self.make_args(),
        )
        self.assertNotIn("--eta-sigma", command)
        self.assertEqual(command[command.index("-cfg") + 1],
                         "/tmp/config_no_srp.txt")

    def test_good_result_is_packed_and_published_atomically(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            project = root / "project"
            project.mkdir()
            experiment = root / "20261005SRP"
            samples = root / "samples.csv"
            samples.write_text("fixed\n")
            opm = root / "DRO-1_init.opm.json"
            opm.write_text("{}")
            config = root / "config_no_srp.txt"
            config.write_text("use_srp = false\n")
            executable = self.make_fake(root)
            row = campaign.run_one(
                "no_srp", "DRO-1", opm, executable, project, experiment,
                samples, config, self.make_args(), os.environ.copy()
            )
            self.assertEqual(row["status"], "succeeded", row)
            final = experiment / "results_no_srp_initial_only" / "DRO-1"
            self.assertTrue((final / "DRO-1_visualization_points.h5").is_file())
            self.assertTrue((final / "DRO-1_keyframe_points.h5").is_file())
            report = json.loads((final / "DRO-1_events.json").read_text())
            self.assertFalse(report["nominal_srp_enabled"])
            self.assertEqual(report["srp_attitude_mode"], "no_srp")
            self.assertEqual(report["campaign_force_model"], "no_srp")
            valid, reasons = campaign.audit_result(
                final, "no_srp", "DRO-1", self.make_args()
            )
            self.assertTrue(valid, reasons)
            self.assertFalse(any(final.parent.glob("DRO-1.__old__*")))


if __name__ == "__main__":
    unittest.main()
