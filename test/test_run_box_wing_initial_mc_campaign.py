#!/usr/bin/env python3
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
import run_box_wing_initial_mc_campaign as campaign


FAKE_PROPAGATOR = r"""#!/usr/bin/env python3
import csv,json,os,sys
from pathlib import Path
args=sys.argv[1:]
def value(name):
    return args[args.index(name)+1]
prefix=Path(value("-o"))
opm=Path(value("-opm"))
count=int(value("--mc-count"))
visual=int(value("--visualization-count"))
config=Path(value("-cfg"))
if os.environ.get("FAKE_COMMAND_LOG"):
    with Path(os.environ["FAKE_COMMAND_LOG"]).open("a") as stream:
        stream.write(json.dumps(args)+"\n")
if "DRO-2" in opm.name:
    raise SystemExit(7)
attitude="moon" if "moon" in config.name else "earth"
hourly=Path(str(prefix)+"_mc_hourly_stats.csv")
fields=["hour","epoch_et","group","sample_count","actual_sample_count",
        "active_dimensions","patch_count","depth_limited_patches",
        *[f"split_dim{i}" for i in range(1,8)],
        "validation_position_error_km","validation_velocity_error_kms"]
with hourly.open("w",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader()
    for label in (10000,20000,50000):
        writer.writerow({"hour":0,"epoch_et":0,"group":attitude,
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
    "sample_count":count,"active_dimensions":6,"dace_storage_dimensions":6,
    "srp_eta":0.0,"srp_model":"box_wing","srp_attitude_mode":attitude,
    "saved_keyframe_hours":[0.0],
    "events":{"first_split":{"observed":False,
        "interval_start_hour":None,"hour":None,"full_criterion":None},
        "first_unconfirmed":{"observed":False,"hour":None},
        "first_validation_failure":{"observed":False,"hour":None,
            "previous_passing_hour":None}}}))
Path(str(prefix)+"_patch_s0000000_coeff.csv").write_text("snapshot\n")
"""


class CampaignTest(unittest.TestCase):
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

    def test_command_uses_attitude_config_and_omits_eta(self):
        args = self.make_args()
        command = campaign.build_command(
            Path("/bin/fake"),
            Path("/tmp/DRO-1_init.opm.json"),
            Path("/tmp/out"),
            Path("/tmp/samples.csv"),
            Path("/tmp/config_box_wing_moon.txt"),
            args,
        )
        self.assertNotIn("--eta-sigma", command)
        self.assertEqual(command[command.index("-cfg") + 1],
                         "/tmp/config_box_wing_moon.txt")
        self.assertEqual(command[command.index("--mc-samples") + 1],
                         "/tmp/samples.csv")

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
            config = root / "config_box_wing_earth.txt"
            config.write_text("earth\n")
            executable = self.make_fake(root)
            row = campaign.run_one(
                "earth", "DRO-1", opm, executable, project, experiment,
                samples, config, self.make_args(), os.environ.copy()
            )
            self.assertEqual(row["status"], "succeeded", row)
            final = (experiment / "results_box_wing_initial_only" /
                     "earth" / "DRO-1")
            self.assertTrue((final / "DRO-1_visualization_points.h5").is_file())
            self.assertTrue((final / "DRO-1_keyframe_points.h5").is_file())
            report = json.loads((final / "DRO-1_events.json").read_text())
            self.assertEqual(report["srp_attitude_mode"], "earth")
            self.assertEqual(len(report["opm_sha256"]), 64)
            self.assertFalse(any(final.parent.glob("DRO-1.__old__*")))

    def test_failed_orbit_does_not_prevent_later_orbit(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            project = root / "project"
            project.mkdir()
            experiment = root / "20261005SRP"
            samples = root / "samples.csv"
            samples.write_text("fixed\n")
            config = root / "config_box_wing_moon.txt"
            config.write_text("moon\n")
            executable = self.make_fake(root)
            opm_bad = root / "DRO-2_init.opm.json"
            opm_good = root / "DRO-1_init.opm.json"
            opm_bad.write_text("{}")
            opm_good.write_text("{}")
            rows = campaign.run_jobs(
                [
                    ("moon", "DRO-2", opm_bad, config),
                    ("moon", "DRO-1", opm_good, config),
                ],
                executable, project, experiment, samples, self.make_args(),
                os.environ.copy(), jobs=2,
            )
            status = {(row["attitude"], row["orbit"]): row["status"]
                      for row in rows}
            self.assertEqual(status[("moon", "DRO-2")], "failed")
            self.assertEqual(status[("moon", "DRO-1")], "succeeded")

    def test_work_items_cover_both_attitudes(self):
        opms = {"DRO-1": Path("/tmp/DRO-1_init.opm.json")}
        configs = {
            "earth": Path("/tmp/config_box_wing_earth.txt"),
            "moon": Path("/tmp/config_box_wing_moon.txt"),
        }
        items = campaign.make_work_items(("earth", "moon"), opms, configs)
        self.assertEqual(
            [(attitude, orbit) for attitude, orbit, _, _ in items],
            [("earth", "DRO-1"), ("moon", "DRO-1")],
        )


if __name__ == "__main__":
    unittest.main()
