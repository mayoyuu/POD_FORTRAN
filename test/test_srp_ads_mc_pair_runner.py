#!/usr/bin/env python3
import csv
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PYTHON=Path('/home/songyu/venvs/astro/bin/python')


class PairedRunnerSmokeTest(unittest.TestCase):
    def test_hour_zero_outputs_are_aligned_and_complete(self):
        opm=ROOT.parent.parent/'OPM'/'L1Halo-1'/'L1Halo-1_init.opm.json'
        self.assertTrue(opm.is_file())
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            samples=tmp/'samples.csv'
            manifest=tmp/'samples.json'
            subprocess.run([
                str(PYTHON),str(ROOT/'tools'/'generate_srp_mc_samples.py'),
                '--csv',str(samples),'--manifest',str(manifest),'--count','256'
            ],check=True,cwd=ROOT,stdout=subprocess.DEVNULL)
            prefix=tmp/'orbit'
            env=os.environ.copy()
            lib_paths=[
                str(ROOT/'lib'),
                str(Path.home()/'local_dace'/'lib'),
                str(Path.home()/'local_openblas'/'lib'),
            ]
            env['LIBRARY_PATH']=':'.join(lib_paths+[env.get('LIBRARY_PATH','')])
            env['LD_LIBRARY_PATH']=':'.join(lib_paths+[env.get('LD_LIBRARY_PATH','')])
            command=[
                'fpm','run','run_srp_ads_mc_pair','--',
                '-opm',str(opm),'-o',str(prefix),'--mc-samples',str(samples),
                '--mc-count','256','--hours','0','--visualization-count','256',
            ]
            result=subprocess.run(command,cwd=ROOT,env=env,text=True,
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode,0,result.stdout[-4000:])

            with (tmp/'orbit_mc_hourly_stats.csv').open(newline='') as handle:
                hourly=list(csv.DictReader(handle))
            with (tmp/'orbit_mc_paired_stats.csv').open(newline='') as handle:
                paired=list(csv.DictReader(handle))
            self.assertEqual(len(hourly),6)
            self.assertEqual(len(paired),3)
            self.assertEqual({int(row['sample_count']) for row in hourly},{10000,20000,50000})
            self.assertEqual({row['group'] for row in hourly},{'treatment','control'})
            self.assertEqual({int(row['actual_sample_count']) for row in hourly},{256})
            self.assertEqual({int(row['actual_sample_count']) for row in paired},{256})
            hourly_required={
                'ads_accuracy_valid','precision_confirmed','maximum_patch_depth',
                'shape_mean_dx_km','shape_mean_dy_km','shape_mean_dz_km',
                'shape_mean_dvx_kms','shape_mean_dvy_kms','shape_mean_dvz_kms',
                'shape_position_max_km','shape_velocity_max_kms',
                'shape_position_max_to_axis1','shape_velocity_max_to_axis1',
                'max_abs_skewness','max_abs_skewness_component',
                'max_abs_excess_kurtosis','max_abs_excess_kurtosis_component',
                'cov_position_xx','cov_position_xy','cov_position_xz',
                'cov_position_yy','cov_position_yz','cov_position_zz',
                'cov_velocity_xx','cov_velocity_xy','cov_velocity_xz',
                'cov_velocity_yy','cov_velocity_yz','cov_velocity_zz',
            }
            paired_required={
                'cov_delta_11','cov_delta_66',
                'position_p95_increment_km','position_p99_increment_km',
                'velocity_p95_increment_kms','velocity_p99_increment_kms',
                'position_p95_amplification','position_p99_amplification',
                'velocity_p95_amplification','velocity_p99_amplification',
            }
            self.assertTrue(hourly_required.issubset(hourly[0]))
            self.assertTrue(paired_required.issubset(paired[0]))

            point_lines=(tmp/'orbit_mc_keyframe_points.csv').read_text().splitlines()
            self.assertEqual(len(point_lines),258)
            rows=list(csv.DictReader(point_lines[1:]))
            self.assertEqual([int(row['sample_id']) for row in rows],list(range(1,257)))
            self.assertTrue(all(row['keyframe_kind']=='fixed_start' for row in rows))
            for row in rows:
                for component in ('dx_km','dy_km','dz_km','dvx_kms','dvy_kms','dvz_kms'):
                    treatment=float(row[f'treatment_{component}'])
                    control=float(row[f'control_{component}'])
                    delta=float(row[f'delta_{component}'])
                    self.assertAlmostEqual(delta,treatment-control,places=12)

            report=json.loads((tmp/'orbit_mc_report.json').read_text())
            self.assertEqual(report['output_schema_version'],'srp_ads_mc_pair_v1')
            self.assertEqual(report['treatment_active_dimensions'],7)
            self.assertEqual(report['control_active_dimensions'],6)
            self.assertEqual(report['dace_storage_dimensions'],7)
            self.assertEqual(report['sample_count'],256)
            self.assertIn('first_split',report['treatment_events'])
            self.assertIn('full_criterion',report['treatment_events']['first_split'])
            self.assertIn('first_split',report['control_events'])
            self.assertTrue(report['pairing_complete'])


if __name__=='__main__':
    unittest.main()
