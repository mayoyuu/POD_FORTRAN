#!/usr/bin/env python3
import csv
import json
import math
import tempfile
import unittest
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import summarize_srp_mc_pair as summary


GROUP_BASE={
    'active_dimensions':'7','patch_count':'1','depth_limited_patches':'0',
    'validation_failed':'0','position_axis1_rms_km':'10',
    'velocity_axis1_rms_kms':'0.01','position_p99_km':'20',
    'velocity_p99_kms':'0.02',
}
for name in ('mean_dx_km','mean_dy_km','mean_dz_km','mean_dvx_kms',
             'mean_dvy_kms','mean_dvz_kms'):
    GROUP_BASE[name]='0'
for name in ('skew_dx','skew_dy','skew_dz','skew_dvx','skew_dvy','skew_dvz',
             'excess_kurt_dx','excess_kurt_dy','excess_kurt_dz',
             'excess_kurt_dvx','excess_kurt_dvy','excess_kurt_dvz'):
    GROUP_BASE[name]='0'

PAIR_BASE={
    'center_position_difference_km':'0',
    'center_velocity_difference_kms':'0',
    'delta_position_rms_km':'0','delta_velocity_rms_kms':'0',
    'delta_position_p95_km':'0','delta_velocity_p95_kms':'0',
    'delta_position_p99_km':'0','delta_velocity_p99_kms':'0',
    'delta_position_p999_km':'0','delta_velocity_p999_kms':'0',
    'delta_position_max_km':'0','delta_velocity_max_kms':'0',
}


class SummaryTest(unittest.TestCase):
    def write_inputs(self,directory,hours=(0,1,2),fail_hour=2,nan_hour=None):
        directory=Path(directory)
        hourly=directory/'orbit_mc_hourly_stats.csv'
        paired=directory/'orbit_mc_paired_stats.csv'
        group_fields=['hour','epoch_et','group','sample_count','actual_sample_count',
                      *GROUP_BASE]
        pair_fields=['hour','epoch_et','sample_count','actual_sample_count',*PAIR_BASE]
        with hourly.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=group_fields)
            writer.writeheader()
            for hour in hours:
                for group in ('treatment','control'):
                    for count in (10000,20000,50000):
                        row=dict(GROUP_BASE,hour=hour,epoch_et=1000+3600*hour,
                                 group=group,sample_count=count,actual_sample_count=count)
                        row['active_dimensions']='7' if group=='treatment' else '6'
                        row['validation_failed']='1' if hour==fail_hour and group=='treatment' else '0'
                        row['mean_dx_km']=str(float('nan')) if hour==nan_hour else str(hour)
                        writer.writerow(row)
        with paired.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=pair_fields)
            writer.writeheader()
            for hour in hours:
                for count in (10000,20000,50000):
                    row=dict(PAIR_BASE,hour=hour,epoch_et=1000+3600*hour,
                             sample_count=count,actual_sample_count=count)
                    row['delta_position_rms_km']=2*hour
                    row['delta_position_max_km']=3*hour
                    writer.writerow(row)
        report=directory/'orbit_mc_report.json'
        report.write_text(json.dumps({'output_schema_version':'srp_ads_mc_pair_v1',
                                      'pairing_complete':True}))
        return hourly,paired,report

    def test_valid_and_all_windows_and_convergence(self):
        with tempfile.TemporaryDirectory() as tmp:
            hourly,paired,report=self.write_inputs(tmp)
            result=summary.summarize_files(hourly,paired,report,expected_hours=3)
            metric=result['metrics']['delta_position_rms_km']
            self.assertEqual(metric['all']['endpoint_value'],4.0)
            self.assertEqual(metric['all']['peak_value'],4.0)
            self.assertEqual(metric['all']['peak_hour'],2.0)
            self.assertEqual(metric['all']['time_mean'],2.0)
            self.assertEqual(metric['all']['trapezoid_integral'],4.0)
            self.assertEqual(metric['valid']['endpoint_hour'],1.0)
            self.assertEqual(metric['valid']['endpoint_value'],2.0)
            self.assertEqual(metric['valid']['time_mean'],1.0)
            self.assertEqual(metric['valid']['trapezoid_integral'],1.0)
            self.assertTrue(result['sampling_convergence']['converged'])
            self.assertEqual(result['first_validation_failure_hour'],2.0)
            self.assertTrue(result['data_complete'])
            flat=Path(tmp)/'summary.csv'
            summary.write_flat_csv(result,flat)
            self.assertIn('sampling_converged',flat.read_text())

    def test_non_hourly_cadence_uses_requested_time_grid(self):
        with tempfile.TemporaryDirectory() as tmp:
            hourly,paired,report=self.write_inputs(
                tmp,hours=(0,0.5,1.0),fail_hour=99)
            result=summary.summarize_files(
                hourly,paired,report,expected_hours=3,save_hours=0.5)
            self.assertTrue(result['data_complete'])
            metric=result['metrics']['delta_position_rms_km']['all']
            self.assertEqual(metric['endpoint_hour'],1.0)
            self.assertEqual(metric['endpoint_value'],2.0)
            self.assertEqual(metric['trapezoid_integral'],1.0)

    def test_gap_is_reported_without_interpolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            hourly,paired,report=self.write_inputs(tmp,hours=(0,2),fail_hour=99)
            result=summary.summarize_files(hourly,paired,report,expected_hours=3)
            self.assertFalse(result['data_complete'])
            self.assertEqual(result['first_data_issue_hour'],1.0)
            self.assertFalse(result['metrics']['delta_position_rms_km']['all']['valid'])

    def test_nonfinite_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            hourly,paired,report=self.write_inputs(tmp,fail_hour=99,nan_hour=1)
            result=summary.summarize_files(hourly,paired,report,expected_hours=3)
            self.assertFalse(result['data_complete'])
            self.assertEqual(result['first_data_issue_hour'],1.0)


if __name__=='__main__':
    unittest.main()
