"""Tests for strict, source-backed complete-case comparison (no propagation)."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


def module():
    path=ROOT/'tools/analyze_srp_modes_case.py'
    if not path.exists(): raise AssertionError('complete-case analysis not implemented')
    spec=importlib.util.spec_from_file_location('case_analysis',path)
    result=importlib.util.module_from_spec(spec); spec.loader.exec_module(result)
    return result


class CaseAnalysis(unittest.TestCase):
    def fixture(self,folder):
        m=module(); folder.mkdir()
        inputs={}
        for name in ('opm','config','mode_csv','mode_metadata'):
            path=folder/name; path.write_text(name)
            inputs[name]={'path':str(path),'sha256':m.sha(path)}
        report=dict(schema_version='srp_modes_ads_v1',formal_result=True,n_variables=8,
            parameterization='SRP_MODES',variable_roles=['x','y','z','vx','vy','vz','s','u'],
            epoch0_tdb_s=100.,duration_hours=1.,save_hours=1.,shape_count=2,validation_count=2,
            da_order=4,max_depth=8,domain_sigma=3.,position_tolerance_km=.1,velocity_tolerance_kms=1e-6,
            rel_tol=1e-12,abs_tol=1e-12,dt_min_s=1e-6,dt_max_s=3600.,
            da_handles_baseline=0,da_handles_final=0,validation_failed=False,precision_unconfirmed=False,
            depth_limited_patches=0,ads_numerical_error={'max_position_km':0.,'max_velocity_kms':0.},
            zero_mode_error={'max_position_km':1e-8,'max_velocity_kms':1e-12},
            provenance={'inputs':inputs,'source_files_sha256':{'src/core.f90':'core'}},
            conditional_srp_offset={'final_earth_position_km':10.,'final_moon_position_km':20.})
        history=[]; shape=[]
        for hour in (0,1):
            history.append([hour,100+hour*3600,*([0]*6)])
            for point,dx in ((1,3),(2,-3)): shape.append([hour,point,dx,4,0,0,0,2])
        self.csv(folder/'c2_history.csv',['hour','epoch_et',*m.NOMINAL],history)
        self.csv(folder/'c2_shape.csv',['hour','point_id',*m.ERROR],shape)
        self.csv(folder/'c2_initial_points.csv',['point_id',*[f'u{i}' for i in range(1,9)]],
                 [[1,*([0]*8)],[2,*([1]*8)]])
        self.csv(folder/'c2_probes.csv',['hour','point_id',*[f'ads_{k}' for k in ('x','y','z','vx','vy','vz')],
                 *[f'real_{k}' for k in ('x','y','z','vx','vy','vz')]],
                 [[hour,point,*([0]*12)] for hour in (0,1) for point in (1,2)])
        report['artifacts_sha256']={name:m.sha(folder/name) for name in
            ('c2_history.csv','c2_shape.csv','c2_initial_points.csv','c2_probes.csv')}
        (folder/'c2_report.json').write_text(json.dumps(report))
        return report

    @staticmethod
    def csv(path,fields,rows):
        with path.open('w',newline='') as stream:
            writer=csv.writer(stream); writer.writerow(fields); writer.writerows(rows)

    def test_known_norms_and_common_nominal(self):
        m=module()
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'a'; self.fixture(folder)
            run=m.load_run(folder)
            envelopes=m.envelopes(run,run['nominal'])
            self.assertEqual(envelopes[1],(5.,2.))
            other=dict(run,nominal={0:(0,)*6,1:(1,0,0,0,0,0)})
            result=m.compare_runs(run,other,'same')
            self.assertAlmostEqual(result['max_state_difference_position_km'],1.)
            self.assertAlmostEqual(result['position_envelope_relative_change'],(32**.5-5)/5)
            self.assertFalse(result['convergence_confirmed'])
            self.assertTrue(m.compare_runs(run,run,'same')['convergence_confirmed'])

    def test_reject_hash_missing_duplicate_and_nonfinite(self):
        m=module()
        for fault in ('hash','missing','duplicate','nonfinite'):
            with self.subTest(fault=fault),tempfile.TemporaryDirectory() as tmp:
                folder=Path(tmp)/'a'; report=self.fixture(folder)
                path=folder/'c2_shape.csv'
                with path.open() as stream: rows=list(csv.reader(stream))
                if fault=='missing': rows.pop()
                if fault=='duplicate': rows.append(rows[-1])
                if fault=='nonfinite': rows[-1][2]='NaN'
                with path.open('w',newline='') as stream: csv.writer(stream).writerows(rows)
                if fault=='hash': path.write_text(path.read_text()+'\n')
                else:
                    report['artifacts_sha256'][path.name]=m.sha(path)
                    (folder/'c2_report.json').write_text(json.dumps(report))
                with self.assertRaises(m.CaseError): m.load_run(folder)

    def test_report_numerical_gates_and_identity(self):
        m=module()
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'a'; report=self.fixture(folder)
            baseline=m.load_run(folder)
            for key,value in (('precision_unconfirmed',True),('validation_failed',True),
                              ('n_variables',7),('da_handles_final',1)):
                bad=dict(report); bad[key]=value
                (folder/'c2_report.json').write_text(json.dumps(bad))
                with self.subTest(key=key),self.assertRaises(m.CaseError): m.load_run(folder)
            bad=dict(report); bad['ads_numerical_error']={'max_position_km':.2,'max_velocity_kms':1e-9}
            (folder/'c2_report.json').write_text(json.dumps(bad))
            with self.assertRaises(m.CaseError): m.load_run(folder)
            different=dict(baseline,report=dict(baseline['report'],domain_sigma=4))
            with self.assertRaisesRegex(m.CaseError,'domain_sigma'): m.compare_runs(baseline,different,'same')
            other=dict(baseline,report=dict(report,mode_csv_sha256='different'))
            other['report']['provenance']=dict(report['provenance'],inputs=dict(report['provenance']['inputs']))
            other['report']['provenance']['inputs']['mode_csv']={'path':'unused','sha256':'different'}
            with self.assertRaisesRegex(m.CaseError,'mode_csv'): m.compare_runs(baseline,other,'tight_rkf')
            # Only the mode-table refinement may change table identity.
            self.assertTrue(m.compare_runs(baseline,other,'half_hour')['convergence_confirmed'])

    def test_probe_failure_not_hidden_by_report(self):
        m=module()
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'a'; report=self.fixture(folder)
            path=folder/'c2_probes.csv'
            with path.open() as stream: data=list(csv.reader(stream))
            data[-1][2]='1.0'
            with path.open('w',newline='') as stream: csv.writer(stream).writerows(data)
            report['artifacts_sha256'][path.name]=m.sha(path)
            (folder/'c2_report.json').write_text(json.dumps(report))
            with self.assertRaisesRegex(m.CaseError,'probe'): m.load_run(folder)

    def test_vector_metrics_keep_proxy_error_separate(self):
        m=module()
        result=m.vector_metrics({0:(0,)*6,1:(3,4,0,0,0,2)},{0:(0,)*6,1:(0,)*6})
        self.assertEqual(result['max_position_km'],5)
        self.assertEqual(result['max_velocity_kms'],2)
        self.assertEqual(result['max_position_hour'],1)

    def test_reused_control_snapshots_without_propagation(self):
        m=module()
        base=m.load_run(ROOT/'../../SRP/boxwing_modes_DRO-2_init_c2_15d_20261014')
        result,nominal=m.audit_controls(ROOT/'../../SRP/SRP_UCERTAINTY_PROP_WITHINIT_ETA0/DRO-2',
            ROOT/'../../SRP/SRP_UCERTAINTY_PROP_WITHINIT/DRO-2',base)
        self.assertEqual(len(nominal),361)
        self.assertEqual(result['C0']['identity_status'],'matched_opm_initial_basis_non_srp_config')
        self.assertFalse(result['C0']['pressure_comparison']['exactly_equal'])
        self.assertFalse(result['C1_30_subdomain']['fully_provenanced_matched_campaign'])
        self.assertAlmostEqual(result['C1_30_subdomain']['source_eta_half_span'],.9)
        self.assertAlmostEqual(result['C1_30_subdomain']['restricted_eta_half_span'],.3)
        self.assertLess(result['C1_30_subdomain']['zero_state_eta_plus_03_position_km'],100.)
        audits=m.audit_snapshot_points(base)
        self.assertEqual(len(audits),2)
        self.assertLess(max(x['max_position_km'] for x in audits),1e-8)


if __name__=='__main__': unittest.main()
