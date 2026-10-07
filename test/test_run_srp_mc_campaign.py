#!/usr/bin/env python3
import json
import tempfile
import unittest
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import run_srp_mc_campaign as campaign


def make_opm(path):
    path.parent.mkdir(parents=True,exist_ok=True)
    pos=100.0/3.0
    vel=3.0e-10
    path.write_text(json.dumps({
        'CX_X':pos,'CY_Y':pos,'CZ_Z':pos,
        'CX_DOT_X_DOT':vel,'CY_DOT_Y_DOT':vel,'CZ_DOT_Z_DOT':vel,
    }))


class CampaignTest(unittest.TestCase):
    def test_discovery_is_strict_orbit_whitelist(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for orbit in campaign.ORBITS:
                make_opm(root/orbit/f'{orbit}_init.opm.json')
            make_opm(root/'L2NRHO-2'/'L2NRHO-2_init.opm.json')
            make_opm(root/'supp_batch'/'DRO-1'/'DRO-1_init.opm.json')
            found=campaign.discover_opms(root)
            self.assertEqual(tuple(found),campaign.ORBITS)
            self.assertTrue(all(found[name].parent.name==name for name in found))

    def test_missing_whitelist_orbit_is_an_explicit_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for orbit in campaign.ORBITS[:-1]:
                make_opm(root/orbit/f'{orbit}_init.opm.json')
            with self.assertRaises(campaign.CampaignError):
                campaign.discover_opms(root)

    def test_split_event_audit_checks_count_and_first_event(self):
        hourly=[{'group':'treatment','sample_count':'50000',
                 **{f'split_dim{i}':'0' for i in range(1,8)}}]
        report={'first_split':{'observed':False,'hour':None,
                              'interval_start_hour':None,'full_criterion':None},
                'first_eta_split':{'observed':False,'hour':None,
                                  'interval_start_hour':None,'full_criterion':None}}
        self.assertEqual(campaign.audit_split_group(
            'treatment',report,[],hourly),[])

        hourly[0]['split_dim7']='1'
        reasons=campaign.audit_split_group('treatment',report,[],hourly)
        self.assertIn('treatment_split_event_count',reasons)

        row={'interval_start_hour':'0','interval_end_hour':'1','event_id':'1',
             'step_event_index':'1','bfs_iteration':'2','parent_depth':'0',
             'component':'1','direction':'7'}
        row.update({f'error{i}':'0.2' for i in range(1,7)})
        row.update({f'tolerance{i}':'0.1' for i in range(1,7)})
        row.update({f'excess{i}':'0.1' for i in range(1,7)})
        row.update({f'direction_score{i}':'1.0' for i in range(1,8)})
        row.update({f'parent_center{i}':'0.0' for i in range(1,8)})
        row.update({f'parent_width{i}':'2.0' for i in range(1,8)})
        report['first_split']={'observed':True,'hour':1.0,
            'interval_start_hour':0.0,'full_criterion':{
                'component':1,'direction':7,'parent_depth':0,'bfs_iteration':2,
                'truncation_estimates':[0.2]*6,'component_tolerances':[0.1]*6,
                'positive_excess':[0.1]*6,'direction_scores':[1.0]*7,
                'parent_unit_center':[0.0]*7,'parent_unit_width':[2.0]*7}}
        report['first_eta_split']=report['first_split']
        self.assertEqual(campaign.audit_split_group(
            'treatment',report,[row],hourly),[])

        report['first_split']['full_criterion']['direction']=6
        reasons=campaign.audit_split_group('treatment',report,[row],hourly)
        self.assertIn('treatment_first_split_mismatch',reasons)

    def test_publish_validates_before_replacing_existing_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            final=root/'final'
            staging=root/'staging'
            final.mkdir()
            (final/'marker').write_text('old')
            staging.mkdir()
            (staging/'marker').write_text('bad')
            with self.assertRaises(campaign.CampaignError):
                campaign.publish_directory(staging,final,
                    lambda path:(False,['invalid']))
            self.assertEqual((final/'marker').read_text(),'old')
            self.assertTrue(staging.exists())

            (staging/'marker').write_text('new')
            campaign.publish_directory(staging,final,lambda path:(True,[]))
            self.assertEqual((final/'marker').read_text(),'new')
            self.assertFalse(staging.exists())


if __name__=='__main__':
    unittest.main()
