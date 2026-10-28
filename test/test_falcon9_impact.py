import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
ENTRY=ROOT/'.sh/falcon9_impact.py'


class FalconImpactPreparationTests(unittest.TestCase):
    def test_four_distinct_covariances_share_state_samples_between_scenarios(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'result'
            proc=subprocess.run([sys.executable,str(ENTRY),'--prepare-only','--count','32',
                                 '--output',str(out)],text=True,capture_output=True)
            self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
            campaign=json.loads((out/'campaign.json').read_text())
            self.assertEqual(len(campaign['jobs']),12)
            durations=[]
            for name in ['FAC_2604_1','FAC_2605_1','FAC_2606_1','FAC_2607_1']:
                opm=json.loads((ROOT/'OPM/Falcon9'/(name+'.opm')).read_text())
                meta=[json.loads((out/name/c/'population.json').read_text()) for c in ['C0','C1','T1']]
                self.assertTrue(all(m['epoch']==opm['EPOCH'] for m in meta))
                self.assertEqual(len({m['input_state_sha256'] for m in meta}),1)
                c0=np.genfromtxt(out/name/'C0/initial_points.csv',delimiter=',',names=True)
                self.assertTrue(np.all(c0['eta']==0))
                t1=np.genfromtxt(out/name/'T1/initial_points.csv',delimiter=',',names=True)
                self.assertTrue(np.all(np.abs(t1['eta'])<=1))
                durations.append(campaign['orbits'][name]['utc_duration_hours'])
            self.assertTrue(all(a>b for a,b in zip(durations,durations[1:])))
            self.assertTrue(all(j['saved_points']==0 for j in campaign['jobs']))
            self.assertTrue(all(j['status']=='prepared' for j in campaign['jobs']))

    def test_c1_t1_selection_prepares_only_eight_final_snapshots(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'result'
            proc=subprocess.run([sys.executable,str(ENTRY),'--prepare-only','--count','32',
                '--cases','C1','T1','--final-only','--output',str(out)],text=True,capture_output=True)
            self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
            campaign=json.loads((out/'campaign.json').read_text())
            self.assertEqual(len(campaign['jobs']),8)
            self.assertEqual(campaign['cases'],['C1','T1'])
            self.assertTrue(campaign['final_only'])
            self.assertTrue(all(j['case'] in ['C1','T1'] for j in campaign['jobs']))
            self.assertFalse(any(out.glob('*/C0')))

    def test_output_is_never_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc=subprocess.run([sys.executable,str(ENTRY),'--prepare-only','--count','2',
                                 '--output',tmp],text=True,capture_output=True)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn('exists',proc.stderr)


if __name__=='__main__': unittest.main()
