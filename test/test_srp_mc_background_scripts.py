#!/usr/bin/env python3
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
START=ROOT/'tools'/'start_srp_mc_campaign_background.sh'
STOP=ROOT/'tools'/'stop_srp_mc_campaign.sh'


class BackgroundScriptTest(unittest.TestCase):
    def test_start_records_live_supervisor_and_stop_terminates_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            fake_python=directory/'fake_python'
            fake_python.write_text('#!/usr/bin/env bash\nsleep 30\n')
            fake_python.chmod(0o755)
            experiment=directory/'experiment'
            environment=os.environ.copy()
            environment['SRP_MC_PYTHON']=str(fake_python)
            environment['SRP_MC_EXPERIMENT_ROOT']=str(experiment)
            started=subprocess.run(
                [str(START),'--dry-run'],cwd=ROOT,env=environment,
                text=True,capture_output=True,check=True)
            self.assertIn('PID=',started.stdout)
            pid=int((experiment/'logs'/'campaign.pid').read_text())
            self.assertTrue(Path(f'/proc/{pid}').exists())
            stopped=subprocess.run(
                [str(STOP)],cwd=ROOT,env=environment,
                text=True,capture_output=True,check=True)
            self.assertIn('Sent SIGTERM to campaign process group',stopped.stdout)


if __name__=='__main__':
    unittest.main()
