"""Selected-state and clearance audit regressions, using small saved fixtures."""
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import report_s702_wake_extent as report


class WakeAuditTests(unittest.TestCase):
    def setUp(self):
        self.study=report.ROOT/'validation/results/s702_wake_extent_20260922'
        self.inputs=report.Inputs()
        self.spec=self.inputs.json(self.study/'study.json')
        self.runtime=self.inputs.json(self.study/'runtime.json')
        self.job=self.spec['jobs'][0]

    def test_nonzero_mach_export_uses_corrected_velocity(self):
        state=report.audit(self.study/self.job['name'],self.job,self.spec,self.runtime,self.inputs)
        a=state['arrays']['fore']; surface=state['surfaces']['fore']
        self.assertGreater(np.max(abs(a['U'][3,:len(surface)][::-1]-surface['ue_over_vinf'])),1e-4)
        self.assertTrue(state['result']['feasible'])

    def test_saved_velocity_tamper_is_detected(self):
        # Inputs.add deliberately confines manifests to ROOT; keep temporary fixture there.
        with tempfile.TemporaryDirectory(dir=self.study.parent,prefix='wake_audit_test_') as tmp:
            folder=Path(tmp)/self.job['name']
            shutil.copytree(self.study/self.job['name'],folder)
            path=folder/'fore_state.npz'
            with np.load(path) as saved: values={k:saved[k].copy() for k in saved.files}
            values['U'][3,20]+=.001
            np.savez_compressed(path,**values)
            with self.assertRaisesRegex(ValueError,'corrected velocity'):
                report.audit(folder,self.job,self.spec,self.runtime,report.Inputs())

    def test_open_wake_clearance_and_crossing(self):
        body=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        # Closing this U-shaped open polyline would falsely intersect the body.
        wake=np.array([[-1.,.5],[-1.,2.],[2.,2.],[2.,.5]])
        self.assertAlmostEqual(report.clearance(wake,body),1.)
        self.assertEqual(report.clearance(np.array([[-1.,.5],[2.,.5]]),body),0.)


if __name__=='__main__':unittest.main()
