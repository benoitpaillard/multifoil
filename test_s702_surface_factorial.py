"""Frozen factorial protocol, controls, failure coverage and arithmetic regressions."""
import copy
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0,str(Path(__file__).parent/'validation'))
import run_s702_surface_factorial as run


class SurfaceFactorialTests(unittest.TestCase):
    def test_matrix_reuses_four_endpoints_and_repeats_explicit_c(self):
        spec=run.specification()
        self.assertEqual(len(spec['jobs']),5)
        self.assertEqual(len(spec['references']),4)
        self.assertEqual({(j['flap'],j['cell']) for j in spec['jobs'] if not j['repeat_control']},
                         {(f,c) for f in (5,10) for c in ('F','A')})
        control=next(j for j in spec['jobs'] if j['repeat_control'])
        self.assertEqual((control['flap'],control['panels'],control['stage']),(10,[440,320],'surface'))
        self.assertFalse(control['control'])  # The reused worker selects the explicit solver.
        for j in spec['jobs']:
            self.assertEqual(j['source_panels'],[440,320])
            self.assertEqual((j['wake_count_recipe'],j['wake_spacing_recipe']),('fine','fine'))
            self.assertEqual((j['target_cl'],j['wake_length']),(.8,1.))
        self.assertEqual((spec['ncrit'],spec['mach'],spec['reynolds']),(6.25,.07,700000.))

    def test_mixed_controls_keep_source_and_wake_recipes_identical(self):
        spec=run.specification(); setups=[]
        for j in spec['jobs']:
            if j['flap']==10: setups.append(run.numerical_setup(j,spec))
        for elements,settings,recipes in setups[1:]:
            for (_,a),(_,b) in zip(setups[0][0],elements): np.testing.assert_array_equal(a,b)
            for k in ('wake_nodes','wake_first_spacing'):self.assertEqual(settings[k],setups[0][1][k])
            self.assertEqual(recipes,setups[0][2])
        self.assertEqual([v[1]['surface_nodes'] for v in setups],
                         [dict(fore=441,aft=241),dict(fore=341,aft=321),dict(fore=441,aft=321)])

    def test_modified_audited_intermediate_state_is_rejected(self):
        target=run.ROOT/run.PARENT/'d5_contour/aft_state.npz'; real=run.sha
        with patch.object(run,'sha',side_effect=lambda p:'altered' if Path(p)==target else real(p)):
            with self.assertRaisesRegex(ValueError,'Audited reference changed'):run.specification()

    def test_interaction_and_conditional_effects_are_distinct(self):
        e=run.effects(dict(B=100.,F=102.,A=103.,C=112.))
        self.assertEqual(e,dict(fore_at_low_aft=2.,fore_at_high_aft=9.,aft_at_low_fore=3.,
                                aft_at_high_fore=10.,interaction=7.,total=12.))
        self.assertEqual(e,run.effects(dict(B=0.,F=2.,A=3.,C=12.)))
        for bad in (dict(B=1.,F=2.,C=3.),dict(B=1.,F=None,A=2.,C=3.),dict(B=1.,F=2.,A=np.nan,C=3.)):
            with self.assertRaisesRegex(ValueError,'Complete finite'):run.effects(bad)

    def fixtures(self,out,spec):
        for j in spec['jobs']:
            folder=out/j['name'];folder.mkdir()
            source=run.ROOT/run.PARENT/f"d{j['flap']}_surface"
            for filename in ('result.json','fore_state.npz','aft_state.npz'):
                shutil.copyfile(source/filename,folder/filename)

    def test_missing_pair_cannot_produce_complete_attribution(self):
        spec=run.specification()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);self.fixtures(out,spec)
            path=out/'d10_aft_only/result.json';r=run.read(path)
            r.update(feasible=False,CL=None,CD=None,CM=None,alpha_deg=None);run.write(path,r)
            s=run.summarize(out,spec)
            self.assertFalse(s['coverage_complete']);self.assertFalse(s['execution_checks_pass'])
            self.assertIsNone(s['paired_effects_counts']);self.assertEqual(s['individual_effects_counts'],{})
            self.assertIsNone(next(p for p in s['paired_cells'] if p['cell']=='A')['flap_delta_CD_counts'])
            self.assertTrue(s['reference_context']['mses_is_published_curve'])
            self.assertFalse(s['reference_context']['transition_settings_matched'])

    def test_control_array_change_blocks_attribution_even_with_matching_forces(self):
        spec=run.specification()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);self.fixtures(out,spec)
            good=run.summarize(out,spec);self.assertTrue(good['execution_checks_pass'])
            path=out/'d10_surface_control/aft_state.npz'
            with np.load(path) as z:arrays={k:z[k].copy() for k in z.files}
            arrays['U'][0,20]*=1.01;np.savez_compressed(path,**arrays)
            s=run.summarize(out,spec)
            self.assertEqual(s['control_CD_difference_counts'],0.)
            self.assertFalse(s['control_pass']);self.assertIsNone(s['paired_effects_counts'])

    def test_preparation_rejects_changed_manifest_on_resume(self):
        spec=run.specification()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);run.prepare(out,spec)
            changed=copy.deepcopy(spec);changed['jobs'][0]['panels']=[450,240]
            with self.assertRaisesRegex(ValueError,'Code, inputs or settings changed'):run.prepare(out,changed)


if __name__=='__main__':unittest.main()
