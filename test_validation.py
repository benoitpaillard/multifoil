"""Checks that experimental comparisons cannot hide missing data or bad matching."""
import unittest
from validate import compare_polar, compare_transition, error_stats, sampled_peak, verify_data, at_angle, reference_info


class ValidationIntegrity(unittest.TestCase):
    def test_failed_solve_and_blank_experimental_drag(self):
        refs=[{'alpha_deg':'0','CL':'0','CD':''},{'alpha_deg':'10','CL':'1','CD':'.02'}]
        preds=[{'alpha_deg':'0','converged':'1','CL':'.1','CD':'.007'},
               {'alpha_deg':'10','converged':'0','CL':'99','CD':'99'}]
        rows=compare_polar(refs,preds,[0,10])
        self.assertEqual(len(rows),3)
        self.assertEqual(error_stats(rows)['unavailable'],2)
        self.assertAlmostEqual(error_stats(rows)['mae'],.1)

    def test_no_alpha_interpolation_or_extrapolation(self):
        rows=compare_polar([{'alpha_deg':'1','CL':'.1'}],
            [{'alpha_deg':'0','converged':'1','CL':'0'},
             {'alpha_deg':'2','converged':'1','CL':'.2'}],[1])
        self.assertIsNone(rows[0]['predicted'])
        with self.assertRaises(ValueError): at_angle([{'alpha_deg':'1'},{'alpha_deg':'1'}],1)

    def test_transition_ranges_and_missing_transition(self):
        refs=[{'alpha_deg':'2','element':'main','side':'upper','x_low':'.3',
               'x_high':'.4','reading_allowance':'.01'}]
        pred=[{'alpha_deg':'2','element':'main','side':'upper',
               'transition_x_over_local_chord':'.405'}]
        row=compare_transition(refs,pred,[2])[0]
        self.assertAlmostEqual(row['distance_outside_interval'],.005)
        self.assertEqual(row['distance_outside_with_reading_allowance'],0)
        self.assertIsNone(compare_transition(refs,[],[2])[0]['predicted'])

    def test_stall_needs_reference_neighbors_and_convergence(self):
        refs=[{'alpha_deg':str(a),'CL':str(cl)} for a,cl in [(10,1),(15,1.3),(20,.8)]]
        preds=[{'alpha_deg':str(a),'CL':str(cl),'converged':'1'} for a,cl in [(10,1),(15,1.4),(20,1.5)]]
        self.assertFalse(sampled_peak(refs,preds,[10,15])['reference_peak_bracket_covered'])
        self.assertLess(sampled_peak(refs,preds,[10,15,20])['predicted_drop_across_reference_peak'],0)
        preds[-1]['converged']='0'
        self.assertIsNone(sampled_peak(refs,preds,[10,15,20])['error_at_reference_peak'])

    def test_microphone_observable_uses_reattachment(self):
        refs=[{'alpha_deg':'0','element':'main','side':'lower','x_low':'.49',
               'x_high':'.51','reading_allowance':'.01','observable':'attached_turbulence'}]
        pred=[{'alpha_deg':'0','element':'main','side':'lower',
               'transition_x_over_local_chord':'.45'}]
        region={'start':{'s_over_cref':.4},'end':{'s_over_cref':.55,'x_over_local_chord':.5},
                'reaches_trailing_edge':False}
        diagnostic={'boundary_layers':{'main':{'lower':{
            'transition':{'s_over_cref':.5},'separation_regions':[region]}}}}
        row=compare_transition(refs,pred,[0],[(0,diagnostic)])[0]
        self.assertEqual(row['predicted'],.5)
        self.assertEqual(row['model_transition_onset'],.45)
        region['reaches_trailing_edge']=True
        self.assertIsNone(compare_transition(refs,pred,[0],[(0,diagnostic)])[0]['predicted'])

    def test_distributed_reference_hashes(self):
        self.assertGreater(len(verify_data()['data_sha256']),10)

    def test_mses_is_classified_as_solver_reference(self):
        info=reference_info({'benchmark':{'reference_kind':'published_solver','reference_solver':'MSES'}})
        self.assertEqual(info['reference_kind'],'published_solver')
        self.assertEqual(info['reference_label'],'Published MSES')
        self.assertEqual(reference_info({})['reference_kind'],'experiment')
        with self.assertRaises(ValueError):
            reference_info({'benchmark':{'reference_kind':'published_solver'}})

    def test_mses_extraction_reproduces_bundled_reference(self):
        import importlib.util
        from pathlib import Path
        path=Path(__file__).parent/'validation/extract_mses_reference.py'
        spec=importlib.util.spec_from_file_location('mses_extract',path)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        rows,audit=mod.extract()
        self.assertEqual(mod.csv_text(rows),(mod.DIRECTORY/'polar.csv').read_text())
        self.assertEqual(mod.csv_text(audit),(mod.DIRECTORY/'digitization_pixels.csv').read_text())
        self.assertTrue(all(r['reference_kind']=='published_solver' for r in rows))

    def test_drag_only_reference_has_no_lift_peak(self):
        refs=[{'alpha_deg':'4','CD':'.006','CD_low':'.005','CD_high':'.007'}]
        preds=[{'alpha_deg':'4','converged':'1','CD':'.008','CL':'.4'}]
        rows=compare_polar(refs,preds,[4])
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['quantity'],'CD')
        self.assertAlmostEqual(rows[0]['error'],.002)
        self.assertAlmostEqual(rows[0]['distance_outside_interval'],.001)
        self.assertIsNone(sampled_peak(refs,preds,[4])['reference_sampled_maximum'])

    def test_mses_drag_extraction_reproduces_bundled_reference(self):
        import importlib.util
        from pathlib import Path
        path=Path(__file__).parent/'validation/extract_mses_drag_reference.py'
        spec=importlib.util.spec_from_file_location('mses_drag_extract',path)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        rows=mod.extract()
        self.assertEqual(mod.csv_text(rows),(mod.DIRECTORY/'polar.csv').read_text())
        self.assertEqual([r['alpha_deg'] for r in rows],[0,4,8,10,12,14,16,18])
        self.assertTrue(all('CL' not in r and r['reference_solver']=='MSES' for r in rows))


if __name__=='__main__': unittest.main()
