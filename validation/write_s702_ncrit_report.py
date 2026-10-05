#!/usr/bin/env python3
"""Write the reproducible narrative after the sweep and transfer reports."""
import argparse
import csv
import json
from pathlib import Path

from report_s702_ncrit_sweep import ROOT, RESULTS, read, score


def table(headers, rows):
    def cell(value):
        if value is None:
            return 'not available'
        if isinstance(value,float):
            return f'{value:.2f}'
        return str(value)
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(cell(v) for v in row)+' |' for row in rows])


def metric(value):
    return 'not available (incomplete)' if value is None else f'{value:.2f}'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analysis',type=Path,default=RESULTS/'s702_ncrit_analysis')
    args = ap.parse_args()
    out = args.analysis.resolve()
    summary = read(out/'summary.json')
    transfer = read(out/'transfer_summary.json')
    ns = summary['best']['ncrit']
    scores = list(csv.DictReader((out/'scores.csv').open()))
    points = list(csv.DictReader((out/'points.csv').open()))
    failed_points = [p for p in points if p['feasible']=='False']
    baseline = next(s for s in scores if float(s['ncrit'])==9.)
    folds = summary['coarse_grid_cross_validation']
    pairs = list(csv.DictReader((out/'pairwise.csv').open()))
    chosen_pairs = [p for p in pairs if float(p['ncrit'])==ns and int(p['first_flap'])==5 and int(p['second_flap'])==10]
    matching_lifts = [float(p['target_cl']) for p in chosen_pairs if p['same_sign']=='True']
    mismatching_lifts = [float(p['target_cl']) for p in chosen_pairs if p['same_sign']=='False']
    holdout_scores = transfer['holdout']['scores']
    holdout_selected = next(r for r in holdout_scores if r['ncrit']==ns)
    common = transfer['holdout'].get('common_subset')
    local = transfer['holdout'].get('after_local_retries')
    local_outcome = ''
    local_worst = ''
    if local:
        local_selected = next(r for r in local['scores'] if r['ncrit']==ns)
        local_outcome = (f'A separate local initialization follow-up completes '
            f'{local_selected["evaluated"]} of nine cases with the same frozen Ncrit and '
            'numerical core. '+
            (f'Its nine-point MAPE is {local_selected["mape_percent"]:.2f}%.' if local_selected['complete'] else
             'Its nine-point score remains unavailable.'))
        recovered = [r for r in csv.DictReader((out/'holdout_after_retries.csv').open())
                     if float(r['ncrit'])==ns and r['feasible']=='True']
        if recovered:
            worst = max(recovered,key=lambda r:abs(float(r['error_percent'])))
            local_worst = (f'The largest reserved point error is {abs(float(worst["error_percent"])):.2f}% '
                f'at flap {worst["flap"]} degrees, CL={worst["target_cl"]}: '
                f'CD={float(worst["CD"]):.7f} versus measured {float(worst["CD_experiment"]):.7f}. '
                f'The published MSES N=9 value there is {float(worst["CD_mses_n9"]):.7f}. '
                'The mean improvement therefore coexists with a substantial pointwise discrepancy.')
    holdout_outcome = (f'The frozen setting completes all nine reserved Reynolds-number cases, '
        f'with MAPE {holdout_selected["mape_percent"]:.2f}%.' if holdout_selected['complete'] else
        f'The frozen setting completes {holdout_selected["evaluated"]} of nine reserved Reynolds-number '
        'cases. The full reserved-validation score is unavailable because of numerical '
        'nonconvergence; the completed subset is a separate diagnostic.')
    required = [(f,c) for f in [-2,0,5,10] for c in [.5,.8,1.2]]
    fine_scores = [dict(ncrit=n,**score([r for r in transfer['fine'] if r['ncrit']==n],required))
                   for n in sorted({r['ncrit'] for r in transfer['fine']})]
    resolution_rows = []
    for n in sorted({r['ncrit'] for r in transfer['fine']}):
        checks = [c for r in transfer['fine'] if r['ncrit']==n for c in r.get('resolution',[])]
        if checks:
            resolution_rows.append([n,
                f'{max(c["maximum_turbulent_Hk_drop"] for c in checks):.4f}',
                f'{max(c["maximum_turbulent_ue_step"] for c in checks):.4f}'])
    fine_shifts = [abs(r['fine_minus_base_counts']) for r in transfer['fine'] if r['fine_minus_base_counts'] is not None]
    mach_shifts = [abs(r['mach0_minus_mach01_counts']) for r in transfer['mach_zero'] if r['feasible']]
    control_shifts = [abs(r['protocol_minus_legacy_counts']) for r in transfer.get('control',[])
                      if r['protocol_minus_legacy_counts'] is not None]
    coarse_values = summary['coarse_ncrit_values']
    refined_values = summary['refined_ncrit_values']
    lines = [f'# S702 transition sweep: selected Ncrit = {ns:g}', '',
        f'The lowest complete calibration error occurs at Ncrit={ns:g} on the sampled grid. '
        f'Mean absolute relative drag error changes from {float(baseline["mape_percent"]):.2f}% at Ncrit=9 '
        f'to {summary["best"]["mape_percent"]:.2f}%. This is an effective calibration of MultiFoil '
        'for the tested S702 configurations. It does not identify a measured tunnel disturbance '
        'threshold or establish accurate contour-design gradients.', '',
        holdout_outcome+' '+local_outcome, '',
        local_worst, '',
        '## What was run', '',
        f'The integer sweep covers {coarse_values}; the additional sampled values are '
        f'{refined_values}. Each candidate uses both elements at the same Ncrit, '
        'flaps -2, 0, 5 and 10 degrees, and CL=0.5, 0.8 and 1.2. Calibration conditions '
        'are Re=1,000,000, M=0.1, 260 fore panels and 180 aft panels. Geometry and viscous '
        'closures are unchanged. The main solver default remains Ncrit=9.', '',
        f'There are {summary["scored_points"]} candidate/operating-point entries, '
        f'{summary["feasible_points"]} converged targets and {summary["reused_points"]} reused targets. '
        f'The new calibration runs retain {summary["recorded_requested_evaluations"]} requested flow '
        f'evaluations, including {summary["failed_requested_evaluations"]} failed evaluations. '
        'Those counts exclude internal continuation steps, which are preserved in each flow metadata file. '
        'A failed requested evaluation can be followed by a successful documented retry.', '',
        f'An execution-service interruption also affected {len(summary.get("interrupted_attempts",[]))} '
        'unfinished calibration targets. Their partial records are archived separately under '
        '`interrupted_attempts`; the counts above describe the finalized target attempts. '
        'The resume helper executed the exact archived worker source and verified every '
        'numerical dependency. Completed targets were reused, and unfinished targets restarted '
        'without claiming to restore an unavailable full Newton state.', '',
        'A candidate must cover all 12 points to receive a global score. For an incomplete '
        'candidate, the lower bound below assumes zero error at every missing point; it is '
        'not an imputed drag prediction. Solver failure is not classified as stall.', '',
        table(['Ncrit','Converged / required','MAPE (%)','RMSE (counts)','Worst error (%)','Minimum possible MAPE (%)'],[
            [s['ncrit'],f"{s['evaluated']} / {s['required']}",
             float(s['mape_percent']) if s['mape_percent'] else None,
             float(s['rmse_counts']) if s['rmse_counts'] else None,
             float(s['maximum_error_percent']) if s['maximum_error_percent'] else None,
             float(s['minimum_possible_mape_percent'])] for s in scores]), '',
        'One drag count is 0.0001. MAPE gives each point the same relative-error weight. '
        'RMSE in counts is a separate objective and can prefer a different threshold; '
        'no ranking penalty or experimental-uncertainty weighting is hidden in the calibration.', '',
        table(['Failed calibration Ncrit','Flap (degrees)','Target CL','Final target error'],[
            [float(p['ncrit']),int(p['flap']),float(p['target_cl']),
             read(ROOT/p['result_path'])['error']] for p in failed_points]), '',
        table(['Diagnostic objective','Best sampled Ncrit','Objective value'],[
            [label,summary['metric_optima'][key]['ncrit'],summary['metric_optima'][key][key]]
            for key,label in [('mape_percent','Relative drag MAPE (%)'),
                              ('rmse_counts','Absolute drag RMSE (counts)'),
                              ('pairwise_delta_rmse_counts','Flap-difference RMSE (counts)')]]), '',
        'The flap-difference objective compares all six configuration pairs at each of the '
        'three lift targets. It is relevant to optimization, but its optimum is reported '
        'as a diagnostic and does not replace the predeclared MAPE selection.', '',
        f'The MAPE advantage over Ncrit=6 is only '
        f'{float(next(s["mape_percent"] for s in scores if float(s["ncrit"])==6.))-summary["best"]["mape_percent"]:.3f} '
        'percentage points. Both relative-error weighting and panel resolution therefore '
        'matter when interpreting the precise winning value.', '',
        '![Sweep and drag differences](ncrit_sweep.png)', '',
        '![Every calibration point against experiment and published MSES](ncrit_pointwise.png)', '',
        '## Does the choice transfer across flaps?', '',
        'For each fold, choose Ncrit on three flaps and evaluate the excluded flap. '
        'Only the predeclared integer grid is eligible. Globally selected fractional '
        'refinements are excluded from this check, preventing their adaptive placement '
        'from contaminating the withheld flap. The original flap measurements had already '
        'been inspected in earlier studies, so this remains an internal grouped check.', '',
        table(['Excluded flap','Selected Ncrit','Training MAPE (%)','Excluded-flap MAPE (%)','N=9 excluded-flap MAPE (%)'],[
            [r['heldout_flap'],r['selected_ncrit'],r['training_mape_percent'],r['test_mape_percent'],r['baseline_n9_mape_percent']] for r in folds]), '',
        f'Combined excluded-flap MAPE (%): {metric(summary["cross_validated_score"]["mape_percent"])}.', '',
        'For the specific 5-degree versus 10-degree disagreement, negative delta means '
        'the 10-degree flap has lower drag:', '',
        table(['CL','Selected model: CD(10)-CD(5), counts','Experiment, counts','Same raw sign'],[
            [float(p['target_cl']),float(p['predicted_delta_counts']),float(p['measured_delta_counts']),p['same_sign']] for p in chosen_pairs]), '',
        'Small raw signs must be considered alongside the mesh sensitivities below.', '',
        f'The selected threshold matches the raw experimental 5-degree/10-degree '
        f'ordering at CL={matching_lifts}; it misses the ordering at '
        f'CL={mismatching_lifts}. A better overall drag fit does not establish '
        'accurate differences at every lift target.', '',
        '## Reserved Reynolds-number check', '',
        'The selected threshold was frozen before the selected-parameter transfer run. '
        'The reserved figures cover Re=700,000, M=0.07, flaps 0, 5 and 10 degrees and the '
        'same three lift targets. Their values were extracted and visually audited before '
        'calibration was finalized, but they do not enter the selection loss and were not '
        'used to retune the threshold. A fresh Ncrit=9 run supplies the numerical baseline.', '',
        f'The execution-service interruption affected '
        f'{len(transfer["holdout"].get("interrupted_attempts",[]))} unfinished reserved-study '
        'targets in addition to the calibration interruptions above. Their partial '
        'attempts and exact-source restart provenance are retained separately.', '',
        'This tests transfer across both Reynolds number and tunnel speed. It does not '
        'isolate Reynolds number from possible changes in the incoming disturbance spectrum.', '',
        table(['Source','Converged / required','MAPE (%)','RMSE (counts)','Worst error (%)'],[
            [f'MultiFoil N={r["ncrit"]:g}',f'{r["evaluated"]} / {r["required"]}',r['mape_percent'],r['rmse_counts'],r['maximum_error_percent']] for r in holdout_scores]+[
            ['Published MSES N=9','9 / 9',transfer['holdout']['published_mses_n9']['mape_percent'],
             transfer['holdout']['published_mses_n9']['rmse_counts'],transfer['holdout']['published_mses_n9']['maximum_error_percent']]]), '',
        holdout_outcome, '',
        (f'Diagnostic on the identical completed subset only: {len(common["points"])} of nine '
         f'conditions, excluding (flap, CL)={common["excluded_points"]}. '
         'This does not substitute for the complete reserved check.' if common else ''), '',
        (table(['Common-subset source','Points','MAPE (%)','RMSE (counts)'],[
            [f'MultiFoil N={r["ncrit"]:g}',r['evaluated'],r['mape_percent'],r['rmse_counts']]
            for r in common['scores']]+[
            ['Published MSES N=9',len(common['points']),common['published_mses_n9']['mape_percent'],
             common['published_mses_n9']['rmse_counts']]]) if common else ''), '',
        ('### Separate local initialization follow-up' if local else ''), '',
        local_outcome, '',
        ('The original incidence seeds came from Re=1,000,000. At the reserved '
         'Reynolds number, two bracket evaluations can both fall below the target lift; '
         'symmetric expansion can then traverse a difficult lower-lift region. '
         'The local follow-up estimates incidence from two admissible model lift values '
         'and cold-starts with a 0.05-degree half-width. It uses no experimental drag '
         'and changes no transition threshold, geometry, panel count or convergence '
         'tolerance. This procedure was added after the reserved failures, so the '
         'original protocol results above remain separate. If cold initialization '
         'fails, a warm follow-up replays up to two initial admissible source evaluations '
         'before using the local bracket. The first successful initialization is used; '
         'drag error never chooses between initialization paths.' if local else ''), '',
        (table(['Flap','CL','Initialization','Converged','Used in follow-up score'],[
            [r['flap'],r['target_cl'],r['initialization'],r['feasible'],r['selected_for_score']]
            for r in local['retries']]) if local else ''), '',
        (table(['After local retries','Converged / required','MAPE (%)','RMSE (counts)'],[
            [f'MultiFoil N={r["ncrit"]:g}',f'{r["evaluated"]} / {r["required"]}',
             r['mape_percent'],r['rmse_counts']] for r in local['scores']]) if local else ''), '',
        local_worst, '',
        ('Hatched bars identify successful local retries. `local_retry_points.csv` '
         'and the retry manifests link every new result to its original failed target. '
         'This recovery separates a numerical initialization problem from the '
         'remaining aerodynamic error; it does not establish reliable shape gradients.' if local else ''), '',
        '![Reserved Reynolds-number comparison](ncrit_holdout.png)', '',
        '## Numerical and independent-observable checks', '',
        'The finer mesh uses 340 fore / 240 aft panels. The comparison below is against '
        'a matched-protocol base mesh, with all failures retained. A difference between two '
        'meshes is an observed sensitivity, not a rigorous discretization-error bound.', '',
        ('The separate base-mesh Ncrit=6 continuation control changes drag by at most '
         f'{max(control_shifts):.4f} counts relative to the reused results. '
         'Fine Ncrit=6 drag changes are measured relative to that control, separating '
         'panel resolution from the change in incidence sequence.' if control_shifts else
         'No converged continuation-control comparisons are available.'), '',
        table(['Fine-mesh Ncrit','Converged / required','MAPE (%)','RMSE (counts)'],[
            [r['ncrit'],f'{r["evaluated"]} / {r["required"]}',r['mape_percent'],r['rmse_counts']] for r in fine_scores]), '',
        f'Largest absolute fine-minus-base drag change among converged checked targets: '
        f'{metric(max(fine_shifts,default=None))} counts. The detailed file `fine_points.csv` identifies '
        'which cases changed; `transfer_summary.json` also records local Hk and velocity '
        'steps for transition-region resolution checks.', '',
        table(['Fine-mesh Ncrit','Largest turbulent Hk drop per panel','Largest turbulent Ue/Uinf step per panel'],
              resolution_rows), '',
        'The XFOIL manual gives approximate post-transition resolution guides of '
        'an Hk drop below 1 and a velocity step below 0.05. These local diagnostics '
        'supplement the two-mesh comparison; they are not a proof of grid independence.', '',
        f'The matched-lift Mach-zero checks change drag by at most {metric(max(mach_shifts,default=None))} '
        'counts relative to M=0.1 at the selected threshold. They assess the intended '
        'incompressible use and are not compared as matched-Mach experimental data.', '',
        'The pressure check uses flap 0, CL=1.04, Re=1,000,000, M=0.1 and NASA Figure '
        '22(b). The caption rounds lift to two decimal places. Cp was excluded from the '
        'calibration loss. Repeated plotting joins on the lower-surface curves are excluded '
        'from the tap comparison. See `pressure_errors.csv` for chord-weighted errors.', '',
        (table(['Surface','N=9 Cp RMSE','Selected Cp RMSE'],[
            [element+' '+side,
             f'{next(r["rms_chord_weighted"] for r in transfer["pressure"]["errors"] if r["source"]=="MultiFoil N=9" and r["element"]==element and r["side"]==side):.4f}',
             f'{next(r["rms_chord_weighted"] for r in transfer["pressure"]["errors"] if r["source"].endswith("selected") and r["element"]==element and r["side"]==side):.4f}']
            for element,side in [('fore','upper'),('aft','upper'),('fore','lower'),('aft','lower')]])
         if transfer['pressure']['complete'] else ''), '',
        'The pressure curves retain visible systematic differences from the taps; '
        'drag calibration does not establish complete boundary-layer or pressure accuracy.', '',
        ('![Pressure comparison](ncrit_pressure.png)' if transfer['pressure']['complete'] else
         'The pressure target did not converge: '+transfer['pressure']['error']), '',
        '![Predicted transition locations](ncrit_transition.png)', '',
        '## Consequences for shape optimization', '',
        f'The sampled values within one percentage point of the minimum calibration MAPE '
        f'are {summary["near_optimal_sampled_values"]}. This is an engineering sensitivity '
        'set, not a confidence interval or evidence that every intervening value converges. '
        'The exact winning grid value should not be treated as a measured physical constant.', '',
        'Keep transition assumptions fixed while changing geometry. Evaluate each candidate '
        'shape at multiple lift points and across plausible transition scenarios, and require '
        'its predicted improvement to survive numerical refinement. The earlier Ncrit=6/9 '
        'shape-sensitivity checks remain relevant; fitting these baseline polars does not '
        'justify discarding the broader disturbance uncertainty.', '',
        'A common Ncrit cannot supply missing wake-induced/bypass transition, wake merging, '
        'roughness, the hinge-cover film or as-built geometry. The coordinate files remain '
        'a reconstruction of the published geometry figure. Pressure discrepancies, unresolved '
        'drag differences and numerical failures are evidence about those remaining limits. '
        'No geometry was optimized or accepted in this sweep.', '',
        '## Files and reproduction', '',
        'The archive includes the solver, sweep driver, calibration metrics, extraction scripts, '
        'all requested flow records, states, source snapshots, tests, plots and this report. '
        '`frozen_selection.json` records the exact follow-up selections. The full execution '
        'protocol is in `../../NCRIT_SWEEP_PROTOCOL.md`.', '',
        'Reproduce the complete calculation and reporting sequence in a fresh directory:', '',
        '```bash',
        'python validation/run_s702_calibration.py --out validation/results/reproduced_s702_calibration --workers 8',
        '```', '',
        'The full run can take hours. To run the individual calibration stages:', '',
        '```bash',
        'export OPENBLAS_NUM_THREADS=1', 'export OMP_NUM_THREADS=1',
        'python validation/s702_ncrit_sweep.py --out validation/results/my_ncrit_coarse --workers 8',
        'python validation/s702_ncrit_sweep.py --out validation/results/my_ncrit_refined --ncrit 5.5,5.75,6.25,6.5 --workers 8',
        'python validation/report_s702_ncrit_sweep.py --coarse validation/results/my_ncrit_coarse --refined validation/results/my_ncrit_refined --out validation/results/my_ncrit_analysis',
        'python -m unittest test_s702_sweep test_s702 test_s702_pressure', '```', '',
        'New runs use a corrected legacy-output metadata path and an expanded fallback bracket '
        'when a cold retry fails to bracket lift. The first coarse runner and every later '
        'execution version are archived. These changes affect metadata/continuation, not '
        'the aerodynamic equations. Existing result folders are guarded against changed inputs.', '',
        'The calibrated MultiFoil and published MSES curves use different transition '
        'thresholds. Their differences combine solver effects with that change in input; '
        'this comparison does not establish that MultiFoil is more accurate than a '
        'similarly calibrated MSES run.', '',
        'References: Somers and Maughmer (2022), '
        '[NASA/CR-20220014964](https://ntrs.nasa.gov/citations/20220014964), '
        'Figures 23(c), 24(a), 25(b), 26(c) for calibration; Figures 23(b), 25(a), 26(b) '
        'for reserved Reynolds-number checks; Figure 22(b) for pressure. '
        '[XFOIL documentation](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt) '
        'describes the meaning and limits of Ncrit and transition-region panel resolution. '
        'All MSES curves are published Ncrit=9 results; MSES was not run here.', '']
    (out/'NCRIT_SWEEP_RESULTS.md').write_text('\n'.join(lines))
    print(out/'NCRIT_SWEEP_RESULTS.md')


if __name__ == '__main__':
    main()
