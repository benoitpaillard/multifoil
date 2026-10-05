# Experimental validation cases

The [reference inventory and comparison candidates](REFERENCES.md) consolidate
the sources already used and the September 2026 literature search, including
published MSES comparisons, experimental data, access checks and remaining gaps.

The current development handover is [CODEX_HANDOVER.md](../CODEX_HANDOVER.md).
The latest [flap-angle comparison](results/s702_flap_effect/README.md) shows
absolute drag and changes from zero flap, with both original and calibrated
MultiFoil settings. All plots are reproducible from archived results without
new flow calculations.

This suite compares the prototype with published wind-tunnel measurements and
published MSES results. References are explicitly labeled as experimental or
code-to-code. It
does **not** certify accuracy. All comparisons retain numerical failures and
report coverage. The original benchmark cases retain their stated inputs.
The later S702 transition sweep explicitly calibrates an effective Ncrit to
drag data; aerodynamic closures and the solver's default Ncrit remain unchanged.

The [completed S702 sweep](results/s702_ncrit_analysis/NCRIT_SWEEP_RESULTS.md)
uses a common Ncrit on both elements and finds a broad useful region around
6–7. The selected Ncrit=6.25 reduces calibration MAPE from 8.84% to 4.85%;
Ncrit=6 minimizes absolute drag and flap-difference RMSE. Independent pressure,
finer-panel, Mach-zero and reserved Reynolds-number checks are reported
separately from calibration. No contour-shape gradient is validated by this fit.
The reserved check gives 5.70% mean drag error versus 10.15% at Ncrit=9 after
two documented initialization recoveries; the original run completed seven of
nine targets. The worst reserved point remains 24.04% high. Original failures,
an unsuccessful cold retry, and the successful warm replay are retained.
Run `python validation/run_s702_calibration.py --out validation/results/reproduced_s702_calibration --workers 8`
to reproduce the complete study; use a fresh directory inside the project.

The [S702 discrepancy study](results/s702_discrepancy/DISCREPANCY_RESULTS.md)
tests Ncrit=6 as an explicit sensitivity case on four flap configurations and
three lift targets, isolates changes on each element, and adds a pressure check
from NASA Figure 22(b). It identifies a strong fore-transition sensitivity but
does not establish a physical correction or validated contour gradients. The
new pressure-data checks run with `python -m unittest -v test_s702_pressure`.

## Run

From the project directory, after installing `requirements.txt`:

```bash
python validate.py --list
python validate.py --out validation_run
python validate.py --case nlr7301 --out nlr_run
python validate.py --case nlr7301_mses --out mses_benchmark
python validate.py --case naca0012_mses_drag --out mses_drag_benchmark
python validate.py --compare-only --out validation/results
python -m unittest -v test_validation.py
```

The last comparison command regenerates plots and error tables from the included
solver outputs. It does not rerun the flow solver. A fresh run requires a new or
empty case output directory. Runs take several minutes; higher resolutions and
failed continuation attempts can take substantially longer.

Explicit sensitivity/refinement runs:

```bash
python validate.py --case s809_somers --ncrit 11 --out s809_n11
python validate.py --case nlr7301 --panel-scale 1.5 --out nlr_refined
python validate.py --case nlr7301_full_polar --out nlr_full
```

The archived results use Ncrit=9 and the base meshes. The three commands above
are available studies, **not studies already executed**. The full NLR polar is
an optional, expensive extension through its experimental lift loss. Its default
two-incidence case does not establish multi-element stall accuracy.

`validate.py` saves a checkpoint after each requested angle. Exit code 0 means
every requested solve and comparison is available; 2 means some are unavailable;
1 means an input or provenance error. None of these is an aerodynamic accuracy
pass. `--compare-only` checks the backend, case conditions, element configuration,
and requested angles in the saved metadata. Reference files are SHA-256 checked
against `sources.json` before comparison.

## Included cases

| Case | Conditions | What is compared |
|---|---|---|
| NACA 0012, Ladson | Re=5.97 million, Mach=0.15, free transition; 260 panels | 18 tabulated CL and CM values, 16 CD values; abrupt lift loss at high incidence |
| S809, Somers | Re=2 million, smooth/free transition; Mach=0 low-speed approximation; 240 panels | 15 digitized CL observations; 14 transition/attached-turbulence ranges; usable lift limit and separated-flow stress points |
| NLR7301 main/flap | Re=2.51 million, Mach=0.185, nominal 20° flap and 2.6% gap; 260/180 panels | Two CL and two CD observations around 6° and 13.1°; four exact digitized incidence values |

The default suite also includes **nlr7301_mses**, a ten-point CL comparison to
the published MSES NLR7301 curve. See [MSES.md](MSES.md) for its source, extraction
audit, incomplete setting match and executed results. MSES was not run here.
It also includes **naca0012_mses_drag**, eight CD observations from COBEM2005
Figure3 at Re=6 million, Mach=0.3, with 260 panels and assumed Ncrit=9.
See [DRAG.md](DRAG.md) for the source audit, reading precision and results.
This single-element check does not validate multi-element wake drag.

Every case uses **Ncrit=9 as an explicit assumption**. The reports do not supply
a directly transferable e^N threshold. NACA tunnel turbulence during the tests
was unknown. S809's Delft tunnel reported turbulence levels of 0.02–0.04%; its
Mach number is not explicitly specified. These condition uncertainties matter
especially for transition and drag. Use sensitivity runs before drawing a
stronger conclusion. No forced-transition cases are scored because the solver
does not currently implement trips.

## Reference provenance and qualifications

**NACA 0012.** Charles L. Ladson, NASA TM-4074 (1988), Table I, printed page 6
(PDF page 8), Re=5.97×10⁶, Mach=0.15, **free transition**.
[Original NASA report](https://ntrs.nasa.gov/citations/19880019495).
The 18 rows were transcribed and checked against the rendered table. At 17.65°
and 18.65° the report gives no CD; those entries remain blank. The geometry uses
the standard finite-trailing-edge NACA polynomial, with the −0.1015 x⁴ term.
It is different from the closed-TE demo and the extended/rescaled geometry used
by the NASA Turbulence Modeling Resource. The quarter-chord moment convention
and common reference chord are retained. The report's repeatability discussion
does not provide a universal uncertainty bound for these observations.

**S809.** D. M. Somers, NREL/SR-440-6918 (1997).
[Original report in the UNT federal document collection](https://digital.library.unt.edu/ark:/67531/metadc681921/).
Geometry is Table 2, printed page 15 (PDF 20); the leading-edge origin is inserted
between the nearest upper and lower coordinates. CL is digitized from
Figure 17(c)(i), printed page 68 (PDF 73). Incidences come from Figure 8's
corrected-angle captions. Pixel coordinates and source-image overlays are in
`data/s809/` and `audit/`. CL reading allowances are ±0.015, widened to ±0.025
where two closely spaced sweep branches cannot be reliably distinguished.
These are **digitization allowances, not measured experimental confidence
intervals**. Approximate alpha reading uncertainty is 0.15°.

Somers identifies usable maximum lift near **CL=1.01**, around 9.22°, and
explicitly questions the validity of later, higher lift values because of
extensive separation. The stall comparison uses this first maximum, while
retaining the later observations as stress comparisons. The larger plotted
CL near 15° is not presented as an experimentally accepted usable CLmax.

S809 transition is digitized from Figure 15(c), printed page 60 (PDF 65). There
are two important differences from an e^N onset coordinate:

* The microphone detects **attached turbulent flow**. Where a laminar bubble
  exists, its symbol marks turbulent reattachment. For those observations, the
  comparison uses the end of the model's separated region containing transition,
  or the transition point if the model remains attached. A separated region
  extending to the TE yields no attached-turbulence prediction.
* Bars mark onset-to-completion ranges; point clusters span the two measured
  span stations. Upper-surface observations at 6° and above compare the modeled
  instantaneous transition with those ranges. An additional ±0.015 chord reading
  allowance is kept separate in the CSV.

Figure 15 alone uses **uncorrected tunnel angles**, whereas Figure 8/17 use
wall-corrected angles. We infer the correspondence of nominal sweep settings
0, 4, 6, 8, 10, 12, 15° to the Figure 8 caption angles −0.01, 4.10, 6.16,
8.20, 10.21, 12.22, 15.24°. This is an explicitly documented correspondence,
**not a published calibration for Figure 15**. Both angles remain in the
reference and comparison files. Treat transition errors as qualified diagnostics
until that mapping and transition criterion are confirmed with the original
measurement records. Additional nominal angles were included in the solver sweep
for diagnostic inspection, but are not substituted for failed corrected-angle
comparisons.

**NLR7301.** Original experiment: B. van den Berg and B. Oskam, NLR MP 79034 U
(1979), *Boundary layer measurements on a two-dimensional wing with flap and a
comparison with calculations*.
[Original report metadata](https://repository.tudelft.nl/record/uuid:5f1bfb3e-71d9-45a8-bb55-c4cdc6f25182).
The full original is campus-restricted and was not accessed. We use openly
published ECARP geometry and experimental digitizations from
[Maximilian Ranneberg's Viiflow examples](https://github.com/mranneberg/viiflow-examples/tree/93bef6c6b779db4605f6bcbb2caf1d88a51c9b85/NLR7301%202-Element%20Airfoil),
pinned to commit `93bef6c6b779db4605f6bcbb2caf1d88a51c9b85`.
The repository's MIT notice is included. Its digitizations are attributed to the
original experiment; they are not original machine-readable measurements.
The benchmark conditions and chord definitions are also described in
[Neto et al., ICAS 2006](https://www.icas.org/icas_archive/ICAS2006/PAPERS/512.PDF).

Keep both ECARP contours in their **assembled coordinates**, with nominal cruise
reference chord 1. The base repaneled geometry has a closest main-TE/flap
distance of about 0.02652 reference chords; it is not adjusted to fit an exact
0.026 clearance. Do not renormalize the main or flap separately and do not
rotate the already deflected flap another 20°. Default geometry is nominal;
the Ranneberg example's optional deformation correction (+0.25° CCW about
(0.94, −0.011), plus +0.0025 vertically) is not applied. That difference must be
resolved before a tighter accuracy claim. The source example also explores
Ncrit/Mach/angle choices; none was fitted here.

The CL and CD digitizations use slightly different measured x coordinates, so
four exact incidence values are solved for two nominal operating points. The
number of decimals reflects digitization output, not experimental precision.
All 20 force observations remain in `data/nlr7301/polar.csv` for the optional
full-polar case. CD is compared with the model's sum of independent-wake
Squire–Young estimates; the model does not resolve wake merging or detailed cove
recirculation. These effects can cause genuine drag disagreement.

`NLRCP55.csv` and `NLRCP126.csv` are retained **unscored** for audit. The two
incidences share 5 exact main-element and 8 exact flap pressure-coordinate pairs;
the original pressure and incidence labels have not been reconciled. The source
plots their Y columns as −Cp and shifts flap x for display. The suite makes no
Cp or NLR transition accuracy claim from these files.

## Read the results

`results/REPORT.md` contains the executed results and plots. Each case also has
`force_errors.csv`, `comparison.json`, complete surface/wake outputs, convergence
diagnostics, and the actual analysis geometry. S809 adds
`transition_errors.csv` and `transition_comparison.png`.

Errors are **prediction minus reference**. MultiFoil solution values are not
interpolated or extrapolated, and no fitted zero-lift correction is used. MSES
references sample a published graphical curve, sometimes between its native
markers, as documented separately. Bias, MAE, RMSE and maximum error use
available comparisons; their denominator and missing count are always reported.
One drag count is 0.0001 CD. The |alpha|≤10° and >10° summaries are simple bins,
not attached-flow or stall classifications. For digitized CL, center errors and
distance outside the reading range are both available.

Reference peak comparisons require its adjacent reference incidences to be in
the requested sweep; missing convergence makes the relevant peak comparison
unavailable. Solver nonconvergence is never classified as physical stall.
Agreement at two multi-element incidences does not validate stall, pressure
distributions, merging wakes or a different geometry. Grid and Ncrit sensitivity,
better matched tunnel conditions and independent benchmark replications remain
necessary before using this prototype for design decisions.

## Follow-on optimization-readiness pilot

See [OPTIMIZATION_PILOT.md](OPTIMIZATION_PILOT.md) for the 26-combination fixed-lift
camber/placement sensitivity study. Its new solves are tracked separately from
the five existing benchmark cases. The official MIT NLR MSES point is now
qualified plot context; it does not change the aggregate benchmark scores.
Three adaptive fine-grid checks follow up the observed transition-dependent
camber ranking reversal, with their own study manifest and results.
