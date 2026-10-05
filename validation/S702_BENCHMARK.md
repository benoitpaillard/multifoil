# S702: published MSES and wind-tunnel comparison

This is an executed trial of the S702 two-element natural-laminar-flow airfoil
from Somers and Maughmer (2022), NASA/CR-20220014964. See the
[results](results/s702_comparison/S702_RESULTS.md),
[drag comparison](results/s702_comparison/s702_drag_comparison.png), and
[fixed-lift comparison](results/s702_comparison/s702_fixed_lift.png).

The trial uses the existing MultiFoil physics without force calibration. It
adds a reproducible geometry/reference extraction, polar runner, fixed-lift
runner and comparisons. Agreement with published MSES is distinct from
agreement with experiment or validation of an optimization gradient.

## Primary source and selected conditions

Dan M. Somers and Mark D. Maughmer, *Design and Experimental Results for a
Low-Speed, Flapped, Slotted, Natural-Laminar-Flow Airfoil*, December 2022.
[NASA report](https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf).

| Input or reference | Selected value |
| --- | --- |
| Geometry | S702, Figure 1, PDF page 27, printed page 21 |
| Flap deflections | -2, 0, +5, +10 degrees |
| Hinge | x/c = 0.9001, z/c = -0.0114 |
| Reynolds number | 1,000,000 on the original assembly chord |
| Mach number | 0.10 |
| Transition | Free; computational Ncrit = 9 as specified for MSES |
| Moment reference | Original x/c = 0.25, z/c = 0 |
| Published comparisons | Figures 24(a), 23(c), 25(b), 26(c), respectively |
| PDF pages | 87, 85, 90, 94, respectively |
| Selected fixed-lift targets | CL = 0.5, 0.8, 1.2 |

The small plain flap on the aft element moves. The entire aft element does
not rotate. The source drawing contains four separately drawn aft contours;
the runner uses the corresponding contour directly. The fore element and
upstream part of the aft element remain fixed, apart from tiny PDF rounding
differences. Original chord, axes, element positions and finite trailing edges
are retained. The aft trailing edge is below the reference chord line; no
alignment correction is made.

The published theoretical curves are MSES results, extracted from the report.
No MSES executable or rerun was used. The report describes MSES refinement of
the fore element after initial Eppler design and reports successful primary
design objectives after wind-tunnel testing. Its zero-flap MSES maximum-lift
overprediction is about 10-13% across Reynolds numbers. The extracted Re=1e6
curves give approximately 11.6%, providing an independent calibration check.

## Geometry and reference audit

The report does not tabulate original geometry coordinates and directs readers
to Airfoils, Incorporated. Its vector PDF, however, contains line vertices for
both elements at all four flap settings. `extract_s702.py` recovers these
vertices and the plotted MSES paths and experimental marker centers. This is
a **vector-figure reconstruction**, not the author's original input deck or
the measured model contour.

The source PDF is pinned by SHA-256. Geometry contains 121 fore-element
vertices, 111 aft vertices at zero flap, and 109 aft vertices at each other
setting. The original chord-line path defines the common coordinate transform.
The hinge cross recovers (0.900049, -0.011395), close to the printed hinge.
Minimum raw-contour clearance is approximately 0.01931c for every setting.
The plotted maximum vertical fore-element thickness is about 0.12943c, versus
the report's stated airfoil thickness 0.1291c. No thickness rescaling is applied;
this difference is another reason to keep exact-input validation qualified.

`source.json` records drawing indices, axes and transforms. Polars use vector
paths for MSES and circular marker centers for experiment. Independent lift
and drag marker blocks agree in CL within 0.000103; their tiny differences
are retained. Figures are checked with colored overlays on the source raster.
`audit_s702.py` also re-extracts the source into a temporary directory and
requires every generated file to match the distributed data byte for byte.

| Audit artifact | Purpose |
| --- | --- |
| `data/s702/source.json` | Source identity, selected cases and axis calibration |
| `data/s702/geometry_extraction.json` | Original vector vertices and geometry transform |
| `data/s702/geometry_audit.png` | Assembly, slot and actual flap movement |
| `data/s702/extraction_audit.png` | MSES/experiment selection on the source plots |
| `data/s702/audit.json` | Geometry checks and extraction round trip |

Vector precision is not experimental uncertainty. No total uncertainty bound
is established by this trial. In particular, the tunnel model's plastic-film
hinge cover, reported manufacturing departures, and their influence on
transition are not represented. Reference CD at common lift is obtained by
linear interpolation between recorded bracketing markers/path vertices, with
no extrapolation. CSV outputs preserve the brackets. No smoothing of measured
forces, calibration to measured drag, or subtraction of a fitted bias occurs.

## Reproduce the calculations

Run from the project root after installing `requirements.txt`. The distributed
geometry and reference tables suffice for solving and plotting. PyMuPDF is
optional and needed only for source re-extraction and source-overlay audits.

```bash
python -m pip install -r requirements.txt
python -m pip install PyMuPDF  # optional source re-extraction only
python validation/extract_s702.py
python validation/audit_s702.py
python -m unittest -v test_s702.py test_validation.py

python validation/s702_benchmark.py \
  --out validation/results/s702_polar_reproduction
python validation/s702_benchmark.py --mode fixed --targets=.5,.8,1.2 \
  --seed-polar validation/results/s702_polar_reproduction \
  --out validation/results/s702_fixed_reproduction

python validation/s702_benchmark.py --settings=coarse --angles=-4,-2,0,2,4,8 \
  --out validation/results/s702_coarse_polar_reproduction
python validation/s702_benchmark.py --mode fixed --settings=coarse --targets=.5,.8,1.2 \
  --seed-polar validation/results/s702_coarse_polar_reproduction \
  --out validation/results/s702_coarse_fixed_reproduction

python validation/report_s702.py \
  --polar validation/results/s702_polar_reproduction \
  --polar validation/results/s702_coarse_polar_reproduction \
  --fixed validation/results/s702_fixed_reproduction \
  --fixed validation/results/s702_coarse_fixed_reproduction \
  --out validation/results/s702_comparison_reproduction
```

The main polar requests alpha = -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10,
12 degrees. Downward and upward branches use separate solver instances,
initializing at -1 degree. Internal continuation is at most one degree per
step and is archived in diagnostics. Fixed-lift targets span the bucket and
the higher-drag range common to the selected polars. Their incidence brackets
come from adjacent admissible points in the corresponding panel-resolution
polar; each target starts a fresh solver. A failed point cannot provide a
bracket. Ncrit, geometry and flow conditions are unchanged between resolutions.

Base panel counts are 260/180; coarse counts are 180/120. Both use the current
curvature paneling, residual tolerance 1e-7, at most 100 Newton iterations,
wake extension 1c, and fixed-lift tolerance 2e-5 with 20 requested evaluations
per target. Every requested evaluation, including any failure, is archived.
CSV uses NaN and JSON uses null for failed forces. Numerical admissibility
requires convergence, finite positive drag and no detected local supersonic
flow. It is not an experimental validation flag.

Each study pins numerical code, input geometry and source metadata hashes.
Rerunning the same unchanged command resumes completed jobs; a changed source
or setting requires a new output directory. The fixed-lift manifest embeds its
seed polar records and their hashes, avoiding dependence on machine-specific
file paths. The report checks current source hashes before comparing results.
It writes plots, CSV metrics, rankings and an execution summary; the narrative
`S702_RESULTS.md` documents this particular execution.

## What this case can establish

The four flap settings provide an independent check of configuration ranking
and drag-bucket movement. They do not validate local contour-shape gradients,
arbitrary flap interpolation, gap optimization, other Reynolds numbers, forced
transition, roughness robustness, or an unconstrained optimizer. Two panel
levels check sensitivity; they do not establish asymptotic grid convergence.
The run stops at 12 degrees and does not sample the measured stall peaks.
No new MultiFoil maximum-lift or post-stall accuracy claim is made.

S702 remains a useful additional challenge alongside TN808 and NLR7301.
Its smoother slot and laminar-flow design do not remove the existing limitations
on wake merging, extensive separation and transition-sensitive optimization.
