# Evidence for using MultiFoil in shape optimization

The consolidated [14 September handover](../CODEX_HANDOVER.md) and
[development priorities](../DEVELOPMENT_PLAN.md) summarize the current state.
The latest [flap-angle plots](results/s702_flap_effect/README.md) make the
remaining high-lift ranking discrepancy explicit. Older literature-stage
statements below are historical; the executed studies above their references
and in the handover define what has subsequently been completed.

Research checked 9 September 2026. This note identifies evidence and a proposed
validation programme. An [incompressible shape-search workflow](INCOMPRESSIBLE_DESIGN.md)
now adds bounded contour modes, paired transition scenarios, fixed-lift
objectives and mandatory finer-grid screening. It provides numerical design
screening; independent validation of contour improvements remains outstanding.
The [executed S702 incompressible trial](results/s702_incompressible_design/INCOMPRESSIBLE_RESULTS.md)
retained the baseline: the best small camber change improved mean drag by only
0.19 counts at Ncrit 9 and 0.47 counts at Ncrit 6, below the predeclared 0.5-count
threshold. Mach zero also failed to resolve the published flap-ranking
disagreement. This is progress in the design workflow, not acceptance of an
experimentally demonstrated shape improvement.

The [new discrepancy investigation](results/s702_discrepancy/DISCREPANCY_RESULTS.md)
isolates the 5°/10° drag disagreement. At CL=0.8, changing only the fore-element
transition threshold from 9 to 6 adds about 30 drag counts to the 5° configuration,
bringing its drag close to experiment. Changing both elements to Ncrit=6 recovers
the experimental ranking, but still leaves errors at other lift/flap combinations.
The −2° case at CL=0.5 actually becomes worse, and an independent Cp comparison
shows remaining pressure differences. This is an identified model sensitivity,
not a validated physical correction. The default transition model and the
status of contour-design derivatives remain unchanged.

The [completed transition sweep](results/s702_ncrit_analysis/NCRIT_SWEEP_RESULTS.md)
now samples 14 common Ncrit values over 3–12, including fractional refinement.
On 12 S702 calibration points, mean absolute relative drag error falls from
8.84% at Ncrit=9 to 4.85% at 6.25. Absolute drag and flap-difference RMSE both
prefer 6; its relative error is only 0.10 percentage points higher. Finer panels
give 4.75% at 6.25 and 4.81% at 6. The precise minimum depends on the objective;
the more useful result is a region around 6–7 that warrants sensitivity checks.
Grouped checks select 6 or 7, with combined excluded-flap error of 5.63%, but
the −2° flap becomes worse. The selected threshold is frozen before reserved
Reynolds-number scoring, and the report retains all numerical failures.
Baseline calibration does not establish correct contour-design derivatives.

With Ncrit=6.25 frozen, the reserved Re=700,000/M=0.07 comparison reduces mean
relative drag error from 10.15% to 5.70%. The original incidence protocol
completed seven of nine points; local incidence estimation and a warm replay
recovered the other two without changing the aerodynamic equations or
convergence tolerances. The warm replay exactly reproduced its two source
lift/drag states before reaching the target. One reserved point (10° flap,
CL=0.8) still overpredicts drag by 24.04%. At CL=1.2, the model still reverses
the measured 5°/10° drag ranking at both tested Reynolds numbers. These limits
prevent treating the fitted threshold as sufficient validation for shape
optimization; maintain transition sensitivity and independent geometry checks.

The earlier numerical pilot is also preserved; see
[OPTIMIZATION_PILOT.md](OPTIMIZATION_PILOT.md) and its
[executed results](results/optimization_pilot/REPORT.md). It checks fixed-lift
geometry sensitivity, without independent validation of rankings.

The pilot identifies an actual caution: the preferred camber direction reverses
between Ncrit 9 and Ncrit 6. A fine-grid follow-up checks this sensitivity.
Nominal numerical consistency therefore does not establish a transition-robust
design improvement; see the executed report for the values and qualifications.

The literature supports using coupled viscous–inviscid methods for aerodynamic
design. Confidence in this particular multi-element implementation requires an
additional test: whether it predicts the performance change between geometries.
Our existing polar comparisons establish only limited agreement for fixed
geometries. They do not establish reliable design gradients or design rankings.

The current recommendation is to develop MultiFoil as a constrained, subcritical
design-screening model, with independent verification of shortlisted geometries.
Maximum-lift, stall, very small-gap and wake-interaction optimization need further
model development and validation before accepting predicted optima.

## The most useful references

### 1. Somers and Maughmer (2022): experimentally tested MSES-assisted design

Dan M. Somers and Mark D. Maughmer, *Design and Experimental Results for a
Low-Speed, Flapped, Slotted, Natural-Laminar-Flow Airfoil*, NASA/CR–20220014964,
December 2022. [NASA report](https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf).

This is the strongest positive design evidence found. Eppler provided initial
element shapes; MSES refined the S702 fore element within the assembled
two-element configuration. Wind-tunnel tests examined four flap deflections,
four Reynolds numbers and free/fixed transition. The authors report achieving
their principal design goals, with qualified overall agreement. MSES still
overpredicted maximum lift by about10–13% at zero flap deflection and missed the
sharpness of stall. See “Execution,” “Tests,” and “Comparison of Theory and
Experiment,” especially Figures23–27. Coordinates are offered through Airfoils,
Inc.; they are not tabulated in this report. This is design-and-test evidence,
not a validated automatic-optimization algorithm.

**Executed follow-up:** the report's vector Figure 1 contains recoverable
contours for all four configurations. The [S702 benchmark](S702_BENCHMARK.md)
now uses those plotted vertices, together with vector-extracted MSES curves
and experimental marker centers, to compare drag-bucket movement and drag at
common lift. See the [executed results](results/s702_comparison/S702_RESULTS.md).
This is a qualified figure reconstruction; exact author coordinates remain
unavailable. The small flap is on the aft element: the benchmark preserves
that distinction and uses the published free-transition conditions. Close
agreement with MSES is insufficient when MSES and experiment favor different
configurations. The report's forced-transition cases are still outside the
current solver's supported inputs.

In the executed S702 matrix, MultiFoil's mean absolute drag difference from
MSES is 1.46% across 12 fixed-lift targets; the difference from experiment is 8.84%.
At CL=0.8, both numerical methods nominally prefer a flap deflection of 5
degrees to 10 degrees, while experiment favors 10 degrees. This reversal
persists at both tested
MultiFoil panel resolutions. The case supports further development of the
analysis method and supplies a concrete optimization-ranking challenge; it
does not establish reliable contour-shape gradients.

### 2. Pomeroy, Williamson and Selig (2012): measured geometry rankings

Brent W. Pomeroy, Gregory A. Williamson and Michael S. Selig, *Experimental
Study of a Multielement Airfoil for Large Wind Turbines*, AIAA2012-2892.
[Author-hosted paper](https://m-selig.ae.illinois.edu/pubs/PomeroyWilliamsonSelig-2012-AIAA-2012-2892.pdf).
[DOI](https://doi.org/10.2514/6.2012-2892).

A four-element system was tested at Re=0.975million while varying gaps,
overhangs and strut settings. Tables1–2 describe the baseline and test matrix;
the figures show resulting lift, drag and efficiency trends. This supplies
experimental design-space information that a single polar cannot provide.
The study's gap findings are specific to its configuration and must not become
a universal minimum-gap rule. Geometry and chord conventions need careful
reconstruction from the cited design source and the test setup.

**Proposed use:** compare the sign and location of measured performance changes
in one-dimensional sweeps first, then selected combined changes. This tests
element placement; it does not alone validate contour-shape derivatives.

### 3. Viken et al. (2017): MSES design with independent RANS assessment

Jeffrey K. Viken, Sally A. Viken, Karen A. Deere and Melissa B. Carter,
*Design of the Cruise and Flap Airfoil for the X-57 Maxwell Distributed Electric
Propulsion Aircraft*, AIAA Aviation2017.
[NASA-hosted paper](https://ntrs.nasa.gov/api/citations/20180003474/downloads/20180003474.pdf).
The2018 NASA accession number is not the conference year.

MSES supported cruise-airfoil and slotted-flap design; USM3D assessed alternative
flap positions using multiple turbulence models. Tables2–4 provide useful
geometry and coefficient comparisons. The paper explicitly reports no planned
developmental wind-tunnel optimization of the flap, and MSES failed to converge
for the40-degree flap cases attempted. This supports using a fast design model
with independent checks, not assuming MSES and RANS agree everywhere. The
powered-wing FUN3D studies are outside our unpowered2D solver's scope.

**Proposed use:** the unpowered flap-position alternatives can test ranking
against RANS once geometry and transition assumptions match. The published
fully turbulent calculations cannot be reproduced merely by selecting an
arbitrary low Ncrit in our natural-transition model.

### 4. Drela (1996): direct instructions for formulating the optimization

Mark Drela, *A User's Guide to LINDOP V2.50*, June1996.
[Official MIT manual](https://web.mit.edu/drela/Public/web/mses/lindop.pdf).
[MSES overview](https://web.mit.edu/drela/Public/web/mses/).

LINDOP uses MSES sensitivities for single- or multipoint optimization with
geometry modes, element-position modes, objectives and constraints. Section14
discusses off-design degradation, excessive shape freedom and optimization that
exploits fragile transition or separation-bubble behavior. This is a direct
precedent for the intended workflow, not experimental validation. Its advice on
transition should be read with the later transition-aware design literature.

**Proposed use:** impose common shape parameters across operating points, solve
incidence to meet lift requirements, constrain thickness and moment, and limit
fine geometric features. Include robustness to operating-condition changes.

### 5. Fidkowski (2022): the foundation actually used by this code

Krzysztof J. Fidkowski, *A Coupled Inviscid–Viscous Airfoil Analysis Solver,
Revisited*, AIAA Journal60(5),2961–2971,2022.
[Author's paper](https://websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf).
[Code and capabilities](https://websites.umich.edu/~kfid/codes.html).

The paper documents mfoil's panel/boundary-layer equations, coupling,
transition treatment and comparisons with XFOIL. Its smooth model functions
are relevant to optimization. MultiFoil reuses the2023 Python implementation,
but adds its own multi-body interactions and wakes. This paper therefore
supports the inherited single-element formulation, not the correctness of our
extension or its shape sensitivities. A Newton state Jacobian is not by itself
a derivative of converged drag with respect to geometry.

**Proposed use:** verify single-element shape perturbations against upstream
mfoil/XFOIL before interpreting disagreement caused by multi-element coupling.
The present MATLAB release and our bundled Python version are different
versions; improvements should not be assumed to transfer automatically.

### 6. Pomeroy and Selig (2017): contour redesign and wake limitations

Brent W. Pomeroy and Michael S. Selig, *Design of Airfoils to Mitigate Wake
Bursting*, AIAA2017-1210.
[Author-hosted paper](https://m-selig.ae.illinois.edu/pubs/PomeroySelig-2017-AIAA-2017-1210-Design-for-WakeBursting.pdf).

Three redesigned multi-element airfoils were assessed against a baseline;
USM3D predicted improvements in lift-to-drag ratio up to17%. MSES provided
transition locations for that assessment, while MFOIL/PROFOIL supported inverse
design. The reported gains are RANS predictions, not wind-tunnel confirmation
of those new geometries. The paper's MFOIL is a different code from Fidkowski's
mfoil used here. It is particularly relevant to the physics an optimizer could
exploit when MultiFoil retains independent non-merging wakes.

**Proposed use:** a later contour-ranking challenge against the published RANS
design alternatives, after exact geometries and transition settings are secured.
Do not attribute the17% result to this MultiFoil implementation or to MSES alone.

### 7. Halila, Martins and Fidkowski (2020): transition-aware optimization

Gustavo L. O. Halila, Joaquim R. R. A. Martins and Krzysztof J. Fidkowski,
*Adjoint-Based Aerodynamic Shape Optimization Including Transition to
Turbulence Effects*, Aerospace Science and Technology107,106243,2020.
[Author's preprint](https://websites.umich.edu/~kfid/MYPUBS/Halila_Martins_Fidkowski_2020_AST.pdf).
[DOI](https://doi.org/10.1016/j.ast.2020.106243).

This RANS study includes transition in the gradient-based design process and
examines the consequences of losing laminar flow through contamination. It
supports assessing both nominal laminar performance and degraded transition
conditions. It does not validate MultiFoil or its simplified e^N closures.

**Proposed use:** include Ncrit/Re sensitivity now; add and verify a true
forced-transition/fully-turbulent option before claiming roughness robustness.
For independent CFD assessment, distinguish transition-matched RANS from a
deliberately fully turbulent stress test.

### 8. UIUC Low-Speed Airfoil Tests: accessible held-out shapes

[UIUC experimental programme, reports and tabulated datasets](https://m-selig.ae.illinois.edu/uiuc_lsat.html).

The programme provides experimental polars for many different low-Reynolds-number
airfoils, with downloadable reports and tabulated data. It offers a practical
source for checking whether a method ranks different single-element geometries
correctly under comparable conditions. It does not automatically establish
that every design was generated using MSES or XFOIL, nor that unrelated shapes
validate local design gradients.

**Proposed use:** select a small set matching the intended Reynolds-number and
lift range, reserve several shapes as holdouts, and compare drag at common lift.
Audit coordinates, trips, tunnel conditions and measurement uncertainty first.

## Additional finding: an official NLR7301 MSES drag reference

[Drela's official NLR7301 example](https://web.mit.edu/drela/Public/web/mses/nlf26.pdf)
shows MSES3.11, Re=2.51million, Mach0.185, Ncrit9 and alpha6degrees, with
CL=2.4840 and CD=0.02264 printed on page1. The other pages show boundary-layer
and wake quantities for that same state. This is one operating point, not a
new drag polar. Its title also includes a configuration suffix that must be
reconciled with the actual coordinates.

Our archived6-degree solution has CL=2.562924 and CD=0.02323949. The nominal
CD difference is+0.00059949, approximately6drag counts or2.65%. This is an
arithmetic comparison of printed and archived values, not an exact-input-deck
validation. No new MSES or MultiFoil run was performed. Exact geometry and
settings must be checked before incorporating this as a scored benchmark.
The follow-on pilot now adds the printed point as qualified plot context and
compares the baseline at fixed lift as well as fixed incidence; it remains unscored.
The earlier three-way figure lacked an NLR MSES drag series because the source
used for that plot supplied only lift; this newly found source closes part of
that evidence gap.

## Proposed acceptance programme for our implementation

The recommendations below are engineering proposals for this project, not
acceptance thresholds claimed by the papers.

1. **Define the design domain.** Specify element count, Reynolds-number range,
   required lift, Mach number and allowed geometry changes. Start in a
   subcritical regime with modest separation. Keep common reference chord and
   force normalization explicit when changing flap placement.
2. **Separate contour and placement tests.** Check small smooth contour changes
   against a single-element reference first. Then test a measured multi-element
   placement sweep. Passing one does not establish the other.
3. **Measure performance differences.** For each geometry g relative to baseline
   g0, compare DeltaCD=CD(g)-CD(g0) at the same required CL and conditions.
   Report bias, error in DeltaCD, pairwise ranking and the selected optimum.
   Do not treat changes smaller than numerical or reference uncertainty as
   resolved improvements. A plot-reading allowance is not the full uncertainty.
4. **Check numerical stability.** Repeat baseline and candidate results on finer
   panel distributions and longer wakes. Check directional differences at more
   than one perturbation size. Require stable signs and rankings. Different
   initial guesses and continuation paths must not silently select different
   branches or manufacture apparent optima.
5. **Use independent holdouts.** Reserve geometry variants and operating points
   that did not influence implementation changes. An optimizer repeatedly using
   the validation set can indirectly fit it even without explicit coefficients.
6. **Verify shortlisted improvements.** Evaluate baseline and finalists with the
   same independently checked RANS/transition setup or with experiments. Check
   constraints there too, not just the optimized objective. RANS discretization
   and model uncertainty still need assessment.

A suitable initial objective is weighted drag at several prescribed lift and
Reynolds-number conditions. Geometry variables should describe smooth,
low-dimensional changes. Thickness, curvature, minimum clearance and moment
constraints should be explicit. Failed solves and flagged out-of-domain states
should be infeasible, never assigned attractive artificial drag values.

For gradient-based optimization of R(u,p)=0, the needed derivative includes the
state response: du/dp=-(dR/du)^(-1)(dR/dp). Our existing Newton solver supplies
state linearizations, but geometry, panel redistribution, wakes and transition
must all be accounted for and verified. No shape adjoint or validated geometry
gradient is claimed by the current package.

The most valuable next step is a small ranking benchmark using exact geometry
variants, before a large optimization campaign. The current51% NLR drag
overprediction near13degrees and the missed experimental stall changes remain
unresolved evidence; the new literature does not erase those discrepancies.

## Historical literature-scan scope and subsequent work

Read: official LINDOP manual; mfoil paper and author documentation; S702 and X-57
NASA papers; Pomeroy2012 and2017 papers; Halila et al. preprint; UIUC dataset
catalogue; official MSES NLR plots. Additional literature PDFs used during the
scan are not all redistributed. Available source reports and extraction audits
for the implemented benchmarks are included under `data/`; the S702, TN808 and
TR664 source reports are among them.

The literature-scan phase itself performed no aerodynamic solves. Subsequent
work implemented and executed the numerical pilot, TN808 and S702 figure
reconstruction/comparison, an incompressible contour search, discrepancy
diagnostics and the transition calibration/transfer study summarized above.
Independent contour-ranking validation, RANS runs and new physical experiments
remain unperformed. Exact S702 author coordinates remain unavailable; figure
reconstruction does not remove that uncertainty for fine contour optimization.

## Subsequent independent-source audit

The [TR664 audit and ranking comparator](INDEPENDENT_RANKING.md) now archive a
primary experimental source with coordinate tables and placement effects.
The candidate remains unscored: full assembled geometry and matched settings
are not yet verified. Its published drag precision cannot independently resolve
the approximately0.7-count camber benefit in the nominal NLR7301 pilot.

## Executed TN808 experimental challenge

The [TN808 reconstruction benchmark](TN808_BENCHMARK.md) now compares individual
flap settings at CL1.5, with results in
[the executed report](results/tn808_comparison/TN808_RESULTS.md).
It exposed and fixed a nonmonotone panel-spacing defect. The experimental
comparison still shows a large drag discrepancy under the stated geometry and
transition assumptions. Numerical convergence therefore does not justify
accepting drag-optimal shapes from this case. The report retains all failed
attempts and successful initialization retries.

The model's larger flap-angle drag ordering can be checked against these data;
the smaller10–20degree difference is unresolved within the declared experimental
comparison allowance. This does not validate the earlier NLR7301 camber
gradient. The cove case also reveals numerical path sensitivity of about1 drag
count at near-identical incidence, which is not a universal noise bound for
other configurations.
