# Independent design-ranking benchmark: source audit and comparison code

A subsequent [TN808 challenge](TN808_BENCHMARK.md) now provides a reconstructed,
executed fixed-lift comparison across flap angles. See
[its numerical results](results/tn808_comparison/TN808_RESULTS.md).
The original TR664 candidate below remains unscored; it was not silently filled
with the different TN808 geometry or data.

8 September 2026. **Status: implemented comparison checks; independent geometry
rankings have not yet been validated.** No new aerodynamic solutions or numerical
experimental polar samples are claimed in this addition. The archived NLR7301
pilot remains an internal sensitivity study.

## What was added

`compare_design_rankings.py` compares drag differences between geometries at a
common lift, keeping experiments, MSES and RANS distinct. It reports agreement,
disagreement or an unresolved ranking after accounting for supplied absolute
bounds. Missing geometry/condition/bound audits, failed solutions, a local
supersonic flag, missing numbers or an unmatched lift block the comparison.
The command also checks the SHA256 of each supplied geometry file.

The new source pack, `data/naca23012_tr664_ranking/`, contains the original NACA
report, eight inspected page images, a source record, machine-readable table ordinates, published precision limits
and an explicitly blocked candidate definition. The candidate is excluded from
the existing default polar benchmarks. Its current executable outcome is in
`results/independent_ranking/status.json`; an empty pair list means no comparison,
not successful validation.

## Recovered experimental evidence

Wenzinger and Harris (1939), *Wind-Tunnel Investigation of an NACA 23012 Airfoil
with Various Arrangements of Slotted Flaps*, NACA-TR-664.
[NASA primary report](https://ntrs.nasa.gov/api/citations/19930091739/downloads/19930091739.pdf).
This public-domain US government report includes main-airfoil and flap
ordinates, slot geometry diagrams, and measured effects of flap placement.

| Evidence | Location in report | Consequence for the benchmark |
|---|---|---|
| Main NACA23012 and slot ordinates; flap2 ordinates | Tables I–II, printed pp.30–31; PDF34–35 | Recoverable source coordinates, but not a ready-to-run assembled two-element coordinate file |
| Slot2-h radii and slot2-i entry rounding | Figure5, printed p.4; PDF8 | Cove reconstruction and connection to the original surface must be checked |
| Physical Re2.19 million; historical effective Re approximately3.5 million | Printed p.5; PDF9 | Use physical Re; do not substitute effective Re or claim a known Ncrit |
| Drag contours at CL1.5,2.0,2.5 for flap2-h and2-i | Figures25–26, printed p.21; PDF25, interpretation pp.19–20 | Useful design-space evidence; flap-angle selection also has to be recovered before assigning a contour value to one geometry |
| Flap-angle envelope comparison | Figure31, printed p.24; PDF28 | This is an envelope over settings, not the polar of one fixed geometry |
| Individual flap placements and section polars in a second tunnel | Figures32 and35, printed pp.25,27; PDF29,31 | These belong to the variable-density campaign at effective Re≈8.2 million; cannot borrow the first campaign's physical Re |

The report's Figure31 discussion places the slotted2-h envelope below the other
flap types at moderate lift around CL1.0–1.7. This is qualitative experimental
support for the value of slot/flap design. It does not establish the derivative
of our NLR7301 camber perturbation, and we have not scored it as a MultiFoil result.

### What the drag accuracy permits

The authors give the following estimated limits on accidental experimental
errors from repeat tests (printed p.5). These are historical repeatability
estimates, not modern total-uncertainty budgets or statistical confidence bounds.
One drag count is CD=0.0001.

| Operating point | Reported per-configuration CD limit | Sum for a two-configuration difference, absent cancellation evidence |
|---|---:|---:|
| CL0 | ±3 counts | ±6 counts |
| CL1 | ±6 counts | ±12 counts |
| CL2.5 | ±20 counts | ±40 counts |

The sums are conservative arithmetic of the reported limits, not experimentally
established uncertainty of paired measurements. Digitization, geometry
reconstruction, lift matching and systematic effects would also need treatment.
Do not interpolate this three-point precision table into a validated uncertainty
model without justification.

The measured profile drag includes hinge fittings. With flaps neutral and
openings sealed, the additional CD was about0.001 (10 counts). The report does
not measure that contribution for deflected flaps. The authors expected relative
merits to be less affected because fittings were shared, but did not quantify
the residual bias. **Subtracting a universal10-count correction is unsupported.**

For scale only, our nominal NLR7301 positive-camber gain was about0.7 counts, on
a different geometry and at different conditions. TR664 cannot establish that
small improvement from its stated precision alone. Larger placement effects
could still be useful challenges after reconstructing the inputs.

### Issues found during reconstruction

The flap2 table is in percent of original wing chord, with a displaced vertical
datum; the nose is at (0,−1.29)% and its nose-arc centre is at (0.91,−1.29)%.
Do not independently normalize the flap to a horizontal nose-to-tail chord and
then apply the figure's placement without transforming the reference point.
Sparse lower-surface coordinates and the small nose radius also require a
controlled contour construction. OCR is insufficient: for example the upper
surface station between9.03 and15.66 is10.31, as checked visually in PDF35.

The full main-element cove and joins have not been reconstructed and verified.
The flap setting associated with every drag-contour sample is not yet established.
Rotating a generic flap or substituting the published envelope would therefore
test an assumed geometry rather than the measured one. No such substitution was
entered into the benchmark.

A second retrieved report, NACA-TR-677 (NASA accession19930091752), has a metadata
title naming NACA23012, whereas the actual document names **NACA23021**. It is
not interchangeable with TR664 and was not added as a matched23012 dataset.

## Running the comparison

From the project root:

```bash
python validation/compare_design_rankings.py \
  validation/data/naca23012_tr664_ranking/candidate.json \
  --out validation/results/independent_ranking/status.json
python -m unittest -v test_design_rankings.py
```

The first command currently exits2 and writes `status: blocked`, because the
candidate's geometry audit is incomplete. This is the expected archived result.
The tests contain synthetic values solely to exercise rejection, interval
arithmetic and ranking classification; they are not aerodynamic validation.

### Input contract for an acquired matched dataset

Use the candidate JSON as the top-level template. `reference_kind` must be
`experiment`, `mses`, or `rans`. Provide `source`, `target_cl`, a small justified
`cl_tolerance`, and written `geometry_audit`, `conditions_audit`, `bounds_audit`.
The three corresponding `*_audited` flags may become true only when those audits
are complete. These are researcher declarations; code cannot certify the
underlying scientific judgments.

Each entry of `geometries` requires a unique `id`, `geometry_file` (relative to
the input JSON, or absolute), `geometry_sha256`, and `reference` and `model`
objects. For multiple elements, use one hashed bundle containing all element
coordinates and placements. Both point objects contain numeric `cl`, positive
`cd`, positive `cd_abs_bound`, and a `provenance` string identifying the source
page/curve or saved model run. The model also requires `converged: true` and
`local_supersonic: false`.

Audit the physical Reynolds number, Mach, chord/force normalization, transition
treatment, geometry and placement for each pair. State approximations explicitly.
The reference bound must include justified digitization and lift-matching effects
as well as applicable measurement limits. A model bound can describe demonstrated
resolution/initialization/transition sensitivity, but must not be presented as
a guaranteed bound on unknown physical-model error.

The comparator does no interpolation, extrapolation or fitting. If source curves
need interpolation to target lift, do that in a separately audited extraction
step with bracketing data, branch choice and an uncertainty allowance. Never
interpolate through a failed model solve or across stall branches. A fixed-lift
run can use the existing `optimization.solve_fixed_lift` with an admissible,
convergent bracket; record any retries and changes in initialization.

For each pair, the code forms ΔCD=CD_B−CD_A. Its bound is the sum of the two
per-geometry absolute bounds. If either the model or reference interval includes
zero, the direction is unresolved. Otherwise equal directions give `agrees` and
opposite directions give `disagrees`. Agreement concerns ordering; it does not
mean that magnitudes agree. MSES agreement remains distinct from agreement with
experiment. Neither certifies an optimum or a contour-shape gradient.

## Remaining work before an independent optimization claim

Complete and inspect assembled geometries and measured flap settings; extract
fixed-configuration observations with a pixel audit; establish the condition and
uncertainty qualifications; then run the paired fixed-lift cases and numerical
sensitivity checks. No MSES executable, original MSES geometry-ranking run or
new experimental ranking has been supplied by this addition. The current
optimization recommendation remains screening with independent verification of
shortlisted designs.

## Executed software checks

23 tests passed: six new design-ranking contract tests, seven existing fixed-lift
optimization tests and ten existing reference-comparison integrity tests. The
new candidate CLI was executed and exited2 as expected, retaining the blocked
status. Aerodynamic solver files were unchanged; no new flow run is claimed.
