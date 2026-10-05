# MultiFoil references and comparison candidates

Search performed 24 September 2026; inventory completed 1 October 2026 after
an interrupted session. This consolidates the references documented in the
project, then ranks additional sources for experimental, published-MSES and
flow-mechanism comparisons. It is a bounded literature search, not a claim to
have found every relevant publication. Multiple figures or registry entries
from one publication are consolidated into one citation.

**The strongest additions are the S103 and S414 MSES/experiment comparisons,
and Pomeroy's three-element wake measurements with tabulated geometry.**
Nakayama's near-wake dataset supplies downloadable measurements for a simpler
single-element diagnostic. These are candidates, not newly passed benchmarks.
No new aerodynamic solves or MSES runs were performed for this review.

Published MSES curves are third-party predictions. The project has neither an
installed MSES executable nor the exact input decks behind most curves.
Agreement with MSES alone is not experimental validation. The frozen
[benchmark registry](sources.json) and existing source manifests remain
authoritative for already scored data; this bibliography does not replace them.

## What to use next

**S414 trip qualification, 3 October 2026:** six 2010 experimental table curves
and six published MSES drag curves are extracted at Re 0.7m/1m. The report's
trip-effectiveness and excess-displacement warnings, and absence of paired
tripped aft Cp, prevent treating the fore-only-trip contrast as a clean wake
intervention. Four lower-boundary free-transition pressure snapshots are also
extracted. See [trip/pressure qualification](results/s414_trip_qualification_20261003/FINDINGS.md).
The next scoped calculation should use this free-transition sequence; no
forced-trip closure calibration is justified. A bounded search additionally
identified Twiss (2019), *High lift configuration of a slotted natural laminar
flow airfoil*, [Illinois thesis PDF](https://www.ideals.illinois.edu/items/112023/bitstreams/366903/data.pdf).
This is a separate experimental campaign, not a recovered 2010 tripped-pressure
dataset. Its geometry, conditions and visualization limitations still need
qualification before use; no measurements from it enter current comparisons.

**S414 update, 2 October 2026:** the original 2010 report (L01) is recovered
from public mirrors. Its modified nominal test contour and two published
envelope-MSES/experiment pressure comparisons are extracted. Exact tested
coordinates remain unavailable; public `.pts` files sharpen finite trailing
edges and have undocumented registration. See the
[geometry/source qualification](results/s414_geometry_qualification_20261002/FINDINGS.md)
for repeated-outline checks, tap inconsistencies and the admission limit.
Keep these 2010 envelope-e^N curves separate from N02's 2014 full-e^N curves.
The matching Figure-13(c) drag/lift curves are now extracted: at CL 0.284 / 1.041,
MSES drag differs from experiment by +5.33 / −13.89 counts. See the
[matching drag evidence](results/s414_drag_reference_20261002/FINDINGS.md).
The [bounded Ncrit-9 reconstruction comparison](S414_RECONSTRUCTION_PROTOCOL.md)
is complete and audited on 3 October: all nine states pass, with an exact repeat.
MultiFoil and published MSES both place the strongest sampled aft pressure rise
one tap interval downstream of experiment, robust to the declared reconstructions.
See the [comparison findings and graphs](results/s414_reconstruction_analysis_20261003/FINDINGS.md).
Exact tested-geometry admission remains open.

**Qualification update, 1 October 2026:** the bounded S702 sampling audit is
complete. Pomeroy's 243 coordinates are extracted and visually checked, but the
2014 absolute placement does not reproduce its relative gap table under the
declared frame. Its 2015 companion has a different placement and sign convention.
See the [geometry and observable audit](results/pomeroy_geometry_qualification_20261001/FINDINGS.md).
Pomeroy remains a useful mechanism source, but cannot yet enter a matched
experimental benchmark. Prioritize S103/S414 geometry qualification for the
three-way comparison. The original literature-review snapshots remain unchanged.

| Priority | Sources | Comparison available | Main qualification before use |
|---|---|---|---|
| 1 | S103, N01 | Experiment versus published MSES: pressure and section characteristics | Reproduce the tested aft-element displacement; exact contour coordinates still need qualification; high-Re test has reported three-dimensional effects |
| 1 | S414, N02 and L01 | Two pressure comparisons and matching 2010 lift/drag curves extracted; nine MultiFoil states audited; 2014 also has OVERFLOW | Aft recovery mismatch survives declared reconstruction sensitivities; exact geometry/edge treatment and tap inconsistencies remain open. 2010 envelope e^N versus 2014 full e^N, both Ncrit 9 |
| 2 | Pomeroy 2014/2015, N03–N04 | Measured three-element wake velocity, off-body pressure and geometry changes; extracted element coordinates | Reconcile absolute/relative placement conflict and paper-specific setups; current model lacks measured wake cross-profile observables; no MSES polar established |
| 2 | Nakayama/TMR, N06 | Downloadable single-element near-wake profiles and coordinates | Match trips and trailing edge; forced transition is not currently a supported MultiFoil input |
| 2 | Adair–Horne, N05; original NLR report, P04 | Boundary-layer/wake interaction and confluence measurements | Complete geometry/condition extraction; original NLR full text remains inaccessible through the tested link |
| 2 | Pomeroy 2012, P16; Coiro 2009, N09 | Measured performance changes with element placement | Recover exact contours and normalization; Coiro institutional full text is restricted |
| 3 | Steiner 2020, N07; Zahle 2012, N08 | Slat design comparisons and thick-airfoil experiments | Source-specific geometry, tunnel and transition audit; do not label every Steiner case a converged MSES comparison |
| 3 | Chin 1993, N10; 30P30N, N11–N12 | High-lift pressures, wake/flowfield measurements | Dataset access and exact configuration audit; unsteady three-dimensional flow may exceed the present model |

For any admitted benchmark, keep original data, separate raw measurements from
digitized curves, record geometry and reference chord, Re, Mach, transition/trips,
tunnel corrections and uncertainty, and compare drag at common lift. Freeze the
protocol before examining a new design's score. Do not retune Ncrit per airfoil
to experimental drag. MSES's transition method/version and our approximate
envelope formulation must also be distinguished.

## References already used in project comparisons or source audits

### P01 — S702 experiment and published MSES

Somers, D. M.; Maughmer, M. D. (2022). *Design and Experimental Results for a
Low-Speed, Flapped, Slotted, Natural-Laminar-Flow Airfoil*.
NASA/CR–20220014964. [NASA report](https://ntrs.nasa.gov/citations/20220014964).

The principal two-element benchmark: reconstructed Figure 1 geometry,
free-transition experimental and MSES polar figures for four flap settings,
reserved Reynolds-number cases, and a separate pressure comparison. Used in
[S702_BENCHMARK.md](S702_BENCHMARK.md), calibration and subsequent diagnostics.
Exact author coordinates and the hinge-film shape remain unresolved. Published
Ncrit 9 MSES and calibrated Ncrit 6.25 MultiFoil are different transition scenarios.

### P02 — NACA 0012 experiment

Ladson, C. L. (1988). *Effects of Independent Variation of Mach and Reynolds
Numbers on the Low-Speed Aerodynamic Characteristics of the NACA 0012 Airfoil
Section*. NASA TM-4074.
[NASA report](https://ntrs.nasa.gov/citations/19880019495).

Table I supplies the Re 5.97 million, Mach 0.15 free-transition force comparison.
The implemented finite-trailing-edge NACA formula is not NASA TMR's modified
closed-trailing-edge geometry. Ncrit 9 is an assumption, not a measured tunnel
condition. See [validation cases](README.md).

### P03 — S809 experiment

Somers, D. M. (1997). *Design and Experimental Results for the S809 Airfoil*.
NREL/SR-440-6918. DOI 10.2172/437668.
[Institutional report](https://digital.library.unt.edu/ark:/67531/metadc681921/).

Table 2 coordinates and experimental figures support Re 2 million lift and
transition comparisons. Mach zero is a low-speed approximation. The source's
transition/reattachment distinction and usable maximum-lift qualification are
retained; massively separated results do not validate stall prediction.

### P04 — Original NLR7301 measurements, metadata only

van den Berg, B.; Oskam, B. (1979). *Boundary Layer Measurements on a
Two-Dimensional Wing with Flap and a Comparison with Calculations*.
NLR MP 79034 U; also AGARD CP-271, paper 18.
[TU Delft record](https://repository.tudelft.nl/record/uuid:5f1bfb3e-71d9-45a8-bb55-c4cdc6f25182),
[NTIS record](https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/N8112036.xhtml).

The experiment underlying our secondary NLR digitizations. The full original
report has not been read. Historical notes describe restricted access; on
24 September the public metadata was accessible but its advertised file URL
returned HTTP 404. That test establishes a failed link, not a current account
requirement. Pressure, boundary-layer and wake measurements make retrieval useful.

### P05 — NLR7301 geometry and experimental digitizations actually used

Ranneberg, M. *viiflow-examples*, commit
`93bef6c6b779db4605f6bcbb2caf1d88a51c9b85`, NLR7301 example.
[Pinned repository](https://github.com/mranneberg/viiflow-examples/tree/93bef6c6b779db4605f6bcbb2caf1d88a51c9b85/NLR7301%202-Element%20Airfoil).

MIT-licensed nominal ECARP contours and digitized experimental lift/drag are in
`data/nlr7301/`. These are not original measurement files. The example's optional
deformation correction is not applied here. Ambiguous pressure files remain
unscored. This data source is distinct from P04 and does not require a Viiflow run.

### P06 — Published NLR7301 MSES lift curve

Neto, J. A. de O.; Junior, C. B.; Cavali, D.; Azevedo, J. L. F. (2006).
*Aerodynamic Study of High-Lift Device Configurations for Take-off and Landing
Conditions*. ICAS 2006, paper 512.
[Conference paper](https://www.icas.org/icas_archive/ICAS2006/PAPERS/512.PDF).

Figure 6 supplies the digitized MSES lift curve for nominal 2.6% gap,
Re 2.51 million and Mach 0.185. Also used for configuration/context checks.
Exact geometry deck, transition settings and grid are unavailable. Registry
entries `neto_2006` and `neto_2006_mses` refer to this same paper.
See [MSES.md](MSES.md).

### P07 — Published NACA 0012 MSES drag curve

Lima e Silva, A. L. F.; Neto, J. A. de O.; Antunes, A. P.; Mendonça, M. T.;
Azevedo, J. L. F.; Silveira Neto, A. (2005). *Numerical Study of Two-Dimensional
High-Lift Configurations Using the MSES Code*. COBEM2005-0277.
[Conference paper](https://abcm.org.br/anais/cobem/2005/PDF/COBEM2005-0277.pdf).

Figure 3 supplies NACA 0012 MSES drag versus incidence at Re 6 million,
Mach 0.3. Its experimental curve is secondary evidence, used visually only.
Figure 10 remains excluded because the configuration/caption is ambiguous.
The 12-count graphical reading allowance is not experimental uncertainty.
See [DRAG.md](DRAG.md) and [THREE_WAY.md](THREE_WAY.md).

### P08 — Official MSES NLR7301 example

Drela, M. *NLR7301, MSES 3.11 example*, undated.
[MIT example PDF](https://web.mit.edu/drela/Public/web/mses/nlf26.pdf).

One printed state: Re 2.51 million, Mach 0.185, Ncrit 9, incidence 6 degrees,
CL 2.4840 and CD 0.02264, with boundary-layer/wake plots. Qualified context,
not a scored polar or an exact-input reproduction. Geometry suffix and moment
reference need reconciliation. See [OPTIMIZATION_READINESS.md](OPTIMIZATION_READINESS.md).

### P09 — NACA 23012 slotted-flap source audit

Wenzinger, C. J.; Harris, T. A. (1939). *Wind-Tunnel Investigation of an NACA
23012 Airfoil with Various Arrangements of Slotted Flaps*. NACA TR-664.
[NASA report](https://ntrs.nasa.gov/citations/19930091739).

Audited as an experimental design-ranking candidate, not scored as a new
validation case. The physical Re 2.19 million must not be replaced by the old
"effective" Re 3.5 million. Some reported envelope polars change geometry;
flap-fitting drag also matters. See [INDEPENDENT_RANKING.md](INDEPENDENT_RANKING.md).

### P10 — NACA 23012 extended-lip reconstruction benchmark

Lowry, J. G. (1941). *Wind-Tunnel Investigation of an NACA 23012 Airfoil with
Several Arrangements of Slotted Flaps with Extended Lips*. NACA TN-808.
[NASA report](https://ntrs.nasa.gov/citations/19930085051).

Table I and Figure 8 underpin the executed, qualified reconstruction comparison.
Uncertain cove/lip reconstruction and transition prevent treating it as exact
geometry validation. Physical Re 2.19 million; Mach 0.105 is an assumption.
See [TN808_BENCHMARK.md](TN808_BENCHMARK.md).

## Method, design and discovery sources already consulted

### P11 — Inherited mfoil formulation

Fidkowski, K. J. (2022). *A Coupled Inviscid–Viscous Airfoil Analysis Solver,
Revisited*. AIAA Journal 60(5), 2961–2971.
[Author paper](https://websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf).

Foundation for the inherited single-element equations, coupling and transition
model. It does not validate MultiFoil's multi-body extension or shape derivatives.

### P12 — mfoil implementation and provenance

Fidkowski, K. J. *mfoil code and documentation*; bundled Python version dated
28 June 2023. [Author page](https://websites.umich.edu/~kfid/codes.html),
[author Python link](https://websites.umich.edu/~kfid/codes/mfoil/mfoil_notex.py),
[actual pinned mirror](https://github.com/SharathSPhD/CINS/blob/514579fb4c65467d26c493c6712490f70ded2afe/vendor/mfoil/mfoil.py).

The bundled source and license are recorded in [THIRD_PARTY.md](../THIRD_PARTY.md).
Current upstream MATLAB and bundled historical Python are distinct versions.

The 3 October 2026 [fore-nose diagnostic](results/s414_fore_nose_20261003/FINDINGS.md)
finds a missing factor of two in the bundled `get_damp` amplification derivative.
The same value/derivative expressions were inspected in the author's
[28 June 2023 Python source](https://public.websites.umich.edu/~kfid/codes/mfoil/mfoil.py);
see the [expression-level receipt](results/s414_fore_nose_20261003/sources/author_source_check.json).
This is not a whole-file comparison or a claim about current MATLAB versions.
The defect's effect on converged aerodynamic predictions remains untested.

### P13 — MSES-based optimization workflow

Drela, M. (1996). *A User's Guide to LINDOP V2.50*.
[MIT manual](https://web.mit.edu/drela/Public/web/mses/lindop.pdf).

Used for multipoint, fixed-lift optimization, geometry/position modes, constraints
and robustness considerations. Method guidance, not experimental validation.

### P14 — Official MSES overview

Drela, M. *MSES overview, documentation and examples*, undated web resource.
[MIT MSES site](https://web.mit.edu/drela/Public/web/mses/).

Source for the MSET/MSES/MPLOT/MPOLAR/LINDOP tool family and public documentation.
Public manuals and plots do not provide the licensed solver executable or an
authorization to claim locally run MSES results.

### P15 — XFOIL operating and resolution guidance

Drela, M.; Youngren, H. *XFOIL documentation*.
[Official text](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt).

Used in the wake/refinement audit for boundary-layer resolution diagnostics.
Local Hk and velocity-step guidance is not a guaranteed drag-error bound.
The consulted copy is preserved in the wake-extent analysis source manifest.

### P16 — Measured multielement placement effects

Pomeroy, B. W.; Williamson, G. A.; Selig, M. S. (2012). *Experimental Study of
a Multielement Airfoil for Large Wind Turbines*. AIAA 2012-2892.
[Author paper](https://m-selig.ae.illinois.edu/pubs/PomeroyWilliamsonSelig-2012-AIAA-2012-2892.pdf).

Already reviewed as a prospective ranking benchmark; not yet run in MultiFoil.
Four-element MFFS-026 tests at baseline Re 0.975 million vary gaps, overhangs
and strut settings. The system/reference chord changes with placement, so
normalization needs explicit handling. See N03 for a related three-element
geometry source, without assuming identical rigging or presence of the strut.

### P17 — X-57 MSES design and RANS assessment

Viken, J. K.; Viken, S. A.; Deere, K. A.; Carter, M. B. (2017). *Design of the
Cruise and Flap Airfoil for the X-57 Maxwell Distributed Electric Propulsion
Aircraft*. AIAA Aviation 2017.
[NASA-hosted paper](https://ntrs.nasa.gov/api/citations/20180003474/downloads/20180003474.pdf).

Reviewed MSES design and USM3D flap-position comparisons, not measured contour
rankings. The NASA accession year 2018 is not the publication year. Powered
wing results exceed the present unpowered two-dimensional scope.

### P18 — Wake-aware contour redesign

Pomeroy, B. W.; Selig, M. S. (2017). *Design of Airfoils to Mitigate Wake
Bursting*. AIAA 2017-1210.
[Author paper](https://m-selig.ae.illinois.edu/pubs/PomeroySelig-2017-AIAA-2017-1210-Design-for-WakeBursting.pdf).

Reviewed contour alternatives assessed with USM3D and MSES transition inputs.
Reported design gains are computational, not experimental validation of the
new contours. Its MFOIL/PROFOIL tools are not Fidkowski's mfoil used here.

### P19 — Transition-aware aerodynamic optimization

Halila, G. L. O.; Martins, J. R. R. A.; Fidkowski, K. J. (2020).
*Adjoint-Based Aerodynamic Shape Optimization Including Transition to Turbulence
Effects*. Aerospace Science and Technology 107, 106243.
[Author paper](https://websites.umich.edu/~kfid/MYPUBS/Halila_Martins_Fidkowski_2020_AST.pdf),
[DOI](https://doi.org/10.1016/j.ast.2020.106243).

Supports transition-aware design and sensitivity scenarios. This RANS study does
not validate our closures. A low Ncrit is not a substitute for an implemented
and verified forced-transition option.

### P20 — UIUC low-speed experimental catalog

Selig and collaborators. *UIUC Low-Speed Airfoil Tests*.
[Program and data](https://m-selig.ae.illinois.edu/uiuc_lsat.html).

Previously reviewed source for future single-element geometry holdouts.
Individual volumes, coordinates and test conditions need selection and audit
before scoring; the whole catalog is not an executed benchmark.

### P21 — Airfoils Incorporated source search

Airfoils, Incorporated. [Website](https://airfoils.com/) and
[publications catalog](https://airfoils.com/publications/).

Consulted in the S702 original-coordinate search. These are discovery sources,
not measured aerodynamic datasets. Saved pages and access results are in
`results/s702_conditions_audit_20260922/`.

### P22 — UIUC coordinate catalog

Selig and collaborators. *UIUC Airfoil Coordinates Database*.
[Coordinate catalog](https://m-selig.ae.illinois.edu/ads/coord_database.html).

Consulted for geometry discovery; a matching airfoil name alone does not prove
identity with the tunnel model or published MSES geometry.

### P23 — S702 journal follow-up, metadata only

Maughmer, M. D.; Metkowski, L. P.; Axten, C. J.; Somers, D. M. (2025, according
to the saved institutional record). *Exploration of a Low-Speed, Flapped,
Slotted, Natural-Laminar-Flow Airfoil*. DOI 10.2514/1.C037802.
[Penn State record](https://pure.psu.edu/en/publications/exploration-of-a-low-speed-flapped-slotted-natural-laminar-flow-a-2/),
[publisher](https://doi.org/10.2514/1.C037802).

The institutional metadata was used in the conditions audit; publisher full
text was not obtained. This is related to the S702 campaign, not an independently
acquired experimental dataset. Preserve the publication-date qualification.

### P24 — Secondary experimental attribution, original not consulted

Seifert and LaTunia (1999), as credited in Figure 3 of P07.
[Source actually consulted: COBEM 2005](https://abcm.org.br/anais/cobem/2005/PDF/COBEM2005-0277.pdf).

Included to expose the complete provenance of the green NACA 0012 curve in our
three-way plot. The original publication and raw records have not been
identified/read here; no original title, DOI or independent uncertainty is inferred.

## Newly found experimental and published-solver candidates

### N01 — S103, experiment and MSES

Somers, D. M. (2012). *An Exploratory Investigation of a Slotted,
Natural-Laminar-Flow Airfoil*. NASA/CR–2012-217560.
[NASA report](https://ntrs.nasa.gov/citations/20120006123).

Full report text inspected and saved. Two elements; low-speed tests at Mach 0.1,
nominal Re 3, 6 and 9 million, free transition; MSES Ncrit 9. Figures 10–11
compare pressure and section characteristics. The tested aft element was rotated
−1 degree and translated relative to the design, with a different vertical
placement at the highest Re. The report identifies three-dimensional effects
there. Use the lower-Re comparisons first. Exact contour availability is not
yet established; a planform picture is not a coordinate deck.

### N02 — S414, experiment, MSES and OVERFLOW

Coder, J. G.; Maughmer, M. D.; Somers, D. M. (2014). *Theoretical and
Experimental Results for the S414, Slotted, Natural-Laminar-Flow Airfoil*.
Journal of Aircraft 51(6), 1883–1890. DOI 10.2514/1.C032566.
[Institutional record](https://pure.psu.edu/en/publications/theoretical-and-experimental-results-for-the-s414-slotted-natural-2),
[university-hosted article](https://www.ae.metu.edu.tr/tuncer/ae546/multi/s414-1.pdf).

Full article saved and inspected. Low-speed experiments span Re 0.5–1.5 million;
the paper compares lift, drag and pressure with MSES and OVERFLOW. MSES uses
full e^N at Ncrit 9; its earlier report uses an approximate envelope method.
The test configuration was modified for low speed, so the rotorcraft design
Mach numbers are not the test conditions. Exact coordinates need qualification.
Reported drag uncertainty is 0.5 count in the low-drag range, rising toward stall.
The 2 October qualification recovers the earlier report and distinguishes its
envelope-method comparisons; it does not add independent experimental data.

### N03 — Three-element wake measurements with coordinate tables

Pomeroy, B. W.; Diebold, J. M.; Ansell, P. J.; Selig, M. S. (2014).
*Study of Burst Wakes in a Multi-Element Airfoil Flowfield*. AIAA Journal
52(4), 821–831. DOI 10.2514/1.J052680.
[Author paper](https://m-selig.ae.illinois.edu/pubs/PomeroyDieboldAnselSelig-2014-JofAC-BurstWakes-FINAL.pdf).

Full article saved and inspected. Split-film velocity and turbulence measurements
at Re 1 million diagnose wake bursting over the three-element MFFS(ns)-026.
Appendix Tables A1–A3 give each element's coordinates; Tables 1–2 give assembled
positions, chord ratios, deflections, gaps and overhangs. This makes a reproducible
geometry plausible without waiting for an author deck. Measurements still need
extraction and condition auditing. It is evidence for a mechanism on this airfoil,
not proof that wake bursting causes the S702 error. Not a newly acquired MSES polar.

### N04 — Wake pressure measurements and placement changes

Pomeroy, B. W.; Selig, M. S. (2015). *Pressure Measurements of Burst Wakes
Over a Three-Element Airfoil*. AIAA 2015-2569.
[Author paper](https://m-selig.ae.illinois.edu/pubs/PomeroySelig-2015-AIAA-2015-2569-ThreeElementAirfoilWakes.pdf).

Full paper saved and inspected. Seven-hole-probe pressure/velocity measurements
at Re 1 million examine incidence, gaps, overhang and flap angles. Coordinates
are explicitly referred to N03; this paper supplies its own configuration tables.
Useful for comparing the direction of wake-thickness and momentum-deficit changes.
Audit reference chord and rigging rather than merging the two papers' setups.

### N05 — Measured confluent boundary layer and downstream wake

Adair, D.; Horne, W. C. (1987). *Characteristics of a Separating Confluent
Boundary Layer and the Downstream Wake*. NASA TM-100046.
[NASA report](https://ntrs.nasa.gov/citations/19880004941).

Full OCR text saved and inspected; publication is December 1987, despite the
1988 accession identifier. Two-element flapped-airfoil pressure, mean-flow and
Reynolds-stress measurements with laser velocimetry/hot-wire methods make it
relevant to wake/boundary-layer merging. Mach 0.09 is stated. Audit the scanned
Reynolds number and geometry against PDF images before entering numeric inputs;
the downloaded OCR corrupts the exponent.

### N06 — Downloadable near-wake experiment and coordinates

Nakayama, A. (1985). *Characteristics of the Flow around Conventional and
Supercritical Airfoils*. Journal of Fluid Mechanics 160, 155–179.
DOI 10.1017/S0022112085003433. Dataset maintained by the Turbulence Modeling
Resource: [2D airfoil near-wake validation](https://tmbwg.github.io/turbmodels/airfoilwake_val.html).

The official dataset page and three numerical files are saved: original
experimental coordinates, a copy with five apparent bad points removed (still
blunt at the trailing edge), and measured wake profiles. DSMA661 Model A,
incidence zero, Re 1.2 million, with trips at 0.16c upper and 0.05c lower.
Profiles span x/c 1.01–2.19; their transverse coordinate is in inches, with a
24-inch chord. This is a single-element diagnostic. The page's separate CFD
geometry has a sharp trailing edge: neither that geometry nor fully turbulent
CFD should be silently interchanged with the tripped experiment.

### N07 — MSES slat design study with experimental benchmark context

Steiner, J.; Viré, A.; Benetti, F.; Timmer, N.; Dwight, R. (2020).
*Parametric Slat Design Study for Thick-Base Airfoils at High Reynolds Numbers*.
Wind Energy Science 5, 1075–1095.
[Open article](https://wes.copernicus.org/articles/5/1075/2020/).

Full HTML saved and inspected. MSES/OpenFOAM slat design study, including a
DU97-W-300/NACA-22 experimental benchmark and an NHLP-90-L1T2 CFD validation.
**MSES failed to converge for NHLP-90-L1T2:** that case is not a successful
three-way comparison. Clean/rough transition scenarios must remain distinct.
The article says code/data are available on request; no corresponding tabular
archive was acquired and no author was contacted.

### N08 — Thick multielement airfoil design and wind-tunnel test

Zahle, F.; Gaunaa, M.; Sørensen, N. N.; Bak, C. (2012). *Design and Wind
Tunnel Testing of a Thick, Multi-Element High-Lift Airfoil*. EWEA 2012.
[DTU paper](https://backend.orbit.dtu.dk/ws/files/7946702/Design_and_Wind_Tunnel_Testing.pdf).

Public institutional paper found: a 40%-thick main element with a slat provides
a more distinct geometry family than S702. Includes measured performance and
RANS assessment. Wall interference and exact contour/placement data need audit.
Do not label the RANS curves as MSES or treat the paper as an executed benchmark.

### N09 — STOL airfoil, MSES and experimental placement optimization

Coiro, D.; Nicolosi, F.; Grasso, F. (2009). *Design and Testing of Multi-Element
Airfoil for Short-Takeoff-and-Landing Ultralight Aircraft*. Journal of Aircraft
46(5), 1795–1807, per the institutional record. DOI 10.2514/1.43429.
[University record](https://www.iris.unina.it/handle/11588/356559).

Institutional abstract explicitly reports XFOIL/MSES/TBVOR analyses, measured
surface/wake pressures and experimental slat/flap-position optimization.
An attractive independent geometry-ranking candidate. Full text has not been
read: the institutional file is explicitly marked private/restricted. Exact
coordinates and conditions therefore remain unqualified.

### N10 — High-Reynolds-number multi-element flowfield measurements

Chin, V. D.; Peters, D. W.; Spaid, F. W.; McGhee, R. J. (1993).
*Flowfield Measurements About a Multi-Element Airfoil at High Reynolds Numbers*.
AIAA 93-3137.
[NASA record](https://ntrs.nasa.gov/citations/19930064303).

Primary metadata inspected and saved. Surface pressures, wake-rake drag and
flowfield surveys across Reynolds-number/flap-gap changes are relevant.
The NASA record provides no downloadable full report; metadata is not a
digitized validation dataset. Recover full text and configuration first.

### N11 — Public 30P30N measurement presentation

Pascioni, K. A.; Cattafesta, L. N.; Choudhari, M. M. (2014).
*Towards a Measurement Database for the Aeroacoustics of the 30P30N Multi-Element
High-Lift Airfoil*. NASA-hosted presentation, accession 20200007770.
[NASA source](https://ntrs.nasa.gov/citations/20200007770).

Presentation text saved and inspected. Useful pressure/flowfield and experimental
setup information for the slat-cove/wake problem. It is a presentation, not the
identically titled full paper in N12. The upload year 2020 is not the experiment's
publication year. It does not by itself establish a matched drag polar.

### N12 — Related 30P30N experimental paper

Pascioni, K. A.; Cattafesta, L. N.; Choudhari, M. M. (2014).
*An Experimental Investigation of the 30P30N Multi-Element High-Lift Airfoil*.
AIAA 2014-3062. [DOI](https://doi.org/10.2514/6.2014-3062).

Bibliographic lead for the fuller experimental account. The NASA file downloaded
under N11 has a different title and format; do not claim it supplies the full
paper or raw data for this entry. Exact geometry, boundary conditions and usable
steady observables must be qualified before benchmarking.

### N13 — A simpler adverse-pressure-gradient wake experiment

Hoffenberg, R.; Sullivan, J. P.; Schneider, S. P. *Wake Measurements in a
Strong Adverse Pressure Gradient*. NASA CR-197272.
[NASA record](https://ntrs.nasa.gov/citations/19950014614).

Public report located; primary metadata describes a flat-plate wake in a diffuser
with measurements through flow reversal. Useful for isolating wake closure
behavior from airfoil geometry. Bibliographic dates conflict: NASA gives 1994,
later papers cite 1995. Verify the report cover before selecting a definitive
year or equating it with a conference version. No new data extracted.

### N14 — Recent S207 experimental context, outside direct 2D scope

Coder, J. G. (2026). *Transonic Wind-Tunnel Testing of a Slotted,
Natural-Laminar-Flow Wing at Full-Scale Conditions*. CEAS Aeronautical Journal;
published 10 March 2026.
[Open article](https://link.springer.com/article/10.1007/s13272-026-00956-0).

Primary article inspected. S207 semi-span wing tests provide forces, pressure
and infrared transition measurements at flight-relevant conditions; low-speed
flap testing is also discussed. This is swept three-dimensional/transonic
evidence, not a direct incompressible two-dimensional benchmark. The article
explicitly states that S207 coordinates are available under license. No license
request is needed for our present literature review.

## Additional methods and numerical verification sources

### N15 — Full MSES operating manual

Drela, M. (2007). *A User's Guide to MSES 3.05*, July 2007.
[MIT manual](https://web.mit.edu/drela/Public/web/mses/mses.pdf).

Full PDF saved and inspected. Useful for force/wake definitions, transition,
fixed-lift operation, grid preparation and comparison protocols. This manual's
version is 3.05; P08's printed example is 3.11. Neither provides a locally
installed solver or the input decks used by the experimental papers.

### N16 — AIRSET geometry manual

Drela, M. *AIRSET user's guide*.
[Official manual](https://web.mit.edu/drela/Public/web/mses/airset.pdf).

Located through the official MSES index; not fully reviewed. Useful follow-up
for element placement and geometry preparation, not experimental evidence.

### N17 — MSES version history

Drela, M. *MSES version notes*.
[Official version notes](https://web.mit.edu/drela/Public/web/mses/version_notes.txt).

Official-index lead for reconciling paper/manual/example versions. Not evidence
that a particular published result used the newest listed version. No solver
capability change has been imported into MultiFoil from these notes.

### N18 — Coupled viscous–inviscid method foundation

Drela, M.; Giles, M. B. (1987). *Viscous-Inviscid Analysis of Transonic and Low
Reynolds Number Airfoils*. AIAA Journal 25(10), 1347–1355.
DOI 10.2514/3.9789.
[Author-hosted source](https://people.maths.ox.ac.uk/gilesm/files/aiaa87b.pdf).

Primary first page/abstract inspected; the web extraction exposed one page,
not the entire article. Foundational ISES coupling, integral boundary layers and
transition treatment. Read the complete formulation before adopting a closure;
its existence does not validate this multi-element implementation.

### N19 — Interactive boundary-layer multielement method

Cebeci, T. (1992). *An Interactive Boundary-Layer Approach to Multielement
Airfoils at High Lift*. Fifth Symposium on Numerical and Physical Aspects of
Aerodynamic Flows.
[NASA source](https://ntrs.nasa.gov/citations/19930018256).

Primary metadata saved; public full text located. Useful alternative treatment
of flap wells and interacting flow regions. Detailed method review and matching
case extraction remain to be done; not a newly established experimental benchmark.

### N20 — TMR two-dimensional multielement CFD verification

Turbulence Modeling Resource. *2D Multi-Element Airfoil Verification Case*.
[Official case](https://tmbwg.github.io/turbmodels/multielementverif.html).

Official page saved and inspected. Geometry/grids and converged CFD results for
a two-dimensional CRM-HL section provide a numerical verification reference.
**The page explicitly supplies no experiment.** Its fully turbulent flow and
high local Mach numbers limit relevance to our incompressible model, even though
freestream Mach is 0.2. Do not turn this into an experimental validation claim.

### N21 — MSES-based wind-turbine airfoil design source

Ragheb, A. M.; Selig, M. S. (2011). *Multi-Element Airfoil Configurations for
Wind Turbines*. AIAA 2011-3971.
[Author paper](https://m-selig.ae.illinois.edu/pubs/RaghebSelig-2011-AIAA-2011-3971-MultiElementAirfoilsWindTurbines.pdf).

Design-source lead connected to the subsequent Pomeroy experiments. Useful for
MSES design choices and configuration provenance; not automatically an
independent measured polar. Detailed geometry/condition extraction remains open.

### N22 — S204 MSES design, without a newly acquired experiment

Somers, D. M. (2012). *Design of a Slotted, Natural-Laminar-Flow Airfoil for
Business-Jet Applications*. NASA/CR–2012-217559.
[NASA report](https://ntrs.nasa.gov/citations/20120013091).

Public report located for a further SNLF design example. Primarily useful for
design/method comparisons; this review did not acquire an associated experimental
dataset. Its higher-speed design conditions are not all within our intended scope.

## Further retrieval leads

These leads were found during the search but lack enough inspected primary
material to make them ready comparison cases. They must not be counted as
full papers already read or as independently acquired experimental datasets.

### L01 — Earlier S414 report, now acquired and inspected

Somers, D. M.; Maughmer, M. D. (2010). *Design and Experimental Results for the
S414 Airfoil*. US Army RDECOM TR-10-D-112, ADA545821.
[University mirror](https://www.ae.metu.edu.tr/tuncer/ae546/multi/ADA545821.pdf),
[DTIC archive mirror](https://archive.org/download/DTIC_ADA545821/DTIC_ADA545821.pdf).
Both 97-page PDFs acquired anonymously on 2 October; archive supplied checksums
and identical Figure-1 vector geometry verified. The report supplies design/test
outlines, pressure-orifice tables, experimental pressure snapshots and MSES
pressure/section comparisons. No original contour table is included; it refers
readers to Airfoils, Incorporated for coordinates. Figure 12(a,b), Re 1 million,
M 0.10, uses the envelope MSES formulation identified in N02, Ncrit 9.
The source-based modified test outline and 188 pressure markers are extracted;
strict geometry admission remains open. See the
[findings and plots](results/s414_geometry_qualification_20261002/FINDINGS.md)
and [access manifest](results/s414_geometry_qualification_20261002/source_manifest.json).
The two public `.pts` files in the teaching catalog are retained as unqualified
geometry candidates, not author coordinate decks. This report overlaps N02's
experiment rather than adding another independent validation dataset.

### L02 — Pomeroy dissertation

Pomeroy, B. W. (2016). *Wake Bursting: A Computational and Experimental
Investigation for Application to High-Lift Multielement Airfoil Design*.
University of Illinois dissertation, cited by
[Pomeroy and Selig 2017](https://m-selig.ae.illinois.edu/pubs/PomeroySelig-2017-AIAA-2017-1210-Design-for-WakeBursting.pdf).
The institutional search returned HTTP 403; no thesis file acquired. That is
an access failure, not proof that an account is required. N03 already supplies
coordinates for the measured three-element baseline.

### L03 — Low-speed S207 coupled-device experiments

*Investigating the Impact of Coupled High-Lift Devices on Natural Laminar Flow
Airfoils*. [Illinois institutional record](https://www.ideals.illinois.edu/items/127506).
Primary search metadata identifies S207 experiments; the direct record returned
HTTP 403. Complete author/date/file qualification remains pending. Potentially
closer to our low-speed scope than N14, but geometry licensing still needs checking.

### L04 — Unpowered subset of a propeller–wing–flap study

Duivenvoorden, R. (2025). *Aerodynamic Phenomena of Propeller–Wing–Flap
Interaction*. TU Delft dissertation.
[Institutional thesis](https://pure.tudelft.nl/ws/portalfiles/portal/243083044/RDuivenvoorden-AerodynamicPhenomenaOfPropellerWingFlapInteraction_20250424_digital.pdf).
A public validation-data lead. Locate the numerical dataset and isolate a
suitable unpowered case before use; propeller flow and spanwise effects cannot
be treated as two-dimensional free-stream airfoil measurements.

### L05 — Further multielement parameter study

*Parameter Analysis of a Multi-Element Airfoil for Application to Airborne Wind
Energy* (2022), Wind Energy Science 7, 1627 onward.
[Open article](https://wes.copernicus.org/articles/7/1627/2022/wes-7-1627-2022.html).
Candidate CFD parameter study for design-space context. Full author/geometry/data
audit remains pending; do not count inherited experimental comparisons as new tests.

### L06 — Direct multielement Newton-method paper

Drela, M. (1990). *Newton Solution of Coupled Viscous/Inviscid Multielement
Airfoil Flows*. AIAA 90-1470.
[Publisher identifier](https://doi.org/10.2514/6.1990-1470).
Bibliographic lead confirmed by later research citations; full primary paper
not acquired through the tested DOI. A high-value method-reading target.

### L07 — Classical high-lift physics

Smith, A. M. O. (1975). *High-Lift Aerodynamics*. Journal of Aircraft
12(6), 501–530. [DOI](https://doi.org/10.2514/3.59830).
Found repeatedly in the inspected papers' reference lists; full original not
read in this review. Background on element interaction, not a substitute for
a matched geometry/measurement dataset.

## Search coverage, saved evidence and checks

The prior inventory was recovered from `sources.json`, `README.md`, `MSES.md`,
`DRAG.md`, `THREE_WAY.md`, `TN808_BENCHMARK.md`, `S702_BENCHMARK.md`,
`INDEPENDENT_RANKING.md`, `OPTIMIZATION_READINESS.md`, `THIRD_PARTY.md`, and
the S702 conditions/wake source audits. This covers documented scientific use;
it does not pretend to recover every incidental browser result from past sessions.

Search families included: S702/S103/S414/S204/S207 MSES experiments; NLR7301
original boundary-layer measurements; multielement wake bursting/confluence;
MSES high-lift wind-tunnel comparisons; thick-airfoil slat optimization;
30P30N/CRM-HL reference data; and Drela/Fidkowski coupled-solver methods.
Primary NASA, MIT, university/author and publisher sources were preferred.
Secondary search hits were used to locate originals, not to manufacture data.

The saved [source manifest](results/reference_review_20260924/source_manifest.json)
records URLs, anonymous access outcomes, timestamps, byte counts and SHA-256
hashes. Twenty successful source downloads and four failed access records were
rechecked after the freeze. Successful acquisitions include metadata, S103/Adair/
30P30N text, S414 and Pomeroy PDFs, the MSES manual, TMR pages and three near-wake
numerical files. Availability is as tested on 24 September, not a guarantee of
future access. Source copyrights/notices remain applicable; downloaded research
papers are not relicensed under MultiFoil's code license.

The [acquisition script](results/reference_review_20260924/acquire_sources.py)
checks saved hashes and fetches only missing items from its bounded selection:

```bash
python validation/results/reference_review_20260924/acquire_sources.py
```

Its selection covers the follow-up downloads, not every discovery-page snapshot
already in the manifest. Failed attempts remain recorded. PDF-derived local
text is a reading aid; coordinate tables and corrupted OCR must be checked
against source page images before numerical extraction.

No solver core, benchmark registry, calibration, acceptance threshold or remote
job changed. Verification for this literature task consists of source/metadata
checks, downloaded-file hashes, bibliography coverage and local-link checks;
no aerodynamic test result is implied.
