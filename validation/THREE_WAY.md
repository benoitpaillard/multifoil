# MultiFoil / MSES / experiment overlays

![Four-panel comparison](results/three_way/multifoil_mses_experiment.png)

The figure uses archived solver results. No new flow solve, coefficient fitting,
angle correction, or interpolation for scoring was performed. Lines connect
available solver samples as visual guides; experimental observations are markers.
The low-angle drag panel repeats a subset of the full NACA0012 case.

| Panel | Observable | Conditions | References |
|---|---|---|---|
| A | NLR7301 CL(alpha) | Re2.51e6, M0.185, nominal2.6%gap | MultiFoil; ICAS2006 MSES curve; official MIT MSES point; NLR experimental digitization |
| B | NLR7301 CD(alpha) | Same nominal configuration and conditions as A | MultiFoil; official MIT MSES point; experimental data |
| C | NACA0012 CD(alpha) | Re6e6, M0.30 | MultiFoil; COBEM2005 MSES and experimental curves |
| D | NACA0012 CD(alpha), 0–10degrees | Same as C | Same three series, enlarged vertical scale |

NLR MultiFoil values are from `results/nlr7301_mses/polar.csv`; the source of
these model outputs does not make their CD an MSES prediction. All seventeen
NLR experimental CL observations and three CD observations are displayed at
their own digitized angles. The experimental drag point near10.14degrees is
visual context and does not add a new exactly matched solve to the validation
statistics. No MSES drag curve was invented or substituted from another case.

The orange diamond is the single printed state from
[MIT's official NLR7301 example](https://web.mit.edu/drela/Public/web/mses/nlf26.pdf):
MSES 3.11, alpha 6 degrees, CL 2.4840, CD 0.02264, Re 2.51 million,
Mach 0.185 and Ncrit 9. It is visually transcribed in
`data/nlr7301_mses_official/`, with the source PDF hash. The exact geometry and
full input deck remain unconfirmed, including an additional configuration suffix
in the title. This point is qualified plot context, not a scored benchmark or a
new drag polar. No physical uncertainty is inferred from printed precision.
The figure and exported table now contain 90 displayed samples including repeated
subsets and both observables for the official point.

For NACA0012, the experimental green curve was extracted from Figure3, page5,
of [Lima e Silva et al., COBEM2005](https://abcm.org.br/anais/cobem/2005/PDF/COBEM2005-0277.pdf).
That paper attributes it to Seifert and LaTunia (1999). The nine visible line
vertices at0,2,4,8,10,12,14,16,18degrees use the same calibrated axes as the red
MSES curve. These are experimental values reproduced in a later publication;
the original raw measurement records were not accessed. The extraction is
stored in `data/naca0012_mses_drag/experiment.csv` and
`experiment_digitization.json`, with hashes in `sources.json`.

The Ladson NACA0012 data at Mach0.15 were not mixed into the Mach0.30 comparison.
Nominal conditions match, but original trips/Ncrit, exact coordinates, tunnel
corrections and grid settings remain unverified. MultiFoil assumes Ncrit9.
Vertical bars show graphical reading allowances, not experimental confidence
intervals. Their incidence allowances are documented separately. The12,14,16,18
degree NACA model points flag local supersonic flow; its compressibility
correction does not resolve shocks. See [DRAG.md](DRAG.md) and [MSES.md](MSES.md)
for fuller limitations and the existing comparison results.

NLR references: [Neto et al., ICAS2006 Figure6](https://www.icas.org/icas_archive/ICAS2006/PAPERS/512.PDF)
and [Ranneberg's published NLR digitization](https://github.com/mranneberg/viiflow-examples/tree/93bef6c6b779db4605f6bcbb2caf1d88a51c9b85/NLR7301%202-Element%20Airfoil).
Their existing provenance and geometry qualifications are retained.

Reproduce from the extracted project directory:

```bash
python validation/plot_three_way.py
# Also check the experimental extraction against the original PDF (PyMuPDF):
python validation/plot_three_way.py --pdf-source /path/to/COBEM2005-0277.pdf
```

Outputs are PNG, vector PDF and `plot_data.csv` in `validation/results/three_way`.
The CSV contains every displayed series with its panel, reading interval and
model flag; panelD repeats the appropriate subset of panelC. The new
experimental series is plot context and is not silently added to the benchmark
scoring or used to change earlier error statistics.
