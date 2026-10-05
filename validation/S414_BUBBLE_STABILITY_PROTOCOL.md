# S414 nose bubble: independent boundary layer and linear-stability amplification (V5, tranche 2), 5 October 2026

Frozen before any S414 calculation with the new solvers. Runner:
`validation/run_s414_bubble_stability.py`. New modules, none sharing a formula
with the inherited closure: `s414_bubble_fd.py` (inverse finite-difference
laminar layer), `s414_stability.py` (Orr–Sommerfeld), `s414_amplification.py`
(e^N). Everything runs locally. No coupled flow solve, force, drag or
transition location is computed.

## Question

Tranche 1 closed the prescribed-pressure replay and cross-checked the attached
laminar closure. The plan's remaining V5 item is the separated layer and the
amplification through the nose recovery. Given the **displacement thickness of
the archived bubble**, do independent calculations reproduce

1. the edge pressure, momentum thickness, shape factor and skin friction of the
   inherited separated-laminar closure, and
2. the amplification N that the inherited envelope correlation accumulates in
   the bubble?

This tests the closure on the model's own displacement body. It cannot show
whether the bubble itself is right, and it says nothing about the measured
flow. The comparison basis is identical displacement thickness downstream of a
declared switch station and identical edge velocity upstream of it.

## Fixed inputs

- Archived nominal states, hash-checked as in tranche 1:
  `s414_lower_bucket_20261003/{cl0056_nominal,cl0175_nominal,cl028_control}/fore_state.npz`.
- Re 1,000,000; fore-lower surface from the stagnation node to the archived
  node at x/c 0.02083 (the last laminar node of the CL 0.056 state and the
  endpoint of the earlier replay). nu = 1e−6; edge speed is the model's
  physical (Karman–Tsien corrected) speed. Ncrit and the closure are untouched.

## Methods

- **Boundary layer.** Finite differences in Falkner–Skan variables, 1201 points
  to eta 30, second-order backward streamwise differences, 8 subintervals per
  archived node interval. Direct mode (edge velocity prescribed) from the
  stagnation similarity start to the switch node: the first archived node with
  Hk ≥ 2.9. Inverse mode beyond it: displacement thickness prescribed, edge
  velocity solved. The prescribed displacement is a PCHIP interpolant in surface
  distance through the solver's own value at the switch node and the archived
  values at the later nodes. Reversed flow uses the FLARE approximation.
- **Stability.** Orr–Sommerfeld, Chebyshev collocation (141 points, domain 60
  displacement thicknesses), spatial eigenvalue at real frequency by Newton
  iteration from a temporal scan. Stations: every marching node from the
  archived node of largest edge speed to the endpoint. 52 frequencies,
  geometrically spaced over omega·theta/ue = 0.004…0.6 at the median station.
  N of each frequency is integrated from its neutral point; N is the envelope.
- **Inherited correlation on the independent state.** The inherited growth
  function evaluated on the finite-difference theta, delta* and edge speed and
  integrated with the same trapezoid rule, at Mach 0. This separates a
  difference in the boundary-layer state from a difference in the growth law.

## Part A — qualification (gate; a failed control stops Parts B and C)

| Control | Reference | Acceptance |
|---|---|---|
| Blasius spatial eigenvalue, Re_δ* 998, ω 0.1122 | 0.3086 − 0.0057i (Jordinson) | ±5e−4 each part |
| Blasius temporal phase speed, Re_δ* 998, α 0.308 | 0.3641 + 0.0080i | ±5e−4 each part |
| Blasius critical Reynolds number | 519.4 | largest growth negative at 515, positive at 525 |
| Collocation 141 against 201 points, same eigenvalue | – | change below 1e−6 |
| Falkner–Skan reversed-flow wall shear, β −0.18, −0.15, −0.10, −0.05 | −0.0976, −0.1334, −0.1405, −0.1082 | ±0.001 |
| Falkner–Skan separation | β −0.1988, H 4.03 | ±0.0005, ±0.01 |
| New solver in direct mode against the tranche-1 solver (cylinder) | identical | theta within 1e−8 |
| Direct solution fed back as prescribed displacement (cylinder) | edge-speed ratio 1 | within 1e−8 |
| Blasius e^N envelope against the inherited correlation | Re_θ at N = 6 | within 10 % |

Reported, not gates: the e^N envelope of similar profiles with H 4.0, 6.4, 8.8
and 11.3 up to Re_θ 152 against the inherited correlation; the change of a
similar reversed profile under FLARE; 26 against 52 frequencies on Blasius.

## Part B — separated-laminar closure on the archived displacement

Three states. At every archived node downstream of the switch report Cp, theta,
H, cf·Re_θ and the strongest reversed velocity, against the archive.

Declared screen, CL 0.056 and 0.175: largest |Cp difference| at those nodes at
most 0.05 (the materiality threshold of stage V4) and theta at the endpoint
within 5 %. CL 0.284 has no bubble and is the control. Reversed velocity above
10 % of the edge speed marks the stations where FLARE is only approximate.

Sensitivities, CL 0.056 and 0.175 only, each changing one thing: switch one
node upstream; linear instead of PCHIP displacement; the archived displacement
scaled by the ratio at the switch node instead of joined; 16 subintervals.

## Part C — amplification

At the endpoint and at every archived node from the switch on: archived N; the
inherited correlation on the independent state; Orr–Sommerfeld N.

Declared screen at the endpoint, CL 0.056 and 0.175: Orr–Sommerfeld N divided
by the inherited correlation on the same independent state within 0.80…1.25.
Inside: the growth law is consistent with linear stability on these profiles.
Outside: a finding about the growth law in this bubble, with the ratio stated.
Whether Orr–Sommerfeld N reaches 9 inside the domain is reported, not
interpreted as a transition location. Amplification sensitivities: the switch
node and the displacement interpolant, both bubble states.

## Interpretation, fixed beforehand

- Both screens inside at both bubble states: the inherited separated-laminar and
  growth closures reproduce independent calculations on the model's own
  displacement body; this check does not support a closure limitation. The
  open question is then whether the bubble exists at all under the true
  pressure input, which this tranche cannot answer.
- A screen outside: a quantified closure difference in this bubble. It is a
  candidate, not a remedy: no model change, no Ncrit change and no drag claim
  follow. A change would need stage V6 with its own justification and
  regression set.
- Every statement is conditional on the archived displacement, on parallel-flow
  stability theory and, where flagged, on FLARE.

## Limits

At most 15 boundary-layer marches and 8 envelope calculations on S414 states;
no coupled solve; no extension; no threshold change after execution. Failed
marches or lost waves are kept and counted. A completed script is not a
physical validation.
