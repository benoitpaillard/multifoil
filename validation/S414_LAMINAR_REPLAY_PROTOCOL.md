# S414 attached-laminar replay controls (V5, tranche 1), 4 October 2026

Frozen before any S414 calculation with the new solver. The runner is
`validation/run_s414_laminar_replay.py`; the independent solver is
`validation/s414_laminar_fd.py`. Everything runs locally in seconds to minutes.
No coupled flow solve, force, drag or transition location is computed.

## Question

The earlier prescribed-pressure replay
(`validation/results/s414_prescribed_pressure_20261004/`) was inconclusive: its
model-at-tap sampling controls did not complete and its upstream seed was not
tested. Stage V5 may not interpret a measured-pressure response until those
controls pass. This tranche asks three things, in this order:

1. Why did the earlier controls fail?
2. Under identical edge velocity, does an independent laminar calculation
   confirm the inherited momentum/shape closure up to laminar separation?
3. Restricted to what direct marching can answer (the attached layer up to
   separation): is the response to tap sampling and to the unmeasured upstream
   history bounded, and what do the measured taps then imply about laminar
   separation on the fore-lower nose?

The amplification that produces the modeled early transition occurs downstream
of separation. It is outside this tranche by construction and needs a separately
qualified interactive or inverse method.

## Fixed inputs

- Archived nominal states (legacy-mode run, hash-checked against the retrieval
  receipt): `s414_lower_bucket_20261003/{cl0056_nominal,cl0175_nominal,cl028_control}/fore_state.npz`.
- Tap positions and experimental Cp: `s414_fore_nose_20261003/pressure_at_taps.csv`
  (hash-checked against that study's receipt). Printed-table coordinates are
  primary; plotted coordinates are one declared sensitivity.
- The eight saved profiles and summary of the earlier replay, unchanged.
- Re 1,000,000; fore-lower surface from the stagnation node to x/c 0.03424 (the
  fourth tap); taps 0–4 (x/c 0.00181, 0.00838, 0.01954, 0.03424, 0.05304) define
  every reconstruction. Ncrit and the closure are untouched.

## Methods

- **FD**: Keller-box solution of the incompressible laminar boundary-layer
  equations, similarity start at the first node (stagnation, m = 1), uniform
  grid of 561 points to eta 14, step halving (12 levels) towards separation,
  separation position from the square-root wall-shear behaviour. Edge velocity
  is the model's physical (Karman–Tsien corrected) speed; nu = 1e−6.
- **Closure**: the inherited integral equations (`mfoil.residual_station`) in
  direct mode with fixed edge velocity, on the same nodes, at Mach 0 with the
  same speed and viscosity. Stagnation start as in the inherited initializer.
  The march ends at the domain end, at the first station with non-positive skin
  friction (separation by linear interpolation of cf·Re_theta), or at a recorded
  nonlinear failure. Nothing downstream of separation is computed.
- **Thwaites** (supporting only, never a gate).
- Marching nodes: native nodes with 8 subintervals; the speed between nodes is
  a PCHIP interpolant in surface distance. Tap positions are inserted as nodes.

## Part A — classification of the earlier failures

Along both archived states, evaluate the determinant of the direct-mode station
Jacobian d(momentum, shape)/d(theta, delta*) and locate its sign change. Compare
with the failing stations recorded in the earlier summary. At the first station
past the sign change, change the prescribed edge velocity by relative
±1e−6…±1e−3 and attempt the unchanged earlier station solver from the saved
state. The classification "direct-mode singularity" requires every recorded
failure to lie at or after the sign-change interval. This is a diagnosis of the
method, not a statement about physical separation.

## Part B — qualification of the independent solver (gate)

| Control | Reference | Acceptance |
|---|---|---|
| Blasius wall shear, shape factor | 0.332057, 2.5911 | ±1e−4, ±1e−3 |
| Hiemenz wall shear | 1.232588 | ±2e−4 |
| Marching ue = xi^m, m = 1, 0, 1/3, −0.05 | similarity profile preserved | wall-shear drift < 1e−8 |
| Howarth ue = 1 − xi, separation | 0.1198 (0.1198–0.1199 in the literature) | ±0.0005 on the two finest of four grids |
| Cylinder ue = 2 sin x, separation | 104.5° | ±0.3° on the two finest of three grids |

A failed control stops Parts C–F. The closure and Thwaites are run on the
Howarth and cylinder flows and reported; they are not gates.

## Part C — closure against FD under the model's own pressure

Three archived states, native edge velocity. Report theta, delta*, H, cf and the
separation position from each method. Declared agreement screen: same topology
(separated or attached to x/c 0.03424), separation positions within 0.0012 of
chord (one native node interval there) and theta at separation within 5 %.
Controls: the closure at Mach 0.10 on native nodes against the archive (reported);
8 against 4 and 16 subintervals and PCHIP against a cubic spline (reported).
Exceeding the screen is a finding about the closure near separation. It does
not stop Parts D–F, which report both methods.

## Reconstructions

For a base model state A (speed q_A, pressure Cp_A) and target tap values Y_k:

| Family | Between and beyond taps (xi ≥ tap 0) |
|---|---|
| F1 linear | Cp linear in surface distance through the taps |
| F2 monotone | Cp by PCHIP through the taps |
| F3 shape, additive | Cp_A + linear interpolation of (Y_k − Cp_A at tap k) |
| F4 shape, ratio | q_A × linear interpolation of (q_Y,k / q_A at tap k) |

Upstream of tap 0 the speed is q_A × g, with g0 the ratio at tap 0:
U1 hold (g = g0 back to stagnation), U2 ramp (g from 1 at stagnation to g0),
U3 late (g = 1 up to the last native node before tap 0, then linear to g0; this
is the earlier replay's seed treatment). F1/F2 contain no between-tap peak;
F3/F4 superpose the base model's between-tap shape. Neither is measured.

## Part D — sampling controls

- Self: target = the base state's own tap values. F3/F4 reduce to the native
  profile (asserted numerically, recorded as trivial). F1/F2 are marched.
- Cross: target = tap values of a *different* archived state B, base A. Pairs
  0.056→0.175, 0.175→0.056, 0.175→0.284, 0.284→0.175. Truth is the response under
  B's native pressure with the same method.

A reconstruction passes a control if topology matches, separation positions
agree within 0.0012 of chord and theta at separation (or at the domain end if
attached) agrees within 5 %. FD is the gate; the closure is reported alongside.

## Part E — upstream sensitivity

The spread over U1–U3 in every cross control and experimental replay, against
the same tolerances. Additionally, closure only: restart at the earlier seed
node (x/c 0.00165) with theta and delta* scaled by 0.8 and 1.25 at constant H
under native pressure at CL 0.056 and 0.175 (four marches).

## Part F — measured taps

Experimental tap values at CL 0.056, 0.175 and 0.284, base = the model state at
the same lift, all 12 reconstructions, plus the plotted-coordinate sensitivity
for U1. A reconstruction is **admitted** at a lift only if it passed every
control with that base state. Outcome per lift, over admitted reconstructions
with FD and closure agreeing in topology:

- *separation retained*: all separate upstream of x/c 0.03424;
- *separation removed*: none separates;
- *not bounded*: mixed, or nothing admitted. Then the result is coverage only
  and the missing measurement is stated.

The full 12-member bracket is always reported, flagged by admission. These are
matched-lift pressure snapshots on a model-dependent stagnation position and
remain labeled conditional.

## Limits

At most 120 FD marches and 130 closure marches (the same set plus the declared
closure-only controls of Parts C and E); no coupled solve; no
extension, no threshold change after execution, no fitting to drag, no closure,
Ncrit or geometry change. Failed marches are kept as failures. A completed
script is not a physical validation.

## Amendment after the first execution (4 October 2026; recorded deviation)

The first execution completed every part (114 FD and 121 closure marches) and
is kept unchanged, with the exact sources that produced it, in
`results/s414_laminar_replay_20261004/initial_run_box_scheme/`. Visual review of
its figures showed a station-to-station sawtooth in the FD wall shear wherever
the edge velocity has a kink (sawtooth index up to 0.15). This is a known
property of the undamped centred box scheme; smooth qualification flows did
not reveal it. The FD results of that execution are therefore not used.

Two numerical changes were made before a complete rerun, and nothing else:

1. FD streamwise differencing: second-order backward differences in ln(xi) at
   each station (first order where a step grows by more than a factor two), in
   place of the centred differencing. The discretization across the layer, the
   grid, the similarity start, the step halving and the separation estimate
   are unchanged. The centred scheme stays available as `scheme='box'`.
2. Closure march: a failed station solve is halved (12 levels) before a
   failure is recorded, as the FD march already did. In the first execution the
   closure ended by a nonlinear failure within one node of the FD separation
   in 73 marches; this resolves how it ends.

One control is added to the Part B gate: a flat plate followed by a 10 % linear
deceleration with a kink at xi 0.5 (attached throughout). The sawtooth index
(median absolute second difference of the wall shear over its median) must be
below 1e−3 on 101 and 201 nodes per segment. The centred scheme gives
0.03–0.05 there.

Every acceptance value of Parts B–F, the reconstructions, the controls, the
march limits and the outcome rules are unchanged. This deviation was decided
during execution and was not approved beforehand; both executions are reported.
