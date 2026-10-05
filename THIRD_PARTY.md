# Third-party source and numerical changes

`vendor/mfoil.py` is the 2023-06-28 Python version of mfoil by
Krzysztof J. Fidkowski. Its MIT license and copyright notice are retained at
the top of that file. The archive's root MIT license covers the new code;
it does not replace the upstream copyright notice.

Author and documentation:
https://websites.umich.edu/~kfid/codes.html

Author's Python download:
https://websites.umich.edu/~kfid/codes/mfoil/mfoil_notex.py

The copy used here was obtained from this repository mirror:
https://github.com/SharathSPhD/CINS/blob/514579fb4c65467d26c493c6712490f70ded2afe/vendor/mfoil/mfoil.py

SHA-256 of the unmodified bundled file:
`c9473e895a7256990c07630eb8bcc7b8310bf3dcb3cc6f62aab2e1bf61acc651`

`viscous.py` reuses upstream residuals and constitutive models. It builds a
multi-body streamfunction system, one Kutta condition per element, separate
boundary layers/wakes, and one Newton system with full displacement-source
coupling. It provides its own geometry setup, wake tracing, dense residual
assembly and scaled linear solve. The dense BL residual/Jacobian are checked
against the original sparse assembly.

A finite zero-gap limit is supplied for TE_info. Actual viscous geometries use
a recorded minimum TE opening because the exactly closed viscous source limit
is singular. Curvature spacing is performed in the local chord frame.

For influence on other bodies, each elementary source streamfunction is
unwrapped along the target contour before assembly. Otherwise upstream branch
cuts can introduce artificial jumps across a different body. The first wake
half-panel source on other bodies includes the terminal surface contributions.
The upstream near-TE approximation on the source body's own streamfunction is
retained. For wake velocities, both weights of the shared midpoint source are
included on both adjacent half panels. This corrects two omitted terms in the
2023 Python velocity assembly. An independent contour integral of analytic
normal velocity verifies the cross-body streamfunction operator.

The simplified e^N, laminar/turbulent, transition, and shear-lag boundary-layer
closures themselves are unchanged. Output Cf includes a density-ratio conversion
at nonzero Mach. mfoil's original geometry constructor and sparse shape-check
issues are bypassed; no installed package is modified.

No Viiflow binary is included. The optional adapter requires the user's separate
licensed installation. No other solver's source is bundled.

## Validation data

The NLR7301 ECARP geometries and experimental digitizations are redistributed
from Maximilian Ranneberg's MIT-licensed `viiflow-examples` repository at commit
`93bef6c6b779db4605f6bcbb2caf1d88a51c9b85`. The original notice is preserved in
`validation/data/nlr7301/LICENSE`. The full original NLR experiment report was
not accessed; provenance and data qualifications are in `validation/README.md`.

NACA 0012 table values and S809 coordinates/figure values were transcribed from
NASA TM-4074 and NREL/SR-440-6918 respectively. Source links, PDF hashes, figure
locations and extraction audits are in `validation/sources.json` and
`validation/audit/`. Digitization allowances are not experimental uncertainties.

The MSES reference is digitized from Figure6 of Neto et al., ICAS2006 paper512
(https://www.icas.org/icas_archive/ICAS2006/PAPERS/512.PDF). A single figure
excerpt, calibration and audit overlay are included for reproducibility, with
attribution in `validation/MSES.md`; copyright remains with the source owner.
No proprietary MSES executable, source code or input deck is included.

The MSES drag reference uses Figure3 of Lima e Silva et al., COBEM2005
(https://abcm.org.br/anais/cobem/2005/PDF/COBEM2005-0277.pdf), copyright2005ABCM.
The figure crop, selected PDF vector coordinates, calibration and audit overlay
support reproducible scientific comparison. Attribution and limitations are in
`validation/DRAG.md`; the original figure is not covered by the code's license.

The fallback curvature measure in `paneling.py` is adapted from the same
MIT-licensed mfoil paneling routine. It changes the inverse spacing interpolator
only when the original routine returns invalid point ordering. The original
vendor file is retained unchanged. NACA TN808 (Lowry,1941) and TR664
(Wenzinger and Harris,1939) are public US government reports; source PDFs and
figure extracts are distributed with citation for reproducible comparison.

The S702 source is Somers and Maughmer (2022), NASA/CR-20220014964,
https://ntrs.nasa.gov/api/citations/20220014964/downloads/NASA-CR-20220014964.pdf.
The report and its original notices are retained. Geometry vertices, MSES
curve paths and experimental marker centers are extracted from its vector
figures for reproducible comparison, with attribution and qualifications in
`validation/S702_BENCHMARK.md`. The code license does not replace source
document notices. No author geometry deck or MSES software is distributed.
