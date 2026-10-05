"""Derivative-only correction to the historical mfoil Jacobian.

The vendored source stays unchanged. CoupledSolver defaults to corrected_af;
explicit legacy mode preserves historical math. A dispatcher reads each model's
copied Param, so legacy and corrected models can coexist in the same process.
Untagged parameters, including direct vendored use, retain legacy behavior.
"""
import numpy as np
from viscous import CoupledSolver, mf

LEGACY = 'legacy'
CORRECTED = 'corrected_af'
MODES = (LEGACY, CORRECTED)
_original_get_damp = mf.get_damp


def get_damp(U, param):
    """Preserve the exact growth value; correct only d[0.1 exp(-20 Hmi)]/dU."""
    mode = getattr(param, 'amplification_jacobian', LEGACY)
    if mode not in MODES:
        raise ValueError('Unknown amplification Jacobian: '+str(mode))
    value, derivative = _original_get_damp(U, param)
    if mode == LEGACY:
        return value, derivative
    hk, hk_u = mf.get_Hk(U, param)
    if hk < 1.05:
        hk, hk_u = 1.05, hk_u*0
    hmi = 1/(hk-1)
    ret = mf.get_Ret(U, param)[0]
    critical = 2.492*hmi**.43+.7*(np.tanh(14*hmi-9.24)+1)
    # Match the original activation arithmetic and branch exactly.
    lr = np.log(ret)/np.log(10)
    if lr < critical-.1:
        return value, derivative
    rn = (lr-(critical-.1))/.2
    rf = 1 if rn >= 1 else 3*rn**2-2*rn**3
    da = .028*(hk-1)-.0345*np.exp(-(3.87*hmi-2.52)**2)
    missing = -np.exp(-20*hmi)*(-hmi**2*hk_u)
    return value, derivative+rf*da/U[0]*missing


def configure(models, mode):
    """Set an explicit per-model option; untagged parameters retain legacy math."""
    if mode not in MODES:
        raise ValueError('Unknown amplification Jacobian: '+str(mode))
    if mf.get_damp not in (_original_get_damp, get_damp):
        raise RuntimeError('Another get_damp replacement is already installed')
    if mode == CORRECTED:
        mf.get_damp = get_damp
    for model in models:
        model.param.amplification_jacobian = mode


class AmplificationJacobianSolver(CoupledSolver):
    """Compatibility name for the correction now integrated into CoupledSolver."""
    def __init__(self, *args, amplification_jacobian=CORRECTED, **kwargs):
        super().__init__(*args, amplification_jacobian=amplification_jacobian, **kwargs)
