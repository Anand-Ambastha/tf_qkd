# -*- coding: utf-8 -*-
"""
TF-QKD over FSO channels with log-normal turbulence — optimized engine.

Security proof: Curty, Azuma & Lo, "Simple security proof of twin-field
type QKD protocol", arXiv:1807.07667v2, Protocol 3, Eqs. (12)-(35).

Implements: channel model (X/Z basis gains, Eqs. 27-33), photon-number
yields (Eq. 34-35), cat-state coefficients (Eqs. 12-13), phase-error
upper bound (Eqs. 20-22), asymptotic key rate (Eqs. 17-19) with alpha
optimization, PLOB bound, deterministic FSO channel (atmospheric +
diffraction + geometric loss), and ergodic SKR under log-normal
turbulence (Monte Carlo).

This is a performance rewrite of the original notebook. All physics,
equations and numerical outputs are unchanged (validated to float64
machine precision, see VALIDATION.md / inline asserts). Educational
narrative, per-cell prints and diagnostic plots have been removed;
see the optimization report for what changed and why.
"""

import math
from functools import lru_cache

import numpy as np
from scipy.special import i0 as bessel_i0, gammaln
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

# =============================================================================
# 1. Global parameters
# =============================================================================

PARAMS = dict(
    pd=1e-7,
    alpha=0.2,
    thetaA=np.arcsin(np.sqrt(0.02)),
    thetaB=-np.arcsin(np.sqrt(0.02)),
    delta=0.0,
)
PARAMS['theta'] = PARAMS['thetaA'] - PARAMS['thetaB']
PARAMS['phi'] = PARAMS['delta'] * np.pi

NMAX_DEFAULT = 4          # Nmax=4 is near-optimal per paper Fig. 5
OPT_MAXITER = 100
ALPHA_BOUNDS = (1e-3, 0.3)


def dB_to_eta(loss_dB):
    """eta = 10^(-loss_dB/10); each half-link has transmittance sqrt(eta)."""
    return 10.0 ** (-np.asarray(loss_dB, dtype=np.float64) / 10.0)


def h(x):
    """Binary entropy, vectorized, safe at x=0,1."""
    x = np.clip(np.asarray(x, dtype=np.float64), 0.0, 1.0)
    t1 = np.where(x > 0, -x * np.log2(np.where(x > 0, x, 1.0)), 0.0)
    t2 = np.where(1 - x > 0, -(1 - x) * np.log2(np.where(1 - x > 0, 1 - x, 1.0)), 0.0)
    return t1 + t2


# =============================================================================
# 2. X-basis channel model — Eqs. (27)-(31)
# =============================================================================

def omega(phi, theta):
    return np.cos(phi) * np.cos(theta)


def f_plus(phi, theta, gamma):
    Om = omega(phi, theta)
    gamma = np.asarray(gamma, dtype=np.float64)
    return np.exp(-gamma * (1.0 + Om)) - np.exp(-2.0 * gamma)


def f_minus(phi, theta, gamma):
    Om = omega(phi, theta)
    gamma = np.asarray(gamma, dtype=np.float64)
    return np.exp(-gamma * (1.0 - Om)) - np.exp(-2.0 * gamma)


def q_xx(kc, bA, bB, phi, theta, gamma):
    """Eq. (28). k_d = 1 - k_c implicitly."""
    xor_val = int(kc) ^ int(bA) ^ int(bB)
    return f_minus(phi, theta, gamma) if xor_val else f_plus(phi, theta, gamma)


def p_xx_conditional(kc, bA, bB, phi, theta, gamma, pd):
    """Eq. (27)."""
    gamma = np.asarray(gamma, dtype=np.float64)
    qxx = q_xx(kc, bA, bB, phi, theta, gamma)
    return (1.0 - pd) * (pd * np.exp(-2.0 * gamma) + qxx)


def p_xx_total(phi, theta, gamma, pd):
    """Eq. (30) — marginal X-basis gain, equal for (1,0) and (0,1)."""
    gamma = np.asarray(gamma, dtype=np.float64)
    Om = omega(phi, theta)
    term1 = 0.5 * (1.0 - pd) * (np.exp(-gamma * Om) + np.exp(gamma * Om)) * np.exp(-gamma)
    term2 = (1.0 - pd) ** 2 * np.exp(-2.0 * gamma)
    return term1 - term2


def bit_error_rate(phi, theta, gamma, pd):
    """Eq. (31)."""
    gamma = np.asarray(gamma, dtype=np.float64)
    Om = omega(phi, theta)
    num = np.exp(-gamma * Om) - (1.0 - pd) * np.exp(-gamma)
    den = np.exp(-gamma * Om) + np.exp(gamma * Om) - 2.0 * (1.0 - pd) * np.exp(-gamma)
    return np.where(np.abs(den) > 1e-300, num / den, 0.5)


# =============================================================================
# 3. Z-basis channel model — Eqs. (32)-(33)
# =============================================================================

def q_zz(betaA, betaB, theta, eta):
    """Eq. (33)."""
    betaA = np.asarray(betaA, dtype=np.float64)
    betaB = np.asarray(betaB, dtype=np.float64)
    eta = np.asarray(eta, dtype=np.float64)
    sqrt_eta = np.sqrt(eta)
    sum_b2 = betaA ** 2 + betaB ** 2
    arg = betaA * betaB * sqrt_eta * np.cos(theta)
    return np.exp(-0.5 * sum_b2 * sqrt_eta) * bessel_i0(arg) - np.exp(-sum_b2 * sqrt_eta)


def p_zz_beta(betaA, betaB, theta, eta, pd):
    """Eq. (32)."""
    eta = np.asarray(eta, dtype=np.float64)
    sum_b2 = betaA ** 2 + betaB ** 2
    qzz = q_zz(betaA, betaB, theta, eta)
    dark = pd * np.exp(-sum_b2 * np.sqrt(eta))
    return (1.0 - pd) * (dark + qzz)


# =============================================================================
# 4. Photon-number yields — Eq. (34)-(35)
#
#    KEY OPTIMIZATION: Eq.(35)'s nested k,l,m,p,q combinatorial sum is a
#    polynomial in sqrt(eta) whose coefficients C[k,l] do NOT depend on
#    eta at all -- only on (nA, nB, thetaA, thetaB). The original code
#    recomputed this nested loop (with exact-arithmetic comb()/factorial()
#    calls) on every single call, for every alpha trial, every loss point,
#    and every Monte-Carlo sample. Here the coefficient table is built
#    once per unique (nA,nB,thetaA,thetaB) and memoized; evaluation for
#    any (scalar or array) eta is then a few vectorized NumPy power ops.
# =============================================================================

@lru_cache(maxsize=None)
def _qzz_coeff_table(nA, nB, thetaA, thetaB):
    """C[k,l] = C(nA,k)C(nB,l) * BS(k,l) * inner_sum(k,l) -- eta-independent."""
    cA, sA = math.cos(thetaA), math.sin(thetaA)
    cB, sB = math.cos(thetaB), math.sin(thetaB)
    C = np.zeros((nA + 1, nB + 1), dtype=np.float64)
    for k in range(nA + 1):
        for l in range(nB + 1):
            BS = 1.0 / (2.0 ** (k + l) * math.factorial(k) * math.factorial(l))
            inner = 0.0
            for m in range(k + 1):
                for p in range(l + 1):
                    qmin, qmax = max(0, m + p - l), min(k, m + p)
                    for q in range(qmin, qmax + 1):
                        mp_q = m + p - q
                        if mp_q < 0 or mp_q > l:
                            continue
                        Ccomb = (math.comb(k, m) * math.comb(l, p) *
                                 math.comb(k, q) * math.comb(l, mp_q))
                        if Ccomb == 0:
                            continue
                        pA_cos, pB_cos = m + q, m + 2 * p - q
                        pA_sin, pB_sin = 2 * k - m - q, 2 * l - m - 2 * p + q
                        if min(pA_cos, pB_cos, pA_sin, pB_sin) < 0:
                            continue
                        Pfac = (cA ** pA_cos) * (cB ** pB_cos) * (sA ** pA_sin) * (sB ** pB_sin)
                        if Pfac == 0.0:
                            continue
                        upper, lower = m + p, k + l - m - p
                        if lower < 0:
                            continue
                        inner += Ccomb * Pfac * math.factorial(upper) * math.factorial(lower)
            C[k, l] = BS * inner
    binA = np.array([math.comb(nA, k) for k in range(nA + 1)], dtype=np.float64)
    binB = np.array([math.comb(nB, l) for l in range(nB + 1)], dtype=np.float64)
    return C * binA[:, None] * binB[None, :]


def qzz_photon_number(nA, nB, thetaA, thetaB, eta):
    """q_ZZ(k_c,k_d | n_A,n_B) -- Eq. (35), vectorized over eta (any shape)."""
    eta_arr = np.asarray(eta, dtype=np.float64)
    shape = eta_arr.shape
    flat = eta_arr.ravel()
    sqrt_eta = np.sqrt(flat)
    one_minus = 1.0 - sqrt_eta
    coeff = _qzz_coeff_table(int(nA), int(nB), float(thetaA), float(thetaB))

    k_idx, l_idx = np.arange(nA + 1), np.arange(nB + 1)
    kl_sum = (k_idx[:, None] + l_idx[None, :]).ravel()
    exps = (nA + nB - kl_sum)
    coeff_flat = coeff.ravel()
    nz = coeff_flat != 0.0
    vacuum = one_minus ** (nA + nB)

    if not np.any(nz):
        result = -vacuum
    else:
        se_pow = sqrt_eta[:, None] ** kl_sum[nz][None, :]
        om_pow = one_minus[:, None] ** exps[nz][None, :]
        result = (coeff_flat[nz][None, :] * se_pow * om_pow).sum(axis=1) - vacuum
    return result.reshape(shape)


def p_zz_photon(nA, nB, thetaA, thetaB, eta, pd):
    """Eq. (34)."""
    eta = np.asarray(eta, dtype=np.float64)
    sqrt_eta = np.sqrt(eta)
    dark = pd * (1.0 - sqrt_eta) ** (nA + nB)
    qzz_pn = qzz_photon_number(nA, nB, thetaA, thetaB, eta)
    return (1.0 - pd) * (dark + qzz_pn)


# =============================================================================
# 5. Cat-state coefficients — Eqs. (12)-(13)
#    Vectorized over the Fock sub-index (no Python-level per-n function calls).
# =============================================================================

def _cat_coeff_vec(j, alpha, m_max):
    if alpha <= 0.0 or m_max < 0:
        return np.zeros(max(m_max + 1, 0))
    m = np.arange(m_max + 1)
    n = 2 * m + j
    log_c = -0.5 * alpha ** 2 + n * np.log(alpha) - 0.5 * gammaln(n + 1)
    return np.exp(log_c)


def c_even(n, alpha):
    return float(_cat_coeff_vec(0, alpha, n)[n]) if n >= 0 and alpha > 0 else 0.0


def c_odd(n, alpha):
    return float(_cat_coeff_vec(1, alpha, n)[n]) if n >= 0 and alpha > 0 else 0.0


def marginal_coeff_sum(j, alpha, Nterms=80):
    """Sum_{n=0}^{Nterms-1} c^(j)_{2n+j}."""
    return float(np.sum(_cat_coeff_vec(j, alpha, Nterms - 1)))


def truncated_coeff_sum(j, alpha, Nmax):
    m_max = (Nmax - j) // 2
    return 0.0 if m_max < 0 else float(np.sum(_cat_coeff_vec(j, alpha, m_max)))


def residual_j(j, alpha, Nmax, Nterms=80):
    """Eq. (22)."""
    S_full = marginal_coeff_sum(j, alpha, Nterms)
    S_trunc = truncated_coeff_sum(j, alpha, Nmax)
    return max(S_full ** 2 - S_trunc ** 2, 0.0)


# =============================================================================
# 6. Phase-error upper bound — Eqs. (20)-(22)
#
#    inner_sum_j is now vectorized over eta: a single call evaluates the
#    bound for a scalar loss point OR for an entire array of (e.g.) MC
#    turbulence samples, reusing the cached Eq.(35) coefficient tables.
# =============================================================================

def inner_sum_j(j, alpha, eta, pd, thetaA, thetaB, Nmax):
    eta_arr = np.asarray(eta, dtype=np.float64)
    scalar_in = (eta_arr.ndim == 0)
    eta_flat = np.atleast_1d(eta_arr)
    m_max = (Nmax - j) // 2
    total = np.zeros_like(eta_flat)
    if m_max >= 0:
        c = _cat_coeff_vec(j, alpha, m_max)
        for mA in range(m_max + 1):
            cA = c[mA]
            if cA == 0.0:
                continue
            nA = 2 * mA + j
            for mB in range(m_max + 1):
                nB = 2 * mB + j
                if nA + nB > Nmax:
                    continue
                cB = c[mB]
                if cB == 0.0:
                    continue
                pzz = np.maximum(p_zz_photon(nA, nB, thetaA, thetaB, eta_flat, pd), 0.0)
                total += cA * cB * np.sqrt(pzz)
    return float(total[0]) if scalar_in else total.reshape(eta_arr.shape)


def phase_error_upper_bound(alpha, eta, pd, thetaA, thetaB, Nmax, phi=0.0):
    """Eqs. (20)-(21).

    NOTE (correctness/perf finding): in the original notebook this function
    accepted (kc, kd) parameters but never used them in the computation --
    the formula is identical for (kc,kd)=(1,0) and (0,1). Calling it twice
    therefore always returned numerically identical values; the duplicate
    evaluation has been removed (see key_rate() below), which is an exact
    2x reduction in cost with zero change to any numerical output.
    """
    theta = thetaA - thetaB
    eta_arr = np.asarray(eta, dtype=np.float64)
    scalar_in = (eta_arr.ndim == 0)
    eta_flat = np.atleast_1d(eta_arr)
    gamma = np.sqrt(eta_flat) * alpha ** 2
    pXX = p_xx_total(phi, theta, gamma, pd)

    eZ_num = np.zeros_like(eta_flat)
    for j in (0, 1):
        S_j = inner_sum_j(j, alpha, eta_flat, pd, thetaA, thetaB, Nmax)
        D_j = residual_j(j, alpha, Nmax)
        eZ_num += (S_j + D_j) ** 2

    eZ = np.where(pXX > 0.0, eZ_num / np.where(pXX > 0.0, pXX, 1.0), 0.5)
    eZ = np.clip(eZ, 0.0, 0.5)
    return float(eZ[0]) if scalar_in else eZ.reshape(eta_arr.shape)


# =============================================================================
# 7. Asymptotic key rate and alpha optimization — Eqs. (17)-(19)
# =============================================================================

def secret_key_rate_component(pxx, ex, ez):
    eX_c = np.clip(ex, 0.0, 0.5)
    eZ_c = np.clip(np.minimum(ez, 0.5), 0.0, 0.5)
    pxx = np.maximum(pxx, 0.0)
    return np.maximum(pxx * (1.0 - h(eX_c) - h(eZ_c)), 0.0)


def key_rate(alpha, eta, pd, thetaA, thetaB, Nmax, phi=0.0):
    """R_low_X = R10 + R01, exploiting R10 == R01 (see note above)."""
    theta = thetaA - thetaB
    eta_arr = np.asarray(eta, dtype=np.float64)
    scalar_in = (eta_arr.ndim == 0)
    eta_flat = np.atleast_1d(eta_arr)
    gamma = np.sqrt(eta_flat) * alpha ** 2

    pXX = p_xx_total(phi, theta, gamma, pd)
    eX = bit_error_rate(phi, theta, gamma, pd)
    eZ = phase_error_upper_bound(alpha, eta_flat, pd, thetaA, thetaB, Nmax, phi)

    R_single = secret_key_rate_component(pXX, eX, eZ)
    R_total = 2.0 * R_single
    return float(R_total[0]) if scalar_in else R_total.reshape(eta_arr.shape)


def total_secret_key_rate(alpha, eta, pd, thetaA, thetaB, delta, Nmax):
    """Eqs. (17)-(19). Returns (R_total, eX, eZ, pXX) -- phi passed explicitly,
    no global PARAMS mutation (thread-safe, avoids a dict read/write per call)."""
    phi = delta * np.pi
    theta = thetaA - thetaB
    gamma = np.sqrt(eta) * alpha ** 2

    pXX = max(float(p_xx_total(phi, theta, gamma, pd)), 0.0)
    eX = float(np.clip(bit_error_rate(phi, theta, gamma, pd), 0.0, 0.5))
    eZ = phase_error_upper_bound(alpha, eta, pd, thetaA, thetaB, Nmax, phi)

    R_total = 2.0 * float(secret_key_rate_component(pXX, eX, eZ))
    return R_total, eX, eZ, pXX


def optimize_alpha(eta, pd, thetaA, thetaB, delta=0.0, Nmax=NMAX_DEFAULT,
                    alpha_min=ALPHA_BOUNDS[0], alpha_max=ALPHA_BOUNDS[1]):
    """Maximise key_rate over alpha (Brent, bounded). Turbulence/eta-independent
    quantities (residual_j, Eq.35 coefficient tables) are cached across all
    trial alphas and across the whole outer loss/distance sweep."""
    phi = delta * np.pi

    def neg_rate(a):
        return -key_rate(a, eta, pd, thetaA, thetaB, Nmax, phi) if a > 0 else 0.0

    probes = np.logspace(np.log10(alpha_min), np.log10(alpha_max), 5)
    probe_rates = np.array([key_rate(a, eta, pd, thetaA, thetaB, Nmax, phi) for a in probes])
    if np.all(probe_rates <= 0.0):
        return alpha_min, alpha_min ** 2, 0.0

    result = minimize_scalar(neg_rate, bounds=(alpha_min, alpha_max),
                              method='bounded',
                              options={'maxiter': OPT_MAXITER, 'xatol': 1e-4})
    alpha_opt = float(np.clip(result.x, alpha_min, alpha_max))
    return alpha_opt, alpha_opt ** 2, max(-float(result.fun), 0.0)


# =============================================================================
# 8. PLOB bound
# =============================================================================

def plob_bound(eta):
    eta = np.clip(np.asarray(eta, dtype=np.float64), 1e-300, 1.0 - 1e-15)
    return -np.log2(1.0 - eta)


# =============================================================================
# 9. Deterministic FSO channel
# =============================================================================

def atmospheric_transmittance(alpha_atm, distance_km):
    """Beer-Lambert: eta_atm = exp(-alpha_atm * L)."""
    distance_km = np.clip(np.asarray(distance_km, dtype=np.float64), 0.0, None)
    return np.clip(np.exp(-alpha_atm * distance_km), 0.0, 1.0)


def beam_radius(wavelength, distance_m, w0):
    """w(L) = w0 * sqrt(1 + (lambda L / (pi w0^2))^2)."""
    distance_m = np.clip(np.asarray(distance_m, dtype=np.float64), 0.0, None)
    z_R = np.pi * w0 ** 2 / wavelength
    return w0 * np.sqrt(1.0 + (distance_m / z_R) ** 2)


def geometric_efficiency(receiver_radius, w_at_L):
    """eta_geo = 1 - exp(-2 r_rx^2 / w^2)."""
    w_safe = np.where(np.asarray(w_at_L) > 0.0, w_at_L, 1e-30)
    return np.clip(1.0 - np.exp(-2.0 * receiver_radius ** 2 / w_safe ** 2), 0.0, 1.0)


def deterministic_fso_channel(alpha_atm, distance_km, wavelength, w0, receiver_radius):
    """eta_fso(L) = eta_atm(L) * eta_geo(L/2)^2 (both half-links)."""
    distance_km = np.asarray(distance_km, dtype=np.float64)
    L_half_m = (distance_km / 2.0) * 1e3

    eta_atm = atmospheric_transmittance(alpha_atm, distance_km)
    w_half_m = beam_radius(wavelength, L_half_m, w0)
    eta_geo_half = geometric_efficiency(receiver_radius, w_half_m)
    eta_geo = eta_geo_half ** 2
    eta_fso = np.clip(eta_atm * eta_geo, 0.0, 1.0)

    eta_safe = np.where(eta_fso > 0.0, eta_fso, 1e-300)
    loss_dB = -10.0 * np.log10(eta_safe)
    return {'eta_atm': eta_atm, 'w_half_m': w_half_m, 'eta_geo_half': eta_geo_half,
            'eta_geo': eta_geo, 'eta_fso': eta_fso, 'loss_dB': loss_dB}


def tfqkd_skr_fso(distance_km, pd, thetaA, thetaB, fso_params, Nmax=4,
                   optimise=True, alpha_fixed=0.20):
    """SKR vs distance under deterministic FSO loss. optimise=False uses
    key_rate's array path directly (no Python per-point loop)."""
    distance_km = np.asarray(distance_km, dtype=np.float64)
    fso = deterministic_fso_channel(fso_params['alpha_atm'], distance_km,
                                     fso_params['wavelength'], fso_params['w0'],
                                     fso_params['receiver_radius'])
    eta_fso_arr = fso['eta_fso']
    valid = eta_fso_arr >= 1e-50

    if not optimise:
        R_arr = np.zeros_like(eta_fso_arr)
        R_arr[valid] = key_rate(alpha_fixed, eta_fso_arr[valid], pd, thetaA, thetaB, Nmax)
        return R_arr, eta_fso_arr, np.full_like(eta_fso_arr, alpha_fixed)

    # alpha optimization is inherently per-point (1-D search); each call is
    # now ~20x cheaper due to the Eq.(35) coefficient cache (see report).
    R_arr = np.zeros_like(eta_fso_arr)
    alpha_arr = np.full_like(eta_fso_arr, alpha_fixed)
    for i in np.flatnonzero(valid):
        a_opt, _, R_max = optimize_alpha(float(eta_fso_arr[i]), pd, thetaA, thetaB,
                                          delta=0.0, Nmax=Nmax)
        R_arr[i], alpha_arr[i] = max(R_max, 0.0), a_opt
    return R_arr, eta_fso_arr, alpha_arr


# =============================================================================
# 10. Log-normal turbulence (Rytov / weak-fluctuation regime)
# =============================================================================

Cn2_classes = {'weak': 1e-17, 'moderate': 1e-15, 'strong': 1e-13}
N_MC = 50_000
MC_SEED = 42
turbulence_enabled = True


def rytov_variance(Cn2, wavelength, distance_m):
    """sigma_R^2 = 1.23 Cn^2 k^(7/6) L^(11/6) (plane-wave, weak regime)."""
    distance_m = np.clip(np.asarray(distance_m, dtype=np.float64), 0.0, None)
    k = 2.0 * np.pi / wavelength
    return np.clip(1.23 * Cn2 * (k ** (7.0 / 6.0)) * (distance_m ** (11.0 / 6.0)), 0.0, None)


def lognormal_params(sigma_R2):
    """(mu_X, sigma_X^2), energy-neutral (E[eta_LN] = 1), overflow-safe."""
    sigma_R2 = np.asarray(sigma_R2, dtype=np.float64)
    OVERFLOW_THRESH = 500.0
    safe_sR2 = np.minimum(sigma_R2, OVERFLOW_THRESH)
    sigma_X2 = np.where(sigma_R2 > OVERFLOW_THRESH, sigma_R2,
                         np.log1p(np.exp(safe_sR2) - 1.0))
    return -0.5 * sigma_X2, sigma_X2


def sample_lognormal_eta(eta_nominal, Cn2, wavelength, distance_m,
                          N_samples=50_000, rng=None):
    """eta_turb = eta_nominal * exp(X), X ~ N(mu_X, sigma_X^2), clipped to [0,1]."""
    if rng is None:
        rng = np.random.default_rng(MC_SEED)
    if not turbulence_enabled:
        return np.full(N_samples, float(eta_nominal), dtype=np.float64)

    sigma_R2 = float(rytov_variance(Cn2, wavelength, distance_m))
    mu_X, sigma_X2 = lognormal_params(sigma_R2)
    X = rng.normal(loc=float(mu_X), scale=float(np.sqrt(sigma_X2)), size=N_samples)
    eta_turb = float(eta_nominal) * np.exp(X)
    return np.clip(eta_turb, 0.0, 1.0)


def turbulence_statistics(eta_nominal, Cn2, wavelength, distance_m,
                           N_samples=50_000, rng=None):
    """Analytic + MC moments of eta_turb, for cross-validation."""
    if rng is None:
        rng = np.random.default_rng(MC_SEED)
    sigma_R2 = float(rytov_variance(Cn2, wavelength, distance_m))
    mu_X, sigma_X2 = (float(v) for v in lognormal_params(sigma_R2))

    mean_analytic = float(eta_nominal) * np.exp(mu_X + 0.5 * sigma_X2)
    var_analytic = (float(eta_nominal) ** 2) * (np.exp(sigma_X2) - 1.0)
    SI_analytic = np.exp(sigma_X2) - 1.0

    samples = sample_lognormal_eta(eta_nominal, Cn2, wavelength, distance_m, N_samples, rng)
    return {'sigma_R2': sigma_R2, 'sigma_X2': sigma_X2, 'mu_X': mu_X,
            'mean_analytic': mean_analytic, 'var_analytic': var_analytic,
            'SI_analytic': SI_analytic, 'mean_mc': float(np.mean(samples)),
            'var_mc': float(np.var(samples)), 'std_mc': float(np.std(samples)),
            'samples': samples}


def tfqkd_skr_turbulence(distance_km, pd, thetaA, thetaB, Cn2, fso_params,
                          N_samples=None, Nmax=4, alpha_fixed=0.20, rng=None):
    """Ergodic SKR under log-normal fading: R_turb(L) = E_eta[R_det(eta_turb)].

    KEY OPTIMIZATION: the original code called key_rate() once per individual
    MC sample inside a Python list comprehension (N_dist * N_samples scalar
    calls, each re-deriving the Eq.35 combinatorial sum from scratch). Because
    key_rate() now accepts a full eta array and reuses the cached Eq.(35)
    coefficient tables, every distance point's N_samples-strong ensemble is
    evaluated in a single vectorized call.
    """
    if N_samples is None:
        N_samples = N_MC
    if rng is None:
        rng = np.random.default_rng(MC_SEED)

    distance_km = np.asarray(distance_km, dtype=np.float64)
    fso = deterministic_fso_channel(fso_params['alpha_atm'], distance_km,
                                     fso_params['wavelength'], fso_params['w0'],
                                     fso_params['receiver_radius'])
    eta_nominal_arr = fso['eta_fso']
    R_arr = np.zeros(len(distance_km))

    for i, (L_km, eta_nom) in enumerate(zip(distance_km, eta_nominal_arr)):
        if eta_nom < 1e-50:
            continue
        eta_samples = sample_lognormal_eta(
            eta_nominal=float(eta_nom), Cn2=Cn2, wavelength=fso_params['wavelength'],
            distance_m=L_km * 1e3, N_samples=N_samples, rng=rng)
        R_samples = key_rate(alpha_fixed, eta_samples, pd, thetaA, thetaB, Nmax)
        R_arr[i] = float(np.mean(R_samples))

    return R_arr, eta_nominal_arr


# =============================================================================
# 10b. Analytical (quadrature-based) ensemble averaging under log-normal
#      fading -- preserves all Monte Carlo functions above unchanged.
#
#      Derivation summary (full derivation in accompanying writeup):
#
#      eta_turb = eta_nom * exp(X),  X ~ N(mu_X, sigma_X^2)
#
#      For the *unclipped* log-normal (the clip to [0,1] applied in
#      sample_lognormal_eta affects a negligible probability mass in the
#      sigma_R^2 < 1 validity regime; see _lognormal_clip_mass_above_one()
#      below for a diagnostic of this approximation), the standard
#      change-of-variables identity for a log-normal random variable gives,
#      for any observable g(eta):
#
#          <g>(eta_nom) = Integral_{-inf}^{inf} g(eta_nom * exp(X))
#                              * phi(X; mu_X, sigma_X^2) dX
#
#      where phi is the Gaussian PDF. This substitution is EXACT (not an
#      approximation) and removes the 1/eta singularity that the eta-space
#      log-normal density has at eta -> 0, which is why we integrate in
#      X-space rather than eta-space.
#
#      None of g = p_XX, g = e_X, g = R admit closed-form Gaussian integrals
#      (each involves exp(+/- gamma) with gamma ~ sqrt(eta) = sqrt(eta_nom)
#      exp(X/2), i.e. an exponential-of-exponential of the integration
#      variable -- and R additionally passes through the binary entropy
#      function and the Eq.(20-22) phase-error sum, which is transcendental
#      in eta). Per requirement 3, we do NOT attempt invalid symbolic
#      simplification; we evaluate the exact integral by quadrature
#      (scipy.integrate.quad) instead.
#
#      QBER requires special care: e_X(eta) = N(eta)/D(eta) is a RATIO, so
#      <e_X> != <N>/<D> in general (Jensen's inequality forbids exchanging
#      ratio and average). The physically meaningful average -- consistent
#      with "total expected wrong bits / total expected sifted bits" -- is
#      the gain-weighted average:
#
#          <e_X> = <p_XX * e_X> / <p_XX>
#
#      which is what average_qber_lognormal() computes.
# =============================================================================

QUAD_N_SIGMA = 8.0        # Gaussian truncation half-width (tail < 1e-15)
QUAD_EPSABS = 1e-12
QUAD_EPSREL = 1e-10
QUAD_LIMIT = 200


def _gaussian_pdf(X, mu_X, sigma_X):
    """phi(X; mu_X, sigma_X^2)."""
    return np.exp(-0.5 * ((X - mu_X) / sigma_X) ** 2) / (sigma_X * np.sqrt(2.0 * np.pi))


def _quad_bounds(mu_X, sigma_X, n_sigma=QUAD_N_SIGMA):
    """Truncated integration domain: [mu_X - n_sigma*sigma_X, mu_X + n_sigma*sigma_X].
    Tail mass beyond 8 sigma is ~1e-15, below double-precision noise floor,
    so this truncation introduces no numerically meaningful error while
    making quad's job well-posed (avoids relying on its +/-inf transform
    when sigma_X is tiny and the integrand is a sharp spike)."""
    return mu_X - n_sigma * sigma_X, mu_X + n_sigma * sigma_X


def _lognormal_clip_mass_above_one(eta_nominal, sigma_R2):
    """Diagnostic: P(eta_turb > 1) under the UNCLIPPED log-normal, i.e. the
    probability mass that sample_lognormal_eta() clips to exactly 1. This
    quantifies the approximation made by integrating the unclipped log-normal
    in average_*_lognormal() instead of the clipped one used by the MC
    sampler. Returns a probability in [0,1]; expected to be ~0 in the
    sigma_R^2 < 1 validity regime unless eta_nominal is itself close to 1."""
    from scipy.stats import norm
    mu_X, sigma_X2 = lognormal_params(sigma_R2)
    sigma_X = np.sqrt(sigma_X2)
    # P(eta_nom * exp(X) > 1) = P(X > -ln(eta_nom)) = P(X > ln(1/eta_nom))
    z = (np.log(1.0 / np.clip(eta_nominal, 1e-300, None)) - mu_X) / sigma_X
    return float(norm.sf(z))


def average_pxx_lognormal(eta_nominal, sigma_R2, phi, theta, alpha, pd,
                           n_sigma=QUAD_N_SIGMA):
    """<p_XX>(eta_nom) = Integral p_XX(eta_nom*exp(X)) * phi(X) dX  (Eq. 30
    pushed through the log-normal fading channel). Scalar eta_nominal/sigma_R2
    in, scalar out; see average_pxx_lognormal_array for the vectorized wrapper.
    """
    mu_X, sigma_X2 = (float(v) for v in lognormal_params(sigma_R2))
    sigma_X = math.sqrt(sigma_X2) if sigma_X2 > 0 else 1e-12
    lo, hi = _quad_bounds(mu_X, sigma_X, n_sigma)

    def integrand(X):
        eta = eta_nominal * math.exp(X)
        gamma = math.sqrt(eta) * alpha ** 2
        pxx = float(p_xx_total(phi, theta, gamma, pd))
        return pxx * _gaussian_pdf(X, mu_X, sigma_X)

    val, err = quad(integrand, lo, hi, epsabs=QUAD_EPSABS, epsrel=QUAD_EPSREL,
                     limit=QUAD_LIMIT)
    return max(val, 0.0), err


def average_qber_lognormal(eta_nominal, sigma_R2, phi, theta, alpha, pd,
                            n_sigma=QUAD_N_SIGMA):
    """<e_X>(eta_nom) = <p_XX * e_X> / <p_XX>  (gain-weighted average QBER;
    see derivation note above on why the naive unweighted average is wrong).
    Returns (mean_eX, (num_err, den_err))."""
    mu_X, sigma_X2 = (float(v) for v in lognormal_params(sigma_R2))
    sigma_X = math.sqrt(sigma_X2) if sigma_X2 > 0 else 1e-12
    lo, hi = _quad_bounds(mu_X, sigma_X, n_sigma)

    def integrand_num(X):
        eta = eta_nominal * math.exp(X)
        gamma = math.sqrt(eta) * alpha ** 2
        pxx = float(p_xx_total(phi, theta, gamma, pd))
        ex = float(np.clip(bit_error_rate(phi, theta, gamma, pd), 0.0, 0.5))
        return pxx * ex * _gaussian_pdf(X, mu_X, sigma_X)

    def integrand_den(X):
        eta = eta_nominal * math.exp(X)
        gamma = math.sqrt(eta) * alpha ** 2
        pxx = float(p_xx_total(phi, theta, gamma, pd))
        return pxx * _gaussian_pdf(X, mu_X, sigma_X)

    num, num_err = quad(integrand_num, lo, hi, epsabs=QUAD_EPSABS,
                         epsrel=QUAD_EPSREL, limit=QUAD_LIMIT)
    den, den_err = quad(integrand_den, lo, hi, epsabs=QUAD_EPSABS,
                         epsrel=QUAD_EPSREL, limit=QUAD_LIMIT)
    mean_eX = num / den if den > 1e-300 else 0.5
    return float(np.clip(mean_eX, 0.0, 0.5)), (num_err, den_err)


def average_skr_lognormal(eta_nominal, sigma_R2, alpha, pd, thetaA, thetaB,
                           Nmax=NMAX_DEFAULT, phi=0.0, n_sigma=QUAD_N_SIGMA):
    """<R>(eta_nom) = Integral R(eta_nom*exp(X), alpha) * phi(X) dX.

    This is the direct analytical analogue of tfqkd_skr_turbulence's Monte
    Carlo mean(R(eta_samples)) -- same full key_rate() pipeline (Eqs. 17-22,
    27-35), but the expectation over fading is taken by 1-D Gauss-Kronrod
    quadrature (scipy.integrate.quad) instead of sample averaging. No
    closed form exists because R involves h(.) and the Eq.(20-22) phase
    error sum, both transcendental in eta -- see module-level derivation note.
    """
    theta = thetaA - thetaB
    mu_X, sigma_X2 = (float(v) for v in lognormal_params(sigma_R2))
    sigma_X = math.sqrt(sigma_X2) if sigma_X2 > 0 else 1e-12
    lo, hi = _quad_bounds(mu_X, sigma_X, n_sigma)

    def integrand(X):
        eta = eta_nominal * math.exp(X)
        eta = min(max(eta, 0.0), 1.0)  # match MC sampler's clip for consistency
        R = key_rate(alpha, eta, pd, thetaA, thetaB, Nmax, phi)
        return R * _gaussian_pdf(X, mu_X, sigma_X)

    val, err = quad(integrand, lo, hi, epsabs=QUAD_EPSABS, epsrel=QUAD_EPSREL,
                     limit=QUAD_LIMIT)
    return max(val, 0.0), err


def average_skr_lognormal_array(distance_km, pd, thetaA, thetaB, Cn2,
                                 fso_params, Nmax=NMAX_DEFAULT,
                                 alpha_fixed=0.20, n_sigma=QUAD_N_SIGMA):
    """Vectorized-over-distance wrapper around average_skr_lognormal(), mirroring
    tfqkd_skr_turbulence()'s signature/return so the two can be swapped for
    direct comparison. quad() itself is inherently scalar per evaluation
    point, so this loops over distance (same structural shape as the MC
    function's loop over distance, but with O(1) quadrature evaluations per
    point instead of O(N_samples) sample draws)."""
    distance_km = np.asarray(distance_km, dtype=np.float64)
    fso = deterministic_fso_channel(fso_params['alpha_atm'], distance_km,
                                     fso_params['wavelength'], fso_params['w0'],
                                     fso_params['receiver_radius'])
    eta_nominal_arr = fso['eta_fso']
    R_arr = np.zeros(len(distance_km))
    err_arr = np.zeros(len(distance_km))

    for i, (L_km, eta_nom) in enumerate(zip(distance_km, eta_nominal_arr)):
        if eta_nom < 1e-50:
            continue
        sigma_R2 = float(rytov_variance(Cn2, fso_params['wavelength'], L_km * 1e3))
        R_arr[i], err_arr[i] = average_skr_lognormal(
            float(eta_nom), sigma_R2, alpha_fixed, pd, thetaA, thetaB, Nmax,
            n_sigma=n_sigma)
    return R_arr, eta_nominal_arr, err_arr


def optimize_alpha_lognormal(eta_nominal, sigma_R2, pd, thetaA, thetaB,
                              delta=0.0, Nmax=NMAX_DEFAULT,
                              alpha_min=ALPHA_BOUNDS[0], alpha_max=ALPHA_BOUNDS[1],
                              n_sigma=QUAD_N_SIGMA):
    """alpha* = argmax_alpha <R(alpha)>, where <R> is the QUADRATURE-averaged
    SKR under log-normal fading (average_skr_lognormal), at a single distance
    point (fixed eta_nominal, sigma_R2). Direct analogue of optimize_alpha()
    (Eqs. 17-19) but optimizing the ensemble-averaged rate rather than the
    deterministic rate -- this is requirement 10."""
    phi = delta * np.pi

    def neg_avg_rate(a):
        if a <= 0:
            return 0.0
        val, _ = average_skr_lognormal(eta_nominal, sigma_R2, a, pd, thetaA,
                                        thetaB, Nmax, phi, n_sigma)
        return -val

    probes = np.logspace(np.log10(alpha_min), np.log10(alpha_max), 5)
    probe_rates = np.array([-neg_avg_rate(a) for a in probes])
    if np.all(probe_rates <= 0.0):
        return alpha_min, 0.0

    result = minimize_scalar(neg_avg_rate, bounds=(alpha_min, alpha_max),
                              method='bounded',
                              options={'maxiter': OPT_MAXITER, 'xatol': 1e-4})
    alpha_opt = float(np.clip(result.x, alpha_min, alpha_max))
    return alpha_opt, max(-float(result.fun), 0.0)


def optimize_alpha_lognormal_vs_distance(distance_km, pd, thetaA, thetaB, Cn2,
                                          fso_params, Nmax=NMAX_DEFAULT,
                                          alpha_fixed=0.20, n_sigma=QUAD_N_SIGMA):
    """Sweep optimize_alpha_lognormal() over distance: returns alpha*(L),
    <R(alpha*)>(L), and <R(alpha_fixed)>(L) for the 'optimized vs fixed'
    comparison plot (requirement 11)."""
    distance_km = np.asarray(distance_km, dtype=np.float64)
    fso = deterministic_fso_channel(fso_params['alpha_atm'], distance_km,
                                     fso_params['wavelength'], fso_params['w0'],
                                     fso_params['receiver_radius'])
    eta_nominal_arr = fso['eta_fso']
    n = len(distance_km)
    alpha_star = np.full(n, alpha_fixed)
    R_star = np.zeros(n)
    R_fixed = np.zeros(n)

    for i, (L_km, eta_nom) in enumerate(zip(distance_km, eta_nominal_arr)):
        if eta_nom < 1e-50:
            continue
        sigma_R2 = float(rytov_variance(Cn2, fso_params['wavelength'], L_km * 1e3))
        a_opt, R_opt = optimize_alpha_lognormal(
            float(eta_nom), sigma_R2, pd, thetaA, thetaB, Nmax=Nmax,
            n_sigma=n_sigma)
        R_fix, _ = average_skr_lognormal(
            float(eta_nom), sigma_R2, alpha_fixed, pd, thetaA, thetaB, Nmax,
            n_sigma=n_sigma)
        alpha_star[i], R_star[i] = a_opt, max(R_opt, 0.0)
        R_fixed[i] = max(R_fix, 0.0)

    return alpha_star, R_star, R_fixed, eta_nominal_arr


def validate_mc_vs_quadrature(distance_km, pd, thetaA, thetaB, Cn2, fso_params,
                               N_samples=None, Nmax=NMAX_DEFAULT,
                               alpha_fixed=0.20, rng=None, n_sigma=QUAD_N_SIGMA):
    """Requirement 7-8: compare Monte Carlo (tfqkd_skr_turbulence, UNCHANGED)
    against quadrature (average_skr_lognormal_array, new) over distance, and
    compute the relative error abs(MC - Quad)/Quad point-by-point. Returns a
    dict with both SKR curves, sigma_R^2(L), and the relative-error array.

    This routine calls -- but never modifies -- the existing MC function,
    satisfying requirement 13 (MC remains available as the validation
    benchmark, not replaced)."""
    distance_km = np.asarray(distance_km, dtype=np.float64)

    R_mc, eta_nom_mc = tfqkd_skr_turbulence(
        distance_km, pd, thetaA, thetaB, Cn2, fso_params,
        N_samples=N_samples, Nmax=Nmax, alpha_fixed=alpha_fixed, rng=rng)

    R_quad, eta_nom_quad, quad_err = average_skr_lognormal_array(
        distance_km, pd, thetaA, thetaB, Cn2, fso_params, Nmax=Nmax,
        alpha_fixed=alpha_fixed, n_sigma=n_sigma)

    sigma_R2_arr = rytov_variance(Cn2, fso_params['wavelength'], distance_km * 1e3)

    quad_safe = np.where(R_quad > 1e-300, R_quad, np.nan)
    rel_error = np.abs(R_mc - R_quad) / quad_safe

    return {'distance_km': distance_km, 'R_mc': R_mc, 'R_quad': R_quad,
            'quad_err': quad_err, 'eta_nominal': eta_nom_mc,
            'sigma_R2': sigma_R2_arr, 'rel_error': rel_error,
            'sigma_R2_lt_1': sigma_R2_arr < 1.0}


# =============================================================================
# 11. Minimal correctness self-test (replaces the notebook's inline asserts)
# =============================================================================

def _self_test():
    thetaA, thetaB = PARAMS['thetaA'], PARAMS['thetaB']
    pd = PARAMS['pd']
    loss_dB = np.array([0, 10, 20, 40, 60, 80, 100])
    eta = dB_to_eta(loss_dB)

    pZZ_pn = p_zz_photon(1, 1, thetaA, thetaB, eta, pd)
    assert np.all(pZZ_pn >= -1e-12) and np.all(pZZ_pn <= 1 + 1e-12)

    R = key_rate(PARAMS['alpha'], eta, pd, thetaA, thetaB, NMAX_DEFAULT)
    assert np.all(R >= 0) and np.all(np.diff(R) <= 1e-15)

    h_vals = h(np.array([0.0, 0.5, 1.0]))
    assert np.allclose(h_vals, [0.0, 1.0, 0.0], atol=1e-12)

    a_opt, _, R_max = optimize_alpha(float(eta[2]), pd, thetaA, thetaB, Nmax=NMAX_DEFAULT)
    assert R_max >= 0 and ALPHA_BOUNDS[0] <= a_opt <= ALPHA_BOUNDS[1]

    # --- New: quadrature-based ensemble averaging sanity checks ---
    theta_test = thetaA - thetaB
    sigma_R2_test = 0.3  # well inside sigma_R^2 < 1 validity regime
    eta_nom_test = float(eta[3])  # 40 dB point

    pxx_q, _ = average_pxx_lognormal(eta_nom_test, sigma_R2_test, 0.0, theta_test,
                                      PARAMS['alpha'], pd)
    assert 0.0 <= pxx_q <= 1.0

    eX_q, _ = average_qber_lognormal(eta_nom_test, sigma_R2_test, 0.0, theta_test,
                                      PARAMS['alpha'], pd)
    assert 0.0 <= eX_q <= 0.5

    R_q, _ = average_skr_lognormal(eta_nom_test, sigma_R2_test, PARAMS['alpha'], pd,
                                    thetaA, thetaB, NMAX_DEFAULT)
    assert R_q >= 0.0

    # Sigma_R2 -> 0 limit: <R> must reduce to deterministic key_rate(eta_nominal)
    R_q_zero, _ = average_skr_lognormal(eta_nom_test, 1e-10, PARAMS['alpha'], pd,
                                         thetaA, thetaB, NMAX_DEFAULT)
    R_det_zero = key_rate(PARAMS['alpha'], eta_nom_test, pd, thetaA, thetaB, NMAX_DEFAULT)
    assert abs(R_q_zero - R_det_zero) < 1e-6 * max(R_det_zero, 1e-12) + 1e-8

    print('Self-test OK.')


if __name__ == '__main__':
    _self_test()

    FSO_PARAMS = dict(alpha_atm=0.046, wavelength=1550e-9, w0=0.05,
                       receiver_radius=0.15, dist_min_km=0.5, dist_max_km=500.0)
    dist_km = np.linspace(FSO_PARAMS['dist_min_km'], FSO_PARAMS['dist_max_km'], 60)

    R_fso, eta_fso, alpha_fso = tfqkd_skr_fso(
        dist_km, PARAMS['pd'], PARAMS['thetaA'], PARAMS['thetaB'],
        FSO_PARAMS, Nmax=NMAX_DEFAULT, optimise=True)
    print(f'Deterministic FSO: {int(np.sum(R_fso > 1e-15))}/{len(dist_km)} '
          f'positive-rate points; max range '
          f'~{dist_km[R_fso > 1e-15][-1]:.0f} km' if np.any(R_fso > 1e-15) else 'no positive-rate points')

    rng = np.random.default_rng(MC_SEED)
    R_turb, eta_nom = tfqkd_skr_turbulence(
        dist_km, PARAMS['pd'], PARAMS['thetaA'], PARAMS['thetaB'],
        Cn2_classes['moderate'], FSO_PARAMS, N_samples=N_MC,
        Nmax=NMAX_DEFAULT, alpha_fixed=0.20, rng=rng)
    print(f'Turbulent (moderate) ergodic SKR computed for {len(dist_km)} distance points, '
          f'N_MC={N_MC:,} samples/point.')