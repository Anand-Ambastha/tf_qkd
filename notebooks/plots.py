# -*- coding: utf-8 -*-
"""Diagnostic plots for the optimized TF-QKD / FSO / turbulence engine.

Reuses tfqkd_fso_turbulence_optimized.py (fast, cached core) to regenerate
the full plot set from the original notebook (channel model, photon yields,
cat states, phase-error convergence, key rate vs loss, alpha optimisation,
deterministic FSO channel, log-normal turbulence). Each figure is saved as
a PNG in the working directory.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import tfqkd_fso_turbulence_optimized as M

plt.rcParams.update({'figure.dpi': 130, 'axes.grid': True, 'grid.alpha': 0.3,
                      'font.size': 10, 'lines.linewidth': 1.8})

OUT = './'
thetaA, thetaB = M.PARAMS['thetaA'], M.PARAMS['thetaB']
theta = thetaA - thetaB
pd = M.PARAMS['pd']
alpha0 = M.PARAMS['alpha']
phi0 = M.PARAMS['phi']
Nmax = M.NMAX_DEFAULT

loss_arr = np.linspace(0, 120, 300)
eta_arr = M.dB_to_eta(loss_arr)

# =============================================================================
# Fig 1 -- X-basis gain, bit-error rate, Z-basis gain
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ax = axes[0]
for pd_v, c in zip([1e-8, 1e-7, 1e-6, 1e-5], plt.cm.viridis(np.linspace(0, 0.85, 4))):
    gamma = np.sqrt(eta_arr) * alpha0**2
    pXX = M.p_xx_total(phi0, theta, gamma, pd_v)
    ax.semilogy(loss_arr, np.clip(pXX, 1e-16, None), color=c, label=f'$p_d=10^{{{int(np.log10(pd_v))}}}$')
ax.semilogy(loss_arr, 0.5*np.sqrt(eta_arr), 'k--', lw=1.2, label=r'$\frac{1}{2}\sqrt{\eta}$ (ref.)')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel(r'$p_{XX}$'); ax.legend(fontsize=8)
ax.set_title('X-basis total gain (Eq. 30)')

ax = axes[1]
for pd_v, c in zip([1e-8, 1e-7, 1e-6, 1e-5], plt.cm.viridis(np.linspace(0, 0.85, 4))):
    gamma = np.sqrt(eta_arr) * alpha0**2
    eX = M.bit_error_rate(phi0, theta, gamma, pd_v)
    ax.semilogy(loss_arr, np.clip(eX, 1e-9, None), color=c, label=f'$p_d=10^{{{int(np.log10(pd_v))}}}$')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel(r'$e_X$'); ax.legend(fontsize=8)
ax.set_title('X-basis bit-error rate (Eq. 31)')

ax = axes[2]
for b2, c in zip([0.001, 0.01, 0.1, 0.3], plt.cm.plasma(np.linspace(0, 0.85, 4))):
    b = np.sqrt(b2)
    pZZ = M.p_zz_beta(b, b, theta, eta_arr, pd)
    ax.semilogy(loss_arr, np.clip(pZZ, 1e-16, None), color=c, label=fr'$\beta^2={b2}$')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel(r'$p_{ZZ}$'); ax.legend(fontsize=8)
ax.set_title('Z-basis gain (Eq. 32)')

fig.suptitle('Channel model', y=1.02)
fig.tight_layout()
fig.savefig(f'{OUT}/fig1_channel_model.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 2 -- Photon-number yields + cat-state coefficients
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
for (nA, nB), c in zip([(0,0),(1,1),(2,2),(1,3),(3,3)], plt.cm.cool(np.linspace(0,0.9,5))):
    pzz = M.p_zz_photon(nA, nB, thetaA, thetaB, eta_arr, pd)
    ax.semilogy(loss_arr, np.clip(pzz, 1e-16, None), color=c, label=f'$(n_A,n_B)=({nA},{nB})$')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel(r'$p_{ZZ}(n_A,n_B)$'); ax.legend(fontsize=8)
ax.set_title('Photon-number yields (Eq. 34-35)')

ax = axes[1]
n_terms = 12
ns = np.arange(n_terms)
c_ev = M._cat_coeff_vec(0, alpha0, n_terms - 1)
c_od = M._cat_coeff_vec(1, alpha0, n_terms - 1)
ax.semilogy(ns, c_ev, 'o-', label=r'$c^{(0)}_{2n}$ (even)')
ax.semilogy(ns, c_od, 's-', label=r'$c^{(1)}_{2n+1}$ (odd)')
ax.set_xlabel('Fock sub-index $n$'); ax.set_ylabel('Coefficient')
ax.set_title(f'Cat-state coefficients (Eq. 12-13), $\\alpha={alpha0}$')
ax.legend()

fig.tight_layout()
fig.savefig(f'{OUT}/fig2_yields_catstates.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 3 -- Phase-error bound + residual convergence vs Nmax
# =============================================================================
loss_fixed = 40.0
eta_fixed = float(M.dB_to_eta(loss_fixed))
nmax_vals = list(range(0, 13, 2))
eZ_vs_Nmax, D0_vs_Nmax, D1_vs_Nmax = [], [], []
for Nm in nmax_vals:
    eZ_vs_Nmax.append(M.phase_error_upper_bound(alpha0, eta_fixed, pd, thetaA, thetaB, Nm))
    D0_vs_Nmax.append(M.residual_j(0, alpha0, Nm))
    D1_vs_Nmax.append(M.residual_j(1, alpha0, Nm))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
ax.plot(nmax_vals, eZ_vs_Nmax, 'o-', ms=6)
ax.set_xlabel('$N_{max}$'); ax.set_ylabel('$e_Z^{upp}$')
ax.set_title(f'Phase-error bound vs $N_{{max}}$ ({loss_fixed:.0f} dB)')

ax = axes[1]
ax.semilogy(nmax_vals, np.maximum(D0_vs_Nmax, 1e-16), 'o-', label=r'$\Delta_0$ (even)')
ax.semilogy(nmax_vals, np.maximum(D1_vs_Nmax, 1e-16), 's-', label=r'$\Delta_1$ (odd)')
ax.set_xlabel('$N_{max}$'); ax.set_ylabel('Residual'); ax.legend()
ax.set_title('Truncation residuals (Eq. 22)')

fig.tight_layout()
fig.savefig(f'{OUT}/fig3_phase_error_convergence.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 4 -- Key rate vs loss (fixed alpha) + PLOB bound
# =============================================================================
R_fixed = M.key_rate(alpha0, eta_arr, pd, thetaA, thetaB, Nmax)
PLOB = M.plob_bound(eta_arr)

fig, ax = plt.subplots(figsize=(8, 5))
pos = R_fixed > 1e-15
ax.semilogy(loss_arr[pos], R_fixed[pos], 'C0-', lw=2, label=f'TF-QKD ($\\alpha={alpha0}$ fixed)')
ax.semilogy(loss_arr, PLOB, 'k--', lw=1.5, label='PLOB bound')
ax.semilogy(loss_arr, np.sqrt(eta_arr), 'k:', lw=1.0, label=r'$\sqrt{\eta}$ (ref.)')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel('Secret key rate [bits/pulse]')
ax.set_ylim(1e-8, 2); ax.legend()
ax.set_title(f'TF-QKD key rate vs loss  ($N_{{max}}={Nmax}$, $p_d={pd:.0e}$)')
fig.tight_layout()
fig.savefig(f'{OUT}/fig4_key_rate_vs_loss.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 5 -- Alpha-optimized key rate vs loss + optimal alpha
# =============================================================================
loss_opt_sweep = np.linspace(0, 100, 50)
eta_opt_sweep = M.dB_to_eta(loss_opt_sweep)
alpha_opt_arr, R_opt_arr = [], []
for e in eta_opt_sweep:
    a, _, R = M.optimize_alpha(float(e), pd, thetaA, thetaB, Nmax=Nmax)
    alpha_opt_arr.append(a); R_opt_arr.append(R)
alpha_opt_arr = np.array(alpha_opt_arr); R_opt_arr = np.array(R_opt_arr)
PLOB_opt = M.plob_bound(eta_opt_sweep)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
pos = R_opt_arr > 1e-15
ax.semilogy(loss_opt_sweep[pos], R_opt_arr[pos], 'C1-', lw=2, label='TF-QKD ($\\alpha$ optimized)')
ax.semilogy(loss_arr[R_fixed > 1e-15], R_fixed[R_fixed > 1e-15], 'C0--', lw=1.3, label=f'$\\alpha={alpha0}$ fixed')
ax.semilogy(loss_opt_sweep, PLOB_opt, 'k--', lw=1.5, label='PLOB bound')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel('Secret key rate [bits/pulse]')
ax.set_ylim(1e-8, 2); ax.legend(fontsize=8)
ax.set_title('Alpha-optimized vs fixed-alpha key rate')

ax = axes[1]
ax.plot(loss_opt_sweep[pos], alpha_opt_arr[pos], 'o-', ms=4, color='C2')
ax.set_xlabel('Loss [dB]'); ax.set_ylabel(r'Optimal $\alpha^*$')
ax.set_title('Optimal signal amplitude vs loss')

fig.tight_layout()
fig.savefig(f'{OUT}/fig5_alpha_optimization.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 6 -- Deterministic FSO channel components + SKR vs distance
# =============================================================================
FSO_PARAMS = dict(alpha_atm=0.046, wavelength=1550e-9, w0=0.05,
                   receiver_radius=0.15, dist_min_km=0.5, dist_max_km=500.0)
dist_km = np.linspace(FSO_PARAMS['dist_min_km'], FSO_PARAMS['dist_max_km'], 120)
fso = M.deterministic_fso_channel(FSO_PARAMS['alpha_atm'], dist_km, FSO_PARAMS['wavelength'],
                                   FSO_PARAMS['w0'], FSO_PARAMS['receiver_radius'])

R_fso, eta_fso, alpha_fso = M.tfqkd_skr_fso(
    dist_km, pd, thetaA, thetaB, FSO_PARAMS, Nmax=Nmax, optimise=True)
PLOB_fso = M.plob_bound(np.clip(eta_fso, 1e-300, None))

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
ax = axes[0]
ax.semilogy(dist_km, np.clip(fso['eta_atm'], 1e-16, 1), 'b-', label=r'$\eta_{atm}$')
ax.semilogy(dist_km, np.clip(fso['eta_geo'], 1e-16, 1), 'g-', label=r'$\eta_{geo}$')
ax.semilogy(dist_km, np.clip(fso['eta_fso'], 1e-16, 1), 'k-', lw=2.2, label=r'$\eta_{fso}$ (total)')
ax.set_xlabel('Distance [km]'); ax.set_ylabel('Transmittance'); ax.legend()
ax.set_title('FSO channel transmittance components')

ax = axes[1]
ax.plot(dist_km, fso['loss_dB'], 'k-', lw=2)
ax.set_xlabel('Distance [km]'); ax.set_ylabel('Total FSO loss [dB]')
ax.set_title('FSO loss vs distance')

ax = axes[2]
pos = R_fso > 1e-15
ax.semilogy(dist_km[pos], R_fso[pos], 'C0-', lw=2.2, label='TF-QKD SKR (FSO)')
ax.semilogy(dist_km, PLOB_fso, 'k--', lw=1.3, label='PLOB bound')
ax.set_xlabel('Distance [km]'); ax.set_ylabel('Secret key rate [bits/pulse]')
ax.set_ylim(1e-8, 2); ax.legend(fontsize=8)
max_range = dist_km[pos][-1] if np.any(pos) else 0
ax.set_title(f'TF-QKD SKR under deterministic FSO\n(max range ~{max_range:.0f} km)')

fig.tight_layout()
fig.savefig(f'{OUT}/fig6_fso_channel_skr.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 7 -- Turbulence: eta_turb distribution (weak/moderate/strong) at 50 km
# =============================================================================
rng = np.random.default_rng(M.MC_SEED)
dist_p = 50e3
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
colors = {'weak': 'C0', 'moderate': 'C1', 'strong': 'C2'}
for cls, Cn2 in M.Cn2_classes.items():
    stats = M.turbulence_statistics(0.3, Cn2, FSO_PARAMS['wavelength'], dist_p,
                                     N_samples=M.N_MC, rng=rng)
    samples = stats['samples']
    axes[0].hist(samples, bins=80, density=True, histtype='step', lw=1.8,
                 color=colors[cls], label=f'{cls} ($\\sigma_R^2$={stats["sigma_R2"]:.3f})')
    sorted_s = np.sort(samples)
    ecdf = np.arange(1, len(sorted_s)+1) / len(sorted_s)
    axes[1].plot(sorted_s[::200], ecdf[::200], color=colors[cls], lw=1.8, label=cls)

axes[0].set_xlabel(r'$\eta_{turb}$'); axes[0].set_ylabel('density'); axes[0].legend(fontsize=8)
axes[0].set_title(r'$\eta_{turb}$ distribution at L=50 km ($\eta_{nom}$=0.3)')
axes[1].set_xlabel(r'$\eta_{turb}$'); axes[1].set_ylabel('CDF'); axes[1].legend(fontsize=8)
axes[1].set_title('Empirical CDF')

fig.tight_layout()
fig.savefig(f'{OUT}/fig7_turbulence_distribution.png', bbox_inches='tight')
plt.close(fig)

# =============================================================================
# Fig 8 -- Ergodic TF-QKD SKR under turbulence vs deterministic FSO
# =============================================================================
dist_km_t = np.linspace(0.5, 300.0, 60)

fso_t = M.deterministic_fso_channel(
    FSO_PARAMS['alpha_atm'],
    dist_km_t,
    FSO_PARAMS['wavelength'],
    FSO_PARAMS['w0'],
    FSO_PARAMS['receiver_radius']
)

# Deterministic FSO
R_det_t, _, _ = M.tfqkd_skr_fso(
    dist_km_t,
    pd,
    thetaA,
    thetaB,
    FSO_PARAMS,
    Nmax=Nmax,
    optimise=False,
    alpha_fixed=0.20
)

# Turbulence cases
R_turb_results = {}
rng = np.random.default_rng(M.MC_SEED)

for cls, Cn2 in M.Cn2_classes.items():
    R_t, _ = M.tfqkd_skr_turbulence(
        dist_km_t,
        pd,
        thetaA,
        thetaB,
        Cn2,
        FSO_PARAMS,
        N_samples=M.N_MC,
        Nmax=Nmax,
        alpha_fixed=0.20,
        rng=rng
    )
    R_turb_results[cls] = R_t

fig, ax = plt.subplots(figsize=(10, 6))

# -------------------------------------------------------------------------
# Deterministic curve
# -------------------------------------------------------------------------
pos = R_det_t > 1e-15
ax.semilogy(
    dist_km_t[pos],
    R_det_t[pos],
    'k-',
    lw=2.5,
    label=r'Deterministic FSO ($\alpha=0.20$)'
)

# -------------------------------------------------------------------------
# Turbulence curves
# -------------------------------------------------------------------------
for cls, c in colors.items():
    R_t = R_turb_results[cls]
    pos_t = R_t > 1e-15

    ax.semilogy(
        dist_km_t[pos_t],
        R_t[pos_t],
        color=c,
        lw=1.8,
        label=f'Log-normal turbulence ({cls})'
    )

# -------------------------------------------------------------------------
# σR² = 1 validity boundaries
# -------------------------------------------------------------------------
validity_lines = {}

for cls, c in colors.items():

    Cn2 = M.Cn2_classes[cls]

    sigmaR2 = M.rytov_variance(
        Cn2,
        FSO_PARAMS['wavelength'],
        dist_km_t * 1e3
    )

    idx = np.where(sigmaR2 >= 1.0)[0]

    if len(idx):

        d_lim = dist_km_t[idx[0]]
        validity_lines[cls] = d_lim

        ax.axvline(
            d_lim,
            color=c,
            linestyle='--',
            linewidth=1.5,
            alpha=0.8,
            label=f'{cls}: $\\sigma_R^2=1$'
        )

        print(
            f"{cls.upper()} log-normal validity limit: "
            f"{d_lim:.2f} km"
        )

# -------------------------------------------------------------------------
# Shade invalid region for weak turbulence
# -------------------------------------------------------------------------
if 'weak' in validity_lines:

    weak_limit = validity_lines['weak']

    ax.axvspan(
        weak_limit,
        dist_km_t.max(),
        alpha=0.10,
        color='grey'
    )

    ax.text(
        weak_limit + 5,
        1e-6,
        'Log-normal model\noutside validity\n($\\sigma_R^2 > 1$)',
        fontsize=8
    )

# -------------------------------------------------------------------------
# Formatting
# -------------------------------------------------------------------------
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Ergodic secret key rate [bits/pulse]')

ax.set_ylim(1e-8, 2)

ax.legend(
    fontsize=8,
    loc='best'
)

ax.set_title(
    'TF-QKD SKR under deterministic FSO and log-normal turbulence\n'
    '(dashed lines indicate $\\sigma_R^2=1$ validity boundaries)'
)

fig.tight_layout()

fig.savefig(
    f'{OUT}/fig8_turbulence_skr.png',
    bbox_inches='tight'
)

plt.close(fig)

print("\n=== Range Check ===")

idx_det = np.where(R_det_t > 1e-15)[0]
if len(idx_det):
    last_det = idx_det[-1]
    print(f"Deterministic cutoff: {dist_km_t[last_det]:.2f} km")
    print(f"Deterministic SKR:    {R_det_t[last_det]:.3e}")

for cls in ["weak", "moderate", "strong"]:
    idx_t = np.where(R_turb_results[cls] > 1e-15)[0]
    if len(idx_t):
        last_t = idx_t[-1]
        print(f"\n{cls.upper()}")
        print(f"Cutoff distance: {dist_km_t[last_t]:.2f} km")
        print(f"SKR at cutoff:   {R_turb_results[cls][last_t]:.3e}")


print('All figures written:')
for fn in ['fig1_channel_model.png', 'fig2_yields_catstates.png',
           'fig3_phase_error_convergence.png', 'fig4_key_rate_vs_loss.png',
           'fig5_alpha_optimization.png', 'fig6_fso_channel_skr.png',
           'fig7_turbulence_distribution.png', 'fig8_turbulence_skr.png']:
    print(' ', fn)



#testing
L = 300.0

fso = M.deterministic_fso_channel(
    FSO_PARAMS['alpha_atm'],
    np.array([L]),
    FSO_PARAMS['wavelength'],
    FSO_PARAMS['w0'],
    FSO_PARAMS['receiver_radius']
)

eta_nom = float(fso['eta_fso'][0])

samples = M.sample_lognormal_eta(
    eta_nom,
    M.Cn2_classes['weak'],
    FSO_PARAMS['wavelength'],
    L*1e3,
    N_samples=500000
)

print("eta_nom =", eta_nom)
print("mean(sampled eta) =", np.mean(samples))
print("max(sampled eta) =", np.max(samples))

print(
    M.rytov_variance(
        M.Cn2_classes["weak"],
        FSO_PARAMS["wavelength"],
        300e3
    )
)