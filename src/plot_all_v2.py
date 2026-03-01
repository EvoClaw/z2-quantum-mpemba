"""Generate all figures for the corrected experiments (publication-quality, English labels)

Fig 1: Exp A  -- DeltaS(t) and Loschmidt rate: ferromagnetic vs paramagnetic phase
Fig 2: Exp B  -- t_M vs initial EA difference DeltaDS (scatter)
Fig 3: Exp B  -- Fractional position test: t_M mod T* histogram
Fig 4: Exp C  -- Finite-size scaling of t_M and t*_DQPT
Fig 0: Bonus  -- Typical DeltaS(t) curves showing Mpemba crossings
"""
import numpy as np, json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import kstest

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'text.usetex': False,
})

# Colorblind-safe 7-color palette (Wong 2011)
COLORS = ['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#F0E442']

os.makedirs('results/figs', exist_ok=True)

# ── Figure 0: Typical DeltaS(t) with Mpemba crossings ────────────────────────
print('Figure 0: typical DeltaS(t) curves...', flush=True)
B = json.load(open('results/exp_b_v2.json'))
t2 = np.array(B['t'])
thetas_all = B['thetas']
thetas_demo = [thetas_all[i] for i in [5, 10, 20, 35]]

# Clean DQPT times: keep only primary peaks (gap > 1.5)
dqpt_raw = B['dqpt_times']
dqpt_clean = []
for x in sorted(dqpt_raw):
    if not dqpt_clean or x - dqpt_clean[-1] > 1.5:
        dqpt_clean.append(x)
dqpt_clean = dqpt_clean[:4]

fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
for ci, theta in enumerate(thetas_demo):
    ds = np.array(B['all_DS'][str(theta)])
    ax.plot(t2, ds, color=COLORS[ci], lw=1.8, label=fr'$\theta={theta:.2f}$')

for td in dqpt_clean:
    ax.axvline(td, color='gray', lw=1, ls='--', alpha=0.5)
if dqpt_clean:
    ax.text(dqpt_clean[0]+0.1, 0.58, r'$t^*_1$', color='gray', fontsize=10)
if len(dqpt_clean) > 1:
    ax.text(dqpt_clean[1]+0.1, 0.58, r'$t^*_2$', color='gray', fontsize=10)

# Mark one Mpemba crossing
pairs_demo = [(p['theta1'], p['theta2'], p['t_M'])
              for p in B['mpemba_pairs']
              if abs(p['theta1']-thetas_demo[0]) < 0.01
              and abs(p['theta2']-thetas_demo[2]) < 0.01]
if pairs_demo:
    t_cross = pairs_demo[0][2]
    ax.axvline(t_cross, color=COLORS[4], lw=1.5, ls=':', alpha=0.9)
    ax.text(t_cross+0.15, 0.30, r'$t_M$', color=COLORS[4], fontsize=11)

ax.set_xlabel('Time $t$')
ax.set_ylabel(r'Entanglement asymmetry $\Delta S$')
ax.set_title(r'$\Delta S(t)$ dynamics in the paramagnetic phase ($g_f=2.0$, $N=12$, $N_A=4$)')
ax.legend(loc='upper right')
ax.set_xlim(0, 15)
ax.set_ylim(bottom=0)
plt.savefig('results/figs/fig0_typical_DS.pdf', bbox_inches='tight')
plt.savefig('results/figs/fig0_typical_DS.png', bbox_inches='tight')
plt.close()
print('  Saved fig0_typical_DS', flush=True)


# ── Figure 1: Phase comparison (Exp A) ───────────────────────────────────────
print('Figure 1: phase comparison...', flush=True)
A = json.load(open('results/exp_a_v2.json'))
t  = np.array(A['t'])
thetas_show = [0.4, 0.6, 0.8, 1.0, 1.2]
N_sys = A['params']['N']

fig, axes = plt.subplots(2, 2, figsize=(10, 7.5), constrained_layout=True)

for col, gf in enumerate([0.2, 2.0]):
    ax_ds  = axes[0, col]
    ax_lsc = axes[1, col]
    phase_label = 'Ferromagnetic phase ($g_f < g_c$)' if gf < 0.5 else 'Paramagnetic phase ($g_f > g_c$)'

    # --- DeltaS(t) ---
    for ci, theta in enumerate(thetas_show):
        ds = np.array(A[f'gf{gf}'][str(theta)]['DS'])
        lw = 2.0 if theta in [0.4, 0.8, 1.2] else 1.0
        ax_ds.plot(t, ds, color=COLORS[ci], lw=lw, label=fr'$\theta={theta:.1f}$')
    ax_ds.set_xlabel('Time $t$')
    ax_ds.set_ylabel(r'$\Delta S(t)$')
    ax_ds.set_title(f'$g_f={gf}$  —  {phase_label}')
    ax_ds.legend(loc='upper right', ncol=2)
    ax_ds.set_xlim(0, 20)
    ax_ds.set_ylim(bottom=0)

    # --- Loschmidt rate function ---
    for ci, theta in enumerate([0.4, 0.8]):
        echo = np.array(A[f'gf{gf}'][str(theta)]['echo'])
        losch = -np.log(np.maximum(echo, 1e-300)) / N_sys
        ax_lsc.plot(t, losch, color=COLORS[ci], lw=1.8, label=fr'$\theta={theta:.1f}$')

    # Mark clean DQPT peaks for g_f=2.0 only
    if gf == 2.0:
        echo04 = np.array(A['gf2.0']['0.4']['echo'])
        losch04 = -np.log(np.maximum(echo04, 1e-300)) / N_sys
        dqpt_a = []
        for k in range(1, len(losch04)-1):
            if (losch04[k] > losch04[k-1] and losch04[k] > losch04[k+1]
                    and losch04[k] > 0.08):
                if not dqpt_a or t[k] - dqpt_a[-1] > 1.5:
                    dqpt_a.append(t[k])
        for td in dqpt_a[:4]:
            ax_lsc.axvline(td, color='gray', lw=1, ls='--', alpha=0.5)
        if dqpt_a:
            ax_lsc.text(dqpt_a[0]+0.1, ax_lsc.get_ylim()[1]*0.7 if ax_lsc.get_ylim()[1] > 0 else 0.05,
                        'DQPT', color='gray', fontsize=9)

    ax_lsc.set_xlabel('Time $t$')
    ax_lsc.set_ylabel(r'Rate function $\lambda(t)$')
    ax_lsc.legend(loc='upper right')
    ax_lsc.set_xlim(0, 20)
    ax_lsc.set_ylim(bottom=0)

for ax, lbl in zip(axes.flatten(), ['(a)', '(b)', '(c)', '(d)']):
    ax.text(0.02, 0.96, lbl, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top')

fig.suptitle(r'Entanglement asymmetry $\Delta S(t)$ and Loschmidt rate across quantum phases'
             f'\n($N={N_sys}$, $N_A=4$, $g_c=0.5$)',
             fontsize=11)
plt.savefig('results/figs/fig1_phase_compare.pdf', bbox_inches='tight')
plt.savefig('results/figs/fig1_phase_compare.png', bbox_inches='tight')
plt.close()
print('  Saved fig1_phase_compare', flush=True)


# ── Figure 2: t_M vs DeltaDS (Exp B) ─────────────────────────────────────────
print('Figure 2: t_M vs DeltaDS...', flush=True)
pairs  = B['mpemba_pairs']
dDS    = np.array([p['delta_DS'] for p in pairs])
tMs    = np.array([p['t_M']      for p in pairs])
DS0_mean = np.array([(p['DS0_1']+p['DS0_2'])/2 for p in pairs])

fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
sc = ax.scatter(dDS, tMs, c=DS0_mean, cmap='viridis', s=15, alpha=0.6)
cb = plt.colorbar(sc, ax=ax)
cb.set_label(r'Mean initial $\Delta S$')

# Polynomial trend line
mask = dDS < 0.6
if mask.sum() > 10:
    coeffs = np.polyfit(dDS[mask], tMs[mask], 2)
    xs = np.linspace(dDS[mask].min(), dDS[mask].max(), 100)
    ax.plot(xs, np.polyval(coeffs, xs), 'r--', lw=1.5,
            label='Quadratic fit', zorder=5)
    ax.legend()

ax.set_xlabel(r'Initial EA difference $\delta(\Delta S) = \Delta S(\theta_2,0) - \Delta S(\theta_1,0)$')
ax.set_ylabel(r'Mpemba crossing time $t_M$')
ax.set_title(r'QME crossing time vs initial EA difference ($g_f=2.0$, paramagnetic phase)')
ax.text(0.02, 0.96, '(a)', transform=ax.transAxes,
        fontsize=13, fontweight='bold', va='top')

plt.savefig('results/figs/fig2_tM_vs_dDS.pdf', bbox_inches='tight')
plt.savefig('results/figs/fig2_tM_vs_dDS.png', bbox_inches='tight')
plt.close()
print('  Saved fig2_tM_vs_dDS', flush=True)


# ── Figure 3: Fractional position test (Exp B) ───────────────────────────────
print('Figure 3: fractional position test...', flush=True)
T_period = B['params']['T_period_theory']
fracs  = np.array([p['frac_T_theory'] for p in pairs])
dists  = np.array([p['dist_to_half']  for p in pairs])
stat, pval = kstest(fracs, 'uniform')

fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

# Left: fractional position histogram
axes[0].hist(fracs, bins=20, color=COLORS[0], edgecolor='white', alpha=0.8,
             density=True, label='Observed')
axes[0].axhline(1.0, color='gray', ls='--', lw=1.5, label='Uniform (null)')
axes[0].axvline(0.5, color='red', ls=':', lw=2, label='DQPT position (0.5)')
axes[0].set_xlabel(r'$(t_M \,\mathrm{mod}\, T^*) / T^*$')
axes[0].set_ylabel('Probability density')
axes[0].set_title(f'Fractional position of $t_M$ relative to DQPT period\n'
                  f'KS test: $p={pval:.2e}$, peak at {fracs.mean():.3f} (DQPT at 0.5)')
axes[0].legend()
axes[0].text(0.02, 0.96, '(a)', transform=axes[0].transAxes,
             fontsize=13, fontweight='bold', va='top')

# Right: distance to nearest DQPT half-period
axes[1].hist(dists, bins=20, color=COLORS[2], edgecolor='white', alpha=0.8,
             density=True, label='Observed')
axes[1].axvline(0.25, color='gray', ls='--', lw=1.5,
                label='Uniform expectation (0.25)')
axes[1].axvline(dists.mean(), color=COLORS[1], ls='-', lw=2,
                label=f'Observed mean ({dists.mean():.3f})')
axes[1].set_xlabel(r'$|t_M \,\mathrm{mod}\, T^* - T^*/2| / T^*$')
axes[1].set_ylabel('Probability density')
axes[1].set_title(r'Normalized distance of $t_M$ to nearest DQPT time')
axes[1].legend()
axes[1].text(0.02, 0.96, '(b)', transform=axes[1].transAxes,
             fontsize=13, fontweight='bold', va='top')

fig.suptitle(f'Statistical test: is $t_M$ locked to DQPT times?  '
             f'($N={B["params"]["N"]}$, $g_f={B["params"]["g_f"]}$, '
             f'$T^*={T_period:.3f}$)',
             fontsize=11)
plt.savefig('results/figs/fig3_fractional_test.pdf', bbox_inches='tight')
plt.savefig('results/figs/fig3_fractional_test.png', bbox_inches='tight')
plt.close()
print('  Saved fig3_fractional_test', flush=True)


# ── Figure 4: Finite-size scaling (Exp C) ────────────────────────────────────
print('Figure 4: finite-size scaling...', flush=True)
C = json.load(open('results/exp_c_v2.json'))
N_list = [8, 10, 12, 14]
theta_pairs = [(0.4, 1.2), (0.6, 0.8), (0.8, 1.4), (1.0, 1.4)]
T_half_theory = np.pi / (2*(2.0-0.5))

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

# Left: t_M vs N
for ci, (th1, th2) in enumerate(theta_pairs):
    key = f'{th1}_{th2}'
    tM_list = [C[f'N{N}'][key]['t_M'] for N in N_list]
    valid = [(N, tm) for N, tm in zip(N_list, tM_list) if tm is not None]
    if valid:
        Ns, tms = zip(*valid)
        axes[0].plot(Ns, tms, 'o-', color=COLORS[ci], lw=1.5, ms=7,
                     label=fr'$\theta={th1}/{th2}$')
axes[0].set_xlabel('System size $N$')
axes[0].set_ylabel(r'Mpemba crossing time $t_M$')
axes[0].set_title(r'Finite-size scaling of $t_M$')
axes[0].legend()
axes[0].set_xticks(N_list)
axes[0].text(0.02, 0.96, '(a)', transform=axes[0].transAxes,
             fontsize=13, fontweight='bold', va='top')

# Right: t_DQPT_1 vs N
for ci, (th1, th2) in enumerate(theta_pairs[:2]):
    key = f'{th1}_{th2}'
    dqpt1_list = []
    for N in N_list:
        dts = (C[f'N{N}'][key].get('dqpt_ref') or
               C[f'N{N}'][key].get('dqpt_times') or [])
        dqpt1_list.append(dts[0] if dts else None)
    valid = [(N, d) for N, d in zip(N_list, dqpt1_list) if d is not None]
    if valid:
        Ns, dts = zip(*valid)
        axes[1].plot(Ns, dts, 'o-', color=COLORS[ci], lw=1.5, ms=7,
                     label=fr'$\theta={th1}/{th2}$')
axes[1].axhline(T_half_theory, color='black', ls='--', lw=1.5,
                label=fr'Thermodynamic limit ({T_half_theory:.4f})')
axes[1].set_xlabel('System size $N$')
axes[1].set_ylabel(r'First DQPT time $t^*_1$')
axes[1].set_title(r'Finite-size scaling of $t^*_1$')
axes[1].legend()
axes[1].set_xticks(N_list)
axes[1].text(0.02, 0.96, '(b)', transform=axes[1].transAxes,
             fontsize=13, fontweight='bold', va='top')

fig.suptitle(r'Finite-size scaling: $t_M$ and $t^*_{\rm DQPT}$  '
             f'($g_f=2.0$, $N_A=4$, $g_c=0.5$)',
             fontsize=11)
plt.savefig('results/figs/fig4_finite_size.pdf', bbox_inches='tight')
plt.savefig('results/figs/fig4_finite_size.png', bbox_inches='tight')
plt.close()
print('  Saved fig4_finite_size', flush=True)

print('\nAll figures saved to results/figs/')
for fn in ['fig0_typical_DS', 'fig1_phase_compare',
           'fig2_tM_vs_dDS', 'fig3_fractional_test', 'fig4_finite_size']:
    for ext in ['png', 'pdf']:
        path = f'results/figs/{fn}.{ext}'
        size_kb = os.path.getsize(path)//1024 if os.path.exists(path) else -1
        print(f'  {path}  ({size_kb} KB)')
