"""Generate publication-quality figures for supplement experiments R1, R2, R4."""
import numpy as np, json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({
    'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 12,
    'legend.fontsize': 9, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'lines.linewidth': 1.8, 'figure.dpi': 150,
    'text.usetex': False,
    'font.family': 'DejaVu Sans',
})
os.makedirs('results/figs', exist_ok=True)

# ─── FIGURE S1: Extended FSS (R1) ──────────────────────────────────────────
print('Generating Figure S1: Extended FSS...', flush=True)
r1 = json.load(open('results/supp_r1_fss.json'))
T_theory = r1['params']['T_theory']
N_list   = [8, 10, 12, 14, 16, 18]

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(N_list)))

pair_label = {'0.8_1.4': r'$\theta_1{=}0.8, \theta_2{=}1.4$',
              '1.0_1.4': r'$\theta_1{=}1.0, \theta_2{=}1.4$',
              '0.4_1.2': r'$\theta_1{=}0.4, \theta_2{=}1.2$'}

# Left: t_M vs 1/N for pairs 0.8_1.4 and 1.0_1.4
ax = axes[0]
for pair, marker, ls in [('0.8_1.4','o','-'), ('1.0_1.4','s','--')]:
    tMs, Ns = [], []
    for N in N_list:
        d = r1.get(f'N{N}', {}).get(pair, {})
        tm = d.get('t_M')
        if tm:
            tMs.append(tm); Ns.append(N)
    ax.plot([1/N for N in Ns], tMs, marker=marker, ls=ls, color='navy' if pair=='0.8_1.4' else 'crimson',
            label=pair_label[pair], markersize=6, zorder=5)
    # Extrapolation: linear fit in 1/N
    if len(Ns) >= 3:
        inv_N = np.array([1/N for N in Ns])
        coeffs = np.polyfit(inv_N, tMs, 1)
        t_inf = coeffs[1]
        x_ext = np.linspace(0, max(inv_N)*1.1, 100)
        ax.plot(x_ext, np.polyval(coeffs, x_ext), '--',
                color='navy' if pair=='0.8_1.4' else 'crimson', alpha=0.4, linewidth=1)
        ax.axhline(t_inf, color='navy' if pair=='0.8_1.4' else 'crimson',
                   linestyle=':', alpha=0.7, linewidth=1)
        ax.annotate(f'$t_M(N\\to\\infty)={t_inf:.3f}$',
                    xy=(0.02, t_inf), xycoords=('axes fraction','data'),
                    fontsize=8, color='navy' if pair=='0.8_1.4' else 'crimson',
                    va='bottom' if pair=='0.8_1.4' else 'top')

ax.set_xlabel(r'$1/N$')
ax.set_ylabel(r'$t_M$')
ax.set_title('(a) Convergence of $t_M$ with system size')
ax.legend(loc='upper right')
ax.set_xlim(left=0)

# Right: t_M vs N (semi-log)
ax = axes[1]
for pair, marker, ls in [('0.8_1.4','o','-'), ('1.0_1.4','s','--'), ('0.4_1.2','^',':')]:
    tMs, Ns = [], []
    for N in N_list:
        d = r1.get(f'N{N}', {}).get(pair, {})
        tm = d.get('t_M')
        if tm:
            tMs.append(tm); Ns.append(N)
    color = {'0.8_1.4':'navy', '1.0_1.4':'crimson', '0.4_1.2':'darkgreen'}[pair]
    ax.plot(Ns, tMs, marker=marker, ls=ls, color=color, label=pair_label[pair], markersize=6)

ax.set_xlabel('$N$')
ax.set_ylabel('$t_M$')
ax.set_title('(b) $t_M$ vs. $N$')
ax.legend(loc='right')
ax.set_xticks(N_list)

fig.suptitle('Supplementary: Finite-size scaling to $N=18$  ($g_f=2.0$, $g_c=0.5$)',
             fontsize=11)
plt.tight_layout()
plt.savefig('results/figs/figS1_fss_extended.pdf', bbox_inches='tight')
plt.savefig('results/figs/figS1_fss_extended.png', bbox_inches='tight')
plt.close()
print('  Saved figS1_fss_extended', flush=True)

# ─── FIGURE S2: QME occurrence map vs g_f (R2) ──────────────────────────────
print('Generating Figure S2: QME phase map...', flush=True)
r2 = json.load(open('results/supp_r2_gf_scan.json'))
g_c = r2['params']['g_c']

gf_arr  = [d['g_f']         for d in r2['gf_data']]
tM_arr  = [d['t_M_04_12']   for d in r2['gf_data']]
dqpt_arr= [d['has_dqpt']    for d in r2['gf_data']]

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# Left: t_M(g_f) for theta=0.4/1.2
ax = axes[0]
# Color by phase
gf_FM  = [g for g,t,d in zip(gf_arr,tM_arr,dqpt_arr) if g < g_c]
tM_FM  = [t for g,t,d in zip(gf_arr,tM_arr,dqpt_arr) if g < g_c]
gf_PM  = [g for g,t,d in zip(gf_arr,tM_arr,dqpt_arr) if g > g_c]
tM_PM  = [t for g,t,d in zip(gf_arr,tM_arr,dqpt_arr) if g > g_c]
gf_cr  = [g for g,t,d in zip(gf_arr,tM_arr,dqpt_arr) if abs(g-g_c)<0.01]
tM_cr  = [t for g,t,d in zip(gf_arr,tM_arr,dqpt_arr) if abs(g-g_c)<0.01]

ax.plot(gf_FM, tM_FM, 'o-', color='steelblue', label='FM phase ($g_f < g_c$)', markersize=5)
ax.plot(gf_PM, tM_PM, 's-', color='tomato', label='PM phase ($g_f > g_c$)', markersize=5)
ax.plot(gf_cr, tM_cr, 'D', color='purple', label='Critical', markersize=8, zorder=6)
ax.axvline(g_c, color='gray', linestyle='--', linewidth=1.2, label=f'$g_c={g_c}$')
ax.set_xlabel(r'$g_f$ (post-quench field)')
ax.set_ylabel(r'$t_M$ (Mpemba crossing time)')
ax.set_title(r'(a) $t_M(g_f)$ for $\theta_1=0.4, \theta_2=1.2$')
ax.legend(fontsize=8)
ax.set_yscale('log')

# Right: QME occurrence (bar) vs DQPT occurrence
ax = axes[1]
gf_vals = np.array(gf_arr)
qme_vals = np.ones(len(gf_arr), dtype=float)  # QME=1 for all
dqpt_vals = np.array([float(d) for d in dqpt_arr])

bar_width = np.diff(gf_vals).min() * 0.4
ax.bar(gf_vals - bar_width/2, qme_vals, width=bar_width, color='steelblue',
       alpha=0.8, label='QME (all $\\theta$ pairs)')
ax.bar(gf_vals + bar_width/2, dqpt_vals, width=bar_width, color='tomato',
       alpha=0.8, label='DQPT')
ax.axvline(g_c, color='gray', linestyle='--', linewidth=1.2)
ax.set_xlabel(r'$g_f$')
ax.set_ylabel('Occurrence (1=yes, 0=no)')
ax.set_title('(b) QME vs. DQPT occurrence')
ax.set_yticks([0, 1])
ax.set_yticklabels(['No', 'Yes'])
ax.legend(fontsize=9)
ax.text(g_c+0.05, 0.5, r'$g_c$', ha='left', va='center', color='gray', fontsize=10)

fig.suptitle('Supplementary: QME occurrence map ($N=14$, $g_c=0.5$)', fontsize=11)
plt.tight_layout()
plt.savefig('results/figs/figS2_qme_map.pdf', bbox_inches='tight')
plt.savefig('results/figs/figS2_qme_map.png', bbox_inches='tight')
plt.close()
print('  Saved figS2_qme_map', flush=True)

# ─── FIGURE S3: Anisotropic XY chain (R4) ───────────────────────────────────
print('Generating Figure S3: XY chain...', flush=True)
r4 = json.load(open('results/supp_r4_xy.json'))

gammas_str = ['0.0', '0.25', '0.5', '0.75', '1.0']
gamma_labels = ['TFIM ($\\gamma{=}0$)', '$\\gamma{=}0.25$', '$\\gamma{=}0.5$',
                '$\\gamma{=}0.75$', 'XX ($\\gamma{=}1$)']
g_c_list = [0.5*np.sqrt(1-float(g)) for g in gammas_str]
g_f = 2.0

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left: t_M vs gamma for both theta pairs (g_f=2.0)
ax = axes[0]
colors_xy = ['navy', 'steelblue', 'seagreen', 'darkorange', 'crimson']
for pair_key, marker, ls, label in [
        ('qme_0.4_1.2','o','-',r'$\theta_1{=}0.4, \theta_2{=}1.2$'),
        ('qme_0.8_1.4','s','--',r'$\theta_1{=}0.8, \theta_2{=}1.4$')]:
    tMs, gm_vals = [], []
    for gm_s in gammas_str:
        d = r4.get(gm_s, {})
        q = d.get(pair_key, {})
        if q.get('has_qme'):
            tMs.append(q['t_M']); gm_vals.append(float(gm_s))
    ax.plot(gm_vals, tMs, marker=marker, ls=ls, color='navy' if 'qme_0.4' in pair_key else 'crimson',
            label=label, markersize=7)

ax.set_xlabel(r'Anisotropy $\gamma$ ($\gamma{=}0$: TFIM, $\gamma{=}1$: XX)')
ax.set_ylabel(r'$t_M$ ($g_f=2.0$)')
ax.set_title('(a) $t_M$ vs. anisotropy $\\gamma$ (fixed $g_f=2.0$)')
ax.legend(fontsize=9)

# Add secondary x-axis for g_c
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks([float(g) for g in gammas_str])
ax2.set_xticklabels([f'$g_c={gc:.3f}$' for gc in g_c_list], fontsize=7, rotation=30)
ax2.set_xlabel('QPT critical field $g_c(\\gamma)$', fontsize=9)

# Right: Critical test (g_f=0.42 comparison)
ax = axes[1]
crit = r4.get('critical_test', {})
g_f2 = crit.get('g_f', 0.42)
crit_data = crit.get('data', {})

gm_list = [float(g) for g in gammas_str]
dqpt_crit = [int(crit_data.get(gs,{}).get('has_dqpt', 0)) for gs in gammas_str]
tM_crit   = [crit_data.get(gs,{}).get('t_M') for gs in gammas_str]
has_qme_crit = [tm is not None for tm in tM_crit]

bar_w = 0.08
ax.bar(np.array(gm_list) - bar_w/2, [int(d) for d in dqpt_crit], width=bar_w,
       color='tomato', alpha=0.8, label='DQPT')
ax.bar(np.array(gm_list) + bar_w/2, [int(d) for d in has_qme_crit], width=bar_w,
       color='steelblue', alpha=0.8, label='QME')
ax.set_xlabel(r'Anisotropy $\gamma$')
ax.set_ylabel('Occurrence (1=yes, 0=no)')
ax.set_title(f'(b) Critical test: $g_f={g_f2}$, all $g_c(\\gamma) < {g_f2}$')
ax.set_yticks([0, 1]); ax.set_yticklabels(['No', 'Yes'])
ax.set_xticks(gm_list)
ax.set_xticklabels([f'$\\gamma={g:.2f}$' for g in gm_list], fontsize=9)
ax.legend(fontsize=9)

# Annotate g_c values
for gm, gc in zip(gm_list, g_c_list):
    ax.annotate(f'$g_c={gc:.3f}$', xy=(gm, 1.05), ha='center', fontsize=7,
                xycoords=('data','axes fraction'))

fig.suptitle(f'Supplementary: Anisotropic XY chain ($N=12$)', fontsize=11)
plt.tight_layout()
plt.savefig('results/figs/figS3_xy_chain.pdf', bbox_inches='tight')
plt.savefig('results/figs/figS3_xy_chain.png', bbox_inches='tight')
plt.close()
print('  Saved figS3_xy_chain', flush=True)

# ─── FIGURE S4: n_k(t) quasiparticle occupation (R5) ───────────────────────
print('Generating Figure S4: Quasiparticle occupations...', flush=True)
r5 = json.load(open('results/supp_r5_nk.json'))
t_arr5 = np.array(r5['t'])
k_sym5 = np.array(r5['k_sym'])
g_f5 = r5['params']['g_f']

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left: n_k(t) at fixed t for different theta
ax = axes[0]
t_indices = [0, 20, 50, 100]  # t = 0, ~2, ~5, ~10
t_labels_idx = [0, len(t_arr5)//6, len(t_arr5)//3, len(t_arr5)//2]
theta_sel = ['0.40', '0.80', '1.20']
colors_theta = ['steelblue', 'seagreen', 'crimson']

# Plot n_k(t) at t_M for theta=0.4 and 1.2
for theta_s, col in zip(theta_sel, colors_theta):
    key = f'theta{theta_s}'
    if key not in r5: continue
    nk_data = r5[key]['n_k_t']
    k_vals_plot = sorted(nk_data.keys(), key=float)
    k_num = np.array([float(k) for k in k_vals_plot])
    # t=0 values
    nk_t0 = np.array([nk_data[k][0] for k in k_vals_plot])
    # t ≈ t_M
    t_mid_idx = len(t_arr5)//2
    nk_tM = np.array([nk_data[k][t_mid_idx] for k in k_vals_plot])
    ax.plot(k_num/np.pi, nk_t0, 'o', color=col, markersize=5,
            label=f'$\\theta={theta_s}$ ($t=0$)', alpha=0.7)

ax.set_xlabel(r'$k/\pi$')
ax.set_ylabel(r'$n_k(0)$ = $\sin^2(\theta/2)$')
ax.set_title(r'(a) Initial quasiparticle occupation $n_k(0)$')
ax.legend(fontsize=9)
ax.set_xlim([0, 2])
ax.axhline(0.5, color='gray', ls=':', label='$n_k=0.5$ (max uncertainty)')

# Right: time evolution of n_k for representative modes (theta=1.0)
ax = axes[1]
theta_r5 = '1.00'
key_r5 = f'theta{theta_r5}'
if key_r5 in r5:
    nk_data = r5[key_r5]['n_k_t']
    k_vals_r5 = sorted(nk_data.keys(), key=float)
    colors_k = ['navy', 'steelblue', 'seagreen', 'darkorange', 'crimson']
    for ki, (k_s, col) in enumerate(zip(k_vals_r5, colors_k)):
        k_val = float(k_s)
        nk_t = np.array(nk_data[k_s])
        ax.plot(t_arr5, nk_t, color=col, alpha=0.85,
                label=f'$k/\\pi={k_val/np.pi:.3f}$', linewidth=1.5)

t_M_ref = r5.get('t_M_ref', {}).get('0.8_1.4')
if t_M_ref:
    ax.axvline(t_M_ref, color='black', linestyle='--', linewidth=1.5,
               label=f'$t_M={t_M_ref:.2f}$')
ax.axhline(0.5, color='gray', ls=':', linewidth=1, label='$n_k=0.5$')
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$n_k(t)$ (quasiparticle occupation)')
ax.set_title(f'(b) Mode occupations under $H(g_f={g_f5})$,  $\\theta={theta_r5}$')
ax.legend(fontsize=8, ncol=2)

fig.suptitle('Supplementary: Quasiparticle mode occupation analysis', fontsize=11)
plt.tight_layout()
plt.savefig('results/figs/figS4_nk.pdf', bbox_inches='tight')
plt.savefig('results/figs/figS4_nk.png', bbox_inches='tight')
plt.close()
print('  Saved figS4_nk', flush=True)

print('\nAll supplement figures saved to results/figs/')
print('  figS1_fss_extended.{pdf,png}')
print('  figS2_qme_map.{pdf,png}')
print('  figS3_xy_chain.{pdf,png}')
print('  figS4_nk.{pdf,png}')
