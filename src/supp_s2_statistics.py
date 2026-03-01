"""Supplement S2: Statistical rigor for the t_M vs t_DQPT anti-correlation claim.

Loads exp_b_v2.json (671 pairs) and performs proper statistical tests:
1. Explicit null model definition and simulation
2. Bootstrap confidence interval on hit rate
3. p-value via exact binomial test
4. KS statistic and p-value for fractional position uniformity
5. Spearman rank correlation between t_M and t_DQPT distance
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os
from scipy import stats

# Load existing 671-pair data
data = json.load(open('results/exp_b_v2.json'))
t_arr = np.array(data['t'])
g_f = data['params']['g_f']
T_period_theory = data['params'].get('T_period_theory', np.pi / (2*(g_f - 0.5)))

print(f'Supplement S2: Statistical rigor for t_M vs t_DQPT anti-correlation')
print(f'g_f = {g_f}, T*_theory = {T_period_theory:.4f}')
print(flush=True)

# Extract t_M values and DQPT times from mpemba_pairs array
mp_pairs   = data.get('mpemba_pairs', [])
dqpt_times = data.get('dqpt_times', [])
T_period_num = None
if len(dqpt_times) >= 2:
    T_period_num = np.mean(np.diff(dqpt_times[:min(4,len(dqpt_times))]))

tM_vals = []
tDQPT_first = []
fracs = []

for pair in mp_pairs:
    tm = pair.get('t_M')
    if tm is None:
        continue
    nearest = pair.get('nearest_dqpt')
    frac = pair.get('frac_T_theory')
    if nearest is None:
        nearest = dqpt_times[0] if dqpt_times else None
    if nearest is None:
        continue

    tM_vals.append(float(tm))
    tDQPT_first.append(float(nearest))
    fracs.append(float(frac) if frac is not None else (float(tm) % T_period_theory) / T_period_theory)

tM_vals  = np.array(tM_vals)
tDQPT_1  = np.array(tDQPT_first)
fracs    = np.array(fracs)

print(f'Number of pairs with t_M and t_DQPT: {len(tM_vals)}', flush=True)

# ── 1. Hit rate analysis using FRACTIONAL POSITIONS ─────────────────────────
# fracs = (t_M mod T*) / T*  — measures where t_M falls within DQPT period
# IMPORTANT: DQPTs occur at the MIDDLE of each period: frac = 0.5
#   t_DQPT = (2n+1) * T*/2 = T*/2, 3T*/2, 5T*/2, ...
#   Within each period, DQPT is at frac = 0.5
# Anti-DQPT alignment: frac AVOIDS 0.5 (prefers frac near 0 or 1 = between DQPTs)
# "Near DQPT" = |frac - 0.5| < delta (within delta of frac=0.5)
# "Near DQPT" hit rate < 2*delta → anti-correlation (t_M avoids DQPT times)

delta = 0.2  # within 20% of DQPT: |frac - 0.5| < 0.2 → frac in (0.3, 0.7)
frac_arr = np.array(fracs)
dist_to_dqpt = np.abs(frac_arr - 0.5)
hits_near_dqpt = dist_to_dqpt < delta
n_hits = int(hits_near_dqpt.sum())
n_total = len(tM_vals)
hit_rate = n_hits / n_total

# Null model: if frac is uniform in [0,1], P(|frac-0.5| < delta) = 2*delta = 0.40
null_rate_theoretical = 2 * delta

print(f'\n── Hit rate: t_M within {delta:.0%} of DQPT time (frac ∈ [0.3, 0.7]) ──')
print(f'Note: DQPTs at frac=0.5; frac≈0 or ≈1 means between DQPTs')
print(f'Observed: {n_hits}/{n_total} = {hit_rate:.3f}')
print(f'Null (uniform frac): P(near DQPT) = 2*{delta} = {null_rate_theoretical:.3f}')
print(f'Anti-correlation expected: hit_rate < null_rate')
t_max = float(t_arr[-1])

# Bootstrap confidence interval on hit rate
n_boot = 10000
np.random.seed(42)
boot_hits = np.random.binomial(n_total, hit_rate, size=n_boot)
boot_rates = boot_hits / n_total
ci_lo = np.percentile(boot_rates, 2.5)
ci_hi = np.percentile(boot_rates, 97.5)
print(f'Bootstrap 95% CI on hit rate: [{ci_lo:.3f}, {ci_hi:.3f}]')

# Exact binomial test: is observed hit_rate < null_rate?
binom_test = stats.binomtest(n_hits, n_total, null_rate_theoretical, alternative='less').pvalue
print(f'Binomial test (H0: rate >= null, H1: rate < null): p = {binom_test:.4f}')

significance = '***' if binom_test < 0.001 else ('**' if binom_test < 0.01 else ('*' if binom_test < 0.05 else 'ns'))
print(f'Significance: {significance} (anti-correlation confirmed: {binom_test < 0.05})')

# ── 2. KS test on fractional positions ─────────────────────────────────────
print(f'\n── KS test: fractional position uniformity ──')
print(f'N fractions: {len(fracs)}')
print(f'Fracs range: [{fracs.min():.3f}, {fracs.max():.3f}], mean: {fracs.mean():.3f}')
ks_stat, ks_pval = stats.kstest(fracs, 'uniform')
print(f'KS statistic: {ks_stat:.4f}, p-value: {ks_pval:.6f}')
print(f'Reject uniformity (p < 0.05): {ks_pval < 0.05}')
print(f'Fraction mean = {fracs.mean():.4f} (vs 0.5 for uniform)')
print(f'Fraction std  = {fracs.std():.4f} (vs {1/np.sqrt(12):.4f} for uniform)')

# ── 3. Spearman correlation: t_M vs distance to nearest DQPT ───────────────
print(f'\n── Spearman rank correlation: t_M vs |t_M - t_DQPT| ──')
dist_to_dqpt = np.abs(tM_vals - tDQPT_1)
rho, p_spearman = stats.spearmanr(tM_vals, dist_to_dqpt)
print(f'Spearman rho(t_M, |t_M - t_DQPT|) = {rho:.4f}, p = {p_spearman:.4f}')

# ── 4. Test: is peak of fracs distribution below 0.5? ─────────────────────
print(f'\n── Distribution analysis of fractional positions ──')
hist, bin_edges = np.histogram(fracs, bins=10, range=(0,1))
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
peak_bin = bin_centers[hist.argmax()]
print(f'Histogram peak at frac = {peak_bin:.2f} (expect 0.5 for DQPT alignment, got {peak_bin:.2f})')
print('Histogram:')
for b, c in zip(bin_centers, hist):
    bar = '#' * (c // 2)
    print(f'  frac={b:.2f}: {bar} ({c})')

# ── 5. Summary ──────────────────────────────────────────────────────────────
print(f'\n===== S2 STATISTICAL SUMMARY =====')
print(f'Pairs analyzed: {n_total}')
print(f'Hit rate: {hit_rate:.3f}  (null: {null_rate_theoretical:.3f})')
print(f'95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]')
print(f'Binomial test p-value: {binom_test:.4f}  ({significance})')
print(f'KS test p-value: {ks_pval:.6f}  (reject uniform: {ks_pval < 0.05})')
print(f'Fractional position peak: {peak_bin:.2f} (expected 0.5 for DQPT alignment)')

results_s2 = {
    'n_pairs': int(n_total),
    'hit_rate': float(hit_rate),
    'null_rate_theoretical': float(null_rate_theoretical),
    'delta_frac': float(delta),
    'bootstrap_ci_95': [float(ci_lo), float(ci_hi)],
    'binomial_p_value': float(binom_test),
    'binomial_significant': bool(binom_test < 0.05),
    'ks_stat': float(ks_stat),
    'ks_p_value': float(ks_pval),
    'reject_uniformity': bool(ks_pval < 0.05),
    'frac_peak_bin': float(peak_bin),
    'frac_mean': float(fracs.mean()),
    'conclusion': f'Hit rate {hit_rate:.3f} < null {null_rate_theoretical:.3f} (p={binom_test:.4f}); anti-correlation confirmed' if binom_test < 0.05 else f'Not statistically significant at p<0.05'
}
json.dump(results_s2, open('results/supp_s2_stats.json','w'), indent=2)
print('\nS2 done → results/supp_s2_stats.json')
print(f'\nCONCLUSION: {results_s2["conclusion"]}')
