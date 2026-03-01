"""Supplement R2: QME occurrence map in (g_f, theta) space.

For each g_f value, determine whether QME occurs for representative theta pairs.
Also measures t_M(g_f) to look for divergence near g_c.
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

N, NA = 14, 4
g_c   = 0.5
t_max = 30.0
n_t   = 300
t_arr = np.linspace(0, t_max, n_t)

# Scan g_f from FM to PM phase, including near-critical region
gf_vals = np.concatenate([
    np.linspace(0.1, 0.45, 8),   # FM phase
    np.array([0.5]),               # exactly at critical
    np.linspace(0.55, 1.0, 8),    # just above critical (crucial region)
    np.linspace(1.2, 3.0, 6)      # deep PM
])
gf_vals = np.round(gf_vals, 4)

# Test three theta pairs with clear initial EA differences
theta_pairs = [(0.4, 1.2), (0.6, 1.0), (0.8, 1.4)]

print(f'Supplement R2: QME phase map  N={N}, NA={NA}, t_max={t_max}')
print(f'g_c = {g_c},  scanning {len(gf_vals)} g_f values')
print(f'Theta pairs: {theta_pairs}', flush=True)

results = {'params': {'N':N,'NA':NA,'g_c':g_c,'t_max':t_max,'n_t':n_t},
           't': t_arr.tolist(),
           'gf_vals': gf_vals.tolist()}

for gf in gf_vals:
    phase = 'FM' if gf < g_c else ('crit' if abs(gf-g_c)<0.02 else 'PM')
    print(f'\ng_f={gf:.4f} ({phase})', flush=True)
    H = build_TFIM_sparse(N, g=float(gf))
    results[f'gf{gf:.4f}'] = {'g_f': float(gf), 'phase': phase}

    # DQPT times (theta=0.4 reference)
    _,_,_,echo_ref = run_ea(N, NA, 0.4, H, t_arr)
    dqpt_raw = find_dqpt_times(t_arr, echo_ref, N, threshold=0.03)
    dqpt_clean = []
    for x in sorted(dqpt_raw):
        if not dqpt_clean or x - dqpt_clean[-1] > 1.5:
            dqpt_clean.append(x)

    # Check if genuine DQPT: rate function must exceed 0.1
    _,_,_,echo04 = run_ea(N, NA, 0.4, H, t_arr)
    losch04 = -np.log(np.maximum(echo04, 1e-300))/N
    max_losch = float(losch04.max())
    has_dqpt = max_losch > 0.05 and gf > g_c and len(dqpt_clean) >= 1

    results[f'gf{gf:.4f}']['dqpt_times'] = dqpt_clean
    results[f'gf{gf:.4f}']['max_losch_rate'] = max_losch
    results[f'gf{gf:.4f}']['has_dqpt'] = has_dqpt
    print(f'  DQPT: {has_dqpt} (max λ={max_losch:.4f}, times={[f"{x:.3f}" for x in dqpt_clean[:3]]})', flush=True)

    # QME for each theta pair
    qme_found = {}
    for th1, th2 in theta_pairs:
        DS1,_,_,_ = run_ea(N, NA, th1, H, t_arr)
        DS2,_,_,_ = run_ea(N, NA, th2, H, t_arr)
        if DS1[0] > DS2[0]:
            DS1, DS2 = DS2, DS1
        cross = find_crossings(t_arr, DS1, DS2)
        tm = cross[0] if cross else None
        qme_found[f'{th1}_{th2}'] = {
            'has_qme': tm is not None,
            't_M': float(tm) if tm is not None else None,
            'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0])
        }
        status = f't_M={tm:.3f}' if tm else 'NO QME'
        print(f'  θ={th1}/{th2}: {status}', flush=True)

    results[f'gf{gf:.4f}']['qme'] = qme_found

# Summary: QME occurrence map
print('\n===== QME OCCURRENCE MAP =====')
print(f'{"g_f":8} {"Phase":6} {"DQPT":6} ' +
      '  '.join(f'QME({th1}/{th2})' for th1,th2 in theta_pairs))
for gf in gf_vals:
    d = results[f'gf{gf:.4f}']
    dqpt = 'YES' if d['has_dqpt'] else 'no '
    qme_strs = []
    for th1,th2 in theta_pairs:
        q = d['qme'][f'{th1}_{th2}']
        qme_strs.append('YES' if q['has_qme'] else 'no ')
    print(f'  {gf:6.4f}  {d["phase"]:6} {dqpt:6}  ' + '    '.join(qme_strs))

# t_M vs g_f near critical region
print('\n===== t_M vs g_f (theta=0.4/1.2) =====')
print('  g_f      t_M      g_f-g_c   DQPT')
for gf in gf_vals:
    d = results[f'gf{gf:.4f}']
    q = d['qme'].get('0.4_1.2', {})
    tm = q.get('t_M')
    dqpt = 'Y' if d['has_dqpt'] else 'n'
    tm_s = f'{tm:.4f}' if tm else 'none  '
    print(f'  {gf:.4f}   {tm_s}   {gf-g_c:+.4f}    {dqpt}')

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/supp_r2_gf_scan.json','w'), indent=2)
print('\nR2 done → results/supp_r2_gf_scan.json')
