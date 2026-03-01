"""Supplement S1: N_A dependence — robustness of QME crossing under subsystem size variation.

Tests whether the QME crossing persists as N_A varies from 2 to 6 (for fixed N=18).
This addresses the single most critical gap: is the crossing an artifact of N_A=4?
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings

g_f = 2.0
N   = 18
NA_list = [2, 3, 4, 5, 6]
theta_pairs = [(0.4, 1.2), (0.6, 1.0), (0.8, 1.4), (1.0, 1.4)]
t_max, n_t = 15.0, 150
t_arr = np.linspace(0, t_max, n_t)

print(f'Supplement S1: N_A dependence  N={N}, g_f={g_f}')
print(f'Testing N_A in {NA_list} for {len(theta_pairs)} theta pairs')
print(flush=True)

H = build_TFIM_sparse(N, g=g_f)
print(f'Built H (dim={2**N}), computing EA...', flush=True)

results = {'params': {'N':N,'g_f':g_f,'t_max':t_max,'n_t':n_t},
           't': t_arr.tolist()}

for NA in NA_list:
    results[f'NA{NA}'] = {}
    print(f'\nN_A = {NA}:', flush=True)
    for th1, th2 in theta_pairs:
        t0 = time.time()
        DS1,_,_,_ = run_ea(N, NA, th1, H, t_arr)
        DS2,_,_,_ = run_ea(N, NA, th2, H, t_arr)
        if DS1[0] > DS2[0]:
            DS1, DS2 = DS2, DS1
        cross = find_crossings(t_arr, DS1, DS2)
        tm = cross[0] if cross else None
        results[f'NA{NA}'][f'{th1}_{th2}'] = {
            'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0]),
            't_M': float(tm) if tm else None,
            'has_qme': tm is not None
        }
        tm_s = f'{tm:.4f}' if tm else 'NO QME'
        print(f'  θ={th1}/{th2}: ΔS(0)={DS1[0]:.4f}/{DS2[0]:.4f}  t_M={tm_s}  ({time.time()-t0:.1f}s)',
              flush=True)

# Summary table
print('\n===== S1 SUMMARY: t_M vs N_A (N=18, g_f=2.0) =====')
print(f'{"Pair":12}' + ''.join(f'  NA={na}' for na in NA_list))
for th1, th2 in theta_pairs:
    key = f'{th1}_{th2}'
    row = f'{th1}/{th2}   '
    for NA in NA_list:
        d = results.get(f'NA{NA}', {}).get(key, {})
        tm = d.get('t_M')
        row += f'  {tm:.3f}' if tm else '  N/A '
    print(row)

# Check: does QME persist for all NA?
print('\n===== QME OCCURRENCE vs N_A =====')
print(f'{"Pair":12}' + ''.join(f'  NA={na}' for na in NA_list))
for th1, th2 in theta_pairs:
    key = f'{th1}_{th2}'
    row = f'{th1}/{th2}   '
    for NA in NA_list:
        d = results.get(f'NA{NA}', {}).get(key, {})
        has = d.get('has_qme', False)
        row += f'  {"YES":3}' if has else '  no '
    print(row)

all_qme = all(
    results.get(f'NA{NA}',{}).get(f'{th1}_{th2}',{}).get('has_qme',False)
    for NA in NA_list for th1,th2 in theta_pairs
)
print(f'\nQME robust for ALL (N_A, theta_pair) combinations: {all_qme}')

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/supp_s1_na_dep.json','w'), indent=2)
print('\nS1 done → results/supp_s1_na_dep.json')
