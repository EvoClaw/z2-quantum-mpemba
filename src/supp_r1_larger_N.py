"""Supplement R1: Finite-size scaling extended to N=16,18,20.

Product-spin-state initial condition requires exact ED (non-Gaussian state,
one-point functions <c_j> != 0 prevent correlation-matrix shortcut).
We push to N=20 using sparse expm_multiply.
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

g_f   = 2.0
g_c   = 0.5
NA    = 4
t_max = 20.0
n_t   = 200
t_arr = np.linspace(0, t_max, n_t)

# Extend FSS to N=16,18,20 for two theta pairs that showed clear convergence
theta_pairs = [(0.8, 1.4), (1.0, 1.4), (0.4, 1.2)]
N_list = [8, 10, 12, 14, 16, 18, 20]

T_theory = np.pi / (2*(g_f - g_c))  # = 1.0472

print(f'Supplement R1: Extended FSS  g_f={g_f}, g_c={g_c}, NA={NA}')
print(f'Theoretical t*_DQPT = {T_theory:.4f}')
print(f'N range: {N_list}  (N=20 → dim=2^20={2**20})')
print(flush=True)

results = {'params': {'g_f':g_f,'g_c':g_c,'NA':NA,'t_max':t_max,'n_t':n_t,
                      'T_theory':T_theory},
           't': t_arr.tolist()}

for N in N_list:
    t0 = time.time()
    print(f'\nN={N} (dim={2**N})...', flush=True)
    H = build_TFIM_sparse(N, g=g_f)
    results[f'N{N}'] = {}

    # Verify Hermitian
    diff = (H - H.conj().T).toarray()
    herm_err = float(np.abs(diff).max())
    print(f'  ||H-H†||={herm_err:.1e}', flush=True)

    # DQPT reference (use theta=0.4 which gives clearest signal)
    _,_,_,echo_ref = run_ea(N, NA, 0.4, H, t_arr)
    dqpt_ref = find_dqpt_times(t_arr, echo_ref, N, threshold=0.05)
    # Keep only primary peaks (gap > 1.5)
    dqpt_clean = []
    for x in sorted(dqpt_ref):
        if not dqpt_clean or x - dqpt_clean[-1] > 1.5:
            dqpt_clean.append(x)
    print(f'  DQPT times (clean): {[f"{x:.4f}" for x in dqpt_clean[:4]]}', flush=True)

    for th1, th2 in theta_pairs:
        DS1,_,_,echo1 = run_ea(N, NA, th1, H, t_arr)
        DS2,_,_,echo2 = run_ea(N, NA, th2, H, t_arr)
        if DS1[0] > DS2[0]:
            DS1, DS2 = DS2, DS1
        cross = find_crossings(t_arr, DS1, DS2)
        tm = cross[0] if cross else None
        nearest = min(dqpt_clean, key=lambda x: abs(x-tm)) if (tm and dqpt_clean) else None
        results[f'N{N}'][f'{th1}_{th2}'] = {
            'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0]),
            't_M': tm, 'dqpt_times': dqpt_clean,
        }
        tm_s = f'{tm:.4f}' if tm else 'N/A'
        nd_s = f'{nearest:.4f}' if nearest else 'N/A'
        print(f'  θ={th1}/{th2}: t_M={tm_s}  t*_1={dqpt_clean[0]:.4f}  nearest={nd_s}', flush=True)
    print(f'  Elapsed: {time.time()-t0:.1f}s', flush=True)

# Summary table
print('\n===== FSS SUMMARY =====')
print(f'{"Pair":12}' + ''.join(f'  N={N:2d}' for N in N_list))
for th1, th2 in theta_pairs:
    key = f'{th1}_{th2}'
    row = f'{th1}/{th2}   '
    for N in N_list:
        d = results.get(f'N{N}', {}).get(key, {})
        tm = d.get('t_M')
        row += f'  {tm:.3f}' if tm else '  N/A '
    print(row)

print(f'\nTheory t*_DQPT = {T_theory:.4f}  (N→∞ limit)')

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/supp_r1_fss.json','w'), indent=2)
print('\nR1 done → results/supp_r1_fss.json')
