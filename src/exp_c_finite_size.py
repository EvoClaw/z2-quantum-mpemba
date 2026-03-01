"""实验 C：有限尺度缩放 N=8,10,12,14（固定 NA=4，theta 对 = 0.8 vs 1.4）"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

g_f = 2.0
t_max, n_t = 20.0, 80
t_arr = np.linspace(0, t_max, n_t)
NA = 4
# 代表性 theta 对（有明显 Mpemba 效应）
theta_pairs = [(0.4, 1.2), (0.6, 0.8), (0.8, 1.4), (1.0, 1.4)]
N_list = [8, 10, 12, 14]

print(f'实验 C: 有限尺度缩放, g_f={g_f}, NA={NA}', flush=True)

results = {'params': {'g_f':g_f,'NA':NA,'t_max':t_max}, 't': t_arr.tolist()}

for N in N_list:
    t0 = time.time()
    print(f'\nN={N} (dim=2^{N}={2**N})...', flush=True)
    H = build_TFIM_sparse(N, g=g_f)
    results[f'N{N}'] = {}

    # 对每个 theta 对
    for th1, th2 in theta_pairs:
        DS1,_,_,echo1 = run_ea(N, NA, th1, H, t_arr)
        DS2,_,_,echo2 = run_ea(N, NA, th2, H, t_arr)
        cross = find_crossings(t_arr, DS1, DS2)
        dqpt = find_dqpt_times(t_arr, echo1, N)
        tm = cross[0] if cross else None
        nearest = min(dqpt, key=lambda x: abs(x-tm)) if (tm and dqpt) else None
        results[f'N{N}'][f'{th1}_{th2}'] = {
            'DS1': DS1.tolist(), 'DS2': DS2.tolist(),
            'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0]),
            't_M': tm, 'dqpt_times': dqpt,
            'nearest_dqpt': float(nearest) if nearest else None
        }
        tm_str = f'{tm:.3f}' if tm else 'none'
        nd_str = f'{nearest:.3f}' if nearest else 'N/A'
        print(f'  theta={th1} vs {th2}: DS0={DS1[0]:.3f}/{DS2[0]:.3f}, t_M={tm_str}, DQPT≈{nd_str}', flush=True)

    print(f'  N={N} 耗时: {time.time()-t0:.1f}s', flush=True)

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/exp_c_finite_size.json','w'), indent=2)
print('\n实验 C 完成')
