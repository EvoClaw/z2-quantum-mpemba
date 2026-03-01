"""实验 C（修正版）：有限尺度缩放

与原版相同，但标注正确的 g_c=0.5 和物理相。
分析 t_M 和 t_DQPT 随 N 的有限尺度行为。
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

g_f = 2.0
g_c = 0.5
t_max, n_t = 20.0, 200
t_arr = np.linspace(0, t_max, n_t)
NA = 4
theta_pairs = [(0.4, 1.2), (0.6, 0.8), (0.8, 1.4), (1.0, 1.4)]
N_list = [8, 10, 12, 14]

print(f'实验 C（修正版）: g_f={g_f}, g_c={g_c}, NA={NA}')
print(f'理论 DQPT 半周期 = π/(2(g_f-g_c)) = {np.pi/(2*(g_f-g_c)):.4f}', flush=True)

results = {'params': {'g_f':g_f,'g_c':g_c,'NA':NA,'t_max':t_max,'n_t':n_t}, 't': t_arr.tolist()}

for N in N_list:
    t0 = time.time()
    print(f'\nN={N} (dim=2^{N}={2**N})...', flush=True)
    H = build_TFIM_sparse(N, g=g_f)

    # 验证 Hermitian
    diff = (H - H.conj().T).toarray()
    print(f'  ||H-H†||_max = {np.abs(diff).max():.2e}', flush=True)

    results[f'N{N}'] = {}
    theta_ref = 0.4
    _,_,_,echo_ref = run_ea(N, NA, theta_ref, H, t_arr)
    dqpt_ref = find_dqpt_times(t_arr, echo_ref, N, threshold=0.01)
    print(f'  DQPT times (theta={theta_ref}): {[f"{x:.3f}" for x in dqpt_ref]}', flush=True)

    for th1, th2 in theta_pairs:
        DS1,_,_,echo1 = run_ea(N, NA, th1, H, t_arr)
        DS2,_,_,echo2 = run_ea(N, NA, th2, H, t_arr)
        if DS1[0] > DS2[0]:
            DS1, DS2 = DS2, DS1
        cross = find_crossings(t_arr, DS1, DS2)
        dqpt = find_dqpt_times(t_arr, echo1, N, threshold=0.01)
        tm = cross[0] if cross else None
        nearest = min(dqpt, key=lambda x: abs(x-tm)) if (tm and dqpt) else None
        results[f'N{N}'][f'{th1}_{th2}'] = {
            'DS1': DS1.tolist(), 'DS2': DS2.tolist(),
            'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0]),
            't_M': tm, 'dqpt_times': dqpt, 'dqpt_ref': dqpt_ref,
            'nearest_dqpt': float(nearest) if nearest else None
        }
        tm_str = f'{tm:.3f}' if tm else 'none'
        nd_str = f'{nearest:.3f}' if nearest else 'N/A'
        diff_str = f'{abs(tm-nearest):.3f}' if (tm and nearest) else 'N/A'
        print(f'  theta={th1} vs {th2}: DS0={DS1[0]:.3f}/{DS2[0]:.3f}, t_M={tm_str}, nearest_DQPT={nd_str}, |diff|={diff_str}', flush=True)

    print(f'  N={N} 耗时: {time.time()-t0:.1f}s', flush=True)

# 有限尺度分析：t_M 和 t_DQPT_1 随 N 的变化
print('\n===== 有限尺度分析 =====')
T_theory = np.pi / (2*(g_f - g_c))
print(f'理论极限（N→∞）: t_DQPT_1 = {T_theory:.4f}')
for th1, th2 in theta_pairs:
    print(f'\ntheta 对 {th1} vs {th2}:')
    for N in N_list:
        key = f'{th1}_{th2}'
        d = results[f'N{N}'][key]
        tm_str = f'{d["t_M"]:.3f}' if d['t_M'] else 'N/A'
        dqpt1_str = f'{d["dqpt_times"][0]:.3f}' if d['dqpt_times'] else 'N/A'
        print(f'  N={N}: t_M={tm_str}  t_DQPT_1={dqpt1_str}')

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/exp_c_v2.json','w'), indent=2)
print('\n实验 C（修正版）完成，结果保存至 results/exp_c_v2.json')
