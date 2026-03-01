"""实验 B：精细 theta 网格，系统研究 t_M 与初始 DeltaS 差值 Delta(DS) 的关系"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

N, NA = 12, 4
g_f = 2.0
t_max, n_t = 20.0, 100
t_arr = np.linspace(0, t_max, n_t)
# 精细 theta 网格
thetas = np.linspace(0.05, np.pi/2 - 0.05, 40).tolist()

print(f'实验 B: N={N}, NA={NA}, g_f={g_f}')
print(f'theta 网格: {len(thetas)} 点, [{thetas[0]:.3f}, {thetas[-1]:.3f}]', flush=True)

H = build_TFIM_sparse(N, g=g_f)
print(f'H built, nnz={H.nnz}', flush=True)

# DQPT 参考时间（用 theta=0.4 附近计算）
theta_ref = min(thetas, key=lambda x: abs(x-0.4))
DS_ref,_,_,echo_ref = run_ea(N, NA, theta_ref, H, t_arr)
dqpt_times = find_dqpt_times(t_arr, echo_ref, N)
print(f'DQPT times: {[f"{x:.3f}" for x in dqpt_times]}', flush=True)

# 计算所有 theta 的 DeltaS(t)
all_DS = {}
all_DS0 = {}
for theta in thetas:
    DS,_,_,_ = run_ea(N, NA, theta, H, t_arr)
    all_DS[theta] = DS
    all_DS0[theta] = DS[0]
    print(f'  theta={theta:.3f}: DS(0)={DS[0]:.4f} DS({t_max:.0f})={DS[-1]:.4f}', flush=True)

# 系统计算所有 theta 对的 t_M
print('\n计算所有 theta 对的 Mpemba 交叉时间...', flush=True)
mpemba_data = []
for i, th1 in enumerate(thetas):
    for j, th2 in enumerate(thetas):
        if j <= i: continue
        DS1, DS2 = all_DS[th1], all_DS[th2]
        if DS1[0] >= DS2[0]: continue  # 确保 th1 的初始 DS 更小
        cross = find_crossings(t_arr, DS1, DS2)
        dDS = DS2[0] - DS1[0]
        if cross:
            tm = cross[0]
            nearest_dqpt = min(dqpt_times, key=lambda x: abs(x-tm)) if dqpt_times else None
            mpemba_data.append({
                'theta1': float(th1), 'theta2': float(th2),
                'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0]),
                'delta_DS': float(dDS), 't_M': float(tm),
                'nearest_dqpt': float(nearest_dqpt) if nearest_dqpt else None,
                'diff_to_dqpt': float(abs(tm - nearest_dqpt)) if nearest_dqpt else None
            })

print(f'共找到 {len(mpemba_data)} 对 Mpemba 交叉', flush=True)
print('\n按 delta_DS 排序的 t_M（前20对）:')
mpemba_sorted = sorted(mpemba_data, key=lambda x: x['delta_DS'])[:20]
for d in mpemba_sorted:
    nd = f"{d['nearest_dqpt']:.3f}" if d['nearest_dqpt'] else 'N/A'
    print(f"  dDS={d['delta_DS']:.4f}: theta={d['theta1']:.3f} vs {d['theta2']:.3f}, t_M={d['t_M']:.3f}, nearest_DQPT={nd}")

results = {
    'params': {'N':N,'NA':NA,'g_f':g_f,'t_max':t_max,'n_t':n_t},
    't': t_arr.tolist(),
    'thetas': thetas,
    'DS0_values': {str(th): float(ds0) for th,ds0 in all_DS0.items()},
    'dqpt_times': dqpt_times,
    'mpemba_pairs': mpemba_data,
    'all_DS': {str(th): ds.tolist() for th,ds in all_DS.items()}
}
os.makedirs('results', exist_ok=True)
json.dump(results, open('results/exp_b_fine_theta.json','w'), indent=2)
print('\n实验 B 完成')
