"""实验 B（修正版）：精细 theta 网格 + 严格的 t_M-DQPT 关联分析

修正：
1. 用分数位置检验（fractional position test）替代平凡 ±1 测试
2. 检测 ΔS(t) 的局部极小值是否落在 DQPT 时刻
3. Kolmogorov-Smirnov 检验 t_M 分布的非均匀性
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os
from scipy.stats import kstest, uniform
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

N, NA = 12, 4
g_f = 2.0
# g_c=0.5, 能隙 Δ=2(g-g_c)=3.0, 理论 DQPT 半周期 T*/2=π/(2*1.5)=1.047
g_c = 0.5
Delta_gap = 2*(g_f - g_c)
T_half = np.pi / Delta_gap          # 理论 DQPT 半周期 ≈ 1.047
T_period = 2 * T_half               # DQPT 完整周期 ≈ 2.094

t_max, n_t = 20.0, 200
t_arr = np.linspace(0, t_max, n_t)
thetas = np.linspace(0.05, np.pi/2 - 0.05, 40).tolist()

print(f'实验 B（修正版）: N={N}, NA={NA}, g_f={g_f}, g_c={g_c}')
print(f'理论 DQPT 半周期 T*/2 = π/Δ = π/{Delta_gap:.2f} = {T_half:.4f}')
print(f'理论 DQPT 完整周期 T* = {T_period:.4f}')
print(f'theta 网格: {len(thetas)} 点', flush=True)

H = build_TFIM_sparse(N, g=g_f)

# 计算所有 theta 的 DeltaS(t) 和 echo
all_DS = {}; all_echo = {}
for theta in thetas:
    DS,_,_,echo = run_ea(N, NA, theta, H, t_arr)
    all_DS[theta] = DS
    all_echo[theta] = echo
    print(f'  theta={theta:.3f}: DS(0)={DS[0]:.4f} DS(end)={DS[-1]:.4f}', flush=True)

# DQPT 参考时间（theta=0.4 附近）
theta_ref = min(thetas, key=lambda x: abs(x-0.4))
dqpt_times = find_dqpt_times(t_arr, all_echo[theta_ref], N, threshold=0.01)
print(f'\n数值 DQPT 时间（theta_ref={theta_ref:.3f}）: {[f"{x:.3f}" for x in dqpt_times]}', flush=True)

# 估算数值 DQPT 周期（相邻 DQPT 时间间隔的平均值）
T_period_num = None
if len(dqpt_times) >= 2:
    gaps = [dqpt_times[k+1]-dqpt_times[k] for k in range(len(dqpt_times)-1)]
    T_period_num = float(np.mean(gaps))
    print(f'数值 DQPT 周期 T*_num = {T_period_num:.4f}  (理论 {T_period:.4f})', flush=True)

# 计算所有 Mpemba 对
print('\n计算 Mpemba 交叉时间...', flush=True)
mpemba_data = []
for i, th1 in enumerate(thetas):
    for j, th2 in enumerate(thetas):
        if j<=i: continue
        DS1, DS2 = all_DS[th1], all_DS[th2]
        if DS1[0] >= DS2[0]: continue
        cross = find_crossings(t_arr, DS1, DS2)
        if not cross: continue
        tm = cross[0]
        dDS = DS2[0] - DS1[0]
        nearest_dqpt = min(dqpt_times, key=lambda x: abs(x-tm)) if dqpt_times else None
        diff_to_dqpt = abs(tm - nearest_dqpt) if nearest_dqpt is not None else None

        # 分数位置分析：t_M mod T*
        frac_theory  = (tm % T_period) / T_period      # 理论周期归一化
        frac_num     = (tm % T_period_num) / T_period_num if T_period_num else None
        # DQPT 时刻的分数位置应为 0.5（半周期）
        dist_to_half = abs(frac_theory - 0.5)

        mpemba_data.append({
            'theta1': float(th1), 'theta2': float(th2),
            'DS0_1': float(DS1[0]), 'DS0_2': float(DS2[0]),
            'delta_DS': float(dDS), 't_M': float(tm),
            'nearest_dqpt': float(nearest_dqpt) if nearest_dqpt is not None else None,
            'diff_to_dqpt': float(diff_to_dqpt) if diff_to_dqpt is not None else None,
            'frac_T_theory': float(frac_theory),
            'frac_T_num': float(frac_num) if frac_num is not None else None,
            'dist_to_half': float(dist_to_half)
        })

print(f'共找到 {len(mpemba_data)} 对 Mpemba 交叉', flush=True)

# ---- 严格统计检验 ----
if mpemba_data:
    tMs = np.array([d['t_M'] for d in mpemba_data])
    fracs = np.array([d['frac_T_theory'] for d in mpemba_data])
    dists = np.array([d['dist_to_half'] for d in mpemba_data])
    diffs = np.array([d['diff_to_dqpt'] for d in mpemba_data if d['diff_to_dqpt'] is not None])

    print('\n===== 严格统计分析 =====')
    print(f'(1) 平凡测试 |t_M - t_DQPT| < ε:')
    for eps in [0.5, 1.0, 2.0]:
        frac_hit = np.mean(diffs < eps) if len(diffs)>0 else 0
        # 随机基准：t_M 均匀分布在 [0, t_max] 内，最近 DQPT 时间
        n_dqpt = len(dqpt_times)
        if n_dqpt>0 and t_max>0:
            # 期望概率：每个 DQPT 附近 2ε 宽度，占总时间比例
            p_random = min(1.0, n_dqpt * 2*eps / t_max)
        else:
            p_random = 0
        print(f'    ε={eps:.1f}: 命中率={frac_hit:.1%}  随机基准={p_random:.1%}  '
              f'{"可能有意义" if frac_hit > 1.5*p_random else "统计上平凡"}')

    print(f'\n(2) 分数位置检验（t_M mod T*）/ T*：')
    print(f'    DQPT 时刻对应分数位置 = 0.5')
    print(f'    观测均值 = {fracs.mean():.3f} ± {fracs.std():.3f}')
    print(f'    |frac - 0.5| 均值 = {dists.mean():.3f}  (均匀分布期望 = 0.25)')

    # KS 检验（是否偏离均匀分布）
    stat, pval = kstest(fracs, 'uniform')
    print(f'\n(3) KS 检验 t_M mod T* 是否均匀：')
    print(f'    KS 统计量 = {stat:.4f}, p-value = {pval:.4f}')
    if pval < 0.05:
        print(f'    结论：分布显著非均匀 (p={pval:.4f} < 0.05)')
        peak = fracs.mean()
        print(f'    峰值位置 ≈ {peak:.3f} (DQPT 时刻对应 0.5)')
    else:
        print(f'    结论：不能排除均匀分布假设 (p={pval:.4f} ≥ 0.05)')
        print(f'    → t_M 与 DQPT 的对应关系在统计上不显著')

    # 检验 ΔS(t) 的局部极小值是否落在 DQPT 时刻
    print(f'\n(4) ΔS(t) 局部极小值 vs DQPT 时间 (theta=0.4 附近):')
    theta_demo = min(thetas, key=lambda x: abs(x-0.4))
    DS_demo = all_DS[theta_demo]
    ds_minima = []
    for k in range(1, len(DS_demo)-1):
        if DS_demo[k] < DS_demo[k-1] and DS_demo[k] < DS_demo[k+1]:
            ds_minima.append(t_arr[k])
    print(f'    ΔS 极小值时刻: {[f"{x:.3f}" for x in ds_minima[:6]]}')
    print(f'    DQPT 时刻:       {[f"{x:.3f}" for x in dqpt_times[:6]]}')
    if dqpt_times and ds_minima:
        for tm_min in ds_minima[:4]:
            nd = min(dqpt_times, key=lambda x: abs(x-tm_min))
            print(f'    ΔS 极小 {tm_min:.3f} → 最近 DQPT {nd:.3f}, 差 {abs(tm_min-nd):.3f}')

results = {
    'params': {'N':N,'NA':NA,'g_f':g_f,'g_c':g_c,'t_max':t_max,'n_t':n_t,
               'T_half_theory':T_half, 'T_period_theory':T_period,
               'T_period_num': T_period_num},
    't': t_arr.tolist(), 'thetas': thetas,
    'dqpt_times': dqpt_times,
    'mpemba_pairs': mpemba_data,
    'all_DS': {str(th): ds.tolist() for th,ds in all_DS.items()},
    'all_echo': {str(th): e.tolist() for th,e in all_echo.items()}
}
os.makedirs('results', exist_ok=True)
json.dump(results, open('results/exp_b_v2.json','w'), indent=2)
print('\n实验 B（修正版）完成，结果保存至 results/exp_b_v2.json')
