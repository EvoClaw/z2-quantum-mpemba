"""实验 A：对比 g_f=2.0（有 DQPT）和 g_f=0.5（无 DQPT）下的量子 Mpemba 效应"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

N, NA = 12, 4
t_max, n_t = 20.0, 80
t_arr = np.linspace(0, t_max, n_t)
thetas = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
gf_list = [0.5, 2.0]   # 0.5=无DQPT, 2.0=有DQPT

results = {'params': {'N':N,'NA':NA,'t_max':t_max,'n_t':n_t},
           't': t_arr.tolist()}

for g_f in gf_list:
    label = f'gf{g_f}'
    print(f'\n=== g_f={g_f} ===', flush=True)
    H = build_TFIM_sparse(N, g=g_f)
    results[label] = {}
    for theta in thetas:
        DS,Srho,Ppi,echo = run_ea(N, NA, theta, H, t_arr)
        results[label][str(theta)] = {
            'DS':DS.tolist(),'Srho':Srho.tolist(),'Ppi':Ppi.tolist(),'echo':echo.tolist()
        }
        print(f'  theta={theta:.1f}: DS(0)={DS[0]:.4f} DS({t_max:.0f})={DS[-1]:.4f}', flush=True)

    # Mpemba 交叉分析
    print(f'  Mpemba 交叉时间 (g_f={g_f}):')
    for i,th1 in enumerate(thetas):
        for j,th2 in enumerate(thetas):
            if j<=i: continue
            DS1=np.array(results[label][str(th1)]['DS'])
            DS2=np.array(results[label][str(th2)]['DS'])
            if DS1[0]>=DS2[0]: continue
            cross=find_crossings(t_arr,DS1,DS2)
            if cross: print(f'    theta={th1} vs {th2}: t_M={cross[0]:.3f}', flush=True)

    # DQPT 分析（用 theta=0.4 的 echo）
    echo04=np.array(results[label]['0.4']['echo'])
    dqpt=find_dqpt_times(t_arr, echo04, N)
    print(f'  DQPT times (g_f={g_f}): {[f"{x:.3f}" for x in dqpt]}', flush=True)
    results[label]['dqpt_times'] = dqpt

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/exp_a_gf_compare.json','w'), indent=2)
print('\n实验 A 完成，结果保存至 results/exp_a_gf_compare.json')
