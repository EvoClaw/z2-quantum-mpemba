"""实验 A（修正版）：对比 g_f=0.2（铁磁相，g<g_c）和 g_f=2.0（顺磁相，g>g_c）

关键修正：g_c = 0.5（Hamiltonian 包含 J/2 因子），
原 g_f=0.5 恰好在临界点，应改为 g_f=0.2（真正的铁磁相）。
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
from exp_common import build_TFIM_sparse, run_ea, find_crossings, find_dqpt_times

N, NA = 12, 4
t_max, n_t = 20.0, 200
t_arr = np.linspace(0, t_max, n_t)
thetas = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

# g_c = 0.5（J=1，耦合项为 J/2）
# g_f=0.2: 铁磁相（g < g_c）
# g_f=2.0: 顺磁相（g >> g_c）
gf_list = [0.2, 2.0]

print(f'实验 A（修正版）: N={N}, NA={NA}, g_c=0.5')
print(f'g_f 列表: {gf_list}  (0.2=铁磁相, 2.0=顺磁相)')

results = {'params': {'N':N,'NA':NA,'t_max':t_max,'n_t':n_t,'g_c':0.5},
           't': t_arr.tolist()}

for g_f in gf_list:
    label = f'gf{g_f}'
    phase = '铁磁相(g<g_c)' if g_f < 0.5 else '顺磁相(g>g_c)'
    print(f'\n=== g_f={g_f} ({phase}) ===', flush=True)
    H = build_TFIM_sparse(N, g=g_f)

    # 验证 Hermitian
    diff = (H - H.conj().T).toarray()
    print(f'  ||H-H†||_max = {np.abs(diff).max():.2e}  (应≈0)', flush=True)

    results[label] = {'phase': phase, 'g_f': g_f}
    for theta in thetas:
        DS,Srho,Ppi,echo = run_ea(N, NA, theta, H, t_arr)
        # 验证范数保持
        results[label][str(theta)] = {
            'DS':DS.tolist(),'Srho':Srho.tolist(),'Ppi':Ppi.tolist(),'echo':echo.tolist()
        }
        print(f'  theta={theta:.1f}: DS(0)={DS[0]:.4f}  DS_min={DS.min():.4f}  DS_mean={DS.mean():.4f}', flush=True)

    # Mpemba 交叉分析
    print(f'\n  Mpemba 交叉时间 (g_f={g_f}):')
    mpemba_list = []
    for i, th1 in enumerate(thetas):
        for j, th2 in enumerate(thetas):
            if j<=i: continue
            DS1=np.array(results[label][str(th1)]['DS'])
            DS2=np.array(results[label][str(th2)]['DS'])
            if DS1[0]>=DS2[0]: continue
            cross=find_crossings(t_arr,DS1,DS2)
            if cross:
                print(f'    theta={th1} vs {th2}: t_M={cross[0]:.3f}', flush=True)
                mpemba_list.append({'theta1':th1,'theta2':th2,'t_M':cross[0]})
    results[label]['mpemba_pairs'] = mpemba_list

    # DQPT 分析（用多个 theta 的 echo 叠加确认）
    all_dqpt = {}
    for theta in [0.4, 0.8, 1.2]:
        echo=np.array(results[label][str(theta)]['echo'])
        dqpt=find_dqpt_times(t_arr, echo, N, threshold=0.01)
        all_dqpt[str(theta)] = dqpt
        print(f'  DQPT times (g_f={g_f}, theta={theta}): {[f"{x:.3f}" for x in dqpt]}', flush=True)
    results[label]['dqpt_by_theta'] = all_dqpt

    # Loschmidt 速率函数
    echo04=np.array(results[label]['0.4']['echo'])
    losch = -np.log(np.maximum(echo04,1e-300))/N
    results[label]['losch_rate_0.4'] = losch.tolist()

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/exp_a_v2.json','w'), indent=2)
print('\n实验 A（修正版）完成，结果保存至 results/exp_a_v2.json')
