"""
Phase 4a 探索性实验：量子 Mpemba 效应 vs DQPT
N=12, NA=4, 扫描不同初始角度 theta
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eigvalsh
import time, json, os

# ======= 内联精确方法 =======

def build_TFIM_sparse(N, J=1.0, g=2.0, pbc=True):
    dim = 2**N
    rows, cols, vals = [], [], []
    for s in range(dim):
        n = [(s >> (N - 1 - j)) & 1 for j in range(N)]
        diag_val = sum(-g * (1 - 2 * n[j]) for j in range(N))
        rows.append(s); cols.append(s); vals.append(diag_val)
        for j in range(N):
            j1 = (j + 1) % N if pbc else (j + 1 if j + 1 < N else None)
            if j1 is None:
                continue
            def jw_left(nb, site):
                return (-1) ** sum(nb[:site])
            if n[j1] == 1 and n[j] == 0:
                n2 = n.copy(); n2[j1] = 0
                sign = jw_left(n, j1) * jw_left(n2, j)
                new_s = (s ^ (1 << (N-1-j1))) | (1 << (N-1-j))
                rows.append(new_s); cols.append(s); vals.append(-J/2 * sign)
            if n[j] == 1 and n[j1] == 0:
                n2 = n.copy(); n2[j] = 0
                sign = jw_left(n, j) * jw_left(n2, j1)
                new_s = (s ^ (1 << (N-1-j))) | (1 << (N-1-j1))
                rows.append(new_s); cols.append(s); vals.append(-J/2 * sign)
            if n[j] == 1 and n[j1] == 1:
                n2 = n.copy(); n2[j1] = 0  # FIX: 先湮灭j1
                sign = jw_left(n, j1) * jw_left(n2, j)
                new_s = (s ^ (1 << (N-1-j))) ^ (1 << (N-1-j1))
                rows.append(new_s); cols.append(s); vals.append(-J/2 * sign)
            if n[j] == 0 and n[j1] == 0:
                n2 = n.copy(); n2[j] = 1
                sign = jw_left(n, j) * jw_left(n2, j1)
                new_s = s | (1 << (N-1-j)) | (1 << (N-1-j1))
                rows.append(new_s); cols.append(s); vals.append(-J/2 * sign)
    return sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim))

def build_psi0(N, theta):
    dim = 2**N
    c, s = np.cos(theta/2), np.sin(theta/2)
    psi = np.zeros(dim, dtype=complex)
    for state in range(dim):
        amp = 1.0
        for j in range(N):
            n_j = (state >> (N-1-j)) & 1
            amp *= s if n_j else c
        psi[state] = amp
    return psi

def partial_trace(psi, N, NA):
    psi_mat = psi.reshape(2**NA, 2**(N-NA))
    return psi_mat @ psi_mat.conj().T

def ea_from_rho(rho_A, NA):
    Q = np.diag(np.array([(-1)**bin(k).count('1') for k in range(2**NA)], dtype=float))
    rho_t = (rho_A + Q @ rho_A @ Q) / 2
    def S(r):
        ev = eigvalsh(r).real; ev = ev[ev>1e-14]
        return -np.sum(ev*np.log(ev))
    return S(rho_t) - S(rho_A), S(rho_A), float(np.real(np.trace(rho_A @ Q)))

# ======= Phase 4a 核心实验 =======

def run_experiment(N, NA, theta, g_f, t_arr, H_cached=None, verbose=False):
    if H_cached is None:
        H_cached = build_TFIM_sparse(N, g=g_f)
    psi0 = build_psi0(N, theta)
    psi_t = psi0.copy()
    t_prev = 0.0
    DS_arr, Srho_arr, Ppi_arr, echo_arr = [], [], [], []
    for t in t_arr:
        dt = t - t_prev
        if dt > 1e-12:
            psi_t = spla.expm_multiply(-1j * H_cached * dt, psi_t)
        t_prev = t
        rho_A = partial_trace(psi_t, N, NA)
        DS, Srho, Ppi = ea_from_rho(rho_A, NA)
        echo = abs(np.dot(psi0.conj(), psi_t))**2
        DS_arr.append(DS); Srho_arr.append(Srho); Ppi_arr.append(Ppi); echo_arr.append(echo)
        if verbose:
            print(f"  t={t:.2f}: DS={DS:.4f} Srho={Srho:.4f} echo={echo:.4f}")
    return np.array(DS_arr), np.array(Srho_arr), np.array(Ppi_arr), np.array(echo_arr)

# ======= 实验参数 =======
N = 12
NA = 4
g_f = 2.0   # 超临界（g_c=1.0 for TFIM）
t_max = 15.0
n_t = 60
t_arr = np.linspace(0, t_max, n_t)

# 初始角度扫描（对应不同初始 DeltaS(0)）
theta_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

print(f"=== Phase 4a 探索：N={N}, NA={NA}, g_f={g_f} ===")
print(f"时间范围: 0 到 {t_max}，{n_t} 个点")
print(f"theta 列表: {theta_list}")
print()

# 预构建 H（共享）
t0 = time.time()
print("构建稀疏 Hamiltonian...")
H = build_TFIM_sparse(N, g=g_f)
print(f"  耗时 {time.time()-t0:.1f}s，nnz={H.nnz}")
print()

all_results = {}
for theta in theta_list:
    t_start = time.time()
    p0 = np.cos(theta)**NA
    pp = (1+p0)/2
    DS0_theory = -pp*np.log(pp)-(1-pp)*np.log(1-pp) if 0<pp<1 else 0.0
    print(f"theta={theta:.2f}: DS0_theory={DS0_theory:.4f}", flush=True)

    DS, Srho, Ppi, echo = run_experiment(N, NA, theta, g_f, t_arr, H_cached=H)
    all_results[str(theta)] = {
        'DS': DS.tolist(),
        'Srho': Srho.tolist(),
        'Ppi': Ppi.tolist(),
        'echo': echo.tolist(),
        'DS0_theory': DS0_theory,
    }
    print(f"  完成，耗时 {time.time()-t_start:.1f}s，DS[0]={DS[0]:.4f}，DS[-1]={DS[-1]:.4f}")

# 同时保存 Loschmidt echo（已在上面计算）
all_results['t'] = t_arr.tolist()
all_results['params'] = {'N': N, 'NA': NA, 'g_f': g_f, 't_max': t_max, 'n_t': n_t}

out_path = '/home/yanlin/wuli/results/phase4a_raw.json'
os.makedirs('/home/yanlin/wuli/results', exist_ok=True)
with open(out_path, 'w') as f:
    json.dump(all_results, f, indent=2)
print(f"\n结果保存至: {out_path}")

# ======= 找 Mpemba 交叉时间 =======
print("\n=== Mpemba 效应分析 ===")
print("初始 DeltaS(0) 排序:")
for theta in theta_list:
    DS = np.array(all_results[str(theta)]['DS'])
    print(f"  theta={theta:.2f}: DS(0)={DS[0]:.4f}")

# 找两条曲线的交叉
print("\n交叉时间矩阵 t_M(theta_i, theta_j)（DS_i(t_M) = DS_j(t_M)）:")
for i, th1 in enumerate(theta_list):
    for j, th2 in enumerate(theta_list):
        if j <= i:
            continue
        DS1 = np.array(all_results[str(th1)]['DS'])
        DS2 = np.array(all_results[str(th2)]['DS'])
        diff = DS1 - DS2
        # 找符号变化
        crossings = []
        for k in range(len(diff)-1):
            if diff[k] * diff[k+1] < 0:
                t_cross = t_arr[k] + abs(diff[k])/(abs(diff[k])+abs(diff[k+1])) * (t_arr[k+1]-t_arr[k])
                crossings.append(t_cross)
        if crossings:
            print(f"  theta={th1:.2f} vs theta={th2:.2f}: t_M={crossings[0]:.3f}")
        else:
            print(f"  theta={th1:.2f} vs theta={th2:.2f}: 无交叉")

# Loschmidt echo 零点（DQPT 临界时间）
print("\nDQPT 分析（Loschmidt echo 零点）:")
theta_ref = 0.4  # 用参考 theta 分析 DQPT
echo_ref = np.array(all_results[str(theta_ref)]['echo'])
log_echo = -np.log(np.maximum(echo_ref, 1e-300)) / N
print(f"  theta={theta_ref}, -log(L)/N 极大值位置:")
for k in range(1, len(log_echo)-1):
    if log_echo[k] > log_echo[k-1] and log_echo[k] > log_echo[k+1]:
        print(f"    t_DQPT ≈ {t_arr[k]:.3f}（极大值 {log_echo[k]:.4f}）")

print("\n=== Phase 4a 完成 ===")
