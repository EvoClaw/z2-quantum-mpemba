"""
快速精确 Z2 纠缠非对称度计算
稀疏矩阵 + scipy.sparse.linalg.expm_multiply，支持 N 到 16
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eigvalsh
import time


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
                n2 = n.copy(); n2[j1] = 0  # FIX
                sign = jw_left(n, j1) * jw_left(n2, j)
                new_s = (s ^ (1 << (N-1-j))) ^ (1 << (N-1-j1))
                rows.append(new_s); cols.append(s); vals.append(-J/2 * sign)
            if n[j] == 0 and n[j1] == 0:
                n2 = n.copy(); n2[j] = 1
                sign = jw_left(n, j) * jw_left(n2, j1)
                new_s = s | (1 << (N-1-j)) | (1 << (N-1-j1))
                rows.append(new_s); cols.append(s); vals.append(-J/2 * sign)
    return sp.csr_matrix((vals, (rows, cols)), shape=(dim, dim))


def build_initial_state(N, theta):
    dim = 2**N
    cos_t = np.cos(theta / 2)
    sin_t = np.sin(theta / 2)
    psi0 = np.zeros(dim, dtype=complex)
    for state in range(dim):
        amp = 1.0
        for j in range(N):
            n_j = (state >> (N - 1 - j)) & 1
            amp *= sin_t if n_j else cos_t
        psi0[state] = amp
    return psi0


def partial_trace_A(psi, N, NA):
    dim_A = 2**NA
    dim_B = 2**(N - NA)
    psi_mat = psi.reshape(dim_A, dim_B)
    return psi_mat @ psi_mat.conj().T


def compute_EA_from_rho(rho_A, NA):
    dim_A = 2**NA
    Q_diag = np.array([(-1)**bin(k).count("1") for k in range(dim_A)], dtype=float)
    Q_mat = np.diag(Q_diag)
    rho_tilde = (rho_A + Q_mat @ rho_A @ Q_mat) / 2
    def vn_entropy(rho):
        evals = eigvalsh(rho).real
        evals = evals[evals > 1e-14]
        return -np.sum(evals * np.log(evals))
    S_rho = vn_entropy(rho_A)
    S_tilde = vn_entropy(rho_tilde)
    P_pi = float(np.real(np.trace(rho_A @ Q_mat)))
    return S_tilde - S_rho, S_rho, S_tilde, P_pi


def evolve_and_compute_EA(N, NA, theta, g_f, t_list, J=1.0, verbose=True):
    t_wall = time.time()
    if verbose:
        print(f"构建稀疏 H (N={N})...")
    H = build_TFIM_sparse(N, J=J, g=g_f)
    if verbose:
        print(f"  nnz={H.nnz}, 耗时 {time.time()-t_wall:.2f}s")
    psi0 = build_initial_state(N, theta)
    results = {"t": list(t_list), "DeltaS": [], "S_rho": [], "S_tilde": [], "P_pi": []}
    psi_t = psi0.copy().astype(complex)
    t_prev = 0.0
    for t_next in sorted(t_list):
        dt = t_next - t_prev
        if dt > 1e-12:
            psi_t = spla.expm_multiply(-1j * H * dt, psi_t)
        t_prev = t_next
        rho_A = partial_trace_A(psi_t, N, NA)
        DS, S_rho, S_tilde, P_pi = compute_EA_from_rho(rho_A, NA)
        results["DeltaS"].append(DS)
        results["S_rho"].append(S_rho)
        results["S_tilde"].append(S_tilde)
        results["P_pi"].append(P_pi)
        if verbose:
            print(f"  t={t_next:.3f}: DeltaS={DS:.5f}, S_rho={S_rho:.5f}, P_pi={P_pi:.5f}")
    if verbose:
        print(f"总耗时: {time.time()-t_wall:.2f}s")
    return results


def compute_loschmidt(N, theta, g_f, t_list, J=1.0):
    H = build_TFIM_sparse(N, J=J, g=g_f)
    psi0 = build_initial_state(N, theta).astype(complex)
    echos, losch = [], []
    psi_t = psi0.copy()
    t_prev = 0.0
    for t in sorted(t_list):
        dt = t - t_prev
        if dt > 1e-12:
            psi_t = spla.expm_multiply(-1j * H * dt, psi_t)
        t_prev = t
        echo = abs(np.dot(psi0.conj(), psi_t))**2
        echos.append(echo)
        losch.append(-np.log(max(echo, 1e-300)) / N)
    return np.array(echos), np.array(losch)
