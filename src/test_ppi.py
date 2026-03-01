"""
测试正确初始态下 P_pi 和 DeltaS 的计算。

关键发现总结:
- 乘积自旋态 |ψ_θ> 不是 BdG Gaussian 态 (有非零单粒子期望值 <c_j>!=0)
- 对任意 BdG Gaussian 态: Q_A ρ_A Q_A = ρ_A -> DeltaS = 0 恒成立
- 正确的 DeltaS 来自 <c_j>(t) != 0 的贡献

方法: 用精确密度矩阵法 (小 NA) 验证正确的 DeltaS.
"""
import numpy as np
import sys
sys.path.insert(0, '/home/yanlin/wuli/src')
from free_fermion_ea import build_TFIM_BdG, solve_BdG, evolve_correlator, initial_correlator_full_BdG

def build_exact_rho_A(N_small, NA, theta, g_f, t):
    """
    用精确密度矩阵法计算 rho_A(t) (对较小的 N=N_small).
    
    初始态: prod_j (cos(th/2)|0> + sin(th/2)|1>)_j
    动力学: 在 H_f(g_f) 下演化.
    
    Returns: rho_A (2^NA x 2^NA density matrix)
    """
    # Build full many-body Hamiltonian for BdG using Fock basis
    # This is exact for small N
    # Using: H = -sum_j (Jc_j^dag c_{j+1} + h.c. + Jc_j c_{j+1} + h.c.) - g sum_j (1-2n_j)
    
    dim = 2**N_small
    H_full = np.zeros((dim, dim), dtype=complex)
    J = 1.0
    
    def fock_action(state_idx, ops):
        """Apply operators to Fock state, return (sign, new_state) or None."""
        sign = 1
        n = [((state_idx >> (N_small - 1 - j)) & 1) for j in range(N_small)]
        for (op_type, site) in ops:
            if op_type == 'create':
                if n[site] == 1:
                    return None  # Already occupied
                # JW sign: count occupied sites to the left
                jw_sign = (-1) ** sum(n[:site])
                n[site] = 1
                sign *= jw_sign
            elif op_type == 'annihilate':
                if n[site] == 0:
                    return None  # Already empty
                jw_sign = (-1) ** sum(n[:site])
                n[site] = 0
                sign *= jw_sign
        new_idx = sum(n[j] * (2**(N_small - 1 - j)) for j in range(N_small))
        return sign, new_idx
    
    # Build H in Fock basis
    for n_state in range(dim):
        n = [(n_state >> (N_small - 1 - j)) & 1 for j in range(N_small)]
        
        # Onsite: -g*(1-2n_j) = -g + 2g*n_j
        H_full[n_state, n_state] += sum(-g_f * (1 - 2*n[j]) for j in range(N_small))
        
        for j in range(N_small):
            j1 = (j + 1) % N_small  # PBC
            
            # Hopping: -J/2 * c_j^dag c_{j+1} + h.c.
            # c_{j+1} then c_j^dag
            res = fock_action(n_state, [('annihilate', j1), ('create', j)])
            if res is not None:
                sign, new_state = res
                H_full[new_state, n_state] += -J/2 * sign
            res = fock_action(n_state, [('annihilate', j), ('create', j1)])
            if res is not None:
                sign, new_state = res
                H_full[new_state, n_state] += -J/2 * sign
            
            # Pairing: -J/2 * c_j c_{j+1} + h.c.
            # c_{j+1} then c_j (create both)
            res = fock_action(n_state, [('annihilate', j1), ('annihilate', j)])
            if res is not None:
                sign, new_state = res
                H_full[new_state, n_state] += -J/2 * sign
            res = fock_action(n_state, [('create', j), ('create', j1)])
            if res is not None:
                sign, new_state = res
                H_full[new_state, n_state] += -J/2 * sign
    
    # Initial state amplitude
    cos_t = np.cos(theta/2)
    sin_t = np.sin(theta/2)
    psi0 = np.zeros(dim, dtype=complex)
    for state in range(dim):
        amp = 1.0
        for j in range(N_small):
            n_j = (state >> (N_small - 1 - j)) & 1
            if n_j == 0:
                amp *= cos_t
            else:
                amp *= sin_t
        psi0[state] = amp
    
    # Time evolve
    evals, evecs = np.linalg.eigh(H_full)
    psi_t = evecs @ (np.exp(-1j * evals * t) * (evecs.conj().T @ psi0))
    
    # Build full density matrix
    rho_full = np.outer(psi_t, psi_t.conj())
    
    # Partial trace over B (last N-NA sites)
    dim_A = 2**NA
    dim_B = 2**(N_small - NA)
    rho_A = np.zeros((dim_A, dim_A), dtype=complex)
    
    for i_A in range(dim_A):
        for j_A in range(dim_A):
            val = 0.0 + 0j
            for k_B in range(dim_B):
                i_full = i_A * dim_B + k_B
                j_full = j_A * dim_B + k_B
                val += rho_full[i_full, j_full]
            rho_A[i_A, j_A] = val
    
    return rho_A


def compute_exact_EA(rho_A, NA):
    """Compute exact Z2 EA from density matrix."""
    # Q_A = prod_j (1-2n_j) = diagonal operator in Fock basis
    Q_A = np.array([(-1)**bin(k).count('1') for k in range(2**NA)], dtype=float)
    Q_mat = np.diag(Q_A)
    
    # rho_tilde = (rho_A + Q rho Q) / 2
    rho_tilde = (rho_A + Q_mat @ rho_A @ Q_mat) / 2
    
    # Entropies
    def entropy(rho):
        evals = np.linalg.eigvalsh(rho).real
        evals = evals[evals > 1e-12]
        return -np.sum(evals * np.log(evals))
    
    S_rho = entropy(rho_A)
    S_tilde = entropy(rho_tilde)
    
    # P_pi = Tr[rho_A Q_A]
    P_pi = np.real(np.trace(rho_A @ Q_mat))
    
    return S_tilde - S_rho, S_rho, S_tilde, P_pi


# === MAIN TEST ===
NA = 4
N_small = 10  # Small enough for exact computation
theta = 0.4
g_f = 2.0

print(f"=== 精确密度矩阵法验证 ===")
print(f"N={N_small}, NA={NA}, theta={theta}, g_f={g_f}")
print()

# Exact values at t=0
P_pi_exact = np.cos(theta)**NA
p_plus = (1 + P_pi_exact) / 2
p_minus = (1 - P_pi_exact) / 2
DS0_exact = -p_plus*np.log(p_plus) - p_minus*np.log(p_minus)
print(f"理论 t=0: P_pi={P_pi_exact:.5f}, DeltaS={DS0_exact:.5f}")
print()

print("时间演化:")
for t in [0.0, 1.0, 2.0, 5.0, 10.0]:
    rho_A = build_exact_rho_A(N_small, NA, theta, g_f, t)
    DS, S_rho, S_tilde, P_pi = compute_exact_EA(rho_A, NA)
    print(f"  t={t:.1f}: DeltaS={DS:.5f}, S_rho={S_rho:.5f}, P_pi={P_pi:.5f}")
    if t == 0.0:
        print(f"    (理论t=0: DeltaS={DS0_exact:.5f}, P_pi={P_pi_exact:.5f})")
