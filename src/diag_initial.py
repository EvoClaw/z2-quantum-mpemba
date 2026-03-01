import numpy as np
import sys
sys.path.insert(0, '/home/yanlin/wuli/src')
from free_fermion_ea import build_TFIM_BdG, solve_BdG, initial_correlator_full_BdG

N = 100; NA = 4; theta = 0.4
lam = np.sin(theta/2)**2
alpha = np.sin(theta)/2
beta = np.cos(theta)
print("=== 初始态问题诊断 ===")
print(f"theta={theta}, lam={lam:.6f}, cos(theta)={beta:.6f}")
print(f"正确P_pi = cos^{NA}(theta) = {beta**NA:.6f}")
print(f"公式给出 = (1-2lam)^{NA} = {(1-2*lam)**NA:.6f} (应相等)")
print()

# 正确的 C_A 矩阵 (包含 JW 跨格关联)
C_correct = np.zeros((NA, NA))
for i in range(NA):
    C_correct[i,i] = lam
    for j in range(i+1, NA):
        m = j - i
        C_correct[i,j] = alpha**2 * beta**(m-1)
        C_correct[j,i] = alpha**2 * beta**(m-1)
print("C_A(0) 正确矩阵 (含JW跨格关联):")
print(np.round(C_correct, 5))
eig_C = np.linalg.eigvalsh(C_correct)
print(f"C_A eigenvalues: {np.sort(eig_C)}")
P_pi_from_C = np.prod(1 - 2*eig_C)
print(f"P_pi from det(1-2C_A): {P_pi_from_C:.6f}")
print()

# 在 TL 中: 每个 k 模式的 2x2 BdG 密度矩阵
k_vals = np.linspace(0, 2*np.pi, 2000, endpoint=False)
C_k = np.zeros(len(k_vals))
F_k = np.zeros(len(k_vals), dtype=complex)
for k_idx, k in enumerate(k_vals):
    C_k[k_idx] = lam
    F_k_val = 0j
    for m in range(1, 300):
        c_m = alpha**2 * beta**(m-1)
        C_k[k_idx] += 2 * c_m * np.cos(k*m)
        F_k_val += (-c_m) * np.exp(-1j*k*m)
        F_k_val += c_m * np.exp(1j*k*m)
    F_k[k_idx] = F_k_val

print("=== k-空间初始态分布 ===")
print(f"C_k 范围: [{np.min(C_k):.4f}, {np.max(C_k):.4f}]")
print(f"|F_k| 范围: [{np.min(np.abs(F_k)):.4f}, {np.max(np.abs(F_k)):.4f}]")
print(f"F_k 实部最大值 (应为0): {np.max(np.abs(F_k.real)):.2e}")

disc = np.sqrt((C_k - 0.5)**2 + np.abs(F_k)**2)
mu_plus = 0.5 + disc
mu_minus = 0.5 - disc
print(f"Gamma_k eigenvalues: mu+ range [{np.min(mu_plus):.4f},{np.max(mu_plus):.4f}]")
print(f"                     mu- range [{np.min(mu_minus):.4f},{np.max(mu_minus):.4f}]")

# 关键结论
print()
print("=== 关键结论 ===")
print("乘积态在 BdG Nambu 空间有 C_k 依赖 k 且 F_k!=0 (虚数)")
print("当前代码误用 C_k=lam (常数), F_k=0 -> 这对应不同的初始态")
print()
print("但是: 对于 Z2 EA 计算来说")
print("  - P_pi = prod(1-2*eig(C_A)) = (1-2*lam)^NA = cos^NA(theta) ✓ (正确!)")
print("  - 真正的问题: S_rho 计算错误 (应为0, 实际给出NA*h(lam))")
print()
print("根本原因: 乘积自旋态是纯态(S=0), 但 BdG 公式把它当成混合态")
print("在热力学极限下才能用 Toeplitz/Gaussian 近似")
print()
print("解决方案: 从 t=0 的正确 DeltaS 出发")
p_plus = (1 + beta**NA)/2
p_minus = (1 - beta**NA)/2
DS0 = -p_plus*np.log(max(p_plus,1e-15)) - p_minus*np.log(max(p_minus,1e-15))
print(f"t=0 正确 DeltaS = H(p+,p-) = {DS0:.4f}  (S_rho=0 纯态)")
print()
print("在 TL (N,NA大) 下, 正确的 TL 近似:")
h = lambda x: -x*np.log(np.clip(x,1e-15,1-1e-15))-(1-x)*np.log(np.clip(1-x,1e-15,1-1e-15))
S_TL = NA * np.mean(0.5*(h(mu_plus)+h(mu_minus)))
print(f"S_rho(TL, NA={NA}) = {S_TL:.4f}")
