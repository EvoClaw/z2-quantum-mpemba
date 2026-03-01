import sys
import numpy as np
from functools import reduce
sys.path.insert(0, '/home/yanlin/wuli/src')
from free_fermion_ea import entanglement_asymmetry, entropy_from_correlator
from free_fermion_kspace import compute_Ck_trajectory, build_toeplitz_from_Ck


def exact_EA_small(C_A):
    NA = C_A.shape[0]
    dim = 2**NA
    lam, U = np.linalg.eigh(C_A)
    lam = np.clip(lam.real, 1e-12, 1-1e-12)
    rho_modes = reduce(np.kron, [np.diag([1-lam[k], lam[k]]) for k in range(NA)])
    V_full = np.zeros((dim,dim),dtype=complex)
    for alpha in range(dim):
        for beta in range(dim):
            oa = [j for j in range(NA) if (alpha>>j)&1]
            ob = [j for j in range(NA) if (beta>>j)&1]
            if len(oa) != len(ob):
                V_full[alpha,beta] = 0.0
            elif len(oa) == 0:
                V_full[alpha,beta] = 1.0
            else:
                V_full[alpha,beta] = np.linalg.det(U[np.ix_(oa,ob)])
    rho_orig = V_full @ rho_modes @ V_full.conj().T
    Z_eig = np.array([(-1)**bin(a).count('1') for a in range(dim)], dtype=float)
    Z_mat = np.diag(Z_eig)
    rho_tilde = (rho_orig + Z_mat @ rho_orig @ Z_mat) / 2
    eig_t = np.linalg.eigvalsh(rho_tilde).real
    eig_t = eig_t[eig_t > 1e-15]; eig_t /= eig_t.sum()
    S_tilde = -np.sum(eig_t * np.log(eig_t))
    S_rho = entropy_from_correlator(C_A)
    P_pi = np.prod(1-2*lam)
    return S_tilde - S_rho, S_rho, S_tilde, P_pi


print("Test 1: Initial product state (C_A = lambda*I)")
NA = 4
for theta in [0.4, 0.8, 1.2, np.pi/4]:
    lam0 = np.sin(theta/2)**2
    C0 = lam0 * np.eye(NA)
    DS_e, Sr_e, St_e, Pp_e = exact_EA_small(C0)
    DS_m, _, _, Pp_m, _, _ = entanglement_asymmetry(C0)
    print("  theta=%.3f: Exact DS=%.6f S_rho=%.4f S_tilde=%.4f Ppi=%.4f  My DS=%.6f" % (
        theta, DS_e, Sr_e, St_e, Pp_e, DS_m))

print()
print("Test 2: After quench (t=5, g_f=2.0)")
N = 100
t_arr = np.array([0.0, 2.0, 5.0, 10.0])
for theta in [0.4, 0.8]:
    lam0 = np.sin(theta/2)**2
    Ckt, _ = compute_Ck_trajectory('TFIM', {'g': 2.0}, lam0, N, t_arr)
    print("  theta=%.2f:" % theta)
    for idx, t in enumerate(t_arr):
        TA = build_toeplitz_from_Ck(Ckt[idx], NA)
        eig_C = np.sort(np.linalg.eigvalsh(TA))
        DS_e, Sr_e, St_e, Pp_e = exact_EA_small(TA)
        DS_m, _, _, Pp_m, _, _ = entanglement_asymmetry(TA)
        print("    t=%.1f: eig=[%.3f,%.3f]  Exact DS=%.5f S_r=%.4f S_t=%.4f  My DS=%.5f Ppi=%.5f" % (
            t, eig_C[0], eig_C[-1], DS_e, Sr_e, St_e, DS_m, Pp_e))

print()
print("Test 3: Synthetic non-uniform state")
lam_test = np.array([0.01, 0.1, 0.5, 0.9])
CA_diag = np.diag(lam_test)
DS_e, Sr_e, St_e, Pp_e = exact_EA_small(CA_diag)
DS_m, _, _, Pp_m, _, _ = entanglement_asymmetry(CA_diag)
print("  Diagonal lam=[0.01,0.1,0.5,0.9]: Exact DS=%.6f S_r=%.4f S_t=%.4f  My DS=%.6f Ppi=%.4f" % (
    DS_e, Sr_e, St_e, DS_m, Pp_e))
