"""
free_fermion_kspace.py
======================
Fast k-space computation for translationally invariant (product state) initial
conditions on free-fermion chains (TFIM, XY chain).

Complexity: O(N) per time step (vs O(N^3) for full-matrix method).
The subsystem Toeplitz correlator is diagonalized in O(N_A^2).

Key: For uniform initial state with occupation lambda = sin^2(theta/2),
the BdG decouples into independent 2x2 blocks for each momentum k.
C_{ij}(t) = (1/N) sum_k C_k(t) e^{ik(i-j)}  (Toeplitz structure)
"""
import numpy as np
from free_fermion_ea import (
    entanglement_asymmetry, entropy_from_correlator,
    _binary_entropy, _r_k, find_mpemba_crossing_time
)


# ── 2x2 BdG k-space Hamiltonians ─────────────────────────────────────────────


def tfim_bdg_k(k, g, J=1.0):
    """2x2 BdG for TFIM at momentum k (PBC)."""
    A_k = -g - J * np.cos(k)
    B_k = -1j * J * np.sin(k)
    return np.array([[A_k, B_k], [B_k, -A_k]])


def xy_bdg_k(k, h, gamma, J=1.0):
    """2x2 BdG for anisotropic XY chain at momentum k (PBC)."""
    A_k = -h / 2.0 - J * np.cos(k)
    B_k = 1j * gamma * J * np.sin(k)
    return np.array([[A_k, B_k], [B_k, -A_k]])


# ── Core k-space computation ──────────────────────────────────────────────────


def compute_Ck_trajectory(model, params_f, lambda_0, N, t_array):
    """
    Compute C_k(t) = <c_k^dag c_k>(t) for all k-modes and all times.

    Returns
    -------
    C_k_t      : (NT, N) real array  — momentum-space occupation
    G_Loschmidt: (NT,) complex array — Loschmidt amplitude
    """
    NT = len(t_array)
    k_arr = 2.0 * np.pi * np.arange(N) / N    # PBC momenta
    C_k_t = np.zeros((NT, N))
    G_Loschmidt = np.ones(NT, dtype=complex)

    for m, k in enumerate(k_arr):
        if model == 'TFIM':
            H2 = tfim_bdg_k(k, params_f['g'])
        elif model == 'XY':
            H2 = xy_bdg_k(k, params_f['h'], params_f['gamma'])
        else:
            raise ValueError(f"Unknown model: {model}")

        eps2, U2 = np.linalg.eigh(H2)    # eps2[0] < 0, eps2[1] >= 0

        Gamma_0 = np.diag([lambda_0, 1.0 - lambda_0])
        M = U2.conj().T @ Gamma_0 @ U2    # constant 2x2 matrix

        for idx, t in enumerate(t_array):
            phase2 = np.exp(-1j * eps2 * t)
            G_row0 = U2[0, :] * phase2         # (G(t))[0,:] = U2[0,:] * phase
            C_k_t[idx, m] = np.real(G_row0 @ M @ G_row0.conj())

            # Loschmidt amplitude: product state initial condition
            # G_k(t) = (1-lambda_0) + lambda_0 * exp(2i eps_pos t)
            G_Loschmidt[idx] *= (1.0 - lambda_0) + lambda_0 * np.exp(2j * eps2[1] * t)

    return C_k_t, G_Loschmidt


def build_toeplitz_from_Ck(C_k_row, NA):
    """
    Build NA x NA Toeplitz correlator from k-space occupations at one time.
    C_{ij} = (1/N) sum_k C_k e^{ik(i-j)} = IFFT(C_k)[|i-j|]
    Since k = 2pi*q/N and IFFT(X)[m] = (1/N) sum_k X[k] exp(2pi*i*k*m/N).
    """
    c_vec = np.real(np.fft.ifft(C_k_row))    # c_vec[m] = C_{0,m}
    T = np.empty((NA, NA))
    for i in range(NA):
        for j in range(NA):
            idx = (i - j) % len(c_vec)
            T[i, j] = c_vec[idx]
    return T


def ea_from_eigenvalues(lam):
    """Compute Z2 EA and related quantities from eigenvalues of C_A."""
    lam_c = np.clip(lam, 1e-15, 1 - 1e-15)

    # Entropy
    a_k = lam_c * np.log(lam_c) + (1 - lam_c) * np.log(1 - lam_c)
    S_rho = -np.sum(a_k)

    # Parity
    P_pi = np.prod(1 - 2 * lam_c)

    # Symmetry-resolved entropies
    r_k_vals = _r_k(lam_c)
    dn_P0 = np.sum(a_k)        # = -S_rho
    dn_Ppi = P_pi * np.sum(r_k_vals)

    p_plus = (1 + P_pi) / 2
    p_minus = (1 - P_pi) / 2

    if p_plus < 1e-12:
        S_plus = 0.0
    else:
        S_plus = -(dn_P0 + dn_Ppi) / (2 * p_plus) + np.log(p_plus)

    if p_minus < 1e-12:
        S_minus = 0.0
    else:
        S_minus = -(dn_P0 - dn_Ppi) / (2 * p_minus) + np.log(p_minus)

    H_pp = (-p_plus * np.log(max(p_plus, 1e-15))
            - p_minus * np.log(max(p_minus, 1e-15)))
    S_tilde = H_pp + p_plus * S_plus + p_minus * S_minus
    DeltaS = S_tilde - S_rho

    return DeltaS, S_rho, P_pi, p_plus, p_minus


def compute_EA_kspace(model, params_f, theta, N, NA, t_array):
    """
    Compute time-dependent Z2 EA and Loschmidt echo using k-space method.
    O(N * NT) complexity (fast for large N with many time steps).

    Returns
    -------
    dict: t, DeltaS, S_rho, P_pi, p_plus, p_minus, g_Loschmidt, phase_Loschmidt, L
    """
    lambda_0 = np.sin(theta / 2) ** 2
    NT = len(t_array)

    C_k_t, G_Loschmidt = compute_Ck_trajectory(model, params_f, lambda_0, N, t_array)

    DeltaS_arr = np.zeros(NT)
    S_rho_arr = np.zeros(NT)
    P_pi_arr = np.zeros(NT)
    p_plus_arr = np.zeros(NT)
    p_minus_arr = np.zeros(NT)

    for idx in range(NT):
        T = build_toeplitz_from_Ck(C_k_t[idx], NA)
        lam = np.linalg.eigvalsh(T)
        DeltaS, S_rho, P_pi, p_plus, p_minus = ea_from_eigenvalues(lam)
        DeltaS_arr[idx] = DeltaS
        S_rho_arr[idx] = S_rho
        P_pi_arr[idx] = P_pi
        p_plus_arr[idx] = p_plus
        p_minus_arr[idx] = p_minus

    L_arr = np.abs(G_Loschmidt) ** 2
    g_L = np.where(L_arr > 1e-300, -np.log(L_arr) / N, np.inf)
    phase_L = np.angle(G_Loschmidt)

    return {
        't': t_array,
        'DeltaS': DeltaS_arr,
        'S_rho': S_rho_arr,
        'P_pi': P_pi_arr,
        'p_plus': p_plus_arr,
        'p_minus': p_minus_arr,
        'g_Loschmidt': g_L,
        'phase_Loschmidt': phase_L,
        'L': L_arr,
    }


# ── Time array builder ────────────────────────────────────────────────────────


def build_time_array(NA, dt_early=0.05, dt_late=0.1,
                      t_crossover=20.0, J=1.0):
    """
    Build N_A-dependent time array per evaluation-protocol.yaml.
    t_max = max(60/J, 3*NA/(4*J))
    """
    t_max = max(60.0 / J, 3.0 * NA / (4.0 * J))
    t1 = np.arange(0, min(t_crossover, t_max) + dt_early * 0.5, dt_early)
    if t_max > t_crossover:
        t2 = np.arange(t_crossover + dt_late, t_max + dt_late * 0.5, dt_late)
        t_arr = np.concatenate([t1, t2])
    else:
        t_arr = t1
    return np.unique(t_arr[t_arr <= t_max + 1e-8])


# ── Validation against full-matrix method ────────────────────────────────────


def validate_kspace_vs_fullmatrix(N=60, NA=10, theta=np.pi/3,
                                   g_f=2.0, tol=1e-5, verbose=True):
    """
    Cross-validate k-space method against full 2N-matrix method.
    Uses TFIM model. Returns True if all checks pass.
    """
    from free_fermion_ea import (
        build_TFIM_BdG, solve_BdG, initial_correlator_full_BdG,
        evolve_correlator, entanglement_asymmetry
    )
    import time

    t_arr = np.array([0.0, 1.0, 5.0, 10.0])

    # Full-matrix
    t0 = time.time()
    HBdG_f = build_TFIM_BdG(N, g_f)
    eps_f, V_f = solve_BdG(HBdG_f)
    C_full_0 = initial_correlator_full_BdG(N, theta)
    res_full = []
    for t in t_arr:
        Cft = evolve_correlator(C_full_0, V_f, eps_f, t)
        CA = Cft[:NA, :NA].real
        DS, Sr, _, Pp, _, _ = entanglement_asymmetry(CA)
        res_full.append((DS, Sr, Pp))
    t_full = time.time() - t0

    # k-space
    t0 = time.time()
    res_k = compute_EA_kspace('TFIM', {'g': g_f}, theta, N, NA, t_arr)
    t_kspace = time.time() - t0

    if verbose:
        print(f"Cross-validation TFIM (N={N}, NA={NA}, theta={theta:.2f}, g_f={g_f}):")
        print(f"  Full-matrix: {t_full:.3f}s,  k-space: {t_kspace:.3f}s")

    all_ok = True
    for i, t in enumerate(t_arr):
        DS_f, Sr_f, Pp_f = res_full[i]
        DS_k = res_k['DeltaS'][i]
        Sr_k = res_k['S_rho'][i]
        Pp_k = res_k['P_pi'][i]
        ok = (abs(DS_f - DS_k) < tol and
              abs(Sr_f - Sr_k) < tol and
              abs(Pp_f - Pp_k) < tol)
        all_ok = all_ok and ok
        if verbose:
            print(f"  t={t:.1f}: DS={DS_k:.5f}/{DS_f:.5f}, "
                  f"S={Sr_k:.5f}/{Sr_f:.5f}, Ppi={Pp_k:.5f}/{Pp_f:.5f}  OK={ok}")

    if verbose:
        print(f"  All passed: {all_ok}")
    return all_ok
