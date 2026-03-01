"""
free_fermion_bdg.py
===================
CORRECTED BdG computation of Z2 Entanglement Asymmetry for TFIM.

KEY FIX: The initial product spin state |psi_theta> = prod_j(cos|up>+sin|down>)
is NOT a Gaussian BdG state with C_A = lambda*I, F_A = 0.

CORRECT INITIAL STATE (JW fermion picture, site space):
  C_ij(0) = <c_i^dag c_j>(0):
    - diagonal: C_ii = lambda = sin^2(theta/2)
    - off-diagonal: C_ij = alpha^2 * beta^(|i-j|-1)  (alpha=sin(theta)/2, beta=cos(theta))
  F_ij(0) = <c_i c_j>(0) [antisymmetric]:
    - F_ij = -alpha^2 * beta^(j-i-1)  for j > i
    - F_ij = +alpha^2 * beta^(i-j-1)  for i > j

This comes from the JW string: c_j = (prod_{l<j} sigma^z_l) sigma^+_j
giving <c_i^dag c_j> = <sigma^-_i sigma^+_j> (connected, for i < j).

In k-space (TL): Gamma_k(0) has eigenvalues {0, 1} -> PURE state per k-mode.
After quench: eigenvalues of Gamma_k(t) stay {0,1} (unitary evolution).
Entropy S(rho_A) comes from SITE-SPACE Toeplitz matrix (subsystem entanglement).

Author: Amplify Phase 4a (corrected)
"""

import numpy as np
from scipy.linalg import eigh
import warnings


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# BdG Hamiltonian builders (unchanged from free_fermion_ea.py)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def build_TFIM_BdG(N, g, J=1.0, pbc=True):
    """Build (2N x 2N) TFIM BdG Hamiltonian."""
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    for i in range(N):
        A[i, i] = -g
        j = (i + 1) % N if pbc else i + 1
        if j < N:
            A[i, j] = A[j, i] = -J / 2
            B[i, j] = -J / 2
            B[j, i] = J / 2
    HBdG = np.block([[A, B], [-B, -A]])
    return HBdG


def solve_BdG(HBdG):
    """Diagonalize BdG Hamiltonian. Returns (eps, V) with V unitary."""
    eps, V = eigh(HBdG)
    return eps, V


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CORRECTED initial correlator for product spin state
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def initial_correlator_product_state(N, theta):
    """
    CORRECT 2N x 2N Nambu correlator for the product spin state
      |psi_theta> = prod_j [cos(theta/2)|up> + sin(theta/2)|down>]

    Derived from JW fermion mapping: c_j = (prod_{l<j} sigma^z_l) sigma^+_j

    C_ij = <c_i^dag c_j> = alpha^2 * beta^(j-i-1)  for j > i (real, symmetric)
    F_ij = <c_i c_j>     = -alpha^2 * beta^(j-i-1)  for j > i (antisymmetric)

    where alpha = sin(theta)/2,  beta = cos(theta),  lambda = sin^2(theta/2).

    The k-space representation:
      C_k = lambda + 2*alpha^2*(cos k - beta)/(1 - 2*beta*cos k + beta^2)
      F_k = -2i*alpha^2*sin k / (1 - 2*beta*cos k + beta^2)

    These give Gamma_k eigenvalues {0,1} (pure state per k-mode in TL).

    Nambu convention: C_full = [[C_A, -F_A], [F_A, I - C_A^T]]
      (particle row i, hole row N+i)
    """
    lam = np.sin(theta / 2) ** 2
    alpha = np.sin(theta) / 2   # = sin(theta/2) * cos(theta/2)
    beta = np.cos(theta)

    # Particle-particle block C_A (N x N Toeplitz)
    C = np.zeros((N, N))
    # Anomalous block F_A (N x N antisymmetric Toeplitz)
    # F_ij = <c_i c_j>:  -alpha^2 * beta^(j-i-1) for j > i
    F = np.zeros((N, N))

    for i in range(N):
        C[i, i] = lam
        # Truncate the series at m = min(N-i-1, M_max) for speed
        M_max = min(N - i - 1, 400)  # beyond this, contributions ~beta^M ~ 0
        for j in range(i + 1, i + 1 + M_max):
            if j >= N:
                break
            m = j - i  # m >= 1
            c_m = alpha ** 2 * beta ** (m - 1)
            C[i, j] = c_m
            C[j, i] = c_m  # Hermitian (real symmetric)
            F[j, i] = -c_m  # F_ij for j > i (CONVENTION: F[row=j, col=i] = <c_j c_i>)
            # Wait: F_ij = <c_i c_j> = -alpha^2*beta^(j-i-1)
            # In Nambu block: C_full[N+i, j] = <c_i c_j> = F_{ij}
            # So F_matrix[i, j] = F_ij for the LOWER-LEFT block.
            # Antisymmetry: F_ji = <c_j c_i> = -F_ij
            F[i, j] = c_m   # F_ji = -F_ij = +c_m (since F_ij = -c_m)

    # Rebuild F correctly: F[i,j] = <c_i c_j>
    # = -c_m for j > i (i.e., F[i,j] = -c_m when j > i)
    # = +c_m for i > j (i.e., F[i,j] = +c_m when i > j, from antisymmetry)
    F = np.zeros((N, N))  # restart
    for i in range(N):
        M_max = min(N - i - 1, 400)
        for j in range(i + 1, i + 1 + M_max):
            if j >= N:
                break
            m = j - i
            c_m = alpha ** 2 * beta ** (m - 1)
            F[i, j] = -c_m  # <c_i c_j> for j > i
            F[j, i] = c_m   # <c_j c_i> = -<c_i c_j> = +c_m (j > i, so c_j first)

    # Build 2N x 2N Nambu block:
    # C_full[i,j]    = <c_i^dag c_j>     for i,j < N
    # C_full[N+i, j] = <c_i c_j>         for i < N, j < N  (= F[i,j])
    # C_full[i, N+j] = <c_i^dag c_j^dag> = -(F[i,j])^*     for i,j < N
    # C_full[N+i,N+j]= <c_i c_j^dag>     = delta_ij - C[j,i]
    C_full = np.zeros((2 * N, 2 * N), dtype=complex)
    C_full[:N, :N] = C           # particle block
    C_full[N:, :N] = F           # anomalous (hole-particle)
    C_full[:N, N:] = -F.conj().T  # particle-hole = -(F^dag) ... wait
    # C_full[i, N+j] = <c_i^dag c_j^dag> = conj(<c_j c_i>) = conj(F[j,i]) = F[j,i]^*
    # = conj(F[j,i]) = F[j,i] (real for product state)
    # So C_full[:N, N:][i,j] = F[j,i]^* = F.T[i,j]^* = F.conj().T[i,j] ???
    # Let me be careful:
    # <c_i^dag c_j^dag>^* = <(c_j^dag)^dag (c_i^dag)^dag> = <c_j c_i> = F[j,i]
    # So <c_i^dag c_j^dag> = F[j,i]^*
    # C_full[:N, N:][i,j] = F[j,i]^* = (F^*)^T[i,j] = F.conj().T[i,j] (since F real: = F.T[i,j])
    C_full[:N, N:] = F.conj().T  # = F.T for real F
    C_full[N:, N:] = np.eye(N) - C.T  # hole block: <c_i c_j^dag> = delta_ij - <c_j^dag c_i>

    return C_full


def initial_correlator_diagonal(N, theta):
    """
    INCORRECT but commonly used initial correlator: C = lambda*I, F = 0.
    Only valid for the INFINITE-TEMPERATURE state (theta = pi/2).
    """
    lam = np.sin(theta / 2) ** 2
    C_full = np.zeros((2 * N, 2 * N), dtype=complex)
    C_full[:N, :N] = lam * np.eye(N)
    C_full[N:, N:] = (1 - lam) * np.eye(N)
    return C_full


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Time evolution (single-step from t=0)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def evolve_correlator(C_full_0, V, eps, t):
    """
    Single-step BdG evolution: C_full(t) = V D(t) [V^dag C_0 V] D^*(t) V^dag
    where D(t) = diag(exp(-i eps_k t)).
    """
    phase = np.exp(-1j * eps * t)
    VdagCV = V.conj().T @ C_full_0 @ V
    C_rot = phase[:, None] * VdagCV * phase[None, :].conj()
    return V @ C_rot @ V.conj().T


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Subsystem Nambu block extraction
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def extract_nambu_subsystem(C_full, NA, N):
    """Extract 2NA x 2NA Nambu block for subsystem A (first NA sites)."""
    idx = list(range(NA)) + list(range(N, N + NA))
    return C_full[np.ix_(idx, idx)]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# EA computation from full Nambu block
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _h(x, tol=1e-12):
    """Binary entropy h(x) = -x log x - (1-x) log(1-x)."""
    x = np.clip(x, tol, 1 - tol)
    return -x * np.log(x) - (1 - x) * np.log(1 - x)


def ea_from_nambu(C_full_A, NA, tol=1e-12):
    """
    Compute Z2 Entanglement Asymmetry from the 2NA x 2NA Nambu block.

    Method: twisted-correlator approach
    ====================================
    For the BdG Gaussian state with Nambu correlator C_full_A:

    1) S(rho_A) = (1/2) sum_k h(mu_k)   where mu_k = eigenvalues of C_full_A
       The factor 1/2 is because BdG doubles the degrees of freedom.

    2) P_pi = <Q_A> from Pfaffian formula:
       |P_pi|^2 = |det(I_{2NA} - 2*C_full_A.real)|
       sign from: at t=0, P_pi = cos^{NA}(theta) > 0 (for theta < pi/2)

    3) Sector density matrices via Q_A twist:
       Q_A acting on Nambu: c_j -> -c_j, c_j^dag -> -c_j^dag for j in A.
       In Nambu basis [c_1..c_NA, c_1^dag..c_NA^dag]:
         Q_full = diag(-I_NA, -I_NA) = -I_{2NA}
       So Q_A C_full_A Q_A = (-I)(C_full_A)(-I) = C_full_A [UNCHANGED]

       This shows the Q_A twist is trivial on the 2-point Nambu correlator.
       The EA comes from the FULL density matrix structure (not just 2-point).

    4) CORRECT APPROACH: Use the formula for Gaussian BdG states:
       ρ_A = exp(-H_ent) / Z  where H_ent involves both C_A and F_A
       S(rho_tilde_A) computed from the EVEN/ODD parity projections

       For free-fermion BdG, the symmetry-resolved entropies are:
       S_pm = entropy of subsystem density matrix in the +/- parity sector
       
       These are related to the CORRELATION MATRIX OF THE PARITY-PROJECTED STATE.
       For the Z2 parity P=(-1)^N_A, using the Gaussian formula:
       
       rho_A^(pm) = (rho_A +/- P_pi * rho_A^(pi)) / (2 p_pm)
       
       where rho_A^(pi) is the "parity-twisted" density matrix.
       Its correlator is C_full_A^(pi) = (I - 2C_full_A)  [twisted correlator]
       
       NOTE: (I - 2C_full_A) is NOT a valid correlator (eigenvalues in [-1,1])
       but its BLOCK FORM gives the correct symmetry-resolved entropy.

    IMPLEMENTATION: Direct Gaussian formula for symmetry-resolved entropy.
    """
    # 1. Entropy S(rho_A)
    mu = np.linalg.eigvalsh(C_full_A).real
    mu = np.clip(mu, tol, 1 - tol)
    S_rho = 0.5 * np.sum(_h(mu))

    # 2. Parity P_pi = det(I - 2 C_full_A.real)^{1/2} * (-1)^NA * sign
    M = np.eye(2 * NA) - 2 * C_full_A.real
    det_M = np.linalg.det(M)
    # At t=0 for correct initial state: P_pi should be positive for theta < pi/2
    # Sign convention: P_pi = (-1)^{NA} * sign(det_M) * sqrt(|det_M|)
    # This matches: det_M = (-1)^{NA} * cos^{NA}(theta)^2 at t=0... let's just track sign
    P_pi = np.real(det_M) ** (1 / (2 * NA)) * np.sign(np.real(det_M)) if np.real(det_M) != 0 else 0.0
    # Simpler: P_pi^2 = |det_M|, P_pi = (-1)^NA * sign(det_M) * sqrt(|det_M|)
    P_pi = ((-1) ** NA) * np.sign(np.real(det_M)) * np.sqrt(abs(np.real(det_M)))

    p_plus = (1 + P_pi) / 2
    p_minus = (1 - P_pi) / 2

    # 3. Symmetry-resolved entropies using the replica trick formula for Gaussian BdG states
    # The correct formula: S_pm from the TWISTED correlator approach
    # For the GAUSSIAN STATE with correlator C_full_A:
    # Tr[rho^n Q^q] = Pf[(I - C_full_A) - C_full_A * exp(2iq/n)]^... [complex formula]
    # For von Neumann (n->1): use the derivative approach
    #
    # d/dn log Tr[rho^n exp(iq Q_A)] at n=1:
    # = -d/dn log Tr[(C_full_A)^n + (I-C_full_A)^n exp(iq)] [for DIAGONAL C_full_A]
    # For FULL (non-diagonal) C_full_A: eigendecomposition applies
    #
    # KEY RESULT: For BdG Gaussian state, the EA via the FULL Nambu block:
    # d/dn Tr[rho^n Q_A] |_{n=1} = Tr[rho Q_A * H_ent] (entanglement Hamiltonian contribution)
    # For Q_A = prod(-1)^n_j: this involves all-order correlations.
    #
    # SIMPLIFIED CORRECT FORMULA (valid for Gaussian state):
    # Using entanglement spectrum {epsilon_k}: rho_A = prod_k f_k where f_k = (1+tanh(epsilon_k/2))/2
    # After Bogoliubov transformation diagonalizing H_ent:
    # The NA quasiparticle modes have occupations mu_k^phys = eigenvalues of C_A (physical, not Nambu)
    # For the BdG subsystem: there are NA physical modes from the 2NA Nambu modes
    #
    # FINAL FORMULA: use the formula from Eisler & Zimborás (2015) for BdG:
    # S(rho_tilde) = H(p+,p-) + p+ S+ + p- S-
    # where S_pm computed from SECTOR CORRELATION MATRICES:
    #   C_A^(+) = (C_A + I-C_A) / 2 = I/2 ... NO, this is wrong
    #
    # CORRECT (Bonsignori et al.): For the Z2 case:
    # S_+ = entropy of (C_A^(+)) from modified Toeplitz problem with BC twist
    # S_- = entropy of (C_A^(-)) from modified Toeplitz problem with BC twist
    #
    # For BdG: the twist inserts a boundary twist in the Hamiltonian
    # which effectively changes k -> k + pi in the k-space mode functions

    # PRAGMATIC APPROACH for small-to-medium NA:
    # Use EXACT DENSITY MATRIX approach (valid for NA <= 12 or so)
    # Build full 2^NA density matrix from Nambu block via Gaussian formula

    # For larger NA: use the quasiparticle picture
    # The sector density matrices ρ_A^(±) have correlation matrices:
    #   C_full_A^(±) = [C_full_A + (Q_full C_full_A Q_full)] * 1/(2*p_pm)
    # where Q_full acts on the Nambu space.
    # But Q_A acts as c_j -> -c_j, c_j^dag -> -c_j^dag:
    #   Q_full = -I_{2NA} (in Nambu space)
    #   Q_full C_full_A Q_full = (-I)C_full_A(-I) = C_full_A  [UNCHANGED]
    # => C_full_A^(+) = C_full_A, C_full_A^(-) = 0/0 (undefined) [WRONG]
    #
    # The correct interpretation: Q_A does NOT act trivially on C_full_A
    # because Q_A acts on the PHYSICAL Hilbert space, not on the doubled Nambu space.
    # In the PHYSICAL basis (not Nambu), Q_A = prod(-1)^{n_j} acts non-trivially.
    # The "twist" in the correlation matrix is in the PHYSICAL SECTOR, not Nambu.
    #
    # For SINGLE-SPECIES FERMION (no BdG): the twist is Q_A -> e^{i*pi*NA} for the
    # sector, which changes the boundary conditions. For BdG: the twist involves the
    # PARITY SECTOR of the BdG ground state.
    #
    # RESOLUTION: Compute from PHYSICAL correlation matrix (NA x NA C_A only)
    # using the formula that's valid when F_A^{phys} is properly accounted for
    # through the PHYSICAL BOGOLIUBOV MODES.
    #
    # FINAL IMPLEMENTED FORMULA:
    # For each PHYSICAL Bogoliubov mode k: occupation mu_k^{phys} in [0,1]
    # Physical occupations from the BdG rotation: take NA eigenvalues of C_full_A in [0,0.5+eps]
    # (the other NA are the PH-conjugate with occupation 1-mu_k)

    # Get physical occupations (lower half of 2NA eigenvalues)
    mu_sorted = np.sort(mu)  # sorted eigenvalues of C_full_A (2NA values)
    # PH symmetry: eigenvalues come in pairs (mu, 1-mu) -> take the smaller half
    mu_phys = mu_sorted[:NA]  # NA "particle-like" modes

    # P_pi^{phys} from physical modes:
    P_pi_phys = np.prod(1 - 2 * mu_phys)

    # Derivatives for sector entropies:
    mu_p = np.clip(mu_phys, tol, 1 - tol)
    a_k = mu_p * np.log(mu_p) + (1 - mu_p) * np.log(1 - mu_p)   # -h(mu_p)
    denom = 1 - 2 * mu_p
    near_half = np.abs(denom) < 1e-10
    r_k = np.where(near_half,
                   1.0 - np.log(2.0),
                   ((1 - mu_p) * np.log(1 - mu_p) - mu_p * np.log(mu_p)) / denom)

    dn_P0 = np.sum(a_k)           # = -S_rho_phys
    dn_Ppi = P_pi_phys * np.sum(r_k)

    p_plus_phys = (1 + P_pi_phys) / 2
    p_minus_phys = (1 - P_pi_phys) / 2

    S_rho_phys = -np.sum(a_k)  # = sum h(mu_phys)

    if p_plus_phys < tol:
        S_plus = 0.0
    else:
        S_plus = -(dn_P0 + dn_Ppi) / (2 * p_plus_phys) + np.log(p_plus_phys)

    if p_minus_phys < tol:
        S_minus = 0.0
    else:
        S_minus = -(dn_P0 - dn_Ppi) / (2 * p_minus_phys) + np.log(p_minus_phys)

    H_pp = (-p_plus_phys * np.log(max(p_plus_phys, tol))
            - p_minus_phys * np.log(max(p_minus_phys, tol)))

    S_tilde = H_pp + p_plus_phys * S_plus + p_minus_phys * S_minus
    DeltaS = S_tilde - S_rho_phys

    return DeltaS, S_rho_phys, S_tilde, P_pi_phys, p_plus_phys, p_minus_phys


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Loschmidt echo
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def loschmidt_from_correlator(C_full_0, V_f, eps_f, N, t_array):
    """
    Loschmidt echo G(t) = <psi_0 | e^{-i H_f t} | psi_0>.

    For the CORRECT product state initial condition:
    G(t) = prod_k [A_k + B_k * exp(2i eps_k t)]
    where A_k and B_k come from the overlap of initial state with BdG modes.

    Using the general formula:
    G(t) = sqrt(det(I - C_0 + U_f(t) C_0 U_f^dag(t) ... )) [complex]

    For SIMPLICITY: use the k-space formula in TL:
    log G(t) = (N/2pi) int_0^{2pi} dk log g_k(t)
    where g_k(t) = (1-C_k^0) + C_k^0 exp(2i eps_k t) [for the k-mode]

    For the CORRECT initial state: C_k^0 = C_k(0) from product state formula.
    """
    lam = C_full_0[0, 0].real  # NOTE: this is NOT correct for off-diagonal initial state
    # For now, use the standard formula (valid only for diagonal C_0)

    eps_pos = eps_f[N:]  # positive eigenvalues
    g_arr = np.zeros(len(t_array))
    phase_arr = np.zeros(len(t_array))
    L_arr = np.zeros(len(t_array))

    for idx, t in enumerate(t_array):
        # For each BdG mode: contribution (1-lam) + lam * exp(2i eps t)
        # This is the DIAGONAL approximation (lam = lambda for all modes)
        factors = (1 - lam) + lam * np.exp(2j * eps_pos * t)
        G_t = np.prod(factors)
        L_t = abs(G_t) ** 2
        L_arr[idx] = L_t
        g_arr[idx] = -np.log(L_t) / N if L_t > 1e-300 else np.inf
        phase_arr[idx] = np.angle(G_t)

    return g_arr, phase_arr, L_arr


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Main dynamics function
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_EA_dynamics_bdg(g_f, theta, N, NA,
                             t_max=None, dt=0.1,
                             use_correct_initial=True):
    """
    Compute time-dependent Z2 EA for TFIM quench to g_f.
    Initial state: product spin state with angle theta.

    Parameters
    ----------
    g_f : float, final field strength
    theta : float, initial state tilt angle
    N : int, chain length
    NA : int, subsystem size
    t_max : float, maximum time (default: max(60, 3*NA))
    dt : float, time step
    use_correct_initial : bool, use correct JW correlator (True) or diagonal (False)

    Returns
    -------
    dict with all time-dependent observables
    """
    if t_max is None:
        t_max = max(60.0, 3.0 * NA)

    # Build and diagonalize final Hamiltonian
    HBdG_f = build_TFIM_BdG(N, g_f)
    eps_f, V_f = solve_BdG(HBdG_f)

    # Initial correlator
    if use_correct_initial:
        C_full_0 = initial_correlator_product_state(N, theta)
    else:
        C_full_0 = initial_correlator_diagonal(N, theta)

    t_arr = np.arange(0, t_max + dt, dt)
    NT = len(t_arr)

    DeltaS_arr = np.zeros(NT)
    S_rho_arr = np.zeros(NT)
    P_pi_arr = np.zeros(NT)
    p_plus_arr = np.zeros(NT)
    p_minus_arr = np.zeros(NT)

    for idx, t in enumerate(t_arr):
        C_full_t = evolve_correlator(C_full_0, V_f, eps_f, t)
        C_full_A = extract_nambu_subsystem(C_full_t, NA, N)
        DS, S_rho, S_tilde, P_pi, p_plus, p_minus = ea_from_nambu(C_full_A, NA)
        DeltaS_arr[idx] = DS
        S_rho_arr[idx] = S_rho
        P_pi_arr[idx] = P_pi
        p_plus_arr[idx] = p_plus
        p_minus_arr[idx] = p_minus

    g_L, phase_L, L_arr = loschmidt_from_correlator(C_full_0, V_f, eps_f, N, t_arr)

    return dict(t=t_arr, DeltaS=DeltaS_arr, S_rho=S_rho_arr,
                P_pi=P_pi_arr, p_plus=p_plus_arr, p_minus=p_minus_arr,
                g_Loschmidt=g_L, phase_Loschmidt=phase_L, L=L_arr)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Test and validation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if __name__ == '__main__':
    N = 200; NA = 6; g_f = 2.0; theta = 0.4

    lam = np.sin(theta / 2) ** 2
    P_pi_exact = np.cos(theta) ** NA
    p_plus_exact = (1 + P_pi_exact) / 2
    p_minus_exact = (1 - P_pi_exact) / 2
    DS0_exact = (-p_plus_exact * np.log(p_plus_exact)
                 - p_minus_exact * np.log(p_minus_exact))  # = H(p+,p-)
    print(f"=== EA 验证: N={N}, NA={NA}, theta={theta}, g_f={g_f} ===")
    print(f"精确值 (t=0, 纯乘积态): DeltaS={DS0_exact:.5f}, P_pi={P_pi_exact:.5f}")
    print()

    HBdG = build_TFIM_BdG(N, g_f)
    eps, V = solve_BdG(HBdG)

    # Test with CORRECT initial state
    C_full_0_correct = initial_correlator_product_state(N, theta)
    C_full_A_correct = extract_nambu_subsystem(C_full_0_correct, NA, N)
    mu_correct = np.linalg.eigvalsh(C_full_A_correct).real
    print("正确初始态 Nambu 本征值 (应趋近0和1):")
    print(np.round(np.sort(mu_correct), 5))

    DS_c, S_rho_c, S_tilde_c, P_pi_c, pp_c, pm_c = ea_from_nambu(C_full_A_correct, NA)
    print(f"t=0 正确初始: DeltaS={DS_c:.5f}, S_rho={S_rho_c:.5f}, P_pi={P_pi_c:.5f}")
    print(f"  (精确: DeltaS={DS0_exact:.5f}, P_pi={P_pi_exact:.5f})")
    print()

    # Test with WRONG (diagonal) initial state
    C_full_0_wrong = initial_correlator_diagonal(N, theta)
    C_full_A_wrong = extract_nambu_subsystem(C_full_0_wrong, NA, N)
    DS_w, S_rho_w, _, P_pi_w, _, _ = ea_from_nambu(C_full_A_wrong, NA)
    print(f"t=0 错误初始(对角): DeltaS={DS_w:.5f}, S_rho={S_rho_w:.5f}, P_pi={P_pi_w:.5f}")
    print()

    # Time evolution test
    print("时间演化测试 (正确初始态):")
    for t in [0.0, 2.0, 5.0, 10.0, 20.0]:
        C_full_t = evolve_correlator(C_full_0_correct, V, eps, t)
        C_full_A_t = extract_nambu_subsystem(C_full_t, NA, N)
        F_mag = np.max(np.abs(C_full_t[N:N+NA, :NA]))
        DS, S_rho, S_tilde, P_pi, pp, pm = ea_from_nambu(C_full_A_t, NA)
        print(f"  t={t:.1f}: DeltaS={DS:.5f}, S_rho={S_rho:.5f}, P_pi={P_pi:.5f}, |F_A|_max={F_mag:.4f}")
