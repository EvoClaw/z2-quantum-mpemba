"""Supplement R5: Quasiparticle mode occupation n_k(t).

For the TFIM, the BdG (Bogoliubov-de Gennes) Hamiltonian in k-space is:
  h_k = [[g - J/2 cos k,  J/2 i sin k],
          [-J/2 i sin k, -(g - J/2 cos k)]]

The initial product state |psi_theta> = prod_j [cos(theta/2)|0> + sin(theta/2)|1>]
has BdG one-particle density matrix M_k where:
  C_k = <c+_k c_k>  (normal)
  F_k = <c_k c_{-k}>  (anomalous)

For the product state in k-space (each k-mode):
  C_k = sin^2(theta/2) (same for all k from Fourier transform of diagonal C_{jj})
  F_k = 0               (no anomalous correlations initially)

After quench to H(g_f), we evolve C_k and F_k using BdG equations.

The QUASIPARTICLE OCCUPATION (for g_f Bogoliubov modes gamma_k) is:
  <gamma+_k gamma_k>(t)  (this tracks how excited the final-state quasiparticles are)

We compute this to reveal: at t_M, are any k-modes special (occupation = 1/2 etc)?

NOTE on one-point functions:
  <c_k> != 0 for the product state (each k-mode is a superposition).
  This shows up in the BdG as 1-point functions, making the problem non-Gaussian.
  We CANNOT use the correlation matrix shortcut for EA, but we CAN use it
  to track the 2-point density matrix n_k(t) which is physical and tractable.

We separate the analysis:
1. Compute the BdG density matrix (2-point correlators only)
2. Show n_k(t) for different modes
3. Correlate mode occupations with t_M
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os
from scipy.linalg import eigh
from exp_common import build_TFIM_sparse, run_ea, find_crossings

# BdG k-modes for large N (correlation matrix level, even though EA needs full ED)
N_bdg = 200   # large N for BdG (k-space, no Hilbert space issue)
J, g_i_bdg = 1.0, 0.0   # we don't need a pre-quench; product state sets initial C_k
t_max, n_t = 15.0, 150
t_arr = np.linspace(0, t_max, n_t)

# k-modes (for PBC, use Ramond sector k = 2pi n / N)
k_vals = 2*np.pi * np.arange(N_bdg) / N_bdg
# k in (-pi, pi]: shift to symmetric representation
k_sym = k_vals.copy()
k_sym[k_sym > np.pi] -= 2*np.pi

def bdg_hamiltonian(k, g, J=1.0):
    """2x2 BdG Hamiltonian at momentum k."""
    eps = g - J/2 * np.cos(k)
    delta = J/2 * np.sin(k)
    return np.array([[eps, 1j*delta], [-1j*delta, -eps]])

def bdg_eigenvalues(k, g, J=1.0):
    """Eigenvalues ±E_k and eigenvectors of h_k."""
    h = bdg_hamiltonian(k, g, J)
    vals, vecs = eigh(h)
    return vals, vecs  # vals[0]<0, vals[1]>0

def time_evolve_bdg(rho_k_init, k, g_f, t, J=1.0):
    """
    Evolve the 2x2 BdG density matrix rho_k(0) for time t under H_BdG(k, g_f).
    rho_k(t) = U(t) rho_k(0) U(t)†  with U(t) = exp(-i h_k t)
    """
    h = bdg_hamiltonian(k, g_f, J)
    vals, vecs = np.linalg.eigh(h)
    # U(t) = V diag(exp(-i vals t)) V†
    phases = np.exp(-1j * vals * t)
    U = vecs @ np.diag(phases) @ vecs.conj().T
    return U @ rho_k_init @ U.conj().T

def quasiparticle_occupation(rho_k, k, g_f, J=1.0):
    """
    Project BdG density matrix onto quasiparticle basis of H(g_f).
    Returns occupation <gamma+_k gamma_k> = (rho_k)_{11} in Bogoliubov basis.
    """
    _, vecs_f = bdg_eigenvalues(k, g_f, J)
    # Transform rho_k to quasiparticle basis
    rho_qp = vecs_f.conj().T @ rho_k @ vecs_f
    return rho_qp[1, 1].real   # occupation of upper quasiparticle mode

# Initial BdG density matrix for product state
# For site product state: <c+_j c_k> = sin^2(theta/2) delta_{jk}, <c_j c_k> = 0
# In k-space: <c+_k c_k'> = sin^2(theta/2) delta_{kk'} (Fourier transform)
# But <c_k c_{-k}> involves cross-correlations: let's compute from k-space sums.
# For the product state: <c+_k c_k> = sin^2(theta/2) for ALL k
# <c_k c_{-k}> = sin^2(theta/2) * ???
# Actually: <c_j c_l> = <c_j><c_l> = (sin theta / 2)^2 for j != l
# In k-space: <c_k c_k'> = sin^2(theta/2) for k != k' if we include 1-pt functions
# This is the non-Gaussian part. The pure 2-pt part:
# C_{kk'}^normal = <c+_k c_k'> = sin^2(theta/2) delta_{kk'}
# The BdG density matrix rho_k is:
#   rho_k = [[<c+_k c_k>,  <c_{-k} c_k>],
#             [<c+_k c+_{-k}>, <c_{-k} c+_{-k}>]]
#          = [[sin^2(th/2),  0             ],
#             [0,             cos^2(th/2)  ]]
# (upper-left = occupation, lower-right = 1-occupation for the -k mode)
# This ignores the 1-point contribution from <c_k> != 0.
# For a QUALITATIVE analysis of quasiparticle content, this is still informative.

g_f_list = [0.8, 1.0, 2.0]
theta_list = [0.4, 0.6, 0.8, 1.0, 1.2]

# For small-N exact EA to get t_M reference
N_ed, NA_ed = 12, 4
theta_pairs_qme = [(0.4, 1.2), (0.8, 1.4)]
g_f_ed = 2.0
H_ed = build_TFIM_sparse(N_ed, g=g_f_ed)

print(f'Supplement R5: Quasiparticle mode occupation analysis')
print(f'BdG large-N: N_bdg={N_bdg}  (2-pt BdG density matrix)')
print(f'Exact EA reference: N={N_ed}, g_f={g_f_ed}\n', flush=True)

# Get t_M for the two theta pairs
t_m_ref = {}
for th1, th2 in theta_pairs_qme:
    t_arr_ed = np.linspace(0, 15, 150)
    DS1,_,_,_ = run_ea(N_ed, NA_ed, th1, H_ed, t_arr_ed)
    DS2,_,_,_ = run_ea(N_ed, NA_ed, th2, H_ed, t_arr_ed)
    if DS1[0] > DS2[0]:
        DS1, DS2 = DS2, DS1
    cross = find_crossings(t_arr_ed, DS1, DS2)
    t_m_ref[f'{th1}_{th2}'] = cross[0] if cross else None
    tm_str = f'{cross[0]:.4f}' if cross else 'N/A'
    print(f'θ={th1}/{th2}: ΔS(0)={DS1[0]:.4f}/{DS2[0]:.4f}, t_M={tm_str}', flush=True)

# Compute n_k(t) for g_f = 2.0 and different theta values
print('\nComputing n_k(t) for g_f=2.0...', flush=True)
g_f_nk = 2.0
results = {'params': {'N_bdg':N_bdg,'g_f':g_f_nk,'t_max':t_max,'n_t':n_t},
           't': t_arr.tolist(),
           'k_sym': k_sym.tolist(),
           't_M_ref': {k: float(v) if v else None for k,v in t_m_ref.items()}}

for theta in theta_list:
    sin2 = np.sin(theta/2)**2
    cos2 = np.cos(theta/2)**2
    # BdG density matrix for this theta
    rho_k0 = np.array([[sin2, 0.], [0., cos2]], dtype=complex)

    # Track n_k(t) for a few representative k-modes
    k_indices = [0, N_bdg//8, N_bdg//4, 3*N_bdg//8, N_bdg//2]
    nk_t = {ki: [] for ki in k_indices}

    for t in t_arr:
        for ki in k_indices:
            k = k_vals[ki]
            rho_t = time_evolve_bdg(rho_k0, k, g_f_nk, t)
            nk_t[ki].append(quasiparticle_occupation(rho_t, k, g_f_nk))

    results[f'theta{theta:.2f}'] = {
        'theta': float(theta), 'DS0_approx': float(-theta*np.log(np.tan(theta/2)**2 + 1)
                                                    if 0 < theta < np.pi/2 else 0),
        'n_k_t': {str(k_vals[ki])[:8]: np.array(nk_t[ki]).tolist() for ki in k_indices}
    }
    print(f'  θ={theta:.2f}: sin²(θ/2)={sin2:.4f}  ' +
          f'  n_k0(t=0)={quasiparticle_occupation(rho_k0, k_vals[0], g_f_nk):.4f}', flush=True)

# Show n_k(t) at t_M
print('\n===== n_k at t=t_M (theta=0.4, 1.2) =====')
for theta in [0.4, 1.2]:
    sin2 = np.sin(theta/2)**2
    cos2 = np.cos(theta/2)**2
    rho_k0 = np.array([[sin2, 0.], [0., cos2]], dtype=complex)
    t_M = t_m_ref.get('0.4_1.2')
    if t_M:
        print(f'  θ={theta:.2f}  (t_M={t_M:.4f}):')
        for ki in [0, N_bdg//8, N_bdg//4, N_bdg//2]:
            k = k_vals[ki]
            rho_tM = time_evolve_bdg(rho_k0, k, g_f_nk, t_M)
            nk = quasiparticle_occupation(rho_tM, k, g_f_nk)
            print(f'    k/pi={k/np.pi:.4f}: n_k(t_M)={nk:.4f}')

# Check k-mode at which n_k(t_M) ≈ 1/2 (maximum uncertainty)
print('\n===== k-modes with n_k(t_M) ≈ 0.5 (for theta=1.0) =====')
theta = 1.0
sin2 = np.sin(theta/2)**2
cos2 = np.cos(theta/2)**2
rho_k0 = np.array([[sin2, 0.], [0., cos2]], dtype=complex)
t_M = t_m_ref.get('0.8_1.4')
if t_M:
    nk_all = []
    for k in k_vals:
        rho_tM = time_evolve_bdg(rho_k0, k, g_f_nk, t_M)
        nk_all.append(quasiparticle_occupation(rho_tM, k, g_f_nk))
    nk_all = np.array(nk_all)
    near_half = np.where(np.abs(nk_all - 0.5) < 0.05)[0]
    print(f'  k-modes with |n_k - 0.5| < 0.05: {len(near_half)} out of {N_bdg}')
    if len(near_half) < 10:
        for ki in near_half[:5]:
            print(f'    k/pi={k_vals[ki]/np.pi:.4f}: n_k={nk_all[ki]:.4f}')

os.makedirs('results', exist_ok=True)
json.dump(results, open('results/supp_r5_nk.json','w'), indent=2)
print('\nR5 done → results/supp_r5_nk.json')
