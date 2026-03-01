"""Supplement R4: Anisotropic XY chain control experiment.

H_XY(gamma) = J/2 Sum_j [ (1+gamma)(c+_j c_{j+1} + h.c.)
                           + (1-gamma)(c+_j c+_{j+1} + h.c.) ] + g Sum_j(2n_j-1)

gamma=1: TFIM  (g_c = J/2 = 0.5)
gamma=0: XX chain (g_c = 0, always in PM for any g>0)

Key point: for gamma in (0,1), QPT critical line g_c(gamma) changes.
DQPT persists for g_f > g_c(gamma).  QME occurrence is tested against both.

For the k-space BdG the dispersion is:
  eps_k = sqrt( (g - J(1+gamma)/2 cos k)^2 + (J(1-gamma)/2 sin k)^2 )
At k=0: eps_0 = |g - J(1+gamma)/2|  => g_c = J(1+gamma)/2 ??? 
No: for XY with our convention J_hop = J(1+gamma), J_pair = J(1-gamma):
  eps_0 = |g - J_hop/2| = |g - J(1+gamma)/2|
BUT the true QPT (gap closing) occurs at the k that minimizes eps_k,
which is k=0 for our convention -> g_c = J(1+gamma)/2.

Wait, that's only correct for gamma<=1 (J_pair<=J_hop).
Check: for gamma=1, g_c = J = 1? But we know g_c=0.5 for TFIM!

I think the issue is the factor of 2 in our TFIM convention.
Our TFIM: H = J/2 Sum [(c+c + c+c+ + h.c.)] + g Sum(2n-1)
The BdG: eps_k = sqrt((g - J/2*cos k)^2 + (J/2*sin k)^2)
At k=0: |g - J/2|, so g_c = J/2 = 0.5 for J=1 ✓

XY generalization with our convention:
H_XY = J_hop/2 Sum (c+c + h.c.) + J_pair/2 Sum (c+c+ + h.c.) + g Sum(2n-1)

TFIM: J_hop = J_pair = J = 1  -> g_c = J/2 = 0.5
XX:   J_hop = J, J_pair = 0   -> g_c = J/2 = 0.5 still??

Hmm, for XX chain (J_pair=0):
eps_k = sqrt((g - J_hop/2 * cos k)^2 + 0) = |g - J_hop/2 * cos k|

Gap closes at k=0 when g = J_hop/2. So g_c = J_hop/2 = 0.5 for J_hop=J=1.
That's the same as TFIM... but XX chain is gapless at g=0 (Luttinger liquid)!

Issue: the XX chain with J_pair=0 is just a hopping chain, which is gapless
for -J_hop/2 < g < J_hop/2 (band filling). So g_c = J_hop/2 as a "PM" onset.

Actually for DQPT in the XY chain: DQPT occurs when g_f > g_c = J(1+gamma)/2 / 2?
Let me just implement it and let the numerics speak.

We use the STANDARD parameterization where:
J_hop = J * (1 + gamma)
J_pair = J * (1 - gamma)
(so gamma=1 -> J_pair/J_hop = 0, gamma=0 -> J_pair/J_hop = 1, which is opposite to TFIM!)

Let me use a different parameterization to match the literature:
gamma=0: isotropic XY (equal hopping and pairing proportions = our TFIM convention!)
gamma=1: XX model (pairing=0)
So define:
J_hop = J
J_pair = J * (1 - gamma)

gamma=0: TFIM (J_hop=J_pair=J) ✓
gamma=1: XX (J_hop=J, J_pair=0)
g_c = sqrt(J_hop * J_pair)/2 = J*sqrt(1-gamma)/2

For gamma=0: g_c = J/2 = 0.5 ✓
For gamma=0.5: g_c = J*sqrt(0.5)/2 ≈ 0.354
For gamma=1: g_c = 0 (XX model, gapless at g=0)
"""
import sys; sys.path.insert(0,'src')
import numpy as np, json, os, time
import scipy.sparse as sp
from exp_common import run_ea, find_crossings, find_dqpt_times

def build_XY_sparse(N, J=1.0, gamma=0.0, g=2.0, pbc=True):
    """
    H_XY = J/2 Sum_j [(c+_j c_{j+1} + h.c.)   <- hopping (always)
                      + (1-gamma)(c+_j c+_{j+1} + h.c.)]  <- pairing reduced
           + g Sum_j (2n_j - 1)
    gamma=0: TFIM (equal hop and pair)
    gamma=1: XX chain (no pairing)
    """
    J_hop  = J
    J_pair = J * (1.0 - gamma)
    dim = 2**N
    data_diag, data_hop, data_pair = [], [], []
    row_d, col_d, row_h, col_h, row_p, col_p = [], [], [], [], [], []

    # diagonal: g(2n_j - 1)
    for state in range(dim):
        diag_val = 0.0
        for j in range(N):
            n_j = (state >> j) & 1
            diag_val += g * (2*n_j - 1)
        data_diag.append(diag_val)
        row_d.append(state); col_d.append(state)

    # hopping and pairing
    for state in range(dim):
        for j in range(N):
            jp1 = (j+1) % N
            # Boundary sign for PBC (Jordan-Wigner)
            pbc_sign = 1
            if pbc and j == N-1:
                n_total = bin(state).count('1')
                pbc_sign = -(-1)**n_total

            n_j   = (state >> j)   & 1
            n_jp1 = (state >> jp1) & 1

            # hopping: -J_hop/2 (c+_j c_{j+1} + h.c.)
            # c+_j c_{j+1} |..0_j..1_{j+1}..> = sgn |..1_j..0_{j+1}..>
            if n_j == 0 and n_jp1 == 1:
                # c+_j c_{j+1}: create at j, destroy at j+1
                new_state = state - (1 << jp1) + (1 << j)
                # Fermionic sign: count flips between j and j+1
                mask = ((1 << jp1) - 1) - ((1 << (j+1)) - 1)
                sign = (-1)**bin(state & mask).count('1')
                val = -J_hop/2.0 * pbc_sign * sign
                data_hop.append(val);  row_h.append(new_state); col_h.append(state)
                data_hop.append(val.conjugate()); row_h.append(state); col_h.append(new_state)

            # pairing: -J_pair/2 (c+_j c+_{j+1} + h.c.)
            if J_pair != 0:
                if n_j == 0 and n_jp1 == 0:
                    # c+_j c+_{j+1}: create at both j and j+1
                    new_state = state + (1 << j) + (1 << jp1)
                    mask = ((1 << jp1) - 1) - ((1 << (j+1)) - 1)
                    sign = (-1)**bin(state & mask).count('1')
                    val = -J_pair/2.0 * pbc_sign * sign
                    data_pair.append(val);  row_p.append(new_state); col_p.append(state)
                    data_pair.append(val.conjugate()); row_p.append(state); col_p.append(new_state)

    H  = sp.csr_matrix((data_diag, (row_d, col_d)), shape=(dim,dim), dtype=complex)
    Hh = sp.csr_matrix((data_hop,  (row_h, col_h)), shape=(dim,dim), dtype=complex)
    Hp = sp.csr_matrix((data_pair, (row_p, col_p)), shape=(dim,dim), dtype=complex)
    return H + Hh + Hp

N, NA = 12, 4
t_max = 25.0
n_t   = 250
t_arr = np.linspace(0, t_max, n_t)

# Test anisotropy values
gammas     = [0.0, 0.25, 0.5, 0.75, 1.0]
gamma_names= ['TFIM','γ=0.25','γ=0.5','γ=0.75','XX']
# g_c(gamma) = J*sqrt(1-gamma)/2
g_c_list   = [0.5*np.sqrt(1-gm) for gm in gammas]
g_f        = 2.0   # fixed post-quench field (well above all g_c values)
# Also test g_f = 0.42: above g_c for gamma=0.25,0.5,0.75,XX but below TFIM g_c=0.5
# This decouples QPT from DQPT occurrence
g_f_critical_test = 0.42
theta_pairs = [(0.4, 1.2), (0.8, 1.4)]

print(f'Supplement R4: Anisotropic XY chain  N={N}, g_f={g_f}')
print(f'{"Model":10} {"gamma":6} {"g_c":6} {"g_f/g_c":8}')
for name, gm, gc in zip(gamma_names, gammas, g_c_list):
    print(f'  {name:10} {gm:6.2f} {gc:6.4f} {g_f/gc if gc>0 else "inf":>8}')
print(flush=True)

results = {'params': {'N':N,'NA':NA,'t_max':t_max,'n_t':n_t,'g_f':g_f},
           't': t_arr.tolist()}

for name, gm, gc in zip(gamma_names, gammas, g_c_list):
    print(f'\n--- {name} (gamma={gm}, g_c={gc:.4f}) ---', flush=True)
    H = build_XY_sparse(N, J=1.0, gamma=float(gm), g=g_f)

    # Check Hermitian
    herm_err = float(np.abs((H-H.conj().T).toarray()).max())
    print(f'  ||H-H†|| = {herm_err:.1e}', flush=True)

    # DQPT
    _,_,_,echo_ref = run_ea(N, NA, 0.4, H, t_arr)
    rate_ref = -np.log(np.maximum(echo_ref, 1e-300))/N
    dqpt_raw = find_dqpt_times(t_arr, echo_ref, N, threshold=0.03)
    dqpt_clean = []
    for x in sorted(dqpt_raw):
        if not dqpt_clean or x - dqpt_clean[-1] > 1.5:
            dqpt_clean.append(x)
    max_rate = float(rate_ref.max())
    has_dqpt = max_rate > 0.05 and len(dqpt_clean) >= 1
    print(f'  DQPT: {has_dqpt} (max λ={max_rate:.4f})', flush=True)
    if dqpt_clean:
        T_num = np.mean(np.diff(dqpt_clean)) if len(dqpt_clean)>1 else None
        print(f'  DQPT times: {[f"{x:.3f}" for x in dqpt_clean[:4]]}', flush=True)
        print(f'  T* numerical: {T_num:.4f}' if T_num else '  T*: single peak', flush=True)

    results[gm] = {'name':name,'gamma':float(gm),'g_c':gc,
                    'has_dqpt':has_dqpt,'max_losch_rate':max_rate,
                    'dqpt_times':dqpt_clean}

    for th1, th2 in theta_pairs:
        DS1,_,_,_ = run_ea(N, NA, th1, H, t_arr)
        DS2,_,_,_ = run_ea(N, NA, th2, H, t_arr)
        if DS1[0] > DS2[0]:
            DS1, DS2 = DS2, DS1
        cross = find_crossings(t_arr, DS1, DS2)
        tm = cross[0] if cross else None
        results[gm][f'qme_{th1}_{th2}'] = {
            'has_qme': tm is not None, 't_M': float(tm) if tm else None,
            'DS0_1':float(DS1[0]),'DS0_2':float(DS2[0])
        }
        status = f't_M={tm:.3f}' if tm else 'NO QME'
        print(f'  θ={th1}/{th2}: {status}  (ΔS(0): {DS1[0]:.4f}/{DS2[0]:.4f})', flush=True)

print('\n===== XY SUMMARY =====')
print(f'{"Model":10} {"g_c":6} {"DQPT":5} ' +
      '  '.join(f'QME({th1}/{th2})' for th1,th2 in theta_pairs))
for name, gm, gc in zip(gamma_names, gammas, g_c_list):
    d = results[gm]
    dqpt = 'YES' if d['has_dqpt'] else 'no '
    qme_strs = []
    for th1, th2 in theta_pairs:
        q = d.get(f'qme_{th1}_{th2}', {})
        qme_strs.append('YES' if q.get('has_qme') else 'no ')
    print(f'  {name:10} {gc:.4f}  {dqpt}  ' + '    '.join(qme_strs))

print('\n===== CRITICAL TEST: g_f=0.42 (below TFIM g_c=0.5, above XY g_c) =====')
print('  For TFIM (g_c=0.5): g_f < g_c → NO DQPT expected')
print('  For γ>0 (g_c<0.43): g_f > g_c → DQPT expected')
print('  This separates QPT boundary from DQPT occurrence.\n')
g_f2 = g_f_critical_test
results_crit = {}
for name, gm, gc in zip(gamma_names, gammas, g_c_list):
    print(f'  {name} (g_c={gc:.4f}, g_f={g_f2}): ', end='', flush=True)
    H2 = build_XY_sparse(N, J=1.0, gamma=float(gm), g=g_f2)
    herm2 = float(np.abs((H2-H2.conj().T).toarray()).max())
    _,_,_,echo2 = run_ea(N, NA, 0.4, H2, t_arr)
    rate2 = -np.log(np.maximum(echo2, 1e-300))/N
    max_r2 = float(rate2.max())
    dqpt2 = find_dqpt_times(t_arr, echo2, N, threshold=0.03)
    dqpt2c = []
    for x in sorted(dqpt2):
        if not dqpt2c or x - dqpt2c[-1] > 1.5:
            dqpt2c.append(x)
    has_dqpt2 = max_r2 > 0.05 and len(dqpt2c) >= 1
    DS1,_,_,_ = run_ea(N, NA, 0.4, H2, t_arr)
    DS2,_,_,_ = run_ea(N, NA, 1.2, H2, t_arr)
    if DS1[0] > DS2[0]: DS1, DS2 = DS2, DS1
    cross2 = find_crossings(t_arr, DS1, DS2)
    tm2 = cross2[0] if cross2 else None
    tm2_s = f't_M={tm2:.3f}' if tm2 else 'NO QME'
    results_crit[gm] = {'has_dqpt':has_dqpt2,'t_M':float(tm2) if tm2 else None,
                        'max_losch':max_r2,'dqpt_times':dqpt2c,'herm_err':herm2}
    print(f'||H-H†||={herm2:.1e}  DQPT={has_dqpt2}(λ={max_r2:.3f})  QME({tm2_s})',flush=True)

results['critical_test'] = {'g_f': g_f2, 'data': results_crit}

os.makedirs('results', exist_ok=True)
results_ser = {str(k): v for k, v in results.items()}
json.dump(results_ser, open('results/supp_r4_xy.json','w'), indent=2)
print('\nR4 done → results/supp_r4_xy.json')
