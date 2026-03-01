import sys
import numpy as np
import time
sys.path.insert(0, '/home/yanlin/wuli/src')
from scipy.linalg import eigh
from free_fermion_ea import build_TFIM_BdG

for N in [50, 100, 200, 500, 1000]:
    HBdG = build_TFIM_BdG(N, g=2.0)
    t0 = time.time()
    eps, V = eigh(HBdG)
    t1 = time.time()
    C0 = np.zeros((2*N, 2*N))
    C0[:N, :N] = 0.25 * np.eye(N)
    C0[N:, N:] = 0.75 * np.eye(N)
    phase = np.exp(-1j * eps * 3.0)
    VdagC0V = V.T @ C0 @ V
    C_rot = phase[:, None] * VdagC0V * phase[None, :].conj()
    C_t = V @ C_rot @ V.T
    t2 = time.time()
    print("N=%4d: eigh=%.3fs  evolve_1step=%.3fs" % (N, t1-t0, t2-t1))
