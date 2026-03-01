"""共享的精确 EA 计算函数（被三个实验脚本 import）"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eigvalsh

def build_TFIM_sparse(N, J=1.0, g=2.0, pbc=True):
    dim = 2**N
    rows, cols, vals = [], [], []
    for s in range(dim):
        n = [(s >> (N-1-j)) & 1 for j in range(N)]
        diag_val = sum(-g*(1-2*n[j]) for j in range(N))
        rows.append(s); cols.append(s); vals.append(diag_val)
        for j in range(N):
            j1 = (j+1)%N if pbc else (j+1 if j+1<N else None)
            if j1 is None: continue
            jw = lambda nb, s: (-1)**sum(nb[:s])
            if n[j1]==1 and n[j]==0:
                n2=n.copy(); n2[j1]=0; sg=jw(n,j1)*jw(n2,j)
                rows.append((s^(1<<(N-1-j1)))|(1<<(N-1-j))); cols.append(s); vals.append(-J/2*sg)
            if n[j]==1 and n[j1]==0:
                n2=n.copy(); n2[j]=0; sg=jw(n,j)*jw(n2,j1)
                rows.append((s^(1<<(N-1-j)))|(1<<(N-1-j1))); cols.append(s); vals.append(-J/2*sg)
            if n[j]==1 and n[j1]==1:
                n2=n.copy(); n2[j1]=0; sg=jw(n,j1)*jw(n2,j)  # FIX: 先湮灭j1
                rows.append((s^(1<<(N-1-j)))^(1<<(N-1-j1))); cols.append(s); vals.append(-J/2*sg)
            if n[j]==0 and n[j1]==0:
                n2=n.copy(); n2[j]=1; sg=jw(n,j)*jw(n2,j1)
                rows.append(s|(1<<(N-1-j))|(1<<(N-1-j1))); cols.append(s); vals.append(-J/2*sg)
    return sp.csr_matrix((vals,(rows,cols)),shape=(dim,dim))

def build_psi0(N, theta):
    dim=2**N; c,s=np.cos(theta/2),np.sin(theta/2)
    psi=np.zeros(dim,dtype=complex)
    for st in range(dim):
        amp=1.0
        for j in range(N):
            amp *= s if (st>>(N-1-j))&1 else c
        psi[st]=amp
    return psi

def partial_trace(psi, N, NA):
    pm=psi.reshape(2**NA,2**(N-NA))
    return pm @ pm.conj().T

def ea_from_rho(rho_A, NA):
    Q=np.diag(np.array([(-1)**bin(k).count('1') for k in range(2**NA)],dtype=float))
    rt=(rho_A+Q@rho_A@Q)/2
    def S(r):
        ev=eigvalsh(r).real; ev=ev[ev>1e-14]
        return -np.sum(ev*np.log(ev))
    return S(rt)-S(rho_A), S(rho_A), float(np.real(np.trace(rho_A@Q)))

def run_ea(N, NA, theta, H, t_arr):
    psi0=build_psi0(N,theta); psi_t=psi0.copy(); t_prev=0.0
    DS,Srho,Ppi,echo=[],[],[],[]
    for t in t_arr:
        dt=t-t_prev
        if dt>1e-12: psi_t=spla.expm_multiply(-1j*H*dt,psi_t)
        t_prev=t
        rho_A=partial_trace(psi_t,N,NA)
        d,s,p=ea_from_rho(rho_A,NA)
        e=abs(np.dot(psi0.conj(),psi_t))**2
        DS.append(d); Srho.append(s); Ppi.append(p); echo.append(e)
    return np.array(DS),np.array(Srho),np.array(Ppi),np.array(echo)

def find_crossings(t_arr, DS1, DS2):
    """找 DS2 > DS1 变为 DS2 < DS1 的时间（Mpemba 交叉）"""
    diff = DS2 - DS1
    crossings = []
    for k in range(len(diff)-1):
        if diff[k]>0 and diff[k+1]<0:
            tm = t_arr[k] + diff[k]/(diff[k]-diff[k+1])*(t_arr[k+1]-t_arr[k])
            crossings.append(tm)
    return crossings

def find_dqpt_times(t_arr, echo, N, threshold=0.05):
    losch = -np.log(np.maximum(echo,1e-300))/N
    times=[]
    for k in range(1,len(losch)-1):
        if losch[k]>losch[k-1] and losch[k]>losch[k+1] and losch[k]>threshold:
            times.append(t_arr[k])
    return times
