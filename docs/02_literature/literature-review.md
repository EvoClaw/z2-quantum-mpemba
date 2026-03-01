# Literature Review: Quantum 1D Models (Phase 1)

## Deep Reading Summary

### P1: arXiv:2512.03341 (Huang, Lin, Igloi 2025)
Problem: Quench dynamics in dimerized XXZ chain (flat-band limit)
Method: Bell basis expansion, exact analytical computation, IBM-Q quantum simulation
Key Results: Closed-form EE and Loschmidt echo for arbitrary even N; periodicity conditions on Delta; Loschmidt zeros identified; IBM-Q agreement up to N=12
Gaps: Only flat-band limit; only sudden quench (finite-rate explicitly listed as future work); no negativity; no entanglement asymmetry

### P2: arXiv:2510.04326 (Chandra, Rahaman et al. 2025)
Problem: Integrable Floquet time crystals in 1D without disorder
Method: BdG free-fermion, Floquet theory, QuTiP, RANSAC finite-size scaling
Key Results: Novel DTC phase in 1D; phase diagram DTC/FTC/OSL/PM; analytically obtained gapless loci; DTC lifetime diverges exponentially with N
CRITICAL GAP: NO entanglement entropy studied; no entanglement spectrum; no information scrambling

### P3: arXiv:2602.12185 (Feb 2026)
Problem: Symmetry-resolved EE from Ballistic Fluctuation Theory
Method: BFT, twist fields, CFT
Key Results: Analytic charged Renyi entropies for free-fermion at equilibrium and out-of-equilibrium
Gaps: Only free-fermion; Floquet systems not covered

### P4: arXiv:2408.00474 (2024)
Problem: GHD in integrable quantum circuits (Trotterized XXZ)
Method: GHD + Bethe ansatz
Key Results: Different nonequilibrium dynamics for large Trotter steps

## Field Landscape

### Hot Subareas (2024-2026)
1. Entanglement asymmetry and quantum Mpemba effect (QME) - fastest growing
2. Integrable Floquet time crystals - disorder-free DTC mechanisms
3. Symmetry-resolved entanglement entropy - exact results emerging
4. Generalized hydrodynamics - maturing
5. Krylov complexity - new chaos diagnostic
6. Quantum many-body scars - beyond PXP

### Key Controversy
- MBL in 1D: Stability debated (thermal avalanche problem)
- DTC lifetime: algebraic in integrable models vs exponential in MBL
- QME for Z2 symmetry: random circuits show NO QME; U(1) shows QME; integrable models unclear

### Computation Estimates
- Free-fermion correlation matrix: N~10000, seconds — exact
- Full ED: N~28 sites, seconds — exact
- DMRG: N~200 sites, minutes — near-exact
- Bethe ansatz / CFT: thermodynamic limit — exact
